/*========================================================================
 *
 * RICIANDECONV3MX.C  3D TV minimization for Rician deconvolution
 *
 * u = riciandeconv3mx(f,K,sigma,lambda,NumIter,dt) performs deconvolution
 * on a 3D volume f with Rician noise with parameter sigma.  The deblurred
 * image u is found as the minimizer of 
 *
 *         /                      / [ (K*u)^2 + f^2       (K*u)f   ]
 *    min  | |grad u| dx + lambda | [ --------- - log I0( ------ ) ] dx.
 *     u   /                      / [ 2 sigma^2          sigma^2   ]
 *
 * Parameter lambda >= 0 determines the strength of the denoising: smaller
 * lambda implies stronger denoising.  NumIter specifies the number of 
 * iterations and dt specifies the timestep.  Parameter dt must be
 * sufficiently small for stability (typically on the order of 0.0001).
 *
 * Pascal Getreuer 2009
 *
 *======================================================================*/

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "papi.h"

/* 3D Gaussian convolution using the method of Alvarez and Mazorra
 * Pascal Getreuer 2009
 */

#define GAUSSIAN_NUMSTEPS  3

/* GaussianBlur: Implements 3D Gaussian convolution with recursive (IIR)
 filters.  The parameter GAUSSIAN_NUMSTEPS determines the quality of the
 convolution.  More steps is a better approximation of the true Gaussian.*/
void GaussianBlur(double *u, const int Size[3], double Ksigma)
{
	const int PlaneStep = Size[0] * Size[1];
	double *uPtr, *uCopy, *uEnd;
	double lambda = (Ksigma * Ksigma) / (2.0 * GAUSSIAN_NUMSTEPS);
	double nu =
	    (1.0 + 2.0 * lambda - sqrt(1.0 + 4.0 * lambda)) / (2.0 * lambda);
	double BoundaryScale = 1.0 / (1.0 - nu);
	double PostScale = pow(nu / lambda, 3 * GAUSSIAN_NUMSTEPS);
	int Step = GAUSSIAN_NUMSTEPS;
	int n[3];
	uEnd = u + PlaneStep * Size[2];

	do {
		for (n[2] = 0, uPtr = u; n[2] < Size[2];
		     ++n[2], uPtr += PlaneStep) {
			uCopy = uPtr;
			for (n[1] = 0; n[1] < Size[1]; ++n[1], uPtr += Size[0]) {

				/* Filter downwards */
				uPtr[0] *= BoundaryScale;
				++uPtr;
				for (n[0] = 1; n[0] < Size[0]; ++n[0], ++uPtr) {
					uPtr[0] += nu * uPtr[-1];
				}

				/* Filter upwards */
				--uPtr;
				uPtr[0] *= BoundaryScale;
				--uPtr;
				for (n[0] = Size[0] - 2; n[0] >= 0;
				     --n[0], --uPtr) {
					uPtr[0] += nu * uPtr[1];
				}
				++uPtr;
			}
			uPtr = uCopy;

			/* Filter right */
			for (n[0] = 0; n[0] < Size[0]; ++n[0], ++uPtr) {
				uPtr[0] *= BoundaryScale;
			}
			for (n[1] = 1; n[1] < Size[1]; ++n[1]) {
				for (n[0] = 0; n[0] < Size[0]; ++n[0], ++uPtr) {
					uPtr[0] += nu * uPtr[-Size[0]];
				}
			}
			--uPtr;

			/* Filter left */
			for (n[0] = Size[0] - 1; n[0] >= 0; --n[0], --uPtr) {
				uPtr[0] *= BoundaryScale;
			}
			for (n[1] = Size[1] - 2; n[1] >= 0; --n[1]) {
				for (n[0] = Size[0] - 1; n[0] >= 0;
				     --n[0], --uPtr) {
					uPtr[0] += nu * uPtr[Size[0]];
				}
			}
			++uPtr;
		}

		/* Filter out */
		n[0] = PlaneStep;
		uPtr = u;

		do {
			uPtr[0] *= BoundaryScale;
			++uPtr;
		} while (--n[0]);
		for (n[2] = 1; n[2] < Size[2]; ++n[2]) {
			n[0] = PlaneStep;

			do {
				uPtr[0] += nu * uPtr[-PlaneStep];
				++uPtr;
			} while (--n[0]);
		}

		/* Filter in */
		n[0] = PlaneStep;

		do {
			--uPtr;
			uPtr[0] *= BoundaryScale;
		} while (--n[0]);
		for (n[2] = Size[2] - 2; n[2] >= 0; --n[2]) {
			n[0] = PlaneStep;

			do {
				--uPtr;
				uPtr[0] += nu * uPtr[PlaneStep];
			} while (--n[0]);
		}
	} while (--Step);

	do {
		u[0] *= PostScale;
	} while (++u < uEnd);
}

/* Method Parameters */
#define DEFAULT_NUMITER  10
#define DEFAULT_DT       0.0001

#define DEFAULT_KSIGMA   1.8
#define DEFAULT_SIGMA    0.008
#define DEFAULT_LAMBDA   0.085

#define EPSILON          1.0E-10
#define EPSILON2         1.0E-5

/* Macro functions */
#define SQR(x) ((x)*(x))

/* Macros for referring to pixel neighbors */
#define CENTER   (n[0]+Size[0]*(n[1]+Size[1]*n[2]))
#define RIGHT    (n[0]+Size[0]*(n[1]+Size[1]*n[2])+Size[0])
#define LEFT     (n[0]+Size[0]*(n[1]+Size[1]*n[2])-Size[0])
#define DOWN     (n[0]+Size[0]*(n[1]+Size[1]*n[2])+1)
#define UP       (n[0]+Size[0]*(n[1]+Size[1]*n[2])-1)
#define ZOUT     (n[0]+Size[0]*(n[1]+Size[1]*n[2]+Size[1]))
#define ZIN      (n[0]+Size[0]*(n[1]+Size[1]*n[2]-Size[1]))
static void riciandeconv3(double *u, const double *f, const int Size[3],
			  double Ksigma, double sigma, double lambda,
			  int NumIter, double dt)
{
	const int NumEl = Size[0] * Size[1] * Size[2];
	double *g;		/* Array storing 1/|grad u| approximation */
	double *conv;		/* Array storing convolutions */
	double sigma2, gamma, r;
	int n[3];
	int Iter;

	/* Initializations */
	sigma2 = SQR(sigma);
	gamma = lambda / sigma2;

	/* Allocate temporary work arrays */
	g = calloc(NumEl, sizeof(double));
	conv = calloc(NumEl, sizeof(double));

    /*** Main gradient descent loop ***/
	for (Iter = 1; Iter <= NumIter; Iter++) {

		/* Approximate g = 1/|grad u| */
		for (n[2] = 1; n[2] < Size[2] - 1; ++n[2])
			for (n[1] = 1; n[1] < Size[1] - 1; ++n[1])
				for (n[0] = 1; n[0] < Size[0] - 1; ++n[0])
					g[CENTER] = 1.0 / sqrt(EPSILON
							       + SQR(u[CENTER] -
								     u[RIGHT])
							       + SQR(u[CENTER] -
								     u[LEFT])
							       + SQR(u[CENTER] -
								     u[DOWN])
							       + SQR(u[CENTER] -
								     u[UP])
							       + SQR(u[CENTER] -
								     u[ZOUT])
							       + SQR(u[CENTER] -
								     u[ZIN]));
		memcpy(conv, u, NumEl * sizeof(double));
		GaussianBlur(conv, Size, Ksigma);
		for (n[2] = 0; n[2] < Size[2]; ++n[2])
			for (n[1] = 0; n[1] < Size[1]; ++n[1])
				for (n[0] = 0; n[0] < Size[0]; ++n[0]) {

					/* Evaluate r = I1((K*u)f/sigma^2) / I0((K*u)f/sigma^2) with
					   a cubic rational approximation. */
					r = conv[CENTER] * f[CENTER] / sigma2;
					r = (r * (2.38944 + r * (0.950037 + r)))
					    / (4.65314 +
					       r * (2.57541 +
						    r * (1.48937 + r)));
					conv[CENTER] -= f[CENTER] * r;
				}
		GaussianBlur(conv, Size, Ksigma);

		/* Update u by a sem-implict step */
		for (n[2] = 1; n[2] < Size[2] - 1; n[2]++)
			for (n[1] = 1; n[1] < Size[1] - 1; n[1]++)
				for (n[0] = 1; n[0] < Size[0] - 1; n[0]++) {
					u[CENTER] =
					    (u[CENTER] +
					     dt * (u[RIGHT] * g[RIGHT]
						   + u[LEFT] * g[LEFT] +
						   u[DOWN] * g[DOWN] +
						   u[UP] * g[UP] +
						   u[ZOUT] * g[ZOUT] +
						   u[ZIN] * g[ZIN] -
						   gamma * conv[CENTER])) /
					    (1.0 +
					     dt * (g[RIGHT] + g[LEFT] +
						   g[DOWN] + g[UP] + g[ZOUT] +
						   g[ZIN]));
				}
		printf("Iter: %d\n", Iter);
	}

	/* Free temporary arrays */
	free(conv);
	free(g);
	return;
}

int main(int argc, char *argv[])
{
	if (argc < 6) {
		printf("riciandeblur3 M N P inputfile outputfile\r\n");
		exit(0);
	}
	int M = atoi(argv[1]);
	int N = atoi(argv[2]);
	int P = atoi(argv[3]);
	int p, n, m;
	FILE *inputfile = fopen(argv[4], "r");
	FILE *outputfile = fopen(argv[5], "w");
	double *f, *u;
	f = calloc(M * N * P, sizeof(double));
	u = calloc(M * N * P, sizeof(double));
	fread(f, sizeof(double), M * N * P, inputfile);

	/* Initialize u = f */
	memcpy(u, f, sizeof(double) * M * N * P);

	/*
	   for(p = 0; p < P; p++)
	   for(n = 0; n < N; n++)
	   for(m = 0; m < M; m++)
	   f[CENTER]=round(f[CENTER]);
	 */

	/* Set up parameters */
	double Ksigma = DEFAULT_KSIGMA;
	double sigma = DEFAULT_SIGMA;
	double lambda = DEFAULT_LAMBDA;
	int numIter = DEFAULT_NUMITER;
	int dt = DEFAULT_DT;
	int size[3];
	size[0] = M;
	size[1] = N;
	size[2] = P;

	/* Call the main denoising routine */
	unsigned g_batchid = 0;
	if (argc == 7) {
		g_batchid = atoi(argv[6]);
	}
	int Events[5];
	u_long_long papi_values[5];
	util_start_papi(g_batchid, Events);
	riciandeconv3(u, f, size, Ksigma, sigma, lambda, numIter, dt);
	util_stop_papi(g_batchid, papi_values);
	util_print_papi(g_batchid, papi_values, (g_batchid == 0));
	fwrite(u, sizeof(double), M * N * P, outputfile);
	printf("Finished\r\n");
}

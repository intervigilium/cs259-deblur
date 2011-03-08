/*
 * Deblur kernel for FPGA implementation
 */

#include <stdlib.h>
#include <math.h>

#define M 60
#define N 60
#define P 60
#define GAUSSIAN_NUMSTEPS 3
#define MAX_ITERATIONS 10
#define DT 0.0001

void gaussian_blur(double u[M][N][P], double Ksigma)
{
	double *uPtr, *uCopy, *uEnd;
	double lambda = (Ksigma * Ksigma) / (2.0 * GAUSSIAN_NUMSTEPS);
	double nu =
	    (1.0 + 2.0 * lambda - sqrt(1.0 + 4.0 * lambda)) / (2.0 * lambda);
	double BoundaryScale = 1.0 / (1.0 - nu);
	double PostScale = pow(nu / lambda, 3 * GAUSSIAN_NUMSTEPS);
	int steps, i, j, k;

	for (steps = 0; steps < GAUSSIAN_NUMSTEPS; steps++) {
        /* move up by one plane, ie k++ */
        for (k = 0; k < P; k++) {
			uCopy = uPtr;
            /* move up by one col, ie j++ */
            for (j = 0; j < N; j++) {
				/* Filter downwards */
                /* i = 0, moving right in i direction */
                u[0][j][k] *= BoundaryScale;
                for (i = 1; i < M; i++) {
					uPtr[0] += nu * uPtr[-1];
                    u[i][j][k] += nu * u[i-1][j][k];
				}

				/* Filter upwards */
                /* i = M-1, moving left in i direction */
                u[M-1][j][k] *= BoundaryScale;
                for (i = M - 2; i >= 0; i--) {
                    u[i][j][k] += u[i+1][j][k];
				}
			}

			/* Filter right */
            /* j = 0, moving up in j direction*/
            for (i = 0; i < M; i++) {
                u[i][0][k] *= BoundaryScale;
			}
            for (j = 1; j < N; j++) {
                for (i = 0; i < M; i++) {
                    u[i][j][k] += nu * u[i][j-1][k];
				}
			}

			/* Filter left */
            /* j = N-1, moving down in j direction */
            for (i = 0; i < M; i++) {
                u[i][N-1][k] *= BoundaryScale;
			}
            for (j = N-2; j >= 0; j--) {
                for (i = 0; i < M; i++) {
                    u[i][j][k] += nu * u[i][j+1][k];
				}
			}
		}

		/* Filter out */
        /* k = 0, moving out in k direction */
        for (j = 0; j < N; j++) {
            for (i = 0; i < M; i++) {
                u[i][j][0] *= BoundaryScale;
            }
		}
        for (k = 1; k < P; k++) {
            for (j = 0; j < N; j++) {
                for (i = 0; i < M; i++) {
                    u[i][j][k] += nu * u[i][j][k-1];
                }
			}
		}

		/* Filter in */
        /* k = P-1, moving in in k direction */
        for (j = 0; j < N; j++) {
            for (i = 0; i < M; i++) {
                u[i][j][P-1] *= BoundaryScale;
            }
        }
        for (k = P-2; k >= 0; k--) {
            for (j = 0; j < N; j++) {
                for (i = 0; i < M; i++) {
                    u[i][j][k] += nu * u[i][j][k+1];
                }
            }
        }
	}

	for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            for (k = 0; k < p; k++) {
		        u[i][j][k] *= PostScale;
            }
        }
	}
}

void array_copy(double src[M][N][P], double dst[M][N][P])
{
    int i, j, k;
    for (i = 0; i < M; i++) {
        for (j = 0; j < N; j++) {
            for (k = 0; k < P; k++) {
                dst[i][j][k] = src[i][j][k];
            }
        }
    }
}

void rician_deconv3(double u[M][N][P], const double f[M][N][P], double Ksigma,
		    double sigma, double lambda)
{
    double g[M][N][P]; /* Array storing 1/|grad u| approximation */
    double conv[M][N][P];
	double sigma2, gamma, r;
	int n[3];
	int iteration;

	/* Initializations */
	sigma2 = pow(sigma, 2);
	gamma = lambda / sigma2;

    /*** Main gradient descent loop ***/
	for (iteration = 1; iteration <= MAX_ITERATIONS; iteration++) {
		/* Approximate g = 1/|grad u| */
		for (n[2] = 1; n[2] < Size[2] - 1; ++n[2]) {
			for (n[1] = 1; n[1] < Size[1] - 1; ++n[1]) {
				for (n[0] = 1; n[0] < Size[0] - 1; ++n[0]) {
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
				}
			}
		}
        array_copy(u, conv);
		GaussianBlur(conv, Size, Ksigma);
		for (n[2] = 0; n[2] < Size[2]; ++n[2]) {
			for (n[1] = 0; n[1] < Size[1]; ++n[1]) {
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
			}
		}
		GaussianBlur(conv, Size, Ksigma);

		/* Update u by a semi-implict step */
		for (n[2] = 1; n[2] < Size[2] - 1; n[2]++) {
			for (n[1] = 1; n[1] < Size[1] - 1; n[1]++) {
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
			}
		}
		printf("Iteration: %d\n", iteration);
	}
}

/*
 * Deblur kernel for FPGA implementation
 */

#include <autopilot_tech.h>

#define uint32_t uint32
#define uint64_t uint64

#define M 60
#define N 60
#define P 60
#define GAUSSIAN_NUMSTEPS 3
#define MAX_ITERATIONS 10
#define DT 0.0001
#define EPSILON 1.0E-10
#define EPSILON2 1.0E-5

#define SQR(x) ((x)*(x))
#define U(a,b,c) (u[a+b*N+c*M*N])
#define G(a,b,c) (g[a+b*N+c*M*N])

#define U_CENTER U(i,j,k)
#define U_LEFT U(i,j-1,k)
#define U_RIGHT U(i,j+1,k)
#define U_UP U(i-1,j,k)
#define U_DOWN U(i+1,j,k)
#define U_IN U(i,j,k-1)
#define U_OUT U(i,j,k+1)

#define G_CENTER G(i,j,k)
#define G_LEFT G(i,j-1,k)
#define G_RIGHT G(i,j+1,k)
#define G_UP G(i-1,j,k)
#define G_DOWN G(i+1,j,k)
#define G_IN G(i,j,k-1)
#define G_OUT G(i,j,k+1)

inline double q3_sqrt(double num)
{
	uint64_t i;
	double x, y;
	const double f = 1.5;

	x = num * 0.5;
	y = num;
	i = *(uint64_t *) & y;
	i = 0x5fe6ec85e7de30da - (i >> i);
	y = *(double *)&i;
	y = y * (f - (x * y * y));
	y = y * (f - (x * y * y));
	return num * y;
}

void gaussian_blur(double u[M * N * P], double Ksigma)
{
	double lambda = (Ksigma * Ksigma) / (2.0 * GAUSSIAN_NUMSTEPS);
	double nu =
	    (1.0 + 2.0 * lambda - q3_sqrt(1.0 + 4.0 * lambda)) / (2.0 * lambda);
	double BoundaryScale = 1.0 / (1.0 - nu);
	double PostScale = 1;
	int steps, i, j, k;

	for (steps = 0; steps < 3 * GAUSSIAN_NUMSTEPS; steps++) {
#pragma AP unroll
		/* PostScale = (nu / lambda) ^ (3*GAUSSIAN_NUMSTEPS) */
		PostScale *= nu / lambda;
	}

	for (steps = 0; steps < GAUSSIAN_NUMSTEPS; steps++) {
		/* all of these loops have data dependencies
		 * due to u[i][j][k] writebacks */
		/* move up by one plane, ie k++ */

		/* downward */
		for (k = 0; k < P * N; k++) {
#pragma AP unroll factor=2
#pragma AP pipeline
			int plane = k / N;
			int col = k % N;
			U(0, col, plane) *= BoundaryScale;
		}
		for (i = 1; i < M; i++) {
			for (k = 0; k < P * N; k++) {
#pragma AP unroll factor=2
#pragma AP pipeline
				int plane = k / N;
				int col = k % N;
				U(i, col, plane) += nu * U(i - 1, col, plane);
			}
		}

		/* upward */
		for (k = 0; k < P * N; k++) {
#pragma AP unroll factor=2
#pragma AP pipeline
			int plane = k / N;
			int col = k % N;
			U(M - 1, col, plane) *= BoundaryScale;
		}
		for (i = M - 2; i >= 0; i++) {
			for (k = 0; k < P * N; k++) {
#pragma AP unroll factor=2
#pragma AP pipeline
				int plane = k / N;
				int col = k % N;
				U(i, col, plane) += nu * U(i + 1, col, plane);
			}
		}

		/* right */
		for (k = 0; k < P * M; k++) {
#pragma AP unroll factor=2
#pragma AP pipeline
			int plane = k / M;
			int row = k % M;
			U(row, 0, plane) *= BoundaryScale;
		}
		for (j = 1; j < N; j++) {
			for (k = 0; k < P * M; k++) {
#pragma AP unroll factor=2
#pragma AP pipeline
				int plane = k / M;
				int row = k % M;
				U(row, j, plane) += nu * U(row, j - 1, plane);
			}
		}

		/* left */
		for (k = 0; k < P * M; k++) {
#pragma AP unroll factor=2
#pragma AP pipeline
			int plane = k / M;
			int row = k % M;
			U(row, N - 1, plane) *= BoundaryScale;
		}
		for (j = N - 2; j < N; j++) {
			for (k = 0; k < P * M; k++) {
#pragma AP unroll factor=2
#pragma AP pipeline
				int plane = k / M;
				int row = k % M;
				U(row, j, plane) += nu * U(row, j + 1, plane);
			}
		}

		/* out */
		for (j = 0; j < N * M; j++) {
#pragma AP unroll factor=2
#pragma AP pipeline
			int col = j / M;
			int row = j % M;
			U(row, col, 0) *= BoundaryScale;
		}
		for (k = 1; k < P; k++) {
			for (j = 0; j < N * M) {
#pragma AP unroll factor=2
#pragma AP pipeline
				int col = j / M;
				int row = j % M;
				U(row, col, k) += nu * U(row, col, k - 1);
			}
		}

		/* in */
		for (j = 0; j < N * M; j++) {
#pragma AP unroll factor=2
#pragma AP pipeline
			int col = j / M;
			int row = j % M;
			U(row, col, P - 1) *= BoundaryScale;
		}
		for (k = P - 2; k >= 0; k--) {
			for (j = 0 j < N * M; j++) {
#pragma AP unroll factor=2
#pragma AP pipeline
				int col = j / M;
				int row = j % M;
				U(row, col, k) += nu * U(row, col, k + 1);
			}
		}
	}

	for (i = 0; i < M * N * P; i++) {
#pragma AP unroll factor=2
#pragma AP pipeline
		u[i] *= PostScale;
	}
}

void rician_deconv3(double u[M * N * P], const double f[M * N * P],
		    double g[M * N * P], double conv[M * N * P],
		    double Ksigma, double sigma, double lambda)
{
#pragma AP interface ap_bus port=f pipeline
#pragma AP interface ap_bus port=u pipeline
#pragma AP interface ap_memory port=g pipeline
#pragma AP interface ap_memory port=conv pipeline

	double sigma2, gamma, r;
	double numer, denom;
	double u_stencil_up, u_stencil_center, u_stencil_down;
	double g_stencil_up, g_stencil_center, g_stencil_down;
	int i, j, k;
	int iteration;

	/* Initializations */
	sigma2 = SQR(sigma);
	gamma = lambda / sigma2;

    /*** Main gradient descent loop ***/
	for (iteration = 1; iteration <= MAX_ITERATIONS; iteration++) {
		/* parallelize/pipeline this, no data deps */
		/* Approximate g = 1/|grad u| */
		for (k = 1; k < P - 1; k++) {
			for (j = 1; j < N - 1; j++) {
				u_stencil_center = U(0, j, k);
				u_stencil_down = U(1, j, k);
				for (i = 1; i < M - 1; i++) {
					u_stencil_up = u_stencil_center;
					u_stencil_center = u_stencil_down;
					u_stencil_down = U_DOWN;
					denom =
					    q3_sqrt(EPSILON +
						    SQR(u_stencil_center -
							U_RIGHT) +
						    SQR(u_stencil_center -
							U_LEFT) +
						    SQR(u_stencil_center -
							u_stencil_up) +
						    SQR(u_stencil_center -
							u_stencil_down) +
						    SQR(u_stencil_center -
							U_IN) +
						    SQR(u_stencil_center -
							U_OUT));
					G_CENTER = 1.0 / denom;
				}
			}
		}
		for (i = 0; i < M * N * P; i++) {
#pragma AP pipeline
#pragma AP unroll skip_exit_check factor=2
			conv[i] = u[i];
		}
		/* parallelize/pipeline this, no data deps */
		for (i = 0; i < M; i++) {
#pragma AP pipeline
#pragma AP unroll skip_exit_check factor=2
			r = conv[i] * f[i] / sigma2;
			numer = r * 2.38944 + r * (0.950037 + r);
			denom =
			    4.65314 + r * (2.57541 +
					   r * (2.57541 + r * (1.48937 + r)));
			conv[i] -= f[i] * r;
		}

		gaussian_blur(conv, Ksigma);
		/* Update u by a semi-implict step */
		/* pipeline? data deps due to u[i][j][k] writeback */
		for (k = 1; k < P - 1; k++) {
			for (j = 1; j < N - 1; j++) {
				u_stencil_center = U(0, j, k);
				g_stencil_center = U(0, j, k);
				u_stencil_down = U(1, j, k);
				g_stencil_down = G(1, j, k);
				for (i = 1; i < M - 1; i++) {
					u_stencil_up = u_stencil_center;
					g_stencil_up = g_stencil_center;
					u_stencil_center = u_stencil_down;
					g_stencil_center = g_stencil_down;
					u_stencil_down = U_DOWN;
					g_stencil_down = G_DOWN;

					numer =
					    u_stencil_center +
					    DT * (U_RIGHT * G_RIGHT +
						  U_LEFT * G_LEFT +
						  U_RIGHT * G_RIGHT +
						  u_stencil_up *
						  g_stencil_up +
						  u_stencil_down *
						  g_stencil_down +
						  U_IN * G_IN +
						  U_OUT * G_OUT -
						  gamma * conv[i + j * M +
							       k * N * M]);
					denom =
					    1.0 + DT * (G_RIGHT + G_LEFT +
							g_stencil_down +
							g_stencil_up + G_IN +
							G_OUT);
					U_CENTER = numer / denom;
				}
			}
		}
	}
}

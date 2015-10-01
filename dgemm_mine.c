const char* dgemm_desc = "My awesome dgemm.";

#ifndef BLOCK_SIZE
#define BLOCK_SIZE ((int) 512)
#endif

#include<immintrin.h>
#include<stdlib.h>

void basic_dgemm(const int lda, const int M, const int N, const int K,
		const double *A, const double *B, double *C) 
{
	int i, j, k;
	const int col_reduced_4 = M - M % 4;
	const int col_reduced_16 = M - M % 16;
	__m256d ymm0, ymm1, ymm2, ymm3, ymm4, ymm5, ymm6, ymm7, ymm8;
	double *b = _mm_malloc(K * sizeof(double), 64);
	double *c = _mm_malloc(M * sizeof(double), 64);
	for (j = 0; j < N; ++j) {
		const int index1 = j * lda;
		for (k = 0; k < K; k++) {
			b[k] = B[index1 + k];
		}
		for (i = 0; i < M; i++) {
			c[i] = C[index1 + i];
		}
		for (k = 0; k < K; ++k) {
			const int index2 = k * lda;
			ymm0 = _mm256_broadcast_sd(&b[k]);
			for (i = 0; i < col_reduced_16; i += 16) {
				ymm1 = _mm256_loadu_pd(&A[index2 + i]);
				ymm2 = _mm256_loadu_pd(&A[index2 + i + 4]);
				ymm3 = _mm256_loadu_pd(&A[index2 + i + 8]);
				ymm4 = _mm256_loadu_pd(&A[index2 + i + 12]);

				ymm5 = _mm256_load_pd(&c[i]);
				ymm6 = _mm256_load_pd(&c[i + 4]);
				ymm7 = _mm256_load_pd(&c[i + 8]);
				ymm8 = _mm256_load_pd(&c[i + 12]);

				_mm256_store_pd(&c[i], _mm256_fmadd_pd(ymm1, ymm0, ymm5));
				_mm256_store_pd(&c[i + 4], _mm256_fmadd_pd(ymm2, ymm0, ymm6));
				_mm256_store_pd(&c[i + 8], _mm256_fmadd_pd(ymm3, ymm0, ymm7));
				_mm256_store_pd(&c[i + 12], _mm256_fmadd_pd(ymm4, ymm0, ymm8));
			}
			for (i = col_reduced_16; i < M; i++) {
				c[i] += A[index2 + i] * b[k];
			}
		}
		for (i = 0; i < col_reduced_4; i += 4) {
			_mm256_storeu_pd(&C[index1 + i], _mm256_loadu_pd(&c[i]));
		}
		for (i = col_reduced_4; i < M; i++) {
			C[index1 + i] = c[i];
		}
	}   
}

/*


	int i,j,k,s;
	double *D = (double*) malloc(M * M * sizeof(double));
	for (i = 0; i < M; i++) {
		for (j = 0; j < M; j++) {
			D[i * M + j] = A[j * M + i]; 
		}
	}   
	for (i = 0; i < M; i++) {
		const int col_reduced = M - M%32;
		const int col_reduced_16 = M - M%16;
		__m256d ymm0, ymm1, ymm2, ymm3, ymm4, ymm5, ymm6, ymm7,
				ymm8, ymm9, ymm10, ymm11, ymm12, ymm13, ymm14, ymm15;
		const double *x = _mm_malloc(M * sizeof(double), 64); 
		x = &B[i * M];
		for (j = 0; j < M; j ++) {
			double res = 0;
			const int index = j * M;
			for (k = 0; k < col_reduced; k += 32) {
				ymm8 = _mm256_load_pd(&x[k]);
				ymm9 = _mm256_load_pd(&x[k + 4]);
				ymm10 = _mm256_load_pd(&x[k + 8]);
				ymm11 = _mm256_load_pd(&x[k + 12]);
				ymm12 = _mm256_load_pd(&x[k + 16]);
				ymm13 = _mm256_load_pd(&x[k + 20]);
				ymm14 = _mm256_load_pd(&x[k + 24]);
				ymm15 = _mm256_load_pd(&x[k + 28]);

				ymm0 = _mm256_loadu_pd(&D[index + k]);
				ymm1 = _mm256_loadu_pd(&D[index + k + 4]);
				ymm2 = _mm256_loadu_pd(&D[index + k + 8]);
				ymm3 = _mm256_loadu_pd(&D[index + k + 12]);
				ymm4 = _mm256_loadu_pd(&D[index + k + 16]);
				ymm5 = _mm256_loadu_pd(&D[index + k + 20]);
				ymm6 = _mm256_loadu_pd(&D[index + k + 24]);
				ymm7 = _mm256_loadu_pd(&D[index + k + 28]);

				ymm0 = _mm256_mul_pd(ymm0, ymm8 );
				ymm1 = _mm256_mul_pd(ymm1, ymm9 );
				ymm2 = _mm256_mul_pd(ymm2, ymm10);
				ymm3 = _mm256_mul_pd(ymm3, ymm11);
				ymm4 = _mm256_mul_pd(ymm4, ymm12);
				ymm5 = _mm256_mul_pd(ymm5, ymm13);
				ymm6 = _mm256_mul_pd(ymm6, ymm14);
				ymm7 = _mm256_mul_pd(ymm7, ymm15);

				ymm0 = _mm256_add_pd(ymm0, ymm1);
				ymm2 = _mm256_add_pd(ymm2, ymm3);
				ymm4 = _mm256_add_pd(ymm4, ymm5);
				ymm6 = _mm256_add_pd(ymm6, ymm7);
				ymm0 = _mm256_add_pd(ymm0, ymm2);
				ymm4 = _mm256_add_pd(ymm4, ymm6);
				ymm0 = _mm256_add_pd(ymm0, ymm4);

				_mm256_store_pd(scratchpad, ymm0);
				for (s = 0; s < 4; s++)
					res += scratchpad[s];
			}
			for (k = col_reduced; k < col_reduced_16; k += 16) {
				ymm8 = _mm256_load_pd(&x[k]);
				ymm9 = _mm256_load_pd(&x[k + 4]);
				ymm10 = _mm256_load_pd(&x[k + 8]);
				ymm11 = _mm256_load_pd(&x[k + 12]);

				ymm0 = _mm256_loadu_pd(&D[index + k]);
				ymm1 = _mm256_loadu_pd(&D[index + k + 4]);
				ymm2 = _mm256_loadu_pd(&D[index + k + 8]);
				ymm3 = _mm256_loadu_pd(&D[index + k + 12]);

				ymm0 = _mm256_mul_pd(ymm0, ymm8 );
				ymm1 = _mm256_mul_pd(ymm1, ymm9 );
				ymm2 = _mm256_mul_pd(ymm2, ymm10);
				ymm3 = _mm256_mul_pd(ymm3, ymm11);

				ymm0 = _mm256_add_pd(ymm0, ymm1);
				ymm2 = _mm256_add_pd(ymm2, ymm3);
				ymm0 = _mm256_add_pd(ymm0, ymm2);

				_mm256_store_pd(scratchpad, ymm0);
				for (s = 0; s < 4; s++)
					res += scratchpad[s];
			}
			for (k = col_reduced_16; k < M; k++) {
				res += D[index + k] * x[k];
			}
			C[j + i * M] = res;
		}
	}

*/


void do_block(const int lda,
		const double *A, const double *B, double *C, 
		const int i, const int j, const int k)
{
	const int M = (i+BLOCK_SIZE > lda? lda-i : BLOCK_SIZE);
	const int N = (j+BLOCK_SIZE > lda? lda-j : BLOCK_SIZE);
	const int K = (k+BLOCK_SIZE > lda? lda-k : BLOCK_SIZE);
	basic_dgemm(lda, M, N, K,
			A + i + k*lda, B + k + j*lda, C + i + j*lda);
}

void square_dgemm(const int M, const double *A, const double *B, double *C) 
{
	const int n_blocks = M / BLOCK_SIZE + (M%BLOCK_SIZE? 1 : 0); 
	int bi, bj, bk; 
	for (bi = 0; bi < n_blocks; ++bi) {
		const int i = bi * BLOCK_SIZE;
		for (bj = 0; bj < n_blocks; ++bj) {
			const int j = bj * BLOCK_SIZE;
			for (bk = 0; bk < n_blocks; ++bk) {
				const int k = bk * BLOCK_SIZE;
				do_block(M, A, B, C, i, j, k); 
			}
		}
	}   
}

/*
void square_dgemm(const int M, const double *A, const double *B, double *C) {
	double scratchpad[8];
	int i,j,k,s;
	double *D = (double*) malloc(M * M * sizeof(double));
	for (i = 0; i < M; i++) {
		for (j = 0; j < M; j++) {
			D[i * M + j] = A[j * M + i]; 
		}
	}   
	for (i = 0; i < M; i++) {
		const int col_reduced = M - M%32;
		const int col_reduced_16 = M - M%16;
		__m256d ymm0, ymm1, ymm2, ymm3, ymm4, ymm5, ymm6, ymm7,
				ymm8, ymm9, ymm10, ymm11, ymm12, ymm13, ymm14, ymm15;
		double *x; 
		const int index1 = i * M;
		for (j = 0; j < M; j ++) {
			double res = 0;
			const int index2 = j * M;
			for (k = 0; k < col_reduced; k += 32) {
				ymm8 = _mm256_loadu_pd(&B[index1 + k]);
				ymm9 = _mm256_loadu_pd(&B[index1 + k + 4]);
				ymm10 = _mm256_loadu_pd(&B[index1 + k + 8]);
				ymm11 = _mm256_loadu_pd(&B[index1 + k + 12]);
				ymm12 = _mm256_loadu_pd(&B[index1 + k + 16]);
				ymm13 = _mm256_loadu_pd(&B[index1 + k + 20]);
				ymm14 = _mm256_loadu_pd(&B[index1 + k + 24]);
				ymm15 = _mm256_loadu_pd(&B[index1 + k + 28]);

				ymm0 = _mm256_loadu_pd(&D[index2 + k]);
				ymm1 = _mm256_loadu_pd(&D[index2 + k + 4]);
				ymm2 = _mm256_loadu_pd(&D[index2 + k + 8]);
				ymm3 = _mm256_loadu_pd(&D[index2 + k + 12]);
				ymm4 = _mm256_loadu_pd(&D[index2 + k + 16]);
				ymm5 = _mm256_loadu_pd(&D[index2 + k + 20]);
				ymm6 = _mm256_loadu_pd(&D[index2 + k + 24]);
				ymm7 = _mm256_loadu_pd(&D[index2 + k + 28]);

				ymm0 = _mm256_mul_pd(ymm0, ymm8 );
				ymm1 = _mm256_mul_pd(ymm1, ymm9 );
				ymm2 = _mm256_mul_pd(ymm2, ymm10);
				ymm3 = _mm256_mul_pd(ymm3, ymm11);
				ymm4 = _mm256_mul_pd(ymm4, ymm12);
				ymm5 = _mm256_mul_pd(ymm5, ymm13);
				ymm6 = _mm256_mul_pd(ymm6, ymm14);
				ymm7 = _mm256_mul_pd(ymm7, ymm15);

				ymm0 = _mm256_add_pd(ymm0, ymm1);
				ymm2 = _mm256_add_pd(ymm2, ymm3);
				ymm4 = _mm256_add_pd(ymm4, ymm5);
				ymm6 = _mm256_add_pd(ymm6, ymm7);
				ymm0 = _mm256_add_pd(ymm0, ymm2);
				ymm4 = _mm256_add_pd(ymm4, ymm6);
				ymm0 = _mm256_add_pd(ymm0, ymm4);

				_mm256_storeu_pd(scratchpad, ymm0);
				for (s = 0; s < 4; s++)
					res += scratchpad[s];
			}
			for (k = col_reduced; k < col_reduced_16; k += 16) {
				ymm8 = _mm256_loadu_pd(&B[index1 + k]);
				ymm9 = _mm256_loadu_pd(&B[index1 + k + 4]);
				ymm10 = _mm256_loadu_pd(&B[index1 + k + 8]);
				ymm11 = _mm256_loadu_pd(&B[index1 + k + 12]);

				ymm0 = _mm256_loadu_pd(&D[index2 + k]);
				ymm1 = _mm256_loadu_pd(&D[index2 + k + 4]);
				ymm2 = _mm256_loadu_pd(&D[index2 + k + 8]);
				ymm3 = _mm256_loadu_pd(&D[index2 + k + 12]);

				ymm0 = _mm256_mul_pd(ymm0, ymm8 );
				ymm1 = _mm256_mul_pd(ymm1, ymm9 );
				ymm2 = _mm256_mul_pd(ymm2, ymm10);
				ymm3 = _mm256_mul_pd(ymm3, ymm11);

				ymm0 = _mm256_add_pd(ymm0, ymm1);
				ymm2 = _mm256_add_pd(ymm2, ymm3);
				ymm0 = _mm256_add_pd(ymm0, ymm2);

				_mm256_storeu_pd(scratchpad, fndef BLOCK_SIZE
#define BLOCK_SIZE ((int) 16)
#endif
				ymm0);
				for (s = 0; s < 4; s++)
					res += scratchpad[s];
			}
			for (k = col_reduced_16; k < M; k++) {
				res += D[index2 + k] * B[index1 + k];
			}
			C[j + index1] = res;
		}
	}


}
*/

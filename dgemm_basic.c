const char* dgemm_desc = "Basic, three-loop dgemm.";

void square_dgemm(const int M, 
                  const double *A, const double *B, double *C)
{
    int i, j, k;
    
    double *a = (double*)malloc(sizeof(double) * (M*K))
    double *b = (double*)malloc(sizeof(double) * (N*K))
    double *c = (double*)malloc(sizeof(double) * (M*N))
    
    for (i = 0; i < M*K; ++i) {
        a[i]=A[i]
    }
    for (i = 0; i < N*K; ++i) {
        b[i]=B[i]
    }
    for (i = 0; i < M*N; ++i) {
        c[i]=C[i]
    }

    
    for (i = 0; i < M; ++i) {
        for (j = 0; j < M; ++j) {
            double cij = c[j*M+i];
            for (k = 0; k < M; ++k)
                cij += a[k*M+i] * b[j*M+k];
            C[j*M+i] = cij;
        }
    }
}

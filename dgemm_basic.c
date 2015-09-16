const char* dgemm_desc = "Basic, three-loop dgemm.";

void square_dgemm(const int M, 
                  const double *A, const double *B, double *C)
{
    int i, j, k;
    
    double *a = (double*)malloc(sizeof(double) * (M*M));
    double *b = (double*)malloc(sizeof(double) * (M*M));
    double *c = (double*)malloc(sizeof(double) * (M*M));
    
    for (i = 0; i < M*M; ++i) {
        a[i]=A[i];
    }
    for (i = 0; i < M*M; ++i) {
        b[i]=B[i];
    }
    for (i = 0; i < M*M; ++i) {
        c[i]=C[i];
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

#include "jacobi_inverse.h"

int main()
{
    ji_t *inst = NULL;
    int M = 3, i = 0, j = 0;
    double Matrix_real[M*M];
    double Matrix_imag[M*M];
    double inv_real[M*M];
    double inv_imag[M*M];

    /*1.00000000e-08+0.j 0.00000000e+00+0.j 0.00000000e+00+0.j]
     [0.00000000e+00+0.j 2.40000001e+00+0.j 4.80000000e+00+0.j]
     [0.00000000e+00+0.j 4.80000000e+00+0.j 9.60000001e+00+0.j*/
    // row 1
    Matrix_real[0] = 1.00000000e-08;  Matrix_real[1] = 0;  Matrix_real[2] = 0; 
    Matrix_imag[0] = 0;         Matrix_imag[1] = 0;          Matrix_imag[2] = 0;
     
    // row 2
    Matrix_real[0+M] = 0;   Matrix_real[1+M] = 2.40000001;  Matrix_real[2+M] = 4.8;
    Matrix_imag[0+M] = 0;           Matrix_imag[1+M] = 0;         Matrix_imag[2+M] = 0;

    // row 3
    Matrix_real[0+M*2] = 0;  Matrix_real[1+M*2] = 4.8;  Matrix_real[2+M*2] = 9.60000001;
    Matrix_imag[0+M*2] = 0;          Matrix_imag[1+M*2] = 0;         Matrix_imag[2+M*2] = 0;

    printf("src:\n");
    for (i=0; i<M; i++) {
        for (j=0; j<M; j++) {
            printf("%.10f+%.10fj ", Matrix_real[i*M+j], Matrix_imag[i*M+j]);
        }
        printf("\n");
    }

    jacobi_inverse_create(&inst, M);
    jacobi_inverse(inst, Matrix_real, Matrix_imag, inv_real, inv_imag, 0);
    jacobi_inverse_destroy(inst);
    printf("inv:\n");
    for (i=0; i<M; i++) {
        for (j=0; j<M; j++) {
            printf("%.10f+%.10fj, ", inv_real[i*M+j], inv_imag[i*M+j]);
        }
        printf("\n");
    }

    
    return 0;
}
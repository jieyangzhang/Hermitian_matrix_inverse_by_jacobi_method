#include "jacobi_inverse.h"

int main()
{
    ji_t *inst = NULL;
    int M = 5, i = 0, j = 0;
    float Matrix_real[M*M];
    float Matrix_imag[M*M];
    /* [[1,    1+5j,  2+4j, 5,     10],
        [1-5j, 2,     6,    1-2j,  7],
        [2-4j, 6,     3,    3-1j,  1],
        [5,    1+2j,  3+1j, 2,     5-2j],
        [10,   7,     1,    5+2j,  1]] */
    // row 1
    Matrix_real[0] = 1;  Matrix_real[1] = 1;  Matrix_real[2] = 2;  Matrix_real[3] = 5;  Matrix_real[4] = 10;
    Matrix_imag[0] = 0;  Matrix_imag[1] = 5;  Matrix_imag[2] = 4;  Matrix_imag[3] = 0;  Matrix_imag[4] = 0;
     
    // row 2
    Matrix_real[0+5] = 1;   Matrix_real[1+5] = 2;  Matrix_real[2+5] = 6;  Matrix_real[3+5] = 1;  Matrix_real[4+5] = 7;
    Matrix_imag[0+5] = -5;  Matrix_imag[1+5] = 0;  Matrix_imag[2+5] = 0;  Matrix_imag[3+5] = -2; Matrix_imag[4+5] = 0;

    // row 3
    Matrix_real[0+5*2] = 2;  Matrix_real[1+5*2] = 6;  Matrix_real[2+5*2] = 3;  Matrix_real[3+5*2] = 3; Matrix_real[4+5*2] = 1;
    Matrix_imag[0+5*2] = -4; Matrix_imag[1+5*2] = 0;  Matrix_imag[2+5*2] = 0;  Matrix_imag[3+5*2] = -1; Matrix_imag[4+5*2] = 0;
  
    // row 4
    Matrix_real[0+5*3] = 5;Matrix_real[1+5*3] = 1;Matrix_real[2+5*3] = 3;Matrix_real[3+5*3] = 2;Matrix_real[4+5*3] = 5;
    Matrix_imag[0+5*3] = 0;Matrix_imag[1+5*3] = 2;Matrix_imag[2+5*3] = 1;Matrix_imag[3+5*3] = 0;Matrix_imag[4+5*3] = -2;

    // row 5
    Matrix_real[0+5*4] = 10;Matrix_real[1+5*4] = 7;Matrix_real[2+5*4] = 1;Matrix_real[3+5*4] = 5;Matrix_real[4+5*4] = 1;
    Matrix_imag[0+5*4] = 0;Matrix_imag[1+5*4] = 0;Matrix_imag[2+5*4] = 0;Matrix_imag[3+5*4] = 2;Matrix_imag[4+5*4] = 0;

    printf("src:\n");
    for (i=0; i<M; i++) {
        for (j=0; j<M; j++) {
            printf("%f+%fj ", Matrix_real[i*M+j], Matrix_imag[i*M+j]);
        }
        printf("\n");
    }

    jacobi_inverse_create(&inst, M);
    jacobi_inverse(inst, Matrix_real, Matrix_imag);
    jacobi_inverse_destroy(inst);
    printf("inv:\n");
    for (i=0; i<M; i++) {
        for (j=0; j<M; j++) {
            printf("%f+%fj ", Matrix_real[i*M+j], Matrix_imag[i*M+j]);
        }
        printf("\n");
    }

    
    return 0;
}
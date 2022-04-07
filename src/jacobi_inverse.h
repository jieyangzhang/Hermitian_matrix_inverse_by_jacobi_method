#ifndef __JACOBI_INVERSE_H__
#define __JACOBI_INVERSE_H__

#include <stdio.h>

typedef struct ji_s {
    int *m_index; /* 上三角部分中每一行的最大值对应的索引 */
    float *V_real; /* 特征向量矩阵实部， 按行展开 */
    float *V_imag; /* 特征向量矩阵虚部， 按行展开 */
    float *T_real; /* 用作临时变量储存， 无实际意义 */
    float *T_imag; /* 用作临时变量储存， 无实际意义 */
    float *eig_values; /* 储存特征值 */
    int M; /* 方阵阶数 */
} ji_t;

/**
 * @brief 创建求逆实例
 * 
 * @param inst handle
 * @param M 方阵阶数
 * @return int 0 -- success; others -- failed
 */
int jacobi_inverse_create(ji_t **inst, int M);

/**
 * @brief 销毁求逆实例
 * 
 * @param inst handle
 * @return int 0 -- success; others -- failed
 */
int jacobi_inverse_destroy(ji_t *inst);

/**
 * @brief 基于Jacobi迭代计算Hermitian矩阵的近似逆
 * 
 * @param Matrix_real(in) Hermitian矩阵的实部， 方阵按行展开
 * @param Matrix_imag(in) Hermitian矩阵的虚部， 方阵按行展开
 * @param inv_real(out) 得到的逆矩阵的实部，按行展开
 * @param inv_imag(out) 得到的逆矩阵的虚部，按行展开
 * @return int 0 -- success; others -- failed
 */
int jacobi_inverse(ji_t *inst,
                   float *Matrix_real,
                   float *Matrix_imag);

#endif
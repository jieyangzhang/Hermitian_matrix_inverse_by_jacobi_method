#include "jacobi_inverse.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

/** 
 * 参考
 * 1. 维基百科：https://en.wikipedia.org/wiki/Jacobi_method_for_complex_Hermitian_matrices
 * 2. python实现Hermitian矩阵对角化：https://github.com/mariarudenko/eigendecomposition
 **/

#define INV_SQRT2 (0.707106781f) /* 1 / sqrt(2) */
#define MAX_JOCIBI_ITERS 1000 /* 最大迭代次数 */
#define PI 3.141592653589793

/* func declare */
/**
 * @brief 左乘Givns rotation矩阵
 * 
 * @param Matrix_real 矩阵实部， 按行展开
 * @param Matrix_imag 矩阵虚部， 按行展开
 * @param theta 角度
 * @param p 操作维度索引1
 * @param q 操作唯独索引2
 * @param M 方阵阶数
 */
static void __rotation_left(double *Matrix_real, double *Matrix_imag, double theta, int p, int q, int M);

/**
 * @brief 右乘Givns rotation矩阵
 * 
 * @param Matrix_real 矩阵实部， 按行展开
 * @param Matrix_imag 矩阵虚部， 按行展开
 * @param theta 角度
 * @param p 操作维度索引1
 * @param q 操作唯独索引2
 * @param M 方阵阶数
 */
static void __rotation_right(double *Matrix_real, double *Matrix_imag, double theta, int p, int q, int M);

/**
 * @brief 计算矩阵的F范数（所有元素模长平方和的开方）
 * 
 * @param Matrix_real(in) 输入矩阵的实部，按行展开
 * @param Matrix_imag(in) 输入矩阵的虚部，按行展开
 * @param M(in) 方阵阶数
 * @param f_norm(out) 输出F范数
 */
static void __frobenius_norm(double *Matrix_real, double *Matrix_imag, int M, double *f_norm);

/**
 * @brief 计算非对角上三角部分的f范数
 * 
 * @param Matrix_real(in) 输入矩阵的实部，按行展开
 * @param Matrix_imag(in) 输入矩阵的虚部，按行展开
 * @param M(in) 方阵阶数
 * @param off(out) 输出F范数
 */
static void __off_diag(double *Matrix_real, double *Matrix_imag, int M, double *off);

static void __cal_transform_matrix(double *V_real,
                                   double *V_imag,
                                   double theta_1,
                                   double theta_2,
                                   int p,
                                   int q,
                                   int M,
                                   double *T_real,
                                   double *T_imag);
/* end func declare */

/**
 * @brief 创建求逆实例
 * 
 * @param inst handle
 * @param M 方阵阶数
 * @return int 0 -- success; others -- failed
 */
int jacobi_inverse_create(ji_t **inst, int M)
{
    int ret = 0;
    
    ji_t *self = (ji_t*)calloc(1, sizeof(ji_t));
    self->M = M;
    self->m_index = (int*)calloc(M-1, sizeof(int));
    self->V_real = (double*)calloc(M*M, sizeof(double));
    self->V_imag = (double*)calloc(M*M, sizeof(double));
    self->T_real = (double*)calloc(M*2, sizeof(double));
    self->T_imag = (double*)calloc(M*2, sizeof(double));
    self->eig_values = (double*)calloc(M, sizeof(double));

    *inst = self;

    return ret;
}

/**
 * @brief 销毁求逆实例
 * 
 * @param inst handle
 * @return int 0 -- success; others -- failed
 */
int jacobi_inverse_destroy(ji_t *inst)
{
    if (inst->m_index) {
        free(inst->m_index);
        inst->m_index = NULL;
    }
    if (inst->V_real) {
        free(inst->V_real);
        inst->V_real = NULL;
    }
    if (inst->V_imag) {
        free(inst->V_imag);
        inst->V_imag = NULL;
    }

    if (inst->T_real) {
        free(inst->T_real);
        inst->T_real = NULL;
    }

    if (inst->T_imag) {
        free(inst->T_imag);
        inst->T_imag = NULL;
    }

    if (inst->eig_values) {
        free(inst->eig_values);
        inst->eig_values = NULL;
    }
    free(inst);
    inst = NULL;

    return 0;
}

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
                   double *Matrix_real,
                   double *Matrix_imag,
                   double *inv_real,
                   double *inv_imag,
                   double diag_compensation)
{
    double tolerance = 1e-8; /* 收敛目标系数 */
    double eps = 0; /* 误差 */
    double f_norm = 0; /* 矩阵F范数 */
    double off = 0; /* 迭代前后非对角部分上三角范数 */ 
    double *V_real = inst->V_real, *V_imag = inst->V_imag;
    double *T_real = inst->T_real, *T_imag = inst->T_imag;
    double tmp_max = 0, tmp = 0;
    double real = 0, imag = 0, magn = 0;
    double phi_1 = 0, phi_2 = 0;
    double theta_1 = 0, theta_2 = 0;
    double *eig_value = inst->eig_values;
    int iter = 0, i = 0, j = 0, k = 0; /* 迭代索引 */
    int M = inst->M;
    int *m_index = inst->m_index;
    int p = 0, q = 0; /* 非对角部分最大模长索引， p为索引中较小值， q为较大值 */

    // TODO 检查输入方阵为Hermitian矩阵

    /* 将特征向量矩阵初始化为单位阵 */
    memset(V_real, 0, sizeof(double)*M*M);
    memset(V_imag, 0, sizeof(double)*M*M);
    memcpy(inv_real, Matrix_real, sizeof(double)*M*M);
    memcpy(inv_imag, Matrix_imag, sizeof(double)*M*M);
    for (i=0; i<M; i++) {

       *(V_real+i*M+i) = 1;
       *(inv_real+i*M+i) += diag_compensation;
    }
    for (i=0; i<M-1; i++) {
        tmp_max = 0;
        for (j=i+1; j<M; j++) {
            real = *(inv_real+M*i+j);
            imag = *(inv_imag+M*i+j);
            magn = real*real+imag*imag;
            if (magn > tmp_max) {
                tmp_max = magn;
                /* 记录列索引 */
                *(m_index+i) = j;
            }
        }
    }
    __off_diag(inv_real, inv_imag, M, &off);
    __frobenius_norm(inv_real, inv_imag, M, &f_norm);
    eps = f_norm*tolerance;
    while ((off>eps) && (iter<MAX_JOCIBI_ITERS)) {
        /* 确认非对角部分中，模长最大的项的索引 */
        /* m_index中记录当前矩阵非对角部分的每一行对应的最大模长的列索引。因此只需要比较m_index中的M-1个值，而不是全部遍历 */
        tmp_max = 0;
        for (i=0; i<M-1; i++) {
            j = *(m_index+i); /* 第i行最大模长对应的列索引 */
            real = *(inv_real+i*M+j);
            imag = *(inv_imag+i*M+j);
            magn = real*real+imag*imag;
            if (magn > tmp_max) {
                tmp_max = magn;
                p = i; /* i一定比j小 */
                q = j;
            }
        }

        real = *(inv_real+p*M+q);
        imag = *(inv_imag+p*M+q);
        magn = sqrt(real*real+imag*imag);
        if (!real) {
            phi_1 = 0;
        } else {
            phi_1 = atan(imag/real);
        }
        tmp = *(inv_real+p*M+p)-*(inv_real+q*M+q);
        if (!tmp) {
            phi_2 = 0;
        } else {
            phi_2 = atan(2*magn/tmp);
        }
        theta_1 = (2*phi_1-PI)/4;
        theta_2 = phi_2/2;
        __rotation_left(inv_real, inv_imag, theta_1, p, q, M);
        __rotation_left(inv_real, inv_imag, theta_2, p, q, M);
        __rotation_right(inv_real, inv_imag, theta_1, p, q, M);
        __rotation_right(inv_real, inv_imag, theta_2, p, q, M);
        /* 构建本次迭代中的旋转变换矩阵 */
        __cal_transform_matrix(V_real, V_imag, theta_1, theta_2, p, q, M, T_real, T_imag);
        /* 更新变化行的模长最大索引 */
        for (i=0; i<M-1; i++) {
            tmp_max = 0;
            if ((p==i) || (q==i)) {
                for (j=i+1; j<M; j++) {
                    real = *(inv_real+i*M+j);
                    imag = *(inv_imag+i*M+j);
                    magn = real*real+imag+imag;
                    if (magn>tmp_max) {
                        tmp_max = magn;
                        *(m_index+i) = j;
                    }
                }
            } else {
                int row_max = *(m_index+i);
                real = *(inv_real+i*M+row_max);
                imag = *(inv_imag+i*M+row_max);
                tmp_max = real*real+imag*imag;
                if (p>i) {
                    real = *(inv_real+i*M+p);
                    imag = *(inv_imag+i*M+p);
                    magn = real*real+imag*imag;
                    if (magn > tmp_max) {
                        tmp_max = magn;
                        *(m_index+i) = p;
                    }
                }
                if (q>i) {
                    real = *(inv_real+i*M+q);
                    imag = *(inv_imag+i*M+q);
                    magn = real*real+imag*imag;
                    if (magn > tmp_max) {
                        tmp_max = magn;
                        *(m_index+i) = q;
                    }
                }
            }
        }
        /* 计算本次迭代后的非对角部分的范数 */
        __off_diag(inv_real, inv_imag, M, &off);
        iter++;
    }

    /* 将V转置处理， 结果储存在V中 */
    for (i=0; i<M; i++) {
        for (j=i+1; j<M; j++) {
            real = *(V_real+i*M+j);
            imag = *(V_imag+i*M+j);

            *(V_real+i*M+j) = *(V_real+j*M+i);
            *(V_imag+i*M+j) = *(V_imag+j*M+i);

            *(V_real+j*M+i) = real;
            *(V_imag+j*M+i) = imag;
        }
    }

    for (i=0; i<M; i++) {
        *(eig_value+i) = *(inv_real+i*M+i);
    }

    double i_k_real = 0;
    double k_j_real = 0;
    double i_k_imag = 0;
    double k_j_imag = 0;
    for (i=0; i<M; i++) {
        /* 计算对角线上元素 */
        tmp = 0;
        for (k=0; k<M; k++) {
            real = *(V_real+i*M+k);
            imag = *(V_imag+i*M+k);
            tmp += (real*real+imag*imag) / *(eig_value+k);
        }
        *(inv_real+i*M+i) = tmp;
        *(inv_imag+i*M+i) = 0;
        for (j=i+1; j<M; j++) {
            real = 0;
            imag = 0;
            /* 第i行的共厄和第j行的转置的内积 */
            for (k=0; k<M; k++) {
                i_k_real = *(V_real+i*M+k), i_k_imag = -1 * *(V_imag+i*M+k);
                k_j_real = *(V_real+j*M+k), k_j_imag = *(V_imag+j*M+k);

                real += (i_k_real*k_j_real-i_k_imag*k_j_imag)/ *(eig_value+k);
                imag += (i_k_imag*k_j_real+i_k_real*k_j_imag)/ *(eig_value+k);
            }
            *(inv_real+i*M+j) = *(inv_real+j*M+i) = real;
            *(inv_imag+i*M+j) = *(inv_imag+j*M+i) = imag;
        }
    }

    return 0;
}

static void __rotation_left(double *Matrix_real, double *Matrix_imag, double theta, int p, int q, int M)
{
    int i = 0;
    double cos_theta = cos(theta), sin_theta = sin(theta);
    double p_real = 0, p_imag = 0, q_real = 0, q_imag = 0;

    for (i=0; i<M; i++) {
        p_real = *(Matrix_real+p*M+i);
        p_imag = *(Matrix_imag+p*M+i);
        q_real = *(Matrix_real+q*M+i);
        q_imag = *(Matrix_imag+q*M+i);
        /* [p, i] */
        *(Matrix_real+p*M+i) = INV_SQRT2*(cos_theta*p_real+sin_theta*p_imag-cos_theta*q_real+sin_theta*q_imag);
        *(Matrix_imag+p*M+i) = INV_SQRT2*(cos_theta*p_imag-sin_theta*p_real-cos_theta*q_imag-sin_theta*q_real);
        /* [q, i] */
        *(Matrix_real+q*M+i) = INV_SQRT2*(cos_theta*p_real+sin_theta*p_imag+cos_theta*q_real-sin_theta*q_imag);
        *(Matrix_imag+q*M+i) = INV_SQRT2*(cos_theta*p_imag-sin_theta*p_real+cos_theta*q_imag+sin_theta*q_real);
    }
}

static void __rotation_right(double *Matrix_real, double *Matrix_imag, double theta, int p, int q, int M)
{
    int i = 0;
    double cos_theta = cos(theta), sin_theta = sin(theta);
    double p_real = 0, p_imag = 0, q_real = 0, q_imag = 0;

    for (i=0; i<M; i++) {
        p_real = *(Matrix_real+i*M+p);
        p_imag = *(Matrix_imag+i*M+p);
        q_real = *(Matrix_real+i*M+q);
        q_imag = *(Matrix_imag+i*M+q);
        /* [i, p] */
        *(Matrix_real+i*M+p) = INV_SQRT2*(p_real*cos_theta-p_imag*sin_theta-q_real*cos_theta-q_imag*sin_theta);
        *(Matrix_imag+i*M+p) = INV_SQRT2*(p_imag*cos_theta+p_real*sin_theta-q_imag*cos_theta+q_real*sin_theta);
        /* [i, q] */
        *(Matrix_real+i*M+q) = INV_SQRT2*(p_real*cos_theta-p_imag*sin_theta+q_real*cos_theta+q_imag*sin_theta);
        *(Matrix_imag+i*M+q) = INV_SQRT2*(p_imag*cos_theta+p_real*sin_theta+q_imag*cos_theta-q_real*sin_theta);
    }
}

static void __frobenius_norm(double *Matrix_real, double *Matrix_imag, int M, double *f_norm)
{
    int i = 0, j = 0;
    double tmp = 0, sum = 0;
    double real = 0, imag = 0;
    for (i=0; i<M; i++) {
        /* 非对角中的上三角部分， 非对角的下三角部分和上三角平方和一致 */
        for (j=i+1; j<M; j++) {
            real = *(Matrix_real+i*M+j);
            imag = *(Matrix_imag+i*M+j);
            tmp += (real*real+imag*imag);
        }
        /* 对角线部分 */
        real = *(Matrix_real+i*M+i);
        imag = *(Matrix_imag+i*M+i);
        sum += (real*real+imag*imag);
    }
    sum += (2*tmp);
    *f_norm = sqrt(sum);
}

static void __off_diag(double *Matrix_real, double *Matrix_imag, int M, double *off)
{
    int i = 0, j = 0;
    double tmp = 0;
    double real = 0, imag = 0;
    for (i=0; i<M; i++) {
        /* 非对角中的上三角部分 */
        for (j=i+1; j<M; j++) {
            real = *(Matrix_real+i*M+j);
            imag = *(Matrix_imag+i*M+j);
            tmp += (real*real+imag*imag);
        }
    }
    *off = sqrt(tmp);
}

static void __cal_transform_matrix(double *V_real,
                                   double *V_imag,
                                   double theta_1,
                                   double theta_2,
                                   int p,
                                   int q,
                                   int M,
                                   double *T_real,
                                   double *T_imag)
{
    double sin_theta1 = sin(theta_1);
    double cos_theta1 = cos(theta_1);
    double sin_theta2 = sin(theta_2);
    double cos_theta2 = cos(theta_2);

    double R_p_p_real = -1*sin_theta1*sin_theta2, R_p_p_imag = -1*cos_theta1*sin_theta2;
    double R_p_q_real = -1*cos_theta1*cos_theta2, R_p_q_imag = -1*sin_theta1*cos_theta2;
    double R_q_p_real = cos_theta1*cos_theta2, R_q_p_imag = -1*sin_theta1*cos_theta2;
    double R_q_q_real = -1*sin_theta1*sin_theta2, R_q_q_imag = cos_theta1*sin_theta2;

    int i = 0, j = 0;
    for (i=0; i<M; i++) {
        if ((i == p) || (i == q)) {
            for (j=0; j<M; j++) {
                double V_p_j_real = *(V_real+p*M+j), V_p_j_imag = *(V_imag+p*M+j);
                double V_q_j_real = *(V_real+q*M+j), V_q_j_imag = *(V_imag+q*M+j);
                if (p == i) {
                    /* V(p, j) = R(p, p)*V(p, j)+R(p, q)*V(q, j) */
                    *(T_real+0*M+j) = (R_p_p_real*V_p_j_real-R_p_p_imag*V_p_j_imag)+(R_p_q_real*V_q_j_real-R_p_q_imag*V_q_j_imag);
                    *(T_imag+0*M+j) = (R_p_p_real*V_p_j_imag+R_p_p_imag*V_p_j_real)+(R_p_q_real*V_q_j_imag+R_p_q_imag*V_q_j_real);
                } else if (q == i) {
                    /* V(q, j) = R(q, p)*V(p, j)+R(q, q)*V(q, j) */
                    *(T_real+1*M+j) = (R_q_p_real*V_p_j_real-R_q_p_imag*V_p_j_imag)+(R_q_q_real*V_q_j_real-R_q_q_imag*V_q_j_imag);
                    *(T_imag+1*M+j) = (R_q_p_real*V_p_j_imag+R_q_p_imag*V_p_j_real)+(R_q_q_real*V_q_j_imag+R_q_q_imag*V_q_j_real);
                }
            }
        }
    }
    for (i=0; i<M; i++) {
        if (i == p) {
            for (j=0; j<M; j++) {
                *(V_real+i*M+j) = *(T_real+0*M+j);
                *(V_imag+i*M+j) = *(T_imag+0*M+j);
            }
        } else if (i == q) {
            for (j=0; j<M; j++) {
                *(V_real+i*M+j) = *(T_real+1*M+j);
                *(V_imag+i*M+j) = *(T_imag+1*M+j);
            }
        }
    }
}
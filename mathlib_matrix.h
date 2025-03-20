#ifndef MATHLIB_MATRIX_H
#define MATHLIB_MATRIX_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>

#ifdef __cplusplus
extern "C" {  // 告诉 C++ 编译器使用 C 的符号规则
#endif

typedef struct {
    int rows;
    int cols;
    double** data;
} Matrix;

// 矩阵操作函数声明
Matrix* create_matrix(int rows, int cols);
void destroy_matrix(Matrix* matrix);
Matrix* matrix_add(Matrix* a, Matrix* b);
Matrix* matrix_subtract(Matrix* a, Matrix* b);
Matrix* matrix_multiply(Matrix* a, Matrix* b);
Matrix* matrix_transpose(Matrix* matrix);
double matrix_determinant(Matrix* matrix);
double matrix_trace(Matrix* matrix);
Matrix* matrix_power(Matrix* matrix, int exponent);
Matrix* matrix_inverse(Matrix* matrix);
Matrix* matrix_adjoint(Matrix* matrix);
int matrix_rank(Matrix* matrix);
double* matrix_eigenvalues(Matrix* matrix, int* size);

#ifdef __cplusplus
}
#endif

#endif
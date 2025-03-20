#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>

typedef struct {
    int rows;
    int cols;
    double** data;
} Matrix;

// 创建矩阵
Matrix* create_matrix(int rows, int cols) {
    if (rows <= 0 || cols <= 0) {
        errno = EINVAL;
        return NULL;
    }
    Matrix* matrix = (Matrix*)malloc(sizeof(Matrix));
    if (matrix == NULL) {
        errno = ENOMEM;
        return NULL;
    }
    matrix->rows = rows;
    matrix->cols = cols;
    matrix->data = (double**)malloc(rows * sizeof(double*));
    if (matrix->data == NULL) {
        free(matrix);
        errno = ENOMEM;
        return NULL;
    }
    for (int i = 0; i < rows; i++) {
        matrix->data[i] = (double*)malloc(cols * sizeof(double));
        if (matrix->data[i] == NULL) {
            for (int j = 0; j < i; j++) {
                free(matrix->data[j]);
            }
            free(matrix->data);
            free(matrix);
            errno = ENOMEM;
            return NULL;
        }
    }
    return matrix;
}

// 销毁矩阵
void destroy_matrix(Matrix* matrix) {
    if (matrix == NULL) return;
    if (matrix->data != NULL) {
        for (int i = 0; i < matrix->rows; i++) {
            free(matrix->data[i]);
        }
        free(matrix->data);
    }
    free(matrix);
}

// 矩阵加法
Matrix* matrix_add(Matrix* a, Matrix* b) {
    if (a == NULL || b == NULL || a->rows != b->rows || a->cols != b->cols) {
        errno = EINVAL;
        return NULL;
    }
    Matrix* result = create_matrix(a->rows, a->cols);
    if (result == NULL) return NULL;
    for (int i = 0; i < a->rows; i++) {
        for (int j = 0; j < a->cols; j++) {
            result->data[i][j] = a->data[i][j] + b->data[i][j];
        }
    }
    return result;
}

// 矩阵减法
Matrix* matrix_subtract(Matrix* a, Matrix* b) {
    if (a == NULL || b == NULL || a->rows != b->rows || a->cols != b->cols) {
        errno = EINVAL;
        return NULL;
    }
    Matrix* result = create_matrix(a->rows, a->cols);
    if (result == NULL) return NULL;
    for (int i = 0; i < a->rows; i++) {
        for (int j = 0; j < a->cols; j++) {
            result->data[i][j] = a->data[i][j] - b->data[i][j];
        }
    }
    return result;
}

// 矩阵乘法
Matrix* matrix_multiply(Matrix* a, Matrix* b) {
    if (a == NULL || b == NULL || a->cols != b->rows) {
        errno = EINVAL;
        return NULL;
    }
    Matrix* result = create_matrix(a->rows, b->cols);
    if (result == NULL) return NULL;
    for (int i = 0; i < a->rows; i++) {
        for (int j = 0; j < b->cols; j++) {
            double sum = 0.0;
            for (int k = 0; k < a->cols; k++) {
                sum += a->data[i][k] * b->data[k][j];
            }
            result->data[i][j] = sum;
        }
    }
    return result;
}

// 矩阵转置
Matrix* matrix_transpose(Matrix* matrix) {
    if (matrix == NULL) {
        errno = EINVAL;
        return NULL;
    }
    Matrix* result = create_matrix(matrix->cols, matrix->rows);
    if (result == NULL) return NULL;
    for (int i = 0; i < matrix->cols; i++) {
        for (int j = 0; j < matrix->rows; j++) {
            result->data[i][j] = matrix->data[j][i];
        }
    }
    return result;
}

// 计算矩阵行列式（仅支持2x2和3x3矩阵）
double matrix_determinant(Matrix* matrix) {
    if (matrix == NULL || matrix->rows != matrix->cols) {
        errno = EINVAL;
        return 0.0;
    }
    if (matrix->rows == 2) {
        return matrix->data[0][0] * matrix->data[1][1] - matrix->data[0][1] * matrix->data[1][0];
    } else if (matrix->rows == 3) {
        return matrix->data[0][0] * (matrix->data[1][1] * matrix->data[2][2] - matrix->data[1][2] * matrix->data[2][1])
               - matrix->data[0][1] * (matrix->data[1][0] * matrix->data[2][2] - matrix->data[1][2] * matrix->data[2][0])
               + matrix->data[0][2] * (matrix->data[1][0] * matrix->data[2][1] - matrix->data[1][1] * matrix->data[2][0]);
    } else {
        errno = ENOSYS;
        return 0.0;
    }
}

// 计算矩阵的迹
double matrix_trace(Matrix* matrix) {
    if (matrix == NULL || matrix->rows != matrix->cols) {
        errno = EINVAL;
        return 0.0;
    }
    double trace = 0.0;
    for (int i = 0; i < matrix->rows; i++) {
        trace += matrix->data[i][i];
    }
    return trace;
}

// 计算矩阵的幂
Matrix* matrix_power(Matrix* matrix, int exponent) {
    if (matrix == NULL || matrix->rows != matrix->cols || exponent < 0) {
        errno = EINVAL;
        return NULL;
    }
    if (exponent == 0) {
        Matrix* identity = create_matrix(matrix->rows, matrix->cols);
        if (identity == NULL) return NULL;
        for (int i = 0; i < matrix->rows; i++) {
            for (int j = 0; j < matrix->cols; j++) {
                identity->data[i][j] = (i == j) ? 1.0 : 0.0;
            }
        }
        return identity;
    }
    Matrix* result = create_matrix(matrix->rows, matrix->cols);
    if (result == NULL) return NULL;
    for (int i = 0; i < matrix->rows; i++) {
        for (int j = 0; j < matrix->cols; j++) {
            result->data[i][j] = (i == j) ? 1.0 : 0.0;
        }
    }
    Matrix* temp = create_matrix(matrix->rows, matrix->cols);
    if (temp == NULL) {
        destroy_matrix(result);
        return NULL;
    }
    for (int i = 0; i < matrix->rows; i++) {
        for (int j = 0; j < matrix->cols; j++) {
            temp->data[i][j] = matrix->data[i][j];
        }
    }
    while (exponent > 0) {
        if (exponent % 2 == 1) {
            Matrix* product = matrix_multiply(result, temp);
            if (product == NULL) {
                destroy_matrix(result);
                destroy_matrix(temp);
                return NULL;
            }
            destroy_matrix(result);
            result = product;
        }
        Matrix* product = matrix_multiply(temp, temp);
        if (product == NULL) {
            destroy_matrix(result);
            destroy_matrix(temp);
            return NULL;
        }
        destroy_matrix(temp);
        temp = product;
        exponent /= 2;
    }
    destroy_matrix(temp);
    return result;
}

// 计算矩阵的逆（仅支持2x2和3x3矩阵）
Matrix* matrix_inverse(Matrix* matrix) {
    if (matrix == NULL || matrix->rows != matrix->cols) {
        errno = EINVAL;
        return NULL;
    }
    double determinant = matrix_determinant(matrix);
    if (determinant == 0) {
        errno = EDOM;
        return NULL;
    }
    if (matrix->rows == 2) {
        Matrix* inverse = create_matrix(2, 2);
        if (inverse == NULL) return NULL;
        inverse->data[0][0] = matrix->data[1][1] / determinant;
        inverse->data[0][1] = -matrix->data[0][1] / determinant;
        inverse->data[1][0] = -matrix->data[1][0] / determinant;
        inverse->data[1][1] = matrix->data[0][0] / determinant;
        return inverse;
    } else if (matrix->rows == 3) {
        Matrix* inverse = create_matrix(3, 3);
        if (inverse == NULL) return NULL;
        inverse->data[0][0] = (matrix->data[1][1] * matrix->data[2][2] - matrix->data[1][2] * matrix->data[2][1]) / determinant;
        inverse->data[0][1] = -(matrix->data[0][1] * matrix->data[2][2] - matrix->data[0][2] * matrix->data[2][1]) / determinant;
        inverse->data[0][2] = (matrix->data[0][1] * matrix->data[1][2] - matrix->data[0][2] * matrix->data[1][1]) / determinant;
        inverse->data[1][0] = -(matrix->data[1][0] * matrix->data[2][2] - matrix->data[1][2] * matrix->data[2][0]) / determinant;
        inverse->data[1][1] = (matrix->data[0][0] * matrix->data[2][2] - matrix->data[0][2] * matrix->data[2][0]) / determinant;
        inverse->data[1][2] = -(matrix->data[0][0] * matrix->data[1][2] - matrix->data[0][2] * matrix->data[1][0]) / determinant;
        inverse->data[2][0] = (matrix->data[1][0] * matrix->data[2][1] - matrix->data[1][1] * matrix->data[2][0]) / determinant;
        inverse->data[2][1] = -(matrix->data[0][0] * matrix->data[2][1] - matrix->data[0][1] * matrix->data[2][0]) / determinant;
        inverse->data[2][2] = (matrix->data[0][0] * matrix->data[1][1] - matrix->data[0][1] * matrix->data[1][0]) / determinant;
        return inverse;
    } else {
        errno = ENOSYS;
        return NULL;
    }
}

// 计算矩阵的伴随矩阵（仅支持2x2和3x3矩阵）
Matrix* matrix_adjoint(Matrix* matrix) {
    if (matrix == NULL || matrix->rows != matrix->cols) {
        errno = EINVAL;
        return NULL;
    }
    if (matrix->rows == 2) {
        Matrix* adjoint = create_matrix(2, 2);
        if (adjoint == NULL) return NULL;
        adjoint->data[0][0] = matrix->data[1][1];
        adjoint->data[0][1] = -matrix->data[0][1];
        adjoint->data[1][0] = -matrix->data[1][0];
        adjoint->data[1][1] = matrix->data[0][0];
        return adjoint;
    } else if (matrix->rows == 3) {
        Matrix* adjoint = create_matrix(3, 3);
        if (adjoint == NULL) return NULL;
        adjoint->data[0][0] = matrix->data[1][1] * matrix->data[2][2] - matrix->data[1][2] * matrix->data[2][1];
        adjoint->data[0][1] = -(matrix->data[0][1] * matrix->data[2][2] - matrix->data[0][2] * matrix->data[2][1]);
        adjoint->data[0][2] = matrix->data[0][1] * matrix->data[1][2] - matrix->data[0][2] * matrix->data[1][1];
        adjoint->data[1][0] = -(matrix->data[1][0] * matrix->data[2][2] - matrix->data[1][2] * matrix->data[2][0]);
        adjoint->data[1][1] = matrix->data[0][0] * matrix->data[2][2] - matrix->data[0][2] * matrix->data[2][0];
        adjoint->data[1][2] = -(matrix->data[0][0] * matrix->data[1][2] - matrix->data[0][2] * matrix->data[1][0]);
        adjoint->data[2][0] = matrix->data[1][0] * matrix->data[2][1] - matrix->data[1][1] * matrix->data[2][0];
        adjoint->data[2][1] = -(matrix->data[0][0] * matrix->data[2][1] - matrix->data[0][1] * matrix->data[2][0]);
        adjoint->data[2][2] = matrix->data[0][0] * matrix->data[1][1] - matrix->data[0][1] * matrix->data[1][0];
        return adjoint;
    } else {
        errno = ENOSYS;
        return NULL;
    }
}

// 计算矩阵的秩（仅支持2x2和3x3矩阵）
int matrix_rank(Matrix* matrix) {
    if (matrix == NULL) {
        errno = EINVAL;
        return -1;
    }
    if (matrix->rows == 0 || matrix->cols == 0) return 0;
    if (matrix->rows == 1 || matrix->cols == 1) return 1;
    if (matrix->rows == 2 && matrix->cols == 2) {
        if (matrix_determinant(matrix) != 0) return 2;
        else return 1;
    }
    if (matrix->rows == 3 && matrix->cols == 3) {
        if (matrix_determinant(matrix) != 0) return 3;
        else {
            Matrix* adjoint = matrix_adjoint(matrix);
            if (adjoint == NULL) return -1;
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    if (adjoint->data[i][j] != 0) {
                        destroy_matrix(adjoint);
                        return 2;
                    }
                }
            }
            destroy_matrix(adjoint);
            return 1;
        }
    }
    errno = ENOSYS;
    return -1;
}

// 计算矩阵的特征值（仅支持2x2和3x3矩阵）
double* matrix_eigenvalues(Matrix* matrix, int* size) {
    if (matrix == NULL || matrix->rows != matrix->cols) {
        errno = EINVAL;
        return NULL;
    }
    if (matrix->rows == 2) {
        double a = matrix->data[0][0];
        double b = matrix->data[0][1];
        double c = matrix->data[1][0];
        double d = matrix->data[1][1];
        double trace = a + d;
        double determinant = a * d - b * c;
        double discriminant = trace * trace - 4 * determinant;
        if (discriminant < 0) {
            errno = EDOM;
            return NULL;
        }
        double sqrt_discriminant = sqrt(discriminant);
        double eigen1 = (trace + sqrt_discriminant) / 2;
        double eigen2 = (trace - sqrt_discriminant) / 2;
        double* eigenvalues = (double*)malloc(2 * sizeof(double));
        if (eigenvalues == NULL) {
            errno = ENOMEM;
            return NULL;
        }
        eigenvalues[0] = eigen1;
        eigenvalues[1] = eigen2;
        *size = 2;
        return eigenvalues;
    } else if (matrix->rows == 3) {
        double a = matrix->data[0][0];
        double b = matrix->data[0][1];
        double c = matrix->data[0][2];
        double d = matrix->data[1][0];
        double e = matrix->data[1][1];
        double f = matrix->data[1][2];
        double g = matrix->data[2][0];
        double h = matrix->data[2][1];
        double i = matrix->data[2][2];
        double trace = a + e + i;
        double determinant = matrix_determinant(matrix);
        double sum_of_products = a * e + a * i + e * i + b * f + c * g + d * h;
        double characteristic_equation[4];
        characteristic_equation[0] = 1.0;
        characteristic_equation[1] = -trace;
        characteristic_equation[2] = sum_of_products - determinant;
        characteristic_equation[3] = -determinant;
        // 使用牛顿迭代法求解三次方程
        double* eigenvalues = (double*)malloc(3 * sizeof(double));
        if (eigenvalues == NULL) {
            errno = ENOMEM;
            return NULL;
        }
        for (int k = 0; k < 3; k++) {
            double guess = trace / 3.0;
            double prev_guess;
            do {
                prev_guess = guess;
                double f = characteristic_equation[0] * pow(guess, 3) +
                           characteristic_equation[1] * pow(guess, 2) +
                           characteristic_equation[2] * guess +
                           characteristic_equation[3];
                double f_prime = 3 * characteristic_equation[0] * pow(guess, 2) +
                                 2 * characteristic_equation[1] * guess +
                                 characteristic_equation[2];
                guess = guess - f / f_prime;
            } while (fabs(guess - prev_guess) > 1e-10);
            eigenvalues[k] = guess;
            // 将方程除以 (x - guess) 以降低次数
            double new_equation[3];
            new_equation[0] = characteristic_equation[0];
            new_equation[1] = characteristic_equation[1] + new_equation[0] * guess;
            new_equation[2] = characteristic_equation[2] + new_equation[1] * guess;
            characteristic_equation[0] = new_equation[1];
            characteristic_equation[1] = new_equation[2];
            characteristic_equation[2] = 0.0;
        }
        *size = 3;
        return eigenvalues;
    } else {
        errno = ENOSYS;
        return NULL;
    }
}
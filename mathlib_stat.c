#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>

// 声明 compare_doubles 函数
static int compare_doubles(const void* a, const void* b);

// 声明 factorial 函数
unsigned long long factorial(int n);

// 计算均值
double calculate_mean(double* data, int size) {
    if (data == NULL || size <= 0) {
        errno = EINVAL;
        return NAN;
    }
    double sum = 0.0;
    for (int i = 0; i < size; i++) {
        sum += data[i];
    }
    return sum / size;
}

// 计算方差
double calculate_variance(double* data, int size) {
    if (data == NULL || size <= 0) {
        errno = EINVAL;
        return NAN;
    }
    double mean = calculate_mean(data, size);
    double sum = 0.0;
    for (int i = 0; i < size; i++) {
        sum += pow(data[i] - mean, 2);
    }
    return sum / size;
}

// 计算标准差
double calculate_standard_deviation(double* data, int size) {
    if (data == NULL || size <= 0) {
        errno = EINVAL;
        return NAN;
    }
    return sqrt(calculate_variance(data, size));
}

// 计算协方差
double calculate_covariance(double* data1, double* data2, int size) {
    if (data1 == NULL || data2 == NULL || size <= 0) {
        errno = EINVAL;
        return NAN;
    }
    double mean1 = calculate_mean(data1, size);
    double mean2 = calculate_mean(data2, size);
    double sum = 0.0;
    for (int i = 0; i < size; i++) {
        sum += (data1[i] - mean1) * (data2[i] - mean2);
    }
    return sum / size;
}

// 计算相关系数
double calculate_correlation(double* data1, double* data2, int size) {
    if (data1 == NULL || data2 == NULL || size <= 0) {
        errno = EINVAL;
        return NAN;
    }
    double cov = calculate_covariance(data1, data2, size);
    double std_dev1 = calculate_standard_deviation(data1, size);
    double std_dev2 = calculate_standard_deviation(data2, size);
    
    if (std_dev1 == 0 || std_dev2 == 0) {
        errno = EDOM;
        return NAN;
    }
    
    return cov / (std_dev1 * std_dev2);
}

// 计算偏度
double calculate_skewness(double* data, int size) {
    if (data == NULL || size <= 0) {
        errno = EINVAL;
        return NAN;
    }
    double mean = calculate_mean(data, size);
    double std_dev = calculate_standard_deviation(data, size);
    if (std_dev == 0) {
        errno = EDOM;
        return NAN;
    }
    double sum = 0.0;
    for (int i = 0; i < size; i++) {
        sum += pow((data[i] - mean) / std_dev, 3);
    }
    return sum / size;
}

// 计算峰度
double calculate_kurtosis(double* data, int size) {
    if (data == NULL || size <= 0) {
        errno = EINVAL;
        return NAN;
    }
    double mean = calculate_mean(data, size);
    double std_dev = calculate_standard_deviation(data, size);
    if (std_dev == 0) {
        errno = EDOM;
        return NAN;
    }
    double sum = 0.0;
    for (int i = 0; i < size; i++) {
        sum += pow((data[i] - mean) / std_dev, 4);
    }
    return sum / size - 3.0;
}

// 计算中位数
double calculate_median(double* data, int size) {
    if (data == NULL || size <= 0) {
        errno = EINVAL;
        return NAN;
    }
    double* sorted_data = (double*)malloc(size * sizeof(double));
    if (sorted_data == NULL) {
        errno = ENOMEM;
        return NAN;
    }
    for (int i = 0; i < size; i++) {
        sorted_data[i] = data[i];
    }
    qsort(sorted_data, size, sizeof(double), compare_doubles);
    
    double median;
    if (size % 2 == 0) {
        median = (sorted_data[size / 2 - 1] + sorted_data[size / 2]) / 2.0;
    } else {
        median = sorted_data[size / 2];
    }
    
    free(sorted_data);
    return median;
}

// 实现 compare_doubles 函数
static int compare_doubles(const void* a, const void* b) {
    double arg1 = *(const double*)a;
    double arg2 = *(const double*)b;
    if (arg1 < arg2) return -1;
    if (arg1 > arg2) return 1;
    return 0;
}

// 计算众数
double calculate_mode(double* data, int size) {
    if (data == NULL || size <= 0) {
        errno = EINVAL;
        return NAN;
    }
    double* sorted_data = (double*)malloc(size * sizeof(double));
    if (sorted_data == NULL) {
        errno = ENOMEM;
        return NAN;
    }
    for (int i = 0; i < size; i++) {
        sorted_data[i] = data[i];
    }
    qsort(sorted_data, size, sizeof(double), compare_doubles);
    
    double mode = sorted_data[0];
    int max_count = 1;
    int current_count = 1;
    
    for (int i = 1; i < size; i++) {
        if (sorted_data[i] == sorted_data[i - 1]) {
            current_count++;
            if (current_count > max_count) {
                max_count = current_count;
                mode = sorted_data[i];
            }
        } else {
            current_count = 1;
        }
    }
    
    free(sorted_data);
    return mode;
}

// 计算范围（极差）
double calculate_range(double* data, int size) {
    if (data == NULL || size <= 0) {
        errno = EINVAL;
        return NAN;
    }
    double* sorted_data = (double*)malloc(size * sizeof(double));
    if (sorted_data == NULL) {
        errno = ENOMEM;
        return NAN;
    }
    for (int i = 0; i < size; i++) {
        sorted_data[i] = data[i];
    }
    qsort(sorted_data, size, sizeof(double), compare_doubles);
    
    double range = sorted_data[size - 1] - sorted_data[0];
    
    free(sorted_data);
    return range;
}

// 计算四分位数
void calculate_quartiles(double* data, int size, double* q1, double* q2, double* q3) {
    if (data == NULL || size <= 0 || q1 == NULL || q2 == NULL || q3 == NULL) {
        errno = EINVAL;
        return;
    }
    double* sorted_data = (double*)malloc(size * sizeof(double));
    if (sorted_data == NULL) {
        errno = ENOMEM;
        return;
    }
    for (int i = 0; i < size; i++) {
        sorted_data[i] = data[i];
    }
    qsort(sorted_data, size, sizeof(double), compare_doubles);
    
    int mid = size / 2;
    *q2 = sorted_data[mid];
    
    if (size % 2 == 0) {
        *q1 = calculate_median(sorted_data, mid);
        *q3 = calculate_median(sorted_data + mid, mid);
    } else {
        *q1 = calculate_median(sorted_data, mid);
        *q3 = calculate_median(sorted_data + mid + 1, size - mid - 1);
    }
    
    free(sorted_data);
}

// 计算累积分布函数（CDF）对于正态分布
double calculate_normal_cdf(double x, double mean, double std_dev) {
    if (std_dev <= 0) {
        errno = EDOM;
        return NAN;
    }
    double z = (x - mean) / (std_dev * sqrt(2));
    return 0.5 * (1 + erf(z));
}

// 计算 z 分数
double calculate_z_score(double x, double mean, double std_dev) {
    if (std_dev == 0) {
        errno = EDOM;
        return NAN;
    }
    return (x - mean) / std_dev;
}

// 计算二项式系数
unsigned long long calculate_binomial_coefficient(int n, int k) {
    if (n < 0 || k < 0 || k > n) {
        errno = EINVAL;
        return 0;
    }
    if (k == 0 || k == n) return 1;
    if (k > n - k) k = n - k; // 优化计算
    unsigned long long result = 1;
    for (int i = 1; i <= k; i++) {
        result = result * (n - k + i) / i;
    }
    return result;
}

// 计算二项式概率
double calculate_binomial_probability(int k, int n, double p) {
    if (k < 0 || k > n || p < 0 || p > 1) {
        errno = EINVAL;
        return NAN;
    }
    unsigned long long coeff = calculate_binomial_coefficient(n, k);
    return coeff * pow(p, k) * pow(1 - p, n - k);
}

// 计算泊松概率
double calculate_poisson_probability(int k, double lambda) {
    if (k < 0 || lambda <= 0) {
        errno = EINVAL;
        return NAN;
    }
    return pow(lambda, k) * exp(-lambda) / factorial(k);
}

// 实现 factorial 函数
unsigned long long factorial(int n) {
    if (n < 0) {
        errno = EINVAL;
        return 0;
    }
    unsigned long long result = 1;
    for (int i = 2; i <= n; i++) {
        result *= i;
    }
    return result;
}

// 计算几何概率
double calculate_geometric_probability(int k, double p) {
    if (k < 0 || p <= 0 || p >= 1) {
        errno = EINVAL;
        return NAN;
    }
    return (1 - p) * pow(p, k);
}

// 计算超几何概率
double calculate_hypergeometric_probability(int k, int n, int K, int N) {
    if (k < 0 || n < 0 || K < 0 || N < 0 || k > n || K > N || n > N) {
        errno = EINVAL;
        return NAN;
    }
    unsigned long long coeff1 = calculate_binomial_coefficient(K, k);
    unsigned long long coeff2 = calculate_binomial_coefficient(N - K, n - k);
    unsigned long long coeff3 = calculate_binomial_coefficient(N, n);
    return (double)(coeff1 * coeff2) / coeff3;
}

// 计算线性回归系数
void calculate_linear_regression(double* x_data, double* y_data, int size, double* slope, double* intercept) {
    if (x_data == NULL || y_data == NULL || size <= 0 || slope == NULL || intercept == NULL) {
        errno = EINVAL;
        return;
    }
    double x_mean = calculate_mean(x_data, size);
    double y_mean = calculate_mean(y_data, size);
    
    double numerator = 0.0;
    double denominator = 0.0;
    for (int i = 0; i < size; i++) {
        numerator += (x_data[i] - x_mean) * (y_data[i] - y_mean);
        denominator += pow(x_data[i] - x_mean, 2);
    }
    
    if (denominator == 0) {
        errno = EDOM;
        return;
    }
    
    *slope = numerator / denominator;
    *intercept = y_mean - *slope * x_mean;
}

// 计算决定系数（R平方）
double calculate_r_squared(double* y_actual, double* y_predicted, int size) {
    if (y_actual == NULL || y_predicted == NULL || size <= 0) {
        errno = EINVAL;
        return NAN;
    }
    double y_mean = calculate_mean(y_actual, size);
    
    double ss_total = 0.0;
    double ss_residual = 0.0;
    for (int i = 0; i < size; i++) {
        ss_total += pow(y_actual[i] - y_mean, 2);
        ss_residual += pow(y_actual[i] - y_predicted[i], 2);
    }
    
    if (ss_total == 0) {
        errno = EDOM;
        return NAN;
    }
    
    return 1.0 - ss_residual / ss_total;
}

// 计算平均绝对误差
double calculate_mae(double* y_actual, double* y_predicted, int size) {
    if (y_actual == NULL || y_predicted == NULL || size <= 0) {
        errno = EINVAL;
        return NAN;
    }
    double sum = 0.0;
    for (int i = 0; i < size; i++) {
        sum += fabs(y_actual[i] - y_predicted[i]);
    }
    return sum / size;
}

// 计算均方误差
double calculate_mse(double* y_actual, double* y_predicted, int size) {
    if (y_actual == NULL || y_predicted == NULL || size <= 0) {
        errno = EINVAL;
        return NAN;
    }
    double sum = 0.0;
    for (int i = 0; i < size; i++) {
        sum += pow(y_actual[i] - y_predicted[i], 2);
    }
    return sum / size;
}

// 计算均方根误差
double calculate_rmse(double* y_actual, double* y_predicted, int size) {
    if (y_actual == NULL || y_predicted == NULL || size <= 0) {
        errno = EINVAL;
        return NAN;
    }
    return sqrt(calculate_mse(y_actual, y_predicted, size));
}

// 计算平均绝对百分比误差
double calculate_mape(double* y_actual, double* y_predicted, int size) {
    if (y_actual == NULL || y_predicted == NULL || size <= 0) {
        errno = EINVAL;
        return NAN;
    }
    double sum = 0.0;
    for (int i = 0; i < size; i++) {
        if (y_actual[i] == 0) {
            errno = EDOM;
            return NAN;
        }
        sum += fabs((y_actual[i] - y_predicted[i]) / y_actual[i]);
    }
    return sum / size * 100.0;
}
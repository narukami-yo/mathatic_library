#ifndef MATHLIB_STAT_H
#define MATHLIB_STAT_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>

// 计算均值
double calculate_mean(double* data, int size);

// 计算方差
double calculate_variance(double* data, int size);

// 计算标准差
double calculate_standard_deviation(double* data, int size);

// 计算协方差
double calculate_covariance(double* data1, double* data2, int size);

// 计算相关系数
double calculate_correlation(double* data1, double* data2, int size);

// 计算偏度
double calculate_skewness(double* data, int size);

// 计算峰度
double calculate_kurtosis(double* data, int size);

// 计算中位数
double calculate_median(double* data, int size);

// 计算众数
double calculate_mode(double* data, int size);

// 计算范围（极差）
double calculate_range(double* data, int size);

// 计算四分位数
void calculate_quartiles(double* data, int size, double* q1, double* q2, double* q3);

// 计算累积分布函数（CDF）对于正态分布
double calculate_normal_cdf(double x, double mean, double std_dev);

// 计算 z 分数
double calculate_z_score(double x, double mean, double std_dev);

// 计算二项式系数
unsigned long long calculate_binomial_coefficient(int n, int k);

// 计算二项式概率
double calculate_binomial_probability(int k, int n, double p);

// 计算泊松概率
double calculate_poisson_probability(int k, double lambda);

// 计算几何概率
double calculate_geometric_probability(int k, double p);

// 计算超几何概率
double calculate_hypergeometric_probability(int k, int n, int K, int N);

// 计算线性回归系数
void calculate_linear_regression(double* x_data, double* y_data, int size, double* slope, double* intercept);

// 计算决定系数（R平方）
double calculate_r_squared(double* y_actual, double* y_predicted, int size);

// 计算平均绝对误差
double calculate_mae(double* y_actual, double* y_predicted, int size);

// 计算均方误差
double calculate_mse(double* y_actual, double* y_predicted, int size);

// 计算均方根误差
double calculate_rmse(double* y_actual, double* y_predicted, int size);

// 计算平均绝对百分比误差
double calculate_mape(double* y_actual, double* y_predicted, int size);

#ifdef __cplusplus
}
#endif

#endif // MATHLIB_STAT_H
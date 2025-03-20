#ifndef MATHLIB_BASE_H
#define MATHLIB_BASE_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>

// 函数声明保持不变
unsigned long long factorial(int n);
int gcd(int a, int b);
int lcm(int a, int b);
unsigned long long fibonacci(int n);
unsigned long long combination(int n, int k);
unsigned long long permutation(int n, int k);
double power(double base, int exponent);
double sqrt_custom(double num);
double ln(double num);
double log10_custom(double num);
double exp_custom(double num);
double sin_custom(double num);
double cos_custom(double num);
double tan_custom(double num);
double sinh_custom(double num);
double cosh_custom(double num);
double tanh_custom(double num);
double asin_custom(double num);
double acos_custom(double num);
double atan_custom(double num);
double asinh_custom(double num);
double acosh_custom(double num);
double atanh_custom(double num);
double abs_custom(double num);
double floor_custom(double num);
double ceil_custom(double num);
double round_custom(double num);
double square(double num);
double cube(double num);

unsigned long long* factorial_array(int n, int* size);
unsigned long long* fibonacci_array(int n, int* size);
unsigned long long* combination_array(int n, int k, int* size);
unsigned long long* permutation_array(int n, int k, int* size);
double* power_array(double base, int exponent, int* size);
double* sqrt_array(double start, double end, int* size);
double* ln_array(double start, double end, int* size);
double* log10_array(double start, double end, int* size);
double* exp_array(double start, double end, int* size);
double* sin_array(double start, double end, int* size);
double* cos_array(double start, double end, int* size);
double* tan_array(double start, double end, int* size);
double* sinh_array(double start, double end, int* size);
double* cosh_array(double start, double end, int* size);
double* tanh_array(double start, double end, int* size);
double* asin_array(double start, double end, int* size);
double* acos_array(double start, double end, int* size);
double* atan_array(double start, double end, int* size);
double* asinh_array(double start, double end, int* size);
double* acosh_array(double start, double end, int* size);
double* atanh_array(double start, double end, int* size);
double* abs_array(double start, double end, int* size);
double* floor_array(double start, double end, int* size);
double* ceil_array(double start, double end, int* size);
double* round_array(double start, double end, int* size);
double* square_array(double start, double end, int* size);
double* cube_array(double start, double end, int* size);

#ifdef __cplusplus
}
#endif

#endif
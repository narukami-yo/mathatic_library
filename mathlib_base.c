#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>

// 计算阶乘
unsigned long long factorial(int n) {
    if (n < 0) {
        errno = EINVAL;
        return 0;
    }
    unsigned long long result = 1;
    for (int i = 1; i <= n; i++) {
        result *= i;
    }
    return result;
}

// 计算最大公约数
int gcd(int a, int b) {
    if (a < 0) a = -a;
    if (b < 0) b = -b;
    while (b != 0) {
        int temp = b;
        b = a % b;
        a = temp;
    }
    return a;
}

// 计算最小公倍数
int lcm(int a, int b) {
    if (a == 0 || b == 0) {
        errno = EINVAL;
        return 0;
    }
    return (a / gcd(a, b)) * b;
}

// 计算斐波那契数列
unsigned long long fibonacci(int n) {
    if (n < 0) {
        errno = EINVAL;
        return 0;
    }
    if (n == 0) return 0;
    if (n == 1) return 1;
    unsigned long long a = 0, b = 1, c;
    for (int i = 2; i <= n; i++) {
        c = a + b;
        a = b;
        b = c;
    }
    return b;
}

// 计算组合数
unsigned long long combination(int n, int k) {
    if (n < 0 || k < 0 || k > n) {
        errno = EINVAL;
        return 0;
    }
    if (k == 0 || k == n) return 1;
    if (k > n - k) k = n - k;
    unsigned long long result = 1;
    for (int i = 1; i <= k; i++) {
        result = result * (n - k + i) / i;
    }
    return result;
}

// 计算排列数
unsigned long long permutation(int n, int k) {
    if (n < 0 || k < 0 || k > n) {
        errno = EINVAL;
        return 0;
    }
    unsigned long long result = 1;
    for (int i = n - k + 1; i <= n; i++) {
        result *= i;
    }
    return result;
}

// 计算幂
double power(double base, int exponent) {
    if (base == 0 && exponent < 0) {
        errno = EDOM;
        return 0.0;
    }
    double result = 1.0;
    int abs_exponent = exponent;
    if (exponent < 0) abs_exponent = -exponent;
    for (int i = 0; i < abs_exponent; i++) {
        result *= base;
    }
    if (exponent < 0) result = 1.0 / result;
    return result;
}

// 计算平方根
double sqrt_custom(double num) {
    if (num < 0) {
        errno = EDOM;
        return 0.0;
    }
    if (num == 0) return 0.0;
    double guess = num / 2.0;
    double prev_guess;
    do {
        prev_guess = guess;
        guess = (guess + num / guess) / 2.0;
    } while (fabs(guess - prev_guess) > 1e-10);
    return guess;
}

// 计算自然对数
double ln(double num) {
    if (num <= 0) {
        errno = EDOM;
        return 0.0;
    }
    return log(num);
}

// 计算常用对数
double log10_custom(double num) {
    if (num <= 0) {
        errno = EDOM;
        return 0.0;
    }
    return log10(num);
}

// 计算指数函数
double exp_custom(double num) {
    return exp(num);
}

// 计算正弦函数
double sin_custom(double num) {
    return sin(num);
}

// 计算余弦函数
double cos_custom(double num) {
    return cos(num);
}

// 计算正切函数
double tan_custom(double num) {
    return tan(num);
}

// 计算双曲正弦函数
double sinh_custom(double num) {
    return sinh(num);
}

// 计算双曲余弦函数
double cosh_custom(double num) {
    return cosh(num);
}

// 计算双曲正切函数
double tanh_custom(double num) {
    return tanh(num);
}

// 计算反正弦函数
double asin_custom(double num) {
    if (num < -1 || num > 1) {
        errno = EDOM;
        return 0.0;
    }
    return asin(num);
}

// 计算反余弦函数
double acos_custom(double num) {
    if (num < -1 || num > 1) {
        errno = EDOM;
        return 0.0;
    }
    return acos(num);
}

// 计算反正切函数
double atan_custom(double num) {
    return atan(num);
}

// 计算反双曲正弦函数
double asinh_custom(double num) {
    return asinh(num);
}

// 计算反双曲余弦函数
double acosh_custom(double num) {
    if (num < 1) {
        errno = EDOM;
        return 0.0;
    }
    return acosh(num);
}

// 计算反双曲正切函数
double atanh_custom(double num) {
    if (num <= -1 || num >= 1) {
        errno = EDOM;
        return 0.0;
    }
    return atanh(num);
}

// 计算绝对值
double abs_custom(double num) {
    return fabs(num);
}

// 计算取整函数
double floor_custom(double num) {
    return floor(num);
}

// 计算取顶函数
double ceil_custom(double num) {
    return ceil(num);
}

// 计算四舍五入函数
double round_custom(double num) {
    return round(num);
}

// 计算平方函数
double square(double num) {
    return num * num;
}

// 计算立方函数
double cube(double num) {
    return num * num * num;
}

// 计算阶乘数组
unsigned long long* factorial_array(int n, int* size) {
    if (n < 0) {
        errno = EINVAL;
        return NULL;
    }
    *size = n + 1;
    unsigned long long* result = (unsigned long long*)malloc((n + 1) * sizeof(unsigned long long));

    result[0] = 1;
    for (int i = 1; i <= n; i++) {
        result[i] = result[i - 1] * i;
    }
    return result;
}

// 计算斐波那契数列数组
unsigned long long* fibonacci_array(int n, int* size) {
    if (n < 0) {
        errno = EINVAL;
        return NULL;
    }
    *size = n + 1;
    unsigned long long* result = (unsigned long long*)malloc((n + 1) * sizeof(unsigned long long));
    result[0] = 0;
    if (n >= 1) result[1] = 1;
    for (int i = 2; i <= n; i++) {
        result[i] = result[i - 1] + result[i - 2];
    }
    return result;
}

// 计算组合数数组
unsigned long long* combination_array(int n, int k, int* size) {
    if (n < 0 || k < 0 || k > n) {
        errno = EINVAL;
        return NULL;
    }
    *size = k + 1;
    unsigned long long* result = (unsigned long long*)malloc((k + 1) * sizeof(unsigned long long));
    result[0] = 1;
    for (int i = 1; i <= k; i++) {
        result[i] = combination(n, i);
    }
    return result;
}

// 计算排列数数组
unsigned long long* permutation_array(int n, int k, int* size) {
    if (n < 0 || k < 0 || k > n) {
        errno = EINVAL;
        return NULL;
    }
    *size = k + 1;
    unsigned long long* result = (unsigned long long*)malloc((k + 1) * sizeof(unsigned long long));
    result[0] = 1;
    for (int i = 1; i <= k; i++) {
        result[i] = permutation(n, i);
    }
    return result;
}

// 计算幂数组
// 计算幂数组
double* power_array(double base, int exponent, int* size) {
    *size = exponent + 1;
    double* result = (double*)malloc((*size) * sizeof(double));
    result[0] = 1.0;
    for (int i = 1; i <= exponent; i++) {
        result[i] = power(base, i);
    }
    return result;
}
// 计算平方根数组
double* sqrt_array(double start, double end, int* size) {
    if (start < 0 || end < 0 || start > end) {
        errno = EDOM;
        return NULL;
    }
    int count = (int)((end - start) / 0.1) + 1;
    *size = count;
    double* result = (double*)malloc(count * sizeof(double));
    for (int i = 0; i < count; i++) {
        double num = start + i * 0.1;
        result[i] = sqrt_custom(num);
    }
    return result;
}

// 计算自然对数数组
double* ln_array(double start, double end, int* size) {
    if (start <= 0 || end <= 0 || start > end) {
        errno = EDOM;
        return NULL;
    }
    int count = (int)((end - start) / 0.1) + 1;
    *size = count;
    double* result = (double*)malloc(count * sizeof(double));
    for (int i = 0; i < count; i++) {
        double num = start + i * 0.1;
        result[i] = ln(num);
    }
    return result;
}

// 计算常用对数数组
double* log10_array(double start, double end, int* size) {
    if (start <= 0 || end <= 0 || start > end) {
        errno = EDOM;
        return NULL;
    }
    int count = (int)((end - start) / 0.1) + 1;
    *size = count;
    double* result = (double*)malloc(count * sizeof(double));
    for (int i = 0; i < count; i++) {
        double num = start + i * 0.1;
        result[i] = log10_custom(num);
    }
    return result;
}

// 计算指数函数数组
double* exp_array(double start, double end, int* size) {
    int count = (int)((end - start) / 0.1) + 1;
    *size = count;
    double* result = (double*)malloc(count * sizeof(double));

    for (int i = 0; i < count; i++) {
        double num = start + i * 0.1;
        result[i] = exp_custom(num);
    }
    return result;
}

// 计算正弦函数数组
double* sin_array(double start, double end, int* size) {
    int count = (int)((end - start) / 0.1) + 1;
    *size = count;
    double* result = (double*)malloc(count * sizeof(double));

    for (int i = 0; i < count; i++) {
        double num = start + i * 0.1;
        result[i] = sin_custom(num);
    }
    return result;
}

// 计算余弦函数数组
double* cos_array(double start, double end, int* size) {
    int count = (int)((end - start) / 0.1) + 1;
    *size = count;
    double* result = (double*)malloc(count * sizeof(double));

    for (int i = 0; i < count; i++) {
        double num = start + i * 0.1;
        result[i] = cos_custom(num);
    }
    return result;
}

// 计算正切函数数组
double* tan_array(double start, double end, int* size) {
    int count = (int)((end - start) / 0.1) + 1;
    *size = count;
    double* result = (double*)malloc(count * sizeof(double));

    for (int i = 0; i < count; i++) {
        double num = start + i * 0.1;
        result[i] = tan_custom(num);
    }
    return result;
}

// 计算双曲正弦函数数组
double* sinh_array(double start, double end, int* size) {
    int count = (int)((end - start) / 0.1) + 1;
    *size = count;
    double* result = (double*)malloc(count * sizeof(double));
 
    for (int i = 0; i < count; i++) {
        double num = start + i * 0.1;
        result[i] = sinh_custom(num);
    }
    return result;
}

// 计算双曲余弦函数数组
double* cosh_array(double start, double end, int* size) {
    int count = (int)((end - start) / 0.1) + 1;
    *size = count;
    double* result = (double*)malloc(count * sizeof(double));

    for (int i = 0; i < count; i++) {
        double num = start + i * 0.1;
        result[i] = cosh_custom(num);
    }
    return result;
}

// 计算双曲正切函数数组
double* tanh_array(double start, double end, int* size) {
    int count = (int)((end - start) / 0.1) + 1;
    *size = count;
    double* result = (double*)malloc(count * sizeof(double));

    for (int i = 0; i < count; i++) {
        double num = start + i * 0.1;
        result[i] = tanh_custom(num);
    }
    return result;
}

// 计算反正弦函数数组
double* asin_array(double start, double end, int* size) {
    int count = (int)((end - start) / 0.1) + 1;
    *size = count;
    double* result = (double*)malloc(count * sizeof(double));

    for (int i = 0; i < count; i++) {
        double num = start + i * 0.1;
        result[i] = asin_custom(num);
    }
    return result;
}

// 计算反余弦函数数组
double* acos_array(double start, double end, int* size) {
    int count = (int)((end - start) / 0.1) + 1;
    *size = count;
    double* result = (double*)malloc(count * sizeof(double));

    for (int i = 0; i < count; i++) {
        double num = start + i * 0.1;
        result[i] = acos_custom(num);
    }
    return result;
}

// 计算反正切函数数组
double* atan_array(double start, double end, int* size) {
    int count = (int)((end - start) / 0.1) + 1;
    *size = count;
    double* result = (double*)malloc(count * sizeof(double));

    for (int i = 0; i < count; i++) {
        double num = start + i * 0.1;
        result[i] = atan_custom(num);
    }
    return result;
}

// 计算反双曲正弦函数数组
double* asinh_array(double start, double end, int* size) {
    int count = (int)((end - start) / 0.1) + 1;
    *size = count;
    double* result = (double*)malloc(count * sizeof(double));

    for (int i = 0; i < count; i++) {
        double num = start + i * 0.1;
        result[i] = asinh_custom(num);
    }
    return result;
}

// 计算反双曲余弦函数数组
double* acosh_array(double start, double end, int* size) {
    int count = (int)((end - start) / 0.1) + 1;
    *size = count;
    double* result = (double*)malloc(count * sizeof(double));

    for (int i = 0; i < count; i++) {
        double num = start + i * 0.1;
        result[i] = acosh_custom(num);
    }
    return result;
}

// 计算反双曲正切函数数组
double* atanh_array(double start, double end, int* size) {
    int count = (int)((end - start) / 0.1) + 1;
    *size = count;
    double* result = (double*)malloc(count * sizeof(double));
    for (int i = 0; i < count; i++) {
        double num = start + i * 0.1;
        result[i] = atanh_custom(num);
    }
    return result;
}

// 计算绝对值数组
double* abs_array(double start, double end, int* size) {
    int count = (int)((end - start) / 0.1) + 1;
    *size = count;
    double* result = (double*)malloc(count * sizeof(double));
    for (int i = 0; i < count; i++) {
        double num = start + i * 0.1;
        result[i] = abs_custom(num);
    }
    return result;
}

// 计算取整函数数组
double* floor_array(double start, double end, int* size) {
    int count = (int)((end - start) / 0.1) + 1;
    *size = count;
    double* result = (double*)malloc(count * sizeof(double));
    for (int i = 0; i < count; i++) {
        double num = start + i * 0.1;
        result[i] = floor_custom(num);
    }
    return result;
}

// 计算取顶函数数组
double* ceil_array(double start, double end, int* size) {
    int count = (int)((end - start) / 0.1) + 1;
    *size = count;
    double* result = (double*)malloc(count * sizeof(double));
    for (int i = 0; i < count; i++) {
        double num = start + i * 0.1;
        result[i] = ceil_custom(num);
    }
    return result;
}

// 计算四舍五入函数数组
double* round_array(double start, double end, int* size) {
    int count = (int)((end - start) / 0.1) + 1;
    *size = count;
    double* result = (double*)malloc(count * sizeof(double));
    for (int i = 0; i < count; i++) {
        double num = start + i * 0.1;
        result[i] = round_custom(num);
    }
    return result;
}

// 计算平方函数数组
double* square_array(double start, double end, int* size) {
    int count = (int)((end - start) / 0.1) + 1;
    *size = count;
    double* result = (double*)malloc(count * sizeof(double));
    for (int i = 0; i < count; i++) {
        double num = start + i * 0.1;
        result[i] = square(num);
    }
    return result;
}

// 计算立方函数数组
double* cube_array(double start, double end, int* size) {
    int count = (int)((end - start) / 0.1) + 1;
    *size = count;
    double* result = (double*)malloc(count * sizeof(double));
    for (int i = 0; i < count; i++) {
        double num = start + i * 0.1;
        result[i] = cube(num);
    }
    return result;
}
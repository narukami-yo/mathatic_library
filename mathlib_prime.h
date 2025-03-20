#ifdef __cplusplus
extern "C" {
#endif

#ifndef MATHLIB_PRIME_H
#define MATHLIB_PRIME_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>

// 素数工具和排序算法声明
int is_prime(int n);
int prime_count(int n);
int* prime_sieve(int start, int end, int* count);
void quick_sort(int* arr, int left, int right);
void merge_sort(int* arr, int left, int right);
void heapify(int* arr, int n, int i);
void heap_sort(int* arr, int n);
void radix_sort(int* arr, int n);

#endif

#ifdef __cplusplus
}
#endif
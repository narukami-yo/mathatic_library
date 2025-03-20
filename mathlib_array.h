#ifndef MATHLIB_ARRAY_H
#define MATHLIB_ARRAY_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

typedef struct {
    int capacity;
    int size;
    int* data;
} DynamicArray;

// 动态数组函数声明
DynamicArray* create_dynamic_array(int initial_capacity);
void destroy_dynamic_array(DynamicArray* array);
int dynamic_array_add(DynamicArray* array, int value);
int dynamic_array_insert(DynamicArray* array, int index, int value);
int dynamic_array_remove(DynamicArray* array, int index);
int dynamic_array_get(DynamicArray* array, int index, int* value);
int dynamic_array_set(DynamicArray* array, int index, int value);
int dynamic_array_size(DynamicArray* array);
int dynamic_array_capacity(DynamicArray* array);
int dynamic_array_find(DynamicArray* array, int value);

// 排序算法声明
void dynamic_array_bubble_sort(DynamicArray* array);
void dynamic_array_selection_sort(DynamicArray* array);
void dynamic_array_insertion_sort(DynamicArray* array);
void dynamic_array_shell_sort(DynamicArray* array);
void dynamic_array_merge_sort(DynamicArray* array);
void dynamic_array_quick_sort(DynamicArray* array);
void dynamic_array_heap_sort(DynamicArray* array);
void dynamic_array_radix_sort(DynamicArray* array);

#ifdef __cplusplus
}
#endif

#endif
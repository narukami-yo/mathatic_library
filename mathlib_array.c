#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

typedef struct {
    int capacity;
    int size;
    int* data;
} DynamicArray;

// 创建动态数组
DynamicArray* create_dynamic_array(int initial_capacity) {
    if (initial_capacity <= 0) {
        errno = EINVAL;
        return NULL;
    }
    DynamicArray* array = (DynamicArray*)malloc(sizeof(DynamicArray));
    if (array == NULL) {
        errno = ENOMEM;
        return NULL;
    }
    array->capacity = initial_capacity;
    array->size = 0;
    array->data = (int*)malloc(initial_capacity * sizeof(int));
    if (array->data == NULL) {
        free(array);
        errno = ENOMEM;
        return NULL;
    }
    return array;
}

// 销毁动态数组
void destroy_dynamic_array(DynamicArray* array) {
    if (array == NULL) return;
    if (array->data != NULL) {
        free(array->data);
    }
    free(array);
}

// 添加元素到动态数组
int dynamic_array_add(DynamicArray* array, int value) {
    if (array == NULL) {
        errno = EINVAL;
        return 0;
    }
    if (array->size >= array->capacity) {
        int new_capacity = array->capacity * 2;
        int* new_data = (int*)realloc(array->data, new_capacity * sizeof(int));
        if (new_data == NULL) {
            errno = ENOMEM;
            return 0;
        }
        array->data = new_data;
        array->capacity = new_capacity;
    }
    array->data[array->size++] = value;
    return 1;
}

// 在指定位置插入元素
int dynamic_array_insert(DynamicArray* array, int index, int value) {
    if (array == NULL || index < 0 || index > array->size) {
        errno = EINVAL;
        return 0;
    }
    if (array->size >= array->capacity) {
        int new_capacity = array->capacity * 2;
        int* new_data = (int*)realloc(array->data, new_capacity * sizeof(int));
        if (new_data == NULL) {
            errno = ENOMEM;
            return 0;
        }
        array->data = new_data;
        array->capacity = new_capacity;
    }
    for (int i = array->size; i > index; i--) {
        array->data[i] = array->data[i - 1];
    }
    array->data[index] = value;
    array->size++;
    return 1;
}

// 删除指定位置的元素
int dynamic_array_remove(DynamicArray* array, int index) {
    if (array == NULL || index < 0 || index >= array->size) {
        errno = EINVAL;
        return 0;
    }
    for (int i = index; i < array->size - 1; i++) {
        array->data[i] = array->data[i + 1];
    }
    array->size--;
    return 1;
}

// 获取动态数组元素
int dynamic_array_get(DynamicArray* array, int index, int* value) {
    if (array == NULL || index < 0 || index >= array->size) {
        errno = EINVAL;
        return 0;
    }
    *value = array->data[index];
    return 1;
}

// 设置动态数组元素
int dynamic_array_set(DynamicArray* array, int index, int value) {
    if (array == NULL || index < 0 || index >= array->size) {
        errno = EINVAL;
        return 0;
    }
    array->data[index] = value;
    return 1;
}

// 获取动态数组大小
int dynamic_array_size(DynamicArray* array) {
    if (array == NULL) {
        errno = EINVAL;
        return -1;
    }
    return array->size;
}

// 获取动态数组容量
int dynamic_array_capacity(DynamicArray* array) {
    if (array == NULL) {
        errno = EINVAL;
        return -1;
    }
    return array->capacity;
}

// 查找元素在动态数组中的位置
int dynamic_array_find(DynamicArray* array, int value) {
    if (array == NULL) {
        errno = EINVAL;
        return -1;
    }
    for (int i = 0; i < array->size; i++) {
        if (array->data[i] == value) {
            return i;
        }
    }
    return -1;
}

// 动态数组冒泡排序
void dynamic_array_bubble_sort(DynamicArray* array) {
    if (array == NULL) return;
    for (int i = 0; i < array->size - 1; i++) {
        for (int j = 0; j < array->size - i - 1; j++) {
            if (array->data[j] > array->data[j + 1]) {
                int temp = array->data[j];
                array->data[j] = array->data[j + 1];
                array->data[j + 1] = temp;
            }
        }
    }
}

// 动态数组选择排序
void dynamic_array_selection_sort(DynamicArray* array) {
    if (array == NULL) return;
    for (int i = 0; i < array->size - 1; i++) {
        int min_index = i;
        for (int j = i + 1; j < array->size; j++) {
            if (array->data[j] < array->data[min_index]) {
                min_index = j;
            }
        }
        if (min_index != i) {
            int temp = array->data[i];
            array->data[i] = array->data[min_index];
            array->data[min_index] = temp;
        }
    }
}

// 动态数组插入排序
void dynamic_array_insertion_sort(DynamicArray* array) {
    if (array == NULL) return;
    for (int i = 1; i < array->size; i++) {
        int key = array->data[i];
        int j = i - 1;
        while (j >= 0 && array->data[j] > key) {
            array->data[j + 1] = array->data[j];
            j--;
        }
        array->data[j + 1] = key;
    }
}

// 动态数组希尔排序
void dynamic_array_shell_sort(DynamicArray* array) {
    if (array == NULL) return;
    int gap = array->size / 2;
    while (gap > 0) {
        for (int i = gap; i < array->size; i++) {
            int temp = array->data[i];
            int j = i;
            while (j >= gap && array->data[j - gap] > temp) {
                array->data[j] = array->data[j - gap];
                j -= gap;
            }
            array->data[j] = temp;
        }
        gap /= 2;
    }
}

// 动态数组归并排序
void dynamic_array_merge_sort(DynamicArray* array) {
    if (array == NULL || array->size <= 1) return;
    int mid = array->size / 2;
    DynamicArray* left = create_dynamic_array(mid);
    DynamicArray* right = create_dynamic_array(array->size - mid);
    for (int i = 0; i < mid; i++) {
        dynamic_array_add(left, array->data[i]);
    }
    for (int i = mid; i < array->size; i++) {
        dynamic_array_add(right, array->data[i]);
    }
    dynamic_array_merge_sort(left);
    dynamic_array_merge_sort(right);
    int i = 0, j = 0, k = 0;
    while (i < left->size && j < right->size) {
        if (left->data[i] <= right->data[j]) {
            array->data[k++] = left->data[i++];
        } else {
            array->data[k++] = right->data[j++];
        }
    }
    while (i < left->size) {
        array->data[k++] = left->data[i++];
    }
    while (j < right->size) {
        array->data[k++] = right->data[j++];
    }
    destroy_dynamic_array(left);
    destroy_dynamic_array(right);
}

// 动态数组快速排序
// 动态数组堆排序
void dynamic_array_heapify(DynamicArray* array, int n, int i) {
    if (array == NULL || n <= 0 || i >= n) return;
    int largest = i;
    int left = 2 * i + 1;
    int right = 2 * i + 2;
    if (left < n && array->data[left] > array->data[largest]) largest = left;
    if (right < n && array->data[right] > array->data[largest]) largest = right;
    if (largest != i) {
        int temp = array->data[i];
        array->data[i] = array->data[largest];
        array->data[largest] = temp;
        dynamic_array_heapify(array, n, largest);
    }
}

void dynamic_array_heap_sort(DynamicArray* array) {
    if (array == NULL || array->size <= 1) return;
    for (int i = array->size / 2 - 1; i >= 0; i--) {
        dynamic_array_heapify(array, array->size, i);
    }
    for (int i = array->size - 1; i >= 0; i--) {
        int temp = array->data[0];
        array->data[0] = array->data[i];
        array->data[i] = temp;
        dynamic_array_heapify(array, i, 0);
    }
}

// 动态数组基数排序
void dynamic_array_radix_sort(DynamicArray* array) {
    if (array == NULL || array->size <= 1) return;
    int max = array->data[0];
    for (int i = 1; i < array->size; i++) {
        if (array->data[i] > max) max = array->data[i];
    }
    for (int exp = 1; max / exp > 0; exp *= 10) {
        int* output = (int*)malloc(array->size * sizeof(int));
        if (output == NULL) {
            errno = ENOMEM;
            return;
        }
        int count[10] = {0};
        for (int i = 0; i < array->size; i++) {
            count[(array->data[i] / exp) % 10]++;
        }
        for (int i = 1; i < 10; i++) {
            count[i] += count[i - 1];
        }
        for (int i = array->size - 1; i >= 0; i--) {
            output[count[(array->data[i] / exp) % 10] - 1] = array->data[i];
            count[(array->data[i] / exp) % 10]--;
        }
        for (int i = 0; i < array->size; i++) {
            array->data[i] = output[i];
        }
        free(output);
    }
}
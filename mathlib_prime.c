#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>

// 检测素数
int is_prime(int n) {
    if (n <= 1) {
        errno = EINVAL;
        return 0;
    }
    if (n == 2) return 1;
    if (n % 2 == 0) return 0;
    for (int i = 3; i <= sqrt(n); i += 2) {
        if (n % i == 0) return 0;
    }
    return 1;
}

// 计算素数计数（小于等于n的素数个数）
int prime_count(int n) {
    if (n < 2) return 0;
    int count = 0;
    for (int i = 2; i <= n; i++) {
        if (is_prime(i)) count++;
    }
    return count;
}

// 计算素数区间筛选（找出[start, end]之间的所有素数）
int* prime_sieve(int start, int end, int* count) {
    if (start < 2) start = 2;
    if (end < start) {
        errno = EINVAL;
        return NULL;
    }
    int* primes = (int*)malloc((end - start + 1) * sizeof(int));

    int index = 0;
    for (int i = start; i <= end; i++) {
        if (is_prime(i)) {
            primes[index++] = i;
        }
    }
    *count = index;
    return primes;
}

// 快速排序
void quick_sort(int* arr, int left, int right) {
    if (arr == NULL || left >= right) return;
    int pivot = arr[(left + right) / 2];
    int i = left, j = right;
    while (i <= j) {
        while (arr[i] < pivot) i++;
        while (arr[j] > pivot) j--;
        if (i <= j) {
            int temp = arr[i];
            arr[i] = arr[j];
            arr[j] = temp;
            i++;
            j--;
        }
    }
    if (left < j) quick_sort(arr, left, j);
    if (i < right) quick_sort(arr, i, right);
}

// 归并排序
void merge_sort(int* arr, int left, int right) {
    if (arr == NULL || left >= right) return;
    int mid = (left + right) / 2;
    merge_sort(arr, left, mid);
    merge_sort(arr, mid + 1, right);
    int* temp = (int*)malloc((right - left + 1) * sizeof(int));
    int i = left, j = mid + 1, k = 0;
    while (i <= mid && j <= right) {
        if (arr[i] <= arr[j]) {
            temp[k++] = arr[i++];
        } else {
            temp[k++] = arr[j++];
        }
    }
    while (i <= mid) {
        temp[k++] = arr[i++];
    }
    while (j <= right) {
        temp[k++] = arr[j++];
    }
    for (i = left, k = 0; i <= right; i++, k++) {
        arr[i] = temp[k];
    }
    free(temp);
}

// 堆排序
void heapify(int* arr, int n, int i) {
    if (arr == NULL || n <= 0 || i >= n) return;
    int largest = i;
    int left = 2 * i + 1;
    int right = 2 * i + 2;
    if (left < n && arr[left] > arr[largest]) largest = left;
    if (right < n && arr[right] > arr[largest]) largest = right;
    if (largest != i) {
        int temp = arr[i];
        arr[i] = arr[largest];
        arr[largest] = temp;
        heapify(arr, n, largest);
    }
}

void heap_sort(int* arr, int n) {
    if (arr == NULL || n <= 0) return;
    for (int i = n / 2 - 1; i >= 0; i--) {
        heapify(arr, n, i);
    }
    for (int i = n - 1; i >= 0; i--) {
        int temp = arr[0];
        arr[0] = arr[i];
        arr[i] = temp;
        heapify(arr, i, 0);
    }
}

// 基数排序
void radix_sort(int* arr, int n) {
    if (arr == NULL || n <= 0) return;
    int max = arr[0];
    for (int i = 1; i < n; i++) {
        if (arr[i] > max) max = arr[i];
    }
    int* output = (int*)malloc(n * sizeof(int));
    for (int exp = 1; max / exp > 0; exp *= 10) {
        int count[10] = {0};
        for (int i = 0; i < n; i++) {
            count[(arr[i] / exp) % 10]++;
        }
        for (int i = 1; i < 10; i++) {
            count[i] += count[i - 1];
        }
        for (int i = n - 1; i >= 0; i--) {
            output[count[(arr[i] / exp) % 10] - 1] = arr[i];
            count[(arr[i] / exp) % 10]--;
        }
        for (int i = 0; i < n; i++) {
            arr[i] = output[i];
        }
    }
    free(output);
}
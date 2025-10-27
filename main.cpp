#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <chrono>
#include <iomanip>

class Sorter {
private:
    std::vector<int> arr;

    // Вспомогательная функция MERGE для Merge Sort
    void merge(int p, int q, int r) {
        std::vector<int> mas(r - p + 1);

        // for i ← q+1 down to p+1
        int i;
        for (i = q + 1; i >= p + 1; i--) {
            mas[i - 1] = arr[i - 1];
        }
        i = i + 1;

        // for j ← q to r-1
        int j;
        for (j = q; j <= r - 1; j++) {
            mas[r + q - j] = arr[j + 1];
        }

        // for k ← p to r
        for (int k = p; k <= r; k++) {
            if (mas[j] < mas[i]) {
                // then A[k] ← mas[j]
                arr[k] = mas[j];
                j = j - 1;
            } else {
                // else A[k] ← mas[i]
                arr[k] = mas[i];
                i = i + 1;
            }
        }
    }

    // Вспомогательная функция MERGE-SORT
    void mergeSortHelper(int p, int r) {
        // if p < r
        if (p < r) {
            // q ← [(p + r) /2]
            int q = (p + r) / 2;
            // MERGE-SORT (A, p, q)
            mergeSortHelper(p, q);
            // MERGE-SORT (A, q+1, r)
            mergeSortHelper(q + 1, r);
            // MERGE (A, p, q, r)
            merge(p, q, r);
        }
    }

    // Вспомогательная функция RX-SORT для Radix Sort
    void rxSort(int t) {
        int n = arr.size();
        std::vector<int> C(10);
        std::vector<int> B(n);

        // for j ← 0 to 9
        for (int j = 0; j <= 9; j++) {
            // do C[j] ← 0
            C[j] = 0;
        }

        // for j ← 0 to length [A] – 1
        for (int j = 0; j <= n - 1; j++) {
            // do C[(A[j] mod (t * 10)) div t] + 1
            int digit = (arr[j] % (t * 10)) / t;
            C[digit] = C[digit] + 1;
        }

        // for j ← 1 to 9
        for (int j = 1; j <= 9; j++) {
            // do C[j] ← C[j – 1] + C[j]
            C[j] = C[j - 1] + C[j];
        }

        // for j ← length [A] – 1 down to 0
        for (int j = n - 1; j >= 0; j--) {
            // do B[C[(A[j] mod (t * 10)) div t)] ← A[j]
            int digit = (arr[j] % (t * 10)) / t;
            B[C[digit] - 1] = arr[j];
            // do C[(A[j] mod (t * 10)) div t] ← C [(A[j] mod (t * 10)) div t] – 1
            C[digit] = C[digit] - 1;
        }

        // for j ← 0 to length [A] – 1
        for (int j = 0; j <= n - 1; j++) {
            // do A[j] ← B[j]
            arr[j] = B[j];
        }
    }

    // Функция для определения максимального количества цифр
    int getMaxDigits() {
        int maxNum = *std::max_element(arr.begin(), arr.end());
        int digits = 0;
        if (maxNum == 0) return 1;
        while (maxNum > 0) {
            digits++;
            maxNum /= 10;
        }
        return digits;
    }

public:
    // Конструктор
    Sorter(std::vector<int> array) : arr(array) {}

    // Установить новый массив
    void setArray(std::vector<int> array) {
        arr = array;
    }

    // Вывод массива
    void printArray(const std::string& message) {
        std::cout << message;
        for (int num : arr) {
            std::cout << num << " ";
        }
        std::cout << std::endl;
    }

    // INSERTION-SORT
    double insertion_sort() {
        std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();

        int n = arr.size();
        // for j←2 to length[A] (в псевдокоде индексы с 1, в C++ с 0, поэтому j=1 соответствует второму элементу)
        for (int j = 1; j <= n - 1; j++) {
            // key ← A[j]
            int key = arr[j];
            // i ← j - 1
            int i = j - 1;
            // while i > 0 & A[i] > key (в C++ индексы с 0, поэтому i >= 0)
            while (i >= 0 && arr[i] > key) {
                // A[i+1] ← A[i]
                arr[i + 1] = arr[i];
                // i ← i - 1
                i = i - 1;
            }
            // A[i+1] ← key
            arr[i + 1] = key;
        }

        std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
        long long nanoseconds = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
        return nanoseconds / 1000000000.0;
    }

    // MERGE-SORT
    double merge_sort() {
        std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();

        mergeSortHelper(0, arr.size() - 1);

        std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
        long long nanoseconds = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
        return nanoseconds / 1000000000.0;
    }

    // RADIX-SORT
    double radix_sort() {
        std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();

        // d - максимальное количество цифр
        int d = getMaxDigits();
        // t ← 1
        int t = 1;
        // for i ← 1 to d
        for (int i = 1; i <= d; i++) {
            // do RX-SORT (A, t)
            rxSort(t);
            // t ← t * 10
            t = t * 10;
        }

        std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
        long long nanoseconds = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
        return nanoseconds / 1000000000.0;
    }
};

int main() {
    // Исходный массив
    std::vector<int> originalArray = {3, 42, 31, 84, 89, 5, 73, 40, 44, 32};

    Sorter sorter(originalArray);
    sorter.printArray("Original array: ");
    std::cout << std::endl;

    // Установка точности вывода
    std::cout << std::fixed << std::setprecision(9);

    // Insertion Sort
    std::cout << "--- INSERTION SORT ---" << std::endl;
    sorter.setArray(originalArray);
    double time1 = sorter.insertion_sort();
    sorter.printArray("Sorted array: ");
    std::cout << time1 << std::endl;
    std::cout << std::endl;

    // Merge Sort
    std::cout << "--- MERGE SORT ---" << std::endl;
    sorter.setArray(originalArray);
    double time2 = sorter.merge_sort();
    sorter.printArray("Sorted array: ");
    std::cout << time2 << std::endl;
    std::cout << std::endl;

    // Radix Sort
    std::cout << "--- RADIX SORT ---" << std::endl;
    sorter.setArray(originalArray);
    double time3 = sorter.radix_sort();
    sorter.printArray("Sorted array: ");
    std::cout << time3 << std::endl;

    return 0;
}

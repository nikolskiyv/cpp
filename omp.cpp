#include <stdio.h>
#include <iostream>
#include <time.h>
#include <omp.h>

using namespace std;

void parallel_multiplication() {
    
    const int n = 10000;
    int i, j;

    float vector[n];
    for (i = 0; i < n; i++)
        vector[i] = 1;

    float result[n] = {};

    float matrix[n][n];
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            matrix[i][j] = 1;
    
    cout << "Enter num of threads: " << endl;
    int th;
    cin >> th;  

    double t_start = clock();

    #pragma omp parallel private(i) num_threads(th)
    {
        #pragma omp parallel for private(j) shared(result, vector, matrix)
            for (i = 0; i < n; i++) 
                for (j = 0; j < n; j++)
                    result[i] += matrix[i][j] * vector[j];
    }  

    double t_end = clock();
    double time = (t_end - t_start) / CLOCKS_PER_SEC;
    cout << "Time: " << time << endl;
}

int main() {
    parallel_multiplication();
}


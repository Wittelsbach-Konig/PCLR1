#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>

# define M_PI 3.14159265358979323846

int main(int argc, char* argv[]) {
    int i, N;
    struct timeval T1, T2;
    long delta_ms;
    N = atoi(argv[1]); /* N равен первому параметру командной строки */
    gettimeofday(&T1, NULL); /* запомнить текущее время T1 */
    int N2 = N/2; /* N2 равен N/2*/
    double M1[N], M2[N2]; /* Массивы M1 разм N и M2 разм N2*/
    double M2_old[N2];
    double A = 490.0; /* Ф*И*О */
    double min = 1; double max = A; double max_2 = max * 10;
    double key, X, intergal_part;
    int j, k, z, min_s;
    unsigned int seed;
    unsigned int* restrict seedp = &seed;
    unsigned int* restrict seedp1 = &seed;
    for (i=0; i<100; ++i) { /* 100 экспериментов */
        srand(i); /* инициализировать начальное значение ГСЧ */
        seed = i;
        /* Заполнить массив исходных данных размером N */
        // GENERATE
        for(j=0; j < N; ++j) {
            M1[j]=((double)rand_r(seedp) / (RAND_MAX)) * (max - min) + min;
        }
        for (k=0; k < N2; ++k) {
            M2[k]=((double)rand_r(seedp1) / (RAND_MAX)) * (max_2 - max) + max;
        }
        // MAP
        for (j = 0; j < N; ++j) {
            M1[j] = pow((M1[j]/M_PI), 3);
        }
        for (k = 0; k < N2; ++k) {
            M2_old[k] = M2[k];
        }
        for(k = 0; k < N2; ++k){
            if(k > 0) {
                M2[k] = M2[k] + M2_old[k-1];
            }
            else {
                M2[k] = M2[k] + 0;
            }
            M2[k] = abs(tan(M2[k]));
        }
        // MERGE
        for(k=0; k < N2; ++k) {
            M2[k]= M1[k] * M2[k];
        }
        // SORT
        for (k = 0; k < (N2-1); ++k) {
            min_s = k;
            for (j = k + 1; j < N2; ++j) {
                if (M2[min_s] > M2[j]) {
                    min_s = j;
                }
            }
            key = M2[min_s];
            for(z = min_s; z > k; --z) {
                M2[z] = M2[z-1];
            }
            M2[k] = key;
        }
        // REDUCE
        X = 0.0;
        key = M2[0];
        for (k = 1; k < N2; ++k) {
            if(M2[k] != 0) {
                if (key == 0 || M2[k] < key) {
                    key = M2[k];
                }
            }
        }
        for(k = 0; k < N2; k++) {
            double del = (double)M2[k] / key;
            modf(del, &intergal_part);
            if (((int)intergal_part % 2) == 0) {
                X= X + sin(M2[k]);
            }
        }
        printf("X = %f", X);
        printf("\n\n");
        /* Решить поставленную задачу, заполнить массив с результатами*/
        /* Отсортировать массив с результатами указанным методом */
    }
    gettimeofday(&T2, NULL); /* запомнить текущее время T2 */
    delta_ms = 1000*(T2.tv_sec - T1.tv_sec) + (T2.tv_usec - T1.tv_usec) / 1000;
    printf("\nN=%d. Milliseconds passed: %ld\n", N, delta_ms); /* T2 - T1 */
    return 0;
}
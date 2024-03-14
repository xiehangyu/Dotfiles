#include <time.h>
#include <immintrin.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <memory.h>

clock_t gemm_baseline(float *A, float *B, float *C, int N);
clock_t gemm_avx(float *A, float *B, float *C, int N);
void gemm_verify_avx(float *A, float *B, float *C, int N);
void gemm_verify_avx_block(float *A, float *B, float *C, int N, int blocksize);
clock_t gemm_avx_block(float *A, float *B, float *C, int N, int blocksize);
void generate_matrix(float *A, float *B, int N);
clock_t gemm_avx_improved(float *A, float *B, float *C, int N);
void gemm_verify_avx_improved(float *A, float *B, float *C, int N);
void transpose8x8(float *src, float *dst, int N);

void printmatrix(float *A, int N)
{
    int i, j;
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            printf("%f ", A[i * N + j]);
        }
        printf("\n");
    }
}
int main(int argc, char **argv)
{
    int n = atoi(argv[1]);
    unsigned seed = (unsigned)time(NULL);
    int blocksize = atoi(argv[2]);
    int N = 1 << n;
    float *A = (float *)malloc(N * N * sizeof(float));
    float *B = (float *)malloc(N * N * sizeof(float));
    float *C = (float *)malloc(N * N * sizeof(float));
    memset(C, 0, N * N * sizeof(float));
    srand(seed);
    generate_matrix(A, B, N);
    gemm_baseline(A, B, C, N);
    memset(C, 0, N * N * sizeof(float));
    srand(seed);
    generate_matrix(A, B, N);
    gemm_avx(A, B, C, N);
    memset(C, 0, N * N * sizeof(float));
    srand(seed);
    generate_matrix(A, B, N);
    gemm_avx_block(A, B, C, N, blocksize);
    memset(C, 0, N * N * sizeof(float));
    srand(seed);
    generate_matrix(A, B, N);
    gemm_avx_improved(A, B, C, N);
}
void generate_matrix(float *A, float *B, int N)
{
    int i, j;
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            A[i * N + j] = (float)rand() / RAND_MAX;
            B[i * N + j] = (float)rand() / RAND_MAX;
        }
    }
}
clock_t gemm_baseline(float *A, float *B, float *C, int N)
{
    clock_t clock1 = clock();
    clock_t clock2;
    int i, j, k;
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            float sum = 0;
            for (k = 0; k < N; k++)
            {
                sum += A[i * N + k] * B[k * N + j];
            }
            C[i * N + j] = sum;
        }
    }
    clock2 = clock();
    printf("Baseline: %ld\n", clock2 - clock1);
    return(clock2 - clock1);
}
inline void transposematrix(float *A, float *B, int N)
{
    if (N < 8)
    {
        int i, j;
        for (i = 0; i < N; i++)
        {
            for (j = 0; j < N; j++)
            {
                B[i * N + j] = A[j * N + i];
            }
        }
    }
    else
    {
        int i, j;
        for (i = 0; i < N; i += 8)
        {
            for (j = 0; j < N; j += 8)
            {
                transpose8x8(&A[i * N + j], &B[j * N + i], N);
            }
        }
    }
}
void transpose8x8(float *src, float *dst, int N)
{
    __m256 row1 = _mm256_loadu_ps(&src[0]);
    __m256 row2 = _mm256_loadu_ps(&src[N]);
    __m256 row3 = _mm256_loadu_ps(&src[2 * N]);
    __m256 row4 = _mm256_loadu_ps(&src[3 * N]);
    __m256 row5 = _mm256_loadu_ps(&src[4 * N]);
    __m256 row6 = _mm256_loadu_ps(&src[5 * N]);
    __m256 row7 = _mm256_loadu_ps(&src[6 * N]);
    __m256 row8 = _mm256_loadu_ps(&src[7 * N]);

    __m256 t1 = _mm256_unpacklo_ps(row1, row2);
    __m256 t2 = _mm256_unpackhi_ps(row1, row2);
    __m256 t3 = _mm256_unpacklo_ps(row3, row4);
    __m256 t4 = _mm256_unpackhi_ps(row3, row4);
    __m256 t5 = _mm256_unpacklo_ps(row5, row6);
    __m256 t6 = _mm256_unpackhi_ps(row5, row6);
    __m256 t7 = _mm256_unpacklo_ps(row7, row8);
    __m256 t8 = _mm256_unpackhi_ps(row7, row8);

    row1 = _mm256_shuffle_ps(t1, t3, _MM_SHUFFLE(1, 0, 1, 0));
    row2 = _mm256_shuffle_ps(t1, t3, _MM_SHUFFLE(3, 2, 3, 2));
    row3 = _mm256_shuffle_ps(t2, t4, _MM_SHUFFLE(1, 0, 1, 0));
    row4 = _mm256_shuffle_ps(t2, t4, _MM_SHUFFLE(3, 2, 3, 2));
    row5 = _mm256_shuffle_ps(t5, t7, _MM_SHUFFLE(1, 0, 1, 0));
    row6 = _mm256_shuffle_ps(t5, t7, _MM_SHUFFLE(3, 2, 3, 2));
    row7 = _mm256_shuffle_ps(t6, t8, _MM_SHUFFLE(1, 0, 1, 0));
    row8 = _mm256_shuffle_ps(t6, t8, _MM_SHUFFLE(3, 2, 3, 2));

    t1 = _mm256_permute2f128_ps(row1, row5, 0x20);
    t2 = _mm256_permute2f128_ps(row2, row6, 0x20);
    t3 = _mm256_permute2f128_ps(row3, row7, 0x20);
    t4 = _mm256_permute2f128_ps(row4, row8, 0x20);
    t5 = _mm256_permute2f128_ps(row1, row5, 0x31);
    t6 = _mm256_permute2f128_ps(row2, row6, 0x31);
    t7 = _mm256_permute2f128_ps(row3, row7, 0x31);
    t8 = _mm256_permute2f128_ps(row4, row8, 0x31);

    _mm256_storeu_ps(&dst[0], t1);
    _mm256_storeu_ps(&dst[N], t2);
    _mm256_storeu_ps(&dst[2 * N], t3);
    _mm256_storeu_ps(&dst[3 * N], t4);
    _mm256_storeu_ps(&dst[4 * N], t5);
    _mm256_storeu_ps(&dst[5 * N], t6);
    _mm256_storeu_ps(&dst[6 * N], t7);
    _mm256_storeu_ps(&dst[7 * N], t8);
}

void gemm_verify_avx(float *A, float *B, float *C, int N)
{
    float *C_copy = (float *)malloc(N * N * sizeof(float));
    double sum = 0;
    double sum2 = 0;
    gemm_baseline(A, B, C_copy, N);
    gemm_avx(A, B, C, N);
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            sum += (C[i * N + j] - C_copy[i * N + j]) * (C[i * N + j] - C_copy[i * N + j]);
            sum2 += (C[i * N + j]) * C[i * N + j];
        }
    }
    printf("Relative Deviation for basic AVX: %lf\n", sqrt(sum / sum2));
}
clock_t gemm_avx(float *A, float *B, float *C, int N)
{
    clock_t clock1 = clock();
    clock_t clock2;
    // generate the transpose matrix for B using AVX
    float *B_copy = (float *)malloc(N * N * sizeof(float));
    transposematrix(B, B_copy, N);
    int i, j, k;
    __m256 a, b, c;
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            float sum = 0;
            float *temp = (float *)malloc(8 * sizeof(float));
            for (k = 0; k < N; k += 8)
            {
                a = _mm256_loadu_ps(A + i * N + k);
                b = _mm256_loadu_ps(B_copy + j * N + k);
                c = _mm256_mul_ps(a, b);
                c = _mm256_hadd_ps(c, c);
                c = _mm256_hadd_ps(c, c);
                temp = (float *)&c;
                sum += temp[0] + temp[4];
            }
            C[i * N + j] = sum;
        }
    }
    free(B_copy);
    clock2 = clock();
    printf("AVX: %ld\n", clock2 - clock1);
    return(clock2 - clock1);
}
clock_t gemm_avx_improved(float *A, float *B, float *C, int N)
{
    clock_t clock1 = clock();
    clock_t clock2;
    for (int i = 0; i < N; ++i)
    {
        for(int k=0; k < N; k++)
        {
            for(int j=0; j < N; j+=8)
            {
                __m256 a = _mm256_set1_ps(A[i*N+k]);
                __m256 b = _mm256_loadu_ps(B+k*N+j);
                __m256 c = _mm256_loadu_ps(C+i*N+j);
                c = _mm256_fmadd_ps(a,b,c);
                _mm256_storeu_ps(C+i*N+j,c);
            }
        }
    }
    clock2=clock();
    printf("AVX Improved: %ld\n", clock2 - clock1);
    return(clock2 - clock1);
}
void gemm_verify_avx_improved(float *A, float *B, float *C, int N)
{
    float *C_copy = (float *)malloc(N * N * sizeof(float));
    double sum = 0;
    double sum2 = 0;
    gemm_baseline(A, B, C_copy, N);
    gemm_avx_improved(A, B, C, N);
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            sum += (C[i * N + j] - C_copy[i * N + j]) * (C[i * N + j] - C_copy[i * N + j]);
            sum2 += (C[i * N + j]) * C[i * N + j];
        }
    }
    printf("Relative Deviation for AVX Improved: %lf\n", sqrt(sum / sum2));
}
void gemm_verify_avx_block(float *A, float *B, float *C, int N, int blocksize)
{
    float *C_copy = (float *)malloc(N * N * sizeof(float));
    double sum = 0;
    double sum2 = 0;
    gemm_baseline(A, B, C_copy, N);
    gemm_avx_block(A, B, C, N, blocksize);
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            sum += (C[i * N + j] - C_copy[i * N + j]) * (C[i * N + j] - C_copy[i * N + j]);
            sum2 += (C[i * N + j]) * C[i * N + j];
        }
    }
    printf("Relative Deviation for block: %lf\n", sqrt(sum / sum2));
}

clock_t gemm_avx_block(float *A, float *B, float *C, int N, int blocksize)
{
    clock_t clock1 = clock();
    clock_t clock2;
    // generate the transpose matrix for B using AVX
    float *B_copy = (float *)malloc(N * N * sizeof(float));
    transposematrix(B, B_copy, N);
    int i, j, k;
    __m256 a, b, c;
    float *temp;
    for (i = 0; i < N; i += blocksize)
    {
        for (j = 0; j < N; j += blocksize)
        {
            for (k = 0; k < N; k += blocksize)
            {
                for (int ii = i; ii < i + blocksize; ii++)
                {
                    for (int jj = j; jj < j + blocksize; jj++)
                    {
                        float sum = 0;
                        for (int kk = k; kk < k + blocksize; kk += 8)
                        {
                            a = _mm256_loadu_ps(A + ii * N + kk);
                            b = _mm256_loadu_ps(B_copy + jj * N + kk);
                            c = _mm256_mul_ps(a, b);
                            c = _mm256_hadd_ps(c, c);
                            c = _mm256_hadd_ps(c, c);
                            temp = (float *)&c;
                            sum += temp[0] + temp[4];
                        }
                        C[ii * N + jj] += sum;
                    }
                }
            }
        }
    }
    free(B_copy);
    clock2 = clock();
    printf("AVX block: %ld\n", clock2 - clock1);
    return(clock2 - clock1);
}
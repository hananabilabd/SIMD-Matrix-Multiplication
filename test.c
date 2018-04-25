#include <xmmintrin.h> //Header for intrinsic functions
#include <stdio.h>
#include <time.h>
#include <x86intrin.h>

int main()
{
// Variables definition and initialization
__m128 acc = _mm_setzero_ps();

int MAX_DIM = 100;
float a[MAX_DIM][MAX_DIM] __attribute__ ((aligned(16)));
float b[MAX_DIM][MAX_DIM] __attribute__ ((aligned(16)));
float c[MAX_DIM][MAX_DIM] __attribute__ ((aligned(16)));
float d[MAX_DIM][MAX_DIM] __attribute__ ((aligned(16)));
for (int i = 0; i < MAX_DIM; ++i)
{
	for (int j = 0; j < MAX_DIM; ++j)
	{
		a[i][j] = 13; // Arbitrary initializations - Replace to test your multiplication
		b[i][j] = 22; // Arbitrary initializations - Replace to test your multiplication
		c[i][j] = 0; // Necessary Initialization - Don't change
		d[i][j] = 0; // Necessary Initialization - Don't change
	}
}
// Unoptimized Matrix Multiplication
clock_t Time1 = clock();
for (int i = 0; i < MAX_DIM; ++i)
{
	for (int j = 0; j < MAX_DIM; ++j)
	{
		for (int k = 0; k < MAX_DIM; k++)
		{
            c[i][j] += a[i][k] * b[k][j];
        }
    }
}
clock_t Time2 = clock();
//=================================================================================================
clock_t Time3 = clock();
/*
float *p1 ;float *p2;float *p3;
p1 =&a[0];
p2 =&b[0];
p3 =&c[0];
float* m; float* v; float* result;
m =&a[0];
v =&b[0];
result =&c[0];

*/
//int k =0;
 for(int i = 0; i < MAX_DIM; i++)
    {
        for(int j = 0; j < MAX_DIM; j++)
        {
            for(int k = 0; k < MAX_DIM; k+=4)
            {
                
                __m128 aa = _mm_load_ps(&(a[i][k]));

                __m128 bb = _mm_load_ps(&(b[k][j]));
                
                
                __m128 result  = _mm_mul_ps(aa, bb);

                acc = _mm_add_ps(acc,result);

                
                
            }
            
            d[i][j] = (acc[0] + acc[1] + acc[2] + acc[3]);

            acc = _mm_setzero_ps();

        }
    }
// YOUR CODE HERE
clock_t Time4 = clock();
// Calculate and print execution times
double TotalTimeLoop = ((double) Time2 - (double) Time1) / CLOCKS_PER_SEC;
double TotalTimeSIMD = ((double) Time4 - (double) Time3) / CLOCKS_PER_SEC;
printf(" Time taken by loop is %.7f \n", TotalTimeLoop);
printf(" Time taken by SIMD optimized code is %.7f \n", TotalTimeSIMD);
printf("Performance gain is x%0.7f \n\n", TotalTimeLoop/TotalTimeSIMD);
return 0;
}
#include <xmmintrin.h> //Header for intrinsic functions
#include <stdio.h>
#include <time.h>
#include <x86intrin.h>
#define Enable 1 
#define Disable 0
#define Single_Precision Enable  
#define Double_Precision Disable 
#define Integer Disable 

#if Single_Precision == Enable 
int main()
{
// Variables definition and initialization
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
for(int i=0;i<MAX_DIM;i+=1){

      for(int j=0;j<MAX_DIM;j+=4){

        for(int k=0;k<MAX_DIM;k+=4){

          __m128 result = _mm_load_ps(&d[i][j]);

          __m128 a_line  = _mm_load_ps(&a[i][k]);

          __m128 b_line0 = _mm_load_ps(&b[k][j+0]);

          __m128 b_line1 = _mm_loadu_ps(&b[k][j+1]);

          __m128 b_line2 = _mm_loadu_ps(&b[k][j+2]);

          __m128 b_line3 = _mm_loadu_ps(&b[k][j+3]);

         result = _mm_add_ps(result, _mm_mul_ps(_mm_shuffle_ps(a_line, a_line, 0x00), b_line0));
         result = _mm_add_ps(result, _mm_mul_ps(_mm_shuffle_ps(a_line, a_line, 0x55), b_line1));
         result = _mm_add_ps(result, _mm_mul_ps(_mm_shuffle_ps(a_line, a_line, 0xaa), b_line2));
         result = _mm_add_ps(result, _mm_mul_ps(_mm_shuffle_ps(a_line, a_line, 0xff), b_line3));
         _mm_store_ps(&d[i][j],result);
        }
      }
    }

// YOUR CODE HERE
clock_t Time4 = clock();
// Calculate and print execution times
double TotalTimeLoop = ((double) Time2 - (double) Time1) / CLOCKS_PER_SEC;
double TotalTimeSIMD = ((double) Time4 - (double) Time3) / CLOCKS_PER_SEC;
printf("Increase in Performance is %.7f \n", TotalTimeLoop/TotalTimeSIMD);
printf(" Time taken by loop is %.7f \n", TotalTimeLoop);
printf(" Time taken by SIMD optimized code is %.7f \n", TotalTimeSIMD);
return 0;
}


#elif Double_Precision == Enable 

#include <xmmintrin.h> //Header for intrinsic functions
#include <stdio.h>
#include <time.h>
#include <x86intrin.h>

int main()
{
// Variables definition and initialization
int MAX_DIM = 100;
double a[MAX_DIM][MAX_DIM] __attribute__ ((aligned(16)));
double b[MAX_DIM][MAX_DIM] __attribute__ ((aligned(16)));
double c[MAX_DIM][MAX_DIM] __attribute__ ((aligned(16)));
double d[MAX_DIM][MAX_DIM] __attribute__ ((aligned(16)));
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
for(int i=0;i<MAX_DIM;i+=1){

      for(int j=0;j<MAX_DIM;j+=4){

        for(int k=0;k<MAX_DIM;k+=2){

          __m128d result = _mm_load_pd(&d[i][j]);

          __m128d a_line  = _mm_load_pd(&a[i][k]);

          __m128d b_line0 = _mm_load_pd(&b[k][j+0]);

          __m128d b_line1 = _mm_loadu_pd(&b[k][j+1]);

          __m128d b_line2 = _mm_loadu_pd(&b[k][j+2]);

          __m128d b_line3 = _mm_loadu_pd(&b[k][j+3]);

         result = _mm_add_pd(result, _mm_mul_pd(_mm_shuffle_pd(a_line, a_line, 0x00), b_line0));
         result = _mm_add_pd(result, _mm_mul_pd(_mm_shuffle_pd(a_line, a_line, 0x55), b_line1));
         result = _mm_add_pd(result, _mm_mul_pd(_mm_shuffle_pd(a_line, a_line, 0xaa), b_line2));
         result = _mm_add_pd(result, _mm_mul_pd(_mm_shuffle_pd(a_line, a_line, 0xff), b_line3));
         _mm_store_pd(&d[i][j],result);
        }
      }
    }

// YOUR CODE HERE
clock_t Time4 = clock();
// Calculate and print execution times
double TotalTimeLoop = ((double) Time2 - (double) Time1) / CLOCKS_PER_SEC;
double TotalTimeSIMD = ((double) Time4 - (double) Time3) / CLOCKS_PER_SEC;
printf("Performance gain is %.7f \n", TotalTimeLoop/TotalTimeSIMD);
printf(" Time taken by loop is %.7f \n", TotalTimeLoop);
printf(" Time taken by SIMD optimized code is %.7f \n", TotalTimeSIMD);
return 0;
}


#elif Integer == Enable 

#include <xmmintrin.h> //Header for intrinsic functions
#include <stdio.h>
#include <time.h>
#include <x86intrin.h>
#include <smmintrin.h>
int main()
{
// Variables definition and initialization
int MAX_DIM = 100;
int a[MAX_DIM][MAX_DIM] __attribute__ ((aligned(16)));
int b[MAX_DIM][MAX_DIM] __attribute__ ((aligned(16)));
int c[MAX_DIM][MAX_DIM] __attribute__ ((aligned(16)));
int d[MAX_DIM][MAX_DIM] __attribute__ ((aligned(16)));
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
for(int i=0;i<MAX_DIM;i+=1){

      for(int j=0;j<MAX_DIM;j+=4){

        for(int k=0;k<MAX_DIM;k+=2){

          __m128i result = _mm_load_si128(&d[i][j]);
          
          __m128i a_line  = _mm_load_si128(&a[i][k]);
          
          __m128i b_line0 = _mm_load_si128(&b[k][j+0]);
        
          __m128i b_line1 = _mm_load_si128(&b[k][j+1]);
        
          __m128i b_line2 = _mm_load_si128(&b[k][j+2]);
     
          __m128i b_line3 = _mm_load_si128(&b[k][j+3]);

         result = _mm_add_epi16(result, _mm_mul_ps(_mm_shuffle_ps(a_line, a_line, 0x00), b_line0));
         result = _mm_add_epi16(result, _mm_mul_ps(_mm_shuffle_ps(a_line, a_line, 0x55), b_line1));
         result = _mm_add_epi16(result, _mm_mul_ps(_mm_shuffle_ps(a_line, a_line, 0xaa), b_line2));
         result = _mm_add_epi16(result, _mm_mul_ps(_mm_shuffle_ps(a_line, a_line, 0xff), b_line3));
         _mm_store_si128(&d[i][j],result);
        }
      }
    }

// YOUR CODE HERE
clock_t Time4 = clock();
// Calculate and print execution times
double TotalTimeLoop = ((double) Time2 - (double) Time1) / CLOCKS_PER_SEC;
double TotalTimeSIMD = ((double) Time4 - (double) Time3) / CLOCKS_PER_SEC;
printf("Performance gain is %.7f \n", TotalTimeLoop/TotalTimeSIMD);
printf(" Time taken by loop is %.7f \n", TotalTimeLoop);
printf(" Time taken by SIMD optimized code is %.7f \n", TotalTimeSIMD);
return 0;
}
#endif
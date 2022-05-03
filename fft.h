#ifndef __FFT_H__
#define __FFT_H__

typedef struct //複數型別
{
  float real;		//實部
  float imag;		//虛部
}complex;

#define PI 3.1415926535897932384626433832795028841971


///////////////////////////////////////////
void conjugate_complex(int n,complex in[],complex out[]);
float c_value(complex f);								//複數取模
void c_plus(complex a,complex b,complex *c);			//複數加
void c_mul(complex a,complex b,complex *c) ;			//複數乘
void c_sub(complex a,complex b,complex *c);				//複數減法
void c_div(complex a,complex b,complex *c);				//複數除法
void fft(int length,complex f[]);							//傅立葉變換 輸出也存在陣列f中
void ifft(int length,complex f[]); 						// 傅立葉逆變換
void c_abs(complex f[],float out[],int n);				//複數陣列取模
////////////////////////////////////////////
#endif
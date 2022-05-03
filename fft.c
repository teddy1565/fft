#include "math.h"
#include "fft.h"

void conjugate_complex(int n,complex in[],complex out[])
{
	int i = 0;
	for(i=0;i<n;i++)
	{
		out[i].imag = -in[i].imag;
		out[i].real = in[i].real;
	}
}

void c_abs(complex f[],float out[],int n)
{
	int i = 0;
	float t;
	for(i=0;i<n;i++)
	{
		t = f[i].real * f[i].real + f[i].imag * f[i].imag;
		out[i] = sqrt(t);
	}
}

float c_value(complex f)
{
	return f.real * f.real + f.imag * f.imag;
}


void c_plus(complex a,complex b,complex *c)
{
	c->real = a.real + b.real;
	c->imag = a.imag + b.imag;
}

void c_sub(complex a,complex b,complex *c)
{
	c->real = a.real - b.real;
	c->imag = a.imag - b.imag;
}

void c_mul(complex a,complex b,complex *c)
{
	c->real = a.real * b.real - a.imag * b.imag;
	c->imag = a.real * b.imag + a.imag * b.real;
}

void c_div(complex a,complex b,complex *c)
{
	c->real = (a.real * b.real + a.imag * b.imag)/(b.real * b.real +b.imag * b.imag);
	c->imag = (a.imag * b.real - a.real * b.imag)/(b.real * b.real +b.imag * b.imag);
}

#define SWAP(a,b)  tempr=(a);(a)=(b);(b)=tempr

void Wn_i(int n,int i,complex *Wn,char flag)
{
	Wn->real = cos(2*PI*i/n);
	if(flag == 1)
		Wn->imag = -sin(2*PI*i/n);
	else if(flag == 0)
		Wn->imag = -sin(2*PI*i/n);
}

//傅立葉變化
void fft(int length,complex f[])
{
	complex t,wn;//中間變數
	int i,j,k,m,n,l,r,M;
	int la,lb,lc;

	/*----計算分解的級數M=log2(length)----*/
	for(i=length,M=1;(i=i/2)!=1;M++);

	/*----按照倒位序重新排列原訊號----*/
	for(i=1,j=length/2;i<=length-2;i++)
	{
		if(i<j)
		{
			t=f[j];
			f[j]=f[i];
			f[i]=t;
		}
		k=length/2;
		while(k<=j)
		{
			j=j-k;
			k=k/2;
		}
		j=j+k;
	}

	/*----FFT演算法----*/
	for(m=1;m<=M;m++)
	{
		la=pow(2,m); 	//la=2^m代表第m級每個分組所含節點數
		lb=la/2;    	//lb代表第m級每個分組所含碟形單元數

		//同時它也表示每個碟形單元上下節點之間的距離
		/*----碟形運算----*/
		for(l=1;l<=lb;l++)
		{
			r=(l-1)*pow(2,M-m);
			for(n=l-1;n<length-1;n=n+la) //遍歷每個分組，分組總數為N/la
			{
				lc=n+lb;  //n,lc分別代表一個碟形單元的上、下節點編號
				Wn_i(length,r,&wn,1);		//wn=Wnr
				c_mul(f[lc],wn,&t);			//t = f[lc] * wn複數運算
				c_sub(f[n],t,&(f[lc]));		//f[lc] = f[n] - f[lc] * Wnr
				c_plus(f[n],t,&(f[n]));		//f[n] = f[n] + f[lc] * Wnr
			}
		}
	}
}

//傅立葉逆變換
void ifft(int length,complex f[])
{
	int i=0;
	conjugate_complex(length,f,f);
	fft(length,f);
	conjugate_complex(length,f,f);
	for(i=0;i<length;i++)
	{
		f[i].imag = (f[i].imag)/length;
		f[i].real = (f[i].real)/length;
	}
}
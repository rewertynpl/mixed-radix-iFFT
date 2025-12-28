# mixed-radix-iFFT
universal mixed radix fast fourier transform FFT c++ source code radix-2 radix-3 radix-4 radix-5 radix-7 radix-11 
created 29.6.2017

author copyright marcin matysek (r)ewertyn.PL

marcin.rewertyn@gmail.com

open-source
#
http://www.mediafire.com/file/hyz4dbski4w00pb/inverse+fourier+transform+iDFT+ifft+4+methods+in+open+office+++something+extra.ods
#
http://www.mediafire.com/file/59bpnci966ulec9/DFT+FFT+RADIX-2+DIT+algorytm+Transformacja+Fouriera+analitycznie+v3.4.xlsx 
#
    //assumption: if N==signal period in table tab[] then resolution = 1 Hz but N=2^b;
    //when we have signal period 22000 Hz and N=2^15=32768 then MAYBE
    //the fundamental frequency is 0,671 Hz
    //that means that in F(1) is 1*0,671 Hz =0,671 Hz
    //in F(2) is 2*0,671 Hz =1,342 Hz
    //in F(9) is 9*0,671 Hz =6,039 Hz
    //in F(30) is 30*0,671 Hz =20,13 Hz that means in F(30) we will see what value has our signal in 20,13 Hz



//-------------------------------<br />
//when you want to have equal results that are in that false modificator in normal FFT then change this:<br /><br /><br />
/*<br /><br />
 void fun_fourier_transform_FFT_mixed_radix(int N,std::complex<double> tab[])<br />
{<br />
	//new:<br />
    for(int j=0;j<N;j++)<br />
    {<br />
     tab[j].real(tab[j].real()*2/N);<br />
     tab[j].imag(tab[j].imag()*2/N);<br />
    }<br />
}<br />
//and:<br />
void fun_inverse_fourier_transform_FFT_mixed_radix(int N,std::complex<double> tab[])<br />
{<br />
	//new:<br />
    for(int j=0;j<N;j++)<br />
    {<br />
     tab[j].real(tab[j].real()*0.5);<br />
     tab[j].imag(tab[j].imag()*0.5);<br />
    }<br />
}<br />
///<br />
///<br />
//for official modificator that is only in inverse FFT:<br />

 void fun_fourier_transform_FFT_mixed_radix(int N,std::complex<double> tab[])<br />
{ <br />
    for(int j=0;j<N;j++)<br />
    {<br />
      //nothing<br />
    }<br />
}<br />
void fun_inverse_fourier_transform_FFT_mixed_radix(int N,std::complex<double> tab[])<br />
{<br />
    for(int j=0;j<N;j++)<br />
    {<br />
     tab[j].real(tab[j].real()*1/(float)N);//??<br />
     tab[j].imag(tab[j].imag()*1/(float)N);//??<br />
    }<br />
<br />
}<br />
*/<br />

update: December 2025

# Mixed-Radix FFT/iFFT with Phase Shift

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A flexible, open-source C++ implementation of the Mixed-Radix Fast Fourier Transform (FFT) and its inverse (iFFT). Designed for educational purposes, it supports signals of arbitrary composite lengths and includes a unique phase-shift feature.

## 📖 About The Project

This project provides a standalone C++ implementation of the **Mixed-Radix Cooley-Tukey algorithm**. Unlike standard Radix-2 FFTs, which require the signal length to be a power of two ($2^n$), this algorithm supports any length $N$ that can be factorized into a product of small prime numbers.

The code is written to be explicit and verbose, serving as a learning tool to demonstrate how "butterfly" computations work for various radixes.

### Key Features

*   **Mixed-Radix Support:** Works with any signal length $N$ factorizable into implemented radixes.
*   **Wide Radix Base:** Native support for Radix-2, 3, 4, 5, 7, 11, and a generic implementation for primes up to 97.
*   **Phase Shift (`fi`):** A unique global parameter allowing for circular signal shifting directly within the frequency domain transform.
*   **Verification:** Includes a generic $O(N^2)$ DFT implementation for accuracy verification.
*   **License:** MIT License (Free for both Open Source and Commercial use).

## ⚙️ Technical Specifications

### The Algorithm
The implementation utilizes a **Decimation-in-Frequency (DIF)** approach. The signal length $N$ is factorized as $N = R_1 \cdot R_2 \cdot \dots \cdot R_m$. The transform executes in $m$ stages, performing smaller DFTs of size $R_i$ at each step.

### Phase Shift Parameter (`fi`)
The code includes a global variable `double fi`.
*   **Function:** It introduces a constant phase offset to the twiddle factors during the first stage of the FFT.
*   **Effect:** Mathematically, a linear phase shift in the frequency domain corresponds to a **circular time-shift** in the time domain. This is useful for signal alignment or modulation simulation without modifying the input buffer directly.

## 🚀 Getting Started

### Prerequisites
*   A C++ compiler (GCC, Clang, MSVC).
*   Standard C++ libraries (`<complex>`, `<cmath>`, `<iostream>`).

### Usage Example

```cpp
#include <iostream>
#include <complex>
// Include the source file or header containing the FFT functions

int main() {
    // Example: N = 2 * 3 * 5 = 30 (Composite number, not a power of 2)
    int N = 30; 
    std::complex<double>* signal = new std::complex<double>[N];

    // 1. Initialize signal
    for(int i = 0; i < N; i++) {
        signal[i] = { (double)i, 0.0 }; // Simple ramp signal
    }

    // 2. Compute FFT
    std::cout << "Calculating FFT..." << std::endl;
    fun_fourier_transform_FFT_mixed_radix(N, signal);
    
    // 3. Reorder Output (Digit-Reversal Permutation)
    // This step is mandatory to get the spectrum in natural order
    fun_inverse_table_FFT(N, signal);

    // 4. Compute Inverse FFT (iFFT)
    std::cout << "Calculating iFFT..." << std::endl;
    fun_inverse_fourier_transform_FFT_mixed_radix(N, signal);
    
    // 5. Reorder Output again to get time-domain signal
    fun_inverse_table_FFT(N, signal);

    // Clean up
    delete[] signal;
    return 0;
}


In this source code i have my own implementation of inverse_bits that is universal for all radix:

# FFT Inverse Table Permutation

This module provides a C++ implementation for permuting an array of complex numbers. It is a key step in Fast Fourier Transform (FFT) algorithms, specifically performing the index reordering (digit-reversal) necessary for mixed-radix FFT implementations.

## Functionality

The `fun_inverse_table_FFT` function:
1. Accepts the array size `M` and an array of complex numbers `tab`.
2. Calculates the new index order based on the prime factor decomposition (radix).
3. Permutes the elements of `tab` in place using temporary buffers.

## Code

```cpp
#include <complex>
#include <iostream>

// NOTE: This function requires the definition of a helper function 'radix_base'.
// int radix_base(int M, int* stg);

void fun_inverse_table_FFT(int M, std::complex<double> tab[])
{
    int rx5=5,rx4=4,rx3=3,rx2=2,rx7=7,rx11=11;
    int stg[100]={};
    int *tab8 = new int[M];
    int *tab9 = new int[M];
    std::complex<double> *tab11 = new std::complex<double>[M];
    int nb_stages=0;
    int nb1=0;
    int nb2=0;
    int nb3=0;

    nb_stages=6;

    // External function to determine radices
    nb_stages=radix_base(M,stg);

    for(int i=0;i<M;i++)
    {
        //tab9[i]=tab2[i];
        //tab8[i]=tab2[i];
        tab9[i]=i;
        tab8[i]=i;
    }

    nb3=1;
    for(int op=nb_stages;op>=2;op--)
    {
        nb1=stg[op];
        nb3=nb3*stg[op];

        if(op==nb_stages)
        {
            nb2=stg[0];
        }
        else
        {
            nb2=nb2*stg[op+1];
        }

           for(int i=0,n=0,p=0;i<M;i=i+M/nb3,n++)
        {
            if(n>=nb1)
            {
                n=0,p=p+M/nb2;
            }
            for(int j=0,k=0;j<M/nb3;j++,k=k+nb1)
            {
                if(op%2==0)
                {
                    tab8[i+j]=tab9[k+n+p];
                }
                else
                {
                    tab9[i+j]=tab8[k+n+p];
                }
            }
        }
    }

    for(int i=0;i<M;i++)
    {
      tab11[i]=tab[tab8[i]];

    }
    for(int i=0;i<M;i++)
    {
      tab[i]=tab11[i];
    }

    delete [] tab8;
    delete [] tab9;
    delete [] tab11;
}

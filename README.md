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

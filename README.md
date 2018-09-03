# mixed-radix-iFFT
universal mixed radix fast fourier transform FFT c++ source code radix-2 radix-3 radix-4 radix-5 radix-7 radix-11 
created 29.6.2017

author copyright marcin matysek (r)ewertyn.PL

marcinmatysek1@gmail.com

open-source
#
http://www.mediafire.com/file/hyz4dbski4w00pb/inverse+fourier+transform+iDFT+ifft+4+methods+in+open+office+++something+extra.ods
#
http://www.mediafire.com/file/59bpnci966ulec9/DFT+FFT+RADIX-2+DIT+algorytm+Transformacja+Fouriera+analitycznie+v3.4.xlsx 
#
    //assumption: if N==signal period in table tab[] then resolution = 1 Hz but N=2^b;
    //when we have signal period 22000 Hz and N=2^15=32768 then
    //the fundamental frequency is 0,671 Hz
    //that means that in F(1) is 1*0,671 Hz =0,671 Hz
    //in F(2) is 2*0,671 Hz =1,342 Hz
    //in F(9) is 9*0,671 Hz =6,039 Hz
    //in F(30) is 30*0,671 Hz =20,13 Hz that means in F(30) we will see what value has our signal in 20,13 Hz


//-------------------------------
//when you want to have equal results that are in that false modificator in normal FFT then change this:
/*
 void fun_fourier_transform_FFT_mixed_radix(int N,std::complex<double> tab[])
{
	//new:
    for(int j=0;j<N;j++)
    {
     tab[j].real(tab[j].real()*2/N);
     tab[j].imag(tab[j].imag()*2/N);
    }
}
//and:
void fun_inverse_fourier_transform_FFT_mixed_radix(int N,std::complex<double> tab[])
{
	//new:
    for(int j=0;j<N;j++)
    {
     tab[j].real(tab[j].real()*0.5);
     tab[j].imag(tab[j].imag()*0.5);
    }
}
///
///
//for official modificator that is only in inverse FFT:

 void fun_fourier_transform_FFT_mixed_radix(int N,std::complex<double> tab[])
{ 
    for(int j=0;j<N;j++)
    {
      //nothing
    }
}
void fun_inverse_fourier_transform_FFT_mixed_radix(int N,std::complex<double> tab[])
{
    for(int j=0;j<N;j++)
    {
     tab[j].real(tab[j].real()*1/(float)N);//??
     tab[j].imag(tab[j].imag()*1/(float)N);//??
    }

}
*/

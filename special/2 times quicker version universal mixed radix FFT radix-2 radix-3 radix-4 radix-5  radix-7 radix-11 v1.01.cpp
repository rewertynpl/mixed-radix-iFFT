


//mixed radix-11 , radix-7 , radix-5  , radix-4  , radix-3  , radix-2 FFT is 2 times quicker that previous version because it calculates FFT witch tables instead of complex number objects??

//working sound analysis program for windows 8.1, v 1.01

    //source:
    //https://www.google.ch/patents/US6957241
    //http://www.google.ch/patents/US20020083107
    //https://www.beechwood.eu/fft-implementation-r2-dit-r4-dif-r8-dif/
    //http://www.cmlab.csie.ntu.edu.tw/cml/dsp/training/coding/transform/fft.html
    //http://dsp.stackexchange.com/questions/3481/radix-4-fft-implementation

    //https://community.arm.com/graphics/b/blog/posts/speeding-up-fast-fourier-transform-mixed-radix-on-mobile-arm-mali-gpu-by-means-of-opencl---part-1
    //book: "Cyfrowe przetwarzanie sygnalow" - Tomasz Zielinski it has quick version of radix-2 because it calculates less sin() and cos()

    //Rabiner L.R., Gold B. Theory and application of digital signal processing p 378 mixed radix

// author marcin matysek (r)ewertyn.PL
// open source
//http://mixedradixfastfouriertransformifft.blogspot.com/



#include <iostream>
#include "conio.h"
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <time.h>
#include <complex>
#include <fstream>

using namespace std;

//complex number method:

void fun_fourier_transform_FFT_radix_2_stg_rest(std::complex<double> tab[],int nb1,int nb2,int nb3,int nb4,int nb5_stg_previous,double tmp6,double (*z)[2]);
void fun_fourier_transform_FFT_radix_3_stg_rest(std::complex<double> tab[],int nb1,int nb2,int nb3,int nb4,int nb5_stg_previous,double tmp6,double (*z)[2]);
void fun_fourier_transform_FFT_radix_4_stg_rest(std::complex<double> tab[],int nb1,int nb2,int nb3,int nb4,int nb5_stg_previous,double tmp6,double (*z)[2]);
void fun_fourier_transform_FFT_radix_5_stg_rest(std::complex<double> tab[],int nb1,int nb2,int nb3,int nb4,int nb5_stg_previous,double tmp6,double (*z)[2]);
void fun_fourier_transform_FFT_radix_7_stg_rest(std::complex<double> tab[],int nb1,int nb2,int nb3,int nb4,int nb5_stg_previous,double tmp99,double (*z)[2]);
void fun_fourier_transform_FFT_radix_11_stg_rest(std::complex<double> tab[],int nb1,int nb2,int nb3,int nb4,int nb5_stg_previous,double tmp99,double (*z)[2]);

void fun_fourier_transform_FFT_radix_2_stg_first(std::complex<double> tab[],int stg_first,double fi2,double (*z)[2]);
void fun_fourier_transform_FFT_radix_3_stg_first(std::complex<double> tab[],int stg_first,double fi2,double (*z)[2]);
void fun_fourier_transform_FFT_radix_4_stg_first(std::complex<double> tab[],int stg_first,double fi2,double (*z)[2]);
void fun_fourier_transform_FFT_radix_5_stg_first(std::complex<double> tab[],int stg_first,double fi2,double (*z)[2]);
void fun_fourier_transform_FFT_radix_7_stg_first(std::complex<double> tab[],int stg_first,double fi2,double (*z)[2]);
void fun_fourier_transform_FFT_radix_11_stg_first(std::complex<double> tab[],int stg_first,double fi2,double (*z)[2]);

void fun_inverse_table_FFT(int M,std::complex<double> tab[]);
int radix_base(int N,int stg[]);

void fun_fourier_transform_FFT_mixed_radix_for_N_1848_points(int N,std::complex<double> tab[]);
void fun_fourier_transform_DFT(int N,std::complex<double> tab[]);




////////////////////////////////////////////////////
/////////////////////
double fi=0.00;//for first stage FFT and for DFT
/////////////////////
///////////////////////////////////////////////////




int tb[600][3]={};


static double diffclock(clock_t clock2,clock_t clock1)
{
    double diffticks=clock1-clock2;
    double diffms=(diffticks)/(CLOCKS_PER_SEC/1000);
    return diffms;
}
int main()
{
    int N;
    //if N==period of signal in table tab[] then resolution = 1 Hz

srand( time( NULL ) );

ifstream file1("11.txt");
ofstream file3("1.txt");
ofstream file4("2.txt");
int a1=0,b1=0;

int aaa=0;
for(int i=0;i<100;i++)
{

    file1>>aaa;
    file1>>aaa;
    file1>>tb[i][0];
    file1>>tb[i][1];
    file1>>tb[i][2];
  //  cout<<aaa<<" "<<tb[i][0]<<" "<<tb[i][1]<<" "<<tb[i][2]<<endl;system("pause");
}





    /////////////////////
    N=2*2*2*2*2;
    ////////////////////






    double c1=0,c2=0;
    std::complex<double> *tab2=new std::complex<double>[N];
    std::complex<double> *tab3=new std::complex<double>[N];
    std::complex<double> *tab4=new std::complex<double>[N];
    std::complex<double> tmp1=0;

    const double pi=3.141592653589793238462;

    for(int i=0;i<N;i++)
    {
          c1=rand()%1255/4.0;
          //  c1=255;
          //  c2=255;
          c2=rand()%13591/9.0;
          //c1=4*sin(1*2*pi*i/float(N));
          //c2=5435*sin(21*2*pi*i/float(N));
          //c1=c1+33.19*sin(7*2*pi*i/float(N));
          //c2=0;
        tab2[i].real(c1);
        tab2[i].imag(c2);
        tab3[i].real(tab2[i].real());
        tab3[i].imag(tab2[i].imag());
        tab4[i].real(tab2[i].real());
        tab4[i].imag(tab2[i].imag());
    }


    double time2;
    double zmienna=0;

    cout<<"signal 1="<<endl;
    cout<<"tab2[j].real():"<<endl;
    for(int j=0;j<N;j++)
    {
    cout.precision(4);
    cout<<round(tab2[j].real()*1000)/1000<<"  ";
    }
    cout<<endl;
    cout<<endl;
    cout<<"tab2[j].imag():"<<endl;
    for(int j=0;j<N;j++)
    {
    cout.precision(4);
    cout<<round(tab2[j].imag()*1000)/1000<<"  ";
    }
    cout<<endl;
    cout<<endl;

    cout<<"signal 2="<<endl;
    cout<<"tab3[j].real():"<<endl;
    for(int j=0;j<N;j++)
    {
    cout.precision(4);
    cout<<round(tab3[j].real()*1000)/1000<<"  ";
    }
    cout<<endl;
    cout<<endl;
    cout<<"tab3[j].imag():"<<endl;
    for(int j=0;j<N;j++)
    {
    cout.precision(4);
    cout<<round(tab3[j].imag()*1000)/1000<<"  ";
    }
    cout<<endl;
    cout<<endl;
    clock_t start = clock();

      cout<<"calculating first: DFT"<<endl;
    start = clock();
    fun_fourier_transform_DFT(N,tab3);
    time2=diffclock( start, clock() );
    cout<<"time: [s] "<<time2/1000<<endl;
    system("pause");
    //////////////////////////////////////////////////////////
     cout<<"frequency Hz from DFT"<<endl;
    cout<<"tab3[j].real():"<<endl;
    for(int j=0;j<N;j++)
    {
    cout.precision(4);
    cout<<round(tab3[j].real()*1000)/1000<<"  ";
    }
    cout<<endl;
    cout<<endl;//system("pause");
    cout<<"tab3[j].imag():"<<endl;
    for(int j=0;j<N;j++)
    {
    cout.precision(4);
    cout<<round(tab3[j].imag()*1000)/1000<<"  ";
    }
    cout<<endl;
    cout<<endl;//system("pause");
    system("pause");
int mx=0;
int fg=0;

    for(int i=0;i<N;i++)
    {
        tab2[i].real(tab4[i].real());
        tab2[i].imag(tab4[i].imag());
    }

    cout<<"calculating: FFT"<<endl;
    start = clock();
    fun_fourier_transform_FFT_mixed_radix_for_N_1848_points(N,tab2);
    time2=diffclock( start, clock() );
    cout<<"time: [ms] "<<time2<<endl;
    cout<<"inverse table"<<endl;
    start = clock();
    fun_inverse_table_FFT(N,tab2);
    time2=diffclock( start, clock() );
    cout<<"time: [ms] "<<time2<<endl;
    system("pause");
    ///////////////////////////////////////////////////////////

/////////////////////////////////////////////////////
    time2=diffclock( start, clock() );


    int flag=0,fl=0;
       for(int j=0;j<N;j++)
        {
            for(int i=0;i<N;i++)
        {
          //if((fabs(tab3[i].real()-tab2[j].real())<=0.03)&&fabs(tab3[0].real()-tab2[0].real())<0.03)
          if((fabs(tab3[i].real()-tab2[j].real())<=0.03)&&(fabs(tab3[i].imag()-tab2[j].imag())<=0.03))
          {

              flag++;
/*
            for(int k=0;k<N;k++)
              {
                  if(fabs(tab3[k].real()-tab2[j].real())<=0.03)
                  {fl++;}

              }
              if(fl==1){
                flag++;}
                fl=0;
.*/
            cout.precision(4);
            //cout<<"DFT["<<i<<"]==FFT["<<j<<"]= "<<round(tab2[j].real()*1000)/1000<<"  ";//system("pause");
          }
            else {
               // cout<<-1<<" . ";//system("pause");
            }
        }
        }
cout<<"if DFT[i]==FFT[j]"<<endl;
    //cout<<"flag="<<flag<<endl;
    if(flag>=N/2){
cout<<endl<<"flag= "<<flag<<endl;
fg++;


/*
    cout<<"frequency Hz FFT radix-4 radix-2"<<endl;
    for(int j=0;j<N;j++)
    {
    cout.precision(4);
    cout<<round(tab2[j].real()*1000)/1000<<"  ";
    }
    cout<<endl;
    cout<<endl;
    for(int j=0;j<N;j++)
    {
    cout.precision(4);
    cout<<round(tab2[j].imag()*1000)/1000<<"  ";
    }
    cout<<endl;
    cout<<endl;system("pause");

    cout<<"frequency Hz DFT"<<endl;
    for(int j=0;j<N;j++)
    {
    cout.precision(4);
    cout<<round(tab3[j].real()*1000)/1000<<"  ";
    }
    cout<<endl;
    cout<<endl;system("pause");
    for(int j=0;j<N;j++)
    {
    cout.precision(4);
    cout<<round(tab3[j].imag()*1000)/1000<<"  ";
    }
    cout<<endl;
    cout<<endl;
    cout<<"if radix-4 == DFT tab2[j].real(): "<<endl;system("pause");
*/
       for(int j=0;j<N;j++)
        {
            if((j-1)%10000==0){system("pause")  ;}

            for(int i=0;i<N;i++)
        {
          //if((fabs(tab3[i].real()-tab2[j].real())<=0.03)&&fabs(tab3[0].real()-tab2[0].real())<0.03)
         if((fabs(tab3[i].real()-tab2[j].real())<=0.03)&&(fabs(tab3[i].imag()-tab2[j].imag())<=0.03))
          {
              if(i==30){   cout<<" ...   DFT["<<i<<"]==FFT["<<j<<"]= "<<round(tab2[j].real()*1000)/1000<<"  "<<endl;system("pause");}
            cout.precision(4);
            file3<<"DFT["<<i<<"]==FFT["<<j<<"]= "<<round(tab2[j].real()*10000)/10000<<"  "<<endl;;//system("pause");
           cout<<"DFT["<<i<<"]==FFT["<<j<<"]= "<<round(tab2[j].real()*10000)/10000<<"  "<<endl;;//system("pause");
          }
            else {
               // cout<<-1<<" . ";//system("pause");
            }
        }
        }
    system("pause");
    file3<<endl;
               for(int j=0;j<N;j++)
        {
            for(int i=0;i<N;i++)
        {
          if((fabs(tab3[i].real()-tab2[j].real())<=0.03)&&fabs(tab3[0].real()-tab2[0].real())<0.03)
          {
            cout.precision(4);
            file3<<" "<<i;//system("pause");
          }
            else {
               // cout<<-1<<" . ";//system("pause");
            }
        }
        }
        /////////////////////////////////////////
        tmp1=0;
        //cout<<"tmp1"<<tmp1<<endl;
        for(int i=1;i<N;i++)
        {
            tmp1.real(tmp1.real()+round(tab2[i].real()*100000000)/100000000);
            tmp1.imag(tmp1.imag()+round(tab2[i].imag()*100000000)/100000000);
        }
         cout<<"tab[0] "<<tab2[0]<<endl;
         cout<<"Sum "<<tmp1<<endl;
        /////////////////////////////////////////
    file3<<endl<<endl;

    }


    if(flag>mx){mx=flag;}
    flag=0;

if(a1==10000000){cout<<b1<<endl;a1=0;b1++;}
a1++;

  cout<<endl<<"max="<<mx<<endl;    file3.close();file1.close();file4.close();
    system("pause");
    cout<<endl;

    cout<<endl;
    delete [] tab2;
    delete [] tab3;
    delete [] tab4;
    system("pause");
    return 0;
}

void fun_inverse_bits_radix_4(int N,std::complex<double> tab[])
{
//code by Sidney Burrus
//http://dsp.stackexchange.com/questions/3481/radix-4-fft-implementation
    std::complex<double> t;
    //N=4^a;
    // Radix-4 bit-reverse
    double T;
    int j = 0;
    int N2 = N>>2;
    int M=0;
    for (int i=0; i < N-1; i++) {
        if (i < j) {
            t = tab[i];
            tab[i] = tab[j];
            tab[j] = t;
        }
        M = N2;
        while ( j >= 3*M ) {
            j -= 3*M;
            M >>= 2;
        }
        j += M;
    }
}
///////////////
void fun_fourier_transform_DFT(int N,std::complex<double> tab[])
{
    const double pi=3.141592653589793238462;
    std::complex<double> *tab2 = new std::complex<double>[N];    // tab2[]==N
    std::complex<double>  w[1]={{1,1}};


double zmienna1=2*pi/(float)N;
double fi2=fi;

for (int i=0;i<N;i++)
{
    for(int j=0;j<N;j++)
    {
          //complex number method:
          w[0].real(cos(i*j*zmienna1+fi2));
          w[0].imag(-sin(i*j*zmienna1+fi2));
          tab2[i]=tab2[i]+tab[j]*w[0];

    }
}

//new:
    for(int j=0;j<N;j++)
    {
      tab[j].real(tab2[j].real()*2/N);
     tab[j].imag(tab2[j].imag()*2/N);
    }

    delete tab2;
}
//////////////////

void fun_fourier_transform_FFT_mixed_radix_for_N_1848_points(int N,std::complex<double> tab[])
{
// author marcin matysek (r)ewertyn.PL 2marcin56@gmail.com
// open source
    const double pi=3.141592653589793238462;
    std::complex<double> *tab2 = new std::complex<double>[N];    // tab2[]==N

    double tmp2;
    double tmp3;
    double tmp6;
    double tmp7;
    double tmp8;
    double tmp9;
    double tmp10;
    double tmp11;

    double fi2=fi; //shift only for stage nr 1

    int i=0;
    int rx5=5,rx4=4,rx3=3,rx2=2,rx7=7,rx11=11;
    int stg[100]={};
    int nb_stages=0;
    int nb1,nb2,nb3,nb4,nb5_stg_previous,stg_first;

    tmp6=2*pi/(N/1);

    tmp2=2*pi/(2/1);
    tmp3=2*pi/(3/1);
    tmp10=2*pi/(4/1);
    tmp8=2*pi/(5/1);
    tmp7=2*pi/(7/1);
    tmp11=2*pi/(11/1);



    double z_rx2[2][2]={{}};
    double z_rx3[3][2]={{}};
    double z_rx4[2][2]={{}};
    double z_rx5[5][2]={{}};
    double z_rx7[7][2]={{}};
    double z_rx11[11][2]={{}};




//radix 2 fundament
          z_rx2[0][0]=cos(0*tmp2);
          z_rx2[0][1]=-sin(0*tmp2);
          z_rx2[1][0]=cos(1*tmp2);
          z_rx2[1][1]=-sin(1*tmp2);
//radix 3 fundament
          z_rx3[0][0]=cos(0*tmp3);
          z_rx3[0][1]=-sin(0*tmp3);
          z_rx3[1][0]=cos(1*tmp3);
          z_rx3[1][1]=-sin(1*tmp3);
          z_rx3[2][0]=cos(2*tmp3);
          z_rx3[2][1]=-sin(2*tmp3);
//radix 4 fundament
          z_rx4[0][0]=cos(0*tmp10);
          z_rx4[0][1]=-sin(0*tmp10);
          z_rx4[1][0]=cos(1*tmp10);
          z_rx4[1][1]=-sin(1*tmp10);
//radix 5 fundament
          z_rx5[0][0]=cos(0*tmp8);
          z_rx5[0][1]=-sin(0*tmp8);
          z_rx5[1][0]=cos(1*tmp8);
          z_rx5[1][1]=-sin(1*tmp8);
          z_rx5[2][0]=cos(2*tmp8);
          z_rx5[2][1]=-sin(2*tmp8);
          z_rx5[3][0]=cos(3*tmp8);
          z_rx5[3][1]=-sin(3*tmp8);
          z_rx5[4][0]=cos(4*tmp8);
          z_rx5[4][1]=-sin(4*tmp8);
//radix 7 fundament
          z_rx7[0][0]=cos(0*tmp7);
          z_rx7[0][1]=-sin(0*tmp7);
          z_rx7[1][0]=cos(1*tmp7);
          z_rx7[1][1]=-sin(1*tmp7);
          z_rx7[2][0]=cos(2*tmp7);
          z_rx7[2][1]=-sin(2*tmp7);
          z_rx7[3][0]=cos(3*tmp7);
          z_rx7[3][1]=-sin(3*tmp7);
          z_rx7[4][0]=cos(4*tmp7);
          z_rx7[4][1]=-sin(4*tmp7);
          z_rx7[5][0]=cos(5*tmp7);
          z_rx7[5][1]=-sin(5*tmp7);
          z_rx7[6][0]=cos(6*tmp7);
          z_rx7[6][1]=-sin(6*tmp7);
//radix 11 fundament
          z_rx11[0][0]= cos(0*tmp11);
          z_rx11[0][1]=-sin(0*tmp11);
          z_rx11[1][0]= cos(1*tmp11);
          z_rx11[1][1]=-sin(1*tmp11);
          z_rx11[2][0]= cos(2*tmp11);
          z_rx11[2][1]=-sin(2*tmp11);
          z_rx11[3][0]= cos(3*tmp11);
          z_rx11[3][1]=-sin(3*tmp11);
          z_rx11[4][0]= cos(4*tmp11);
          z_rx11[4][1]=-sin(4*tmp11);
          z_rx11[5][0]= cos(5*tmp11);
          z_rx11[5][1]=-sin(5*tmp11);
          z_rx11[6][0]= cos(6*tmp11);
          z_rx11[6][1]=-sin(6*tmp11);
          z_rx11[7][0]= cos(7*tmp11);
          z_rx11[7][1]=-sin(7*tmp11);
          z_rx11[8][0]= cos(8*tmp11);
          z_rx11[8][1]=-sin(8*tmp11);
          z_rx11[9][0]= cos(9*tmp11);
          z_rx11[9][1]=-sin(9*tmp11);
          z_rx11[10][0]= cos(10*tmp11);
          z_rx11[10][1]=-sin(10*tmp11);

/*
for(int j=0;j<102;j++)
{
    for(int i=0;i<102;i++)
    {
      if(((fabs(round(z_rx11[j].imag()*1000)/1000-round(z_rx11[i].imag()*1000)/1000)<0.001)
         &&(fabs(round(z_rx11[j].real()*1000)/1000-round(z_rx11[i].real()*1000)/1000)<0.001))
         ||((fabs(round(z_rx11[j].imag()*1000)/1000+round(z_rx11[i].imag()*1000)/1000)<0.001)
         &&(fabs(round(z_rx11[j].real()*1000)/1000+round(z_rx11[i].real()*1000)/1000)<0.001)))
         {
             cout<<j<<" "<<round(z_rx11[j].real()*1000)/1000<<" "<<round(z_rx11[j].imag()*1000)/1000<<"   ";
             cout<<i<<" "<<round(z_rx11[i].real()*1000)/1000<<" "<<round(z_rx11[i].imag()*1000)/1000<<endl;
         }
             //cout<<endl;
    }
    //system("pause");
}
*/

	nb_stages=radix_base(N,stg);

        if(nb_stages>=1){cout<<"N= "<<N<<endl;}
        for(int m=1;m<=nb_stages;m++)
        {
        cout <<"stage:"<<m<<" = radix-"<<stg[m] << endl;
        }
        cout << endl;



    if(nb_stages>=1)
    {
        stg_first=N/stg[1];
        if(stg[1]==2)
        {
            fun_fourier_transform_FFT_radix_2_stg_first(tab,stg_first,fi2,z_rx2);
        }
        else if(stg[1]==3)
        {
            fun_fourier_transform_FFT_radix_3_stg_first(tab,stg_first,fi2,z_rx3);
        }
        else if(stg[1]==4)
        {
            fun_fourier_transform_FFT_radix_4_stg_first(tab,stg_first,fi2,z_rx4);
        }
        else if(stg[1]==5)
        {
            fun_fourier_transform_FFT_radix_5_stg_first(tab,stg_first,fi2,z_rx5);
        }
        else if(stg[1]==7)
        {
            fun_fourier_transform_FFT_radix_7_stg_first(tab,stg_first,fi2,z_rx7);
        }
        else if(stg[1]==11)
        {
            fun_fourier_transform_FFT_radix_11_stg_first(tab,stg_first,fi2,z_rx11);
        }
        else{}
        nb1=N;
        nb4=1;
        for(int i=0;i<nb_stages-1;i++)
        {
            nb1=nb1/stg[0+i];
            nb2=nb1/stg[1+i];
            nb3=nb2/stg[2+i];
            nb4=nb4*stg[0+i];
            nb5_stg_previous=stg[1+i];

            if(stg[i+2]==2)
            {
                fun_fourier_transform_FFT_radix_2_stg_rest(tab,nb1,nb2,nb3,nb4,nb5_stg_previous,tmp6,z_rx2);
            }
            else if(stg[i+2]==3)
            {
                fun_fourier_transform_FFT_radix_3_stg_rest(tab,nb1,nb2,nb3,nb4,nb5_stg_previous,tmp6,z_rx3);
            }
            else if(stg[i+2]==4)
            {
                fun_fourier_transform_FFT_radix_4_stg_rest(tab,nb1,nb2,nb3,nb4,nb5_stg_previous,tmp6,z_rx4);
            }
            else if(stg[i+2]==5)
            {
                fun_fourier_transform_FFT_radix_5_stg_rest(tab,nb1,nb2,nb3,nb4,nb5_stg_previous,tmp6,z_rx5);
            }
            else if(stg[i+2]==7)
            {
                fun_fourier_transform_FFT_radix_7_stg_rest(tab,nb1,nb2,nb3,nb4,nb5_stg_previous,tmp6,z_rx7);
            }
            else if(stg[i+2]==11)
            {
                fun_fourier_transform_FFT_radix_11_stg_rest(tab,nb1,nb2,nb3,nb4,nb5_stg_previous,tmp6,z_rx11);
            }
            else{}

        }
    }

//new:
    for(int j=0;j<N;j++)
    {
     tab[j].real(tab[j].real()*2/N);
     tab[j].imag(tab[j].imag()*2/N);
    }
    delete [] tab2;
}
///////////////////////////////////////////////////


  void fun_fourier_transform_FFT_radix_2_stg_rest(std::complex<double> tab[],int nb1,int nb2,int nb3,int nb4,int nb5_stg_previous,double tmp6,double (*z)[2])
    {
        double tmp1[2],tmp2[2],tmp20[2];
        double w[2][2]={{}};
        double tmp100=0.0;
        double tmp200=0.0;
        double tmp300=0.0;
        int nb_tmp1=0.0;
        int nb_tmp2=0.0;
        int nb_tmp3=0.0;

        for(int b=0;b<nb5_stg_previous;b=b+1)
        {
            tmp300=nb4*b*tmp6;
            tmp100=nb3*tmp300;
            for(int j=0;j<nb3;j=j+1)
            {
                tmp200=j*tmp300;
                w[0][0]=cos(0*tmp100+tmp200);
                w[0][1]=-sin(0*tmp100+tmp200);
                w[1][0]=cos(1*tmp100+tmp200);
                w[1][1]=-sin(1*tmp100+tmp200);
                //w[1].imag(-sin(nb4*b*(1*nb3+j)*tmp5));

                nb_tmp1=b*nb2+j;
                for(int i=0;i<nb4;i=i+1)
                {
                    nb_tmp2=i*nb1;
                    nb_tmp3=nb_tmp1+nb_tmp2;

                  tmp1[0]=w[0][0]*tab[nb_tmp3+0*nb3].real()-w[0][1]*tab[nb_tmp3+0*nb3].imag();
                  tmp1[1]=w[0][0]*tab[nb_tmp3+0*nb3].imag()+w[0][1]*tab[nb_tmp3+0*nb3].real();
                  tmp2[0]=w[1][0]*tab[nb_tmp3+1*nb3].real()-w[1][1]*tab[nb_tmp3+1*nb3].imag();
                  tmp2[1]=w[1][0]*tab[nb_tmp3+1*nb3].imag()+w[1][1]*tab[nb_tmp3+1*nb3].real();

                  tmp20[0]=z[0][0]*(tmp1[0]+tmp2[0])-z[0][1]*(tmp1[1]+tmp2[1]);
                  tmp20[1]=z[0][0]*(tmp1[1]+tmp2[1])+z[0][1]*(tmp1[0]+tmp2[0]);
                    //radix-2
                  tab[nb_tmp3+0*nb3].real(tmp20[0]);
                  tab[nb_tmp3+0*nb3].imag(tmp20[1]);

                  tab[nb_tmp3+1*nb3].real(z[0][0]*tmp1[0]-z[0][1]*tmp1[1]+z[1][0]*tmp2[0]-z[1][1]*tmp2[1]);
                  tab[nb_tmp3+1*nb3].imag(z[0][0]*tmp1[1]+z[0][1]*tmp1[0]+z[1][0]*tmp2[1]+z[1][1]*tmp2[0]);
                }
            }
        }
    }
///////////////////////////////////////////////////////////////


    void fun_fourier_transform_FFT_radix_3_stg_rest(std::complex<double> tab[],int nb1,int nb2,int nb3,int nb4,int nb5_stg_previous,double tmp6,double (*z)[2])
    {
        double tmp1[2],tmp2[2],tmp3[2],tmp20[2],tmp30[2];
        double w[3][2]={{}};
        double tmp100=0.0;
        double tmp200=0.0;
        double tmp300=0.0;
        int nb_tmp1=0.0;
        int nb_tmp2=0.0;
        int nb_tmp3=0.0;

        for(int b=0;b<nb5_stg_previous;b=b+1)
        {
            tmp300=nb4*b*tmp6;
            tmp100=nb3*tmp300;
            for(int j=0;j<nb3;j=j+1)
            {
                tmp200=j*tmp300;
                w[0][0]=cos(0*tmp100+tmp200);
                w[0][1]=-sin(0*tmp100+tmp200);
                w[1][0]=cos(1*tmp100+tmp200);
                w[1][1]=-sin(1*tmp100+tmp200);
                w[2][0]=cos(2*tmp100+tmp200);
                w[2][1]=-sin(2*tmp100+tmp200);
                //w[2].real(cos(nb4*b*(2*nb3+j)*tmp99));

                nb_tmp1=b*nb2+j;
                for(int i=0;i<nb4;i=i+1)
                {
                    nb_tmp2=i*nb1;
                    nb_tmp3=nb_tmp1+nb_tmp2;

                  tmp1[0]=w[0][0]*tab[nb_tmp3+0*nb3].real()-w[0][1]*tab[nb_tmp3+0*nb3].imag();
                  tmp1[1]=w[0][0]*tab[nb_tmp3+0*nb3].imag()+w[0][1]*tab[nb_tmp3+0*nb3].real();
                  tmp2[0]=w[1][0]*tab[nb_tmp3+1*nb3].real()-w[1][1]*tab[nb_tmp3+1*nb3].imag();
                  tmp2[1]=w[1][0]*tab[nb_tmp3+1*nb3].imag()+w[1][1]*tab[nb_tmp3+1*nb3].real();
                  tmp3[0]=w[2][0]*tab[nb_tmp3+2*nb3].real()-w[2][1]*tab[nb_tmp3+2*nb3].imag();
                  tmp3[1]=w[2][0]*tab[nb_tmp3+2*nb3].imag()+w[2][1]*tab[nb_tmp3+2*nb3].real();

                  tmp30[0]=z[0][0]*tmp1[0]-z[0][1]*tmp1[1];
                  tmp30[1]=z[0][0]*tmp1[1]+z[0][1]*tmp1[0];
                  tmp20[0]=z[0][0]*(tmp1[0]+tmp2[0]+tmp3[0])-z[0][1]*(tmp1[1]+tmp2[1]+tmp3[1]);
                  tmp20[1]=z[0][0]*(tmp1[1]+tmp2[1]+tmp3[1])+z[0][1]*(tmp1[0]+tmp2[0]+tmp3[0]);


                 //radix-3
                  tab[nb_tmp3+0*nb3].real(tmp20[0]);
                  tab[nb_tmp3+0*nb3].imag(tmp20[1]);

                  tab[nb_tmp3+1*nb3].real(tmp30[0]+z[1][0]*tmp2[0]-z[1][1]*tmp2[1]+z[2][0]*tmp3[0]-z[2][1]*tmp3[1]);
                  tab[nb_tmp3+1*nb3].imag(tmp30[1]+z[1][0]*tmp2[1]+z[1][1]*tmp2[0]+z[2][0]*tmp3[1]+z[2][1]*tmp3[0]);

                  tab[nb_tmp3+2*nb3].real(tmp30[0]+z[2][0]*tmp2[0]-z[2][1]*tmp2[1]+z[1][0]*tmp3[0]-z[1][1]*tmp3[1]);
                  tab[nb_tmp3+2*nb3].imag(tmp30[1]+z[2][0]*tmp2[1]+z[2][1]*tmp2[0]+z[1][0]*tmp3[1]+z[1][1]*tmp3[0]);
                }
            }
        }
    }
///////////////////////////////////////////////////////////////


    void fun_fourier_transform_FFT_radix_4_stg_rest(std::complex<double> tab[],int nb1,int nb2,int nb3,int nb4,int nb5_stg_previous,double tmp6,double (*z)[2])
    {
        double tmp1[2],tmp2[2],tmp3[2],tmp4[2];
        double tmp101[2],tmp102[2],tmp103[2],tmp104[2];
        double tmp111[2],tmp112[2],tmp113[2],tmp114[2];
        double  w[4][2]={{1,0}};
        double tmp100=0.0;
        double tmp200=0.0;
        double tmp300=0.0;
        int nb_tmp1=0.0;
        int nb_tmp2=0.0;
        int nb_tmp3=0.0;

        for(int b=0;b<nb5_stg_previous;b=b+1)
        {
            tmp300=nb4*b*tmp6;
            tmp100=nb3*tmp300;
            for(int j=0;j<nb3;j=j+1)
            {
              tmp200=j*tmp300;
              w[0][0]=cos(0*tmp100+tmp200);
              w[0][1]=-sin(0*tmp100+tmp200);
              w[1][0]=cos(1*tmp100+tmp200);
              w[1][1]=-sin(1*tmp100+tmp200);
              w[2][0]=cos(2*tmp100+tmp200);
              w[2][1]=-sin(2*tmp100+tmp200);
              //w[3].real(cos(nb4*b*(3*nb3+j)*tmp5));
              w[3][0]=cos(3*tmp100+tmp200);
              w[3][1]=-sin(3*tmp100+tmp200);

                nb_tmp1=b*nb2+j;
                for(int i=0;i<nb4;i=i+1)
                {
                    nb_tmp2=i*nb1;
                    nb_tmp3=nb_tmp1+nb_tmp2;
                    tmp1[0]=w[0][0]*tab[nb_tmp3+0*nb3].real()-w[0][1]*tab[nb_tmp3+0*nb3].imag();
                    tmp1[1]=w[0][0]*tab[nb_tmp3+0*nb3].imag()+w[0][1]*tab[nb_tmp3+0*nb3].real();
                    tmp2[0]=w[1][0]*tab[nb_tmp3+1*nb3].real()-w[1][1]*tab[nb_tmp3+1*nb3].imag();
                    tmp2[1]=w[1][0]*tab[nb_tmp3+1*nb3].imag()+w[1][1]*tab[nb_tmp3+1*nb3].real();
                    tmp3[0]=w[2][0]*tab[nb_tmp3+2*nb3].real()-w[2][1]*tab[nb_tmp3+2*nb3].imag();
                    tmp3[1]=w[2][0]*tab[nb_tmp3+2*nb3].imag()+w[2][1]*tab[nb_tmp3+2*nb3].real();
                    tmp4[0]=w[3][0]*tab[nb_tmp3+3*nb3].real()-w[3][1]*tab[nb_tmp3+3*nb3].imag();
                    tmp4[1]=w[3][0]*tab[nb_tmp3+3*nb3].imag()+w[3][1]*tab[nb_tmp3+3*nb3].real();


                    tmp101[0]=tmp2[0]-tmp4[0];
                    tmp101[1]=tmp2[1]-tmp4[1];
                    tmp102[0]=tmp2[0]+tmp4[0];
                    tmp102[1]=tmp2[1]+tmp4[1];
                    tmp103[0]=tmp1[0]-tmp3[0];
                    tmp103[1]=tmp1[1]-tmp3[1];
                    tmp104[0]=tmp1[0]+tmp3[0];
                    tmp104[1]=tmp1[1]+tmp3[1];

                    tmp111[0]=z[0][0]*(tmp104[0]+tmp102[0])-z[0][1]*(tmp104[1]+tmp102[1]);
                    tmp111[1]=z[0][0]*(tmp104[1]+tmp102[1])+z[0][1]*(tmp104[0]+tmp102[0]);
                    tmp112[0]=z[0][0]*tmp103[0]-z[0][1]*tmp103[1];
                    tmp112[1]=z[0][0]*tmp103[1]+z[0][1]*tmp103[0];
                    tmp113[0]=z[0][0]*(tmp104[0]-tmp102[0])-z[0][1]*(tmp104[1]-tmp102[1]);
                    tmp113[1]=z[0][0]*(tmp104[1]-tmp102[1])+z[0][1]*(tmp104[0]-tmp102[0]);
                    tmp114[0]=z[0][0]*tmp103[0]-z[0][1]*tmp103[1];
                    tmp114[1]=z[0][0]*tmp103[1]+z[0][1]*tmp103[0];

                    //radix-4
                    tab[nb_tmp3+0*nb3].real(tmp111[0]);
                    tab[nb_tmp3+0*nb3].imag(tmp111[1]);
                    tab[nb_tmp3+1*nb3].real(tmp112[0]+z[1][0]*tmp101[0]-z[1][1]*tmp101[1]);
                    tab[nb_tmp3+1*nb3].imag(tmp112[1]+z[1][0]*tmp101[1]+z[1][1]*tmp101[0]);
                    tab[nb_tmp3+2*nb3].real(tmp113[0]);
                    tab[nb_tmp3+2*nb3].imag(tmp113[1]);
                    tab[nb_tmp3+3*nb3].real(tmp114[0]-(z[1][0]*tmp101[0]-z[1][1]*tmp101[1]));
                    tab[nb_tmp3+3*nb3].imag(tmp114[1]-(z[1][0]*tmp101[1]+z[1][1]*tmp101[0]));
                }
            }
        }
    }
/////////////////////////////////////////////////////////////////

   void fun_fourier_transform_FFT_radix_5_stg_rest(std::complex<double> tab[],int nb1,int nb2,int nb3,int nb4,int nb5_stg_previous,double tmp6,double  (*z)[2])
    {
        double tmp1[2],tmp2[2],tmp3[2],tmp4[2],tmp5[2],tmp20[2],tmp30[2];
        double  w[5][2]={{}};
        double tmp100=0.0;
        double tmp200=0.0;
        double tmp300=0.0;
        int nb_tmp1=0.0;
        int nb_tmp2=0.0;
        int nb_tmp3=0.0;

        for(int b=0;b<nb5_stg_previous;b=b+1)
        {
            tmp300=nb4*b*tmp6;
            tmp100=nb3*tmp300;
            for(int j=0;j<nb3;j=j+1)
            {
                tmp200=j*tmp300;
                w[0][0]=cos(0*tmp100+tmp200);
                w[0][1]=-sin(0*tmp100+tmp200);
                w[1][0]=cos(1*tmp100+tmp200);
                w[1][1]=-sin(1*tmp100+tmp200);
                w[2][0]=cos(2*tmp100+tmp200);
                w[2][1]=-sin(2*tmp100+tmp200);
                //w[2].real(cos(nb4*b*(2*nb3+j)*tmp99));
                w[3][0]=cos(3*tmp100+tmp200);
                w[3][1]=-sin(3*tmp100+tmp200);
                w[4][0]=cos(4*tmp100+tmp200);
                w[4][1]=-sin(4*tmp100+tmp200);

                nb_tmp1=b*nb2+j;
                for(int i=0;i<nb4;i=i+1)
                {
                    nb_tmp2=i*nb1;
                    nb_tmp3=nb_tmp1+nb_tmp2;

                  tmp1[0]=w[0][0]*tab[nb_tmp3+0*nb3].real()-w[0][1]*tab[nb_tmp3+0*nb3].imag();
                  tmp1[1]=w[0][0]*tab[nb_tmp3+0*nb3].imag()+w[0][1]*tab[nb_tmp3+0*nb3].real();
                  tmp2[0]=w[1][0]*tab[nb_tmp3+1*nb3].real()-w[1][1]*tab[nb_tmp3+1*nb3].imag();
                  tmp2[1]=w[1][0]*tab[nb_tmp3+1*nb3].imag()+w[1][1]*tab[nb_tmp3+1*nb3].real();
                  tmp3[0]=w[2][0]*tab[nb_tmp3+2*nb3].real()-w[2][1]*tab[nb_tmp3+2*nb3].imag();
                  tmp3[1]=w[2][0]*tab[nb_tmp3+2*nb3].imag()+w[2][1]*tab[nb_tmp3+2*nb3].real();
                  tmp4[0]=w[3][0]*tab[nb_tmp3+3*nb3].real()-w[3][1]*tab[nb_tmp3+3*nb3].imag();
                  tmp4[1]=w[3][0]*tab[nb_tmp3+3*nb3].imag()+w[3][1]*tab[nb_tmp3+3*nb3].real();
                  tmp5[0]=w[4][0]*tab[nb_tmp3+4*nb3].real()-w[4][1]*tab[nb_tmp3+4*nb3].imag();
                  tmp5[1]=w[4][0]*tab[nb_tmp3+4*nb3].imag()+w[4][1]*tab[nb_tmp3+4*nb3].real();

                  tmp30[0]=z[0][0]*tmp1[0]-z[0][1]*tmp1[1];
                  tmp30[1]=z[0][0]*tmp1[1]+z[0][1]*tmp1[0];
                  tmp20[0]=z[0][0]*(tmp1[0]+tmp2[0]+tmp3[0]+tmp4[0]+tmp5[0])-z[0][1]*(tmp1[1]+tmp2[1]+tmp3[1]+tmp4[1]+tmp5[1]);
                  tmp20[1]=z[0][0]*(tmp1[1]+tmp2[1]+tmp3[1]+tmp4[1]+tmp5[1])+z[0][1]*(tmp1[0]+tmp2[0]+tmp3[0]+tmp4[0]+tmp5[0]);


                 //radix-5
                  tab[nb_tmp3+0*nb3].real(tmp20[0]);
                  tab[nb_tmp3+0*nb3].imag(tmp20[1]);

                  tab[nb_tmp3+1*nb3].real(tmp30[0]+z[1][0]*tmp2[0]-z[1][1]*tmp2[1]+z[2][0]*tmp3[0]-z[2][1]*tmp3[1]+z[3][0]*tmp4[0]-z[3][1]*tmp4[1]+z[4][0]*tmp5[0]-z[4][1]*tmp5[1]);
                  tab[nb_tmp3+1*nb3].imag(tmp30[1]+z[1][0]*tmp2[1]+z[1][1]*tmp2[0]+z[2][0]*tmp3[1]+z[2][1]*tmp3[0]+z[3][0]*tmp4[1]+z[3][1]*tmp4[0]+z[4][0]*tmp5[1]+z[4][1]*tmp5[0]);

                  tab[nb_tmp3+2*nb3].real(tmp30[0]+z[2][0]*tmp2[0]-z[2][1]*tmp2[1]+z[4][0]*tmp3[0]-z[4][1]*tmp3[1]+z[1][0]*tmp4[0]-z[1][1]*tmp4[1]+z[3][0]*tmp5[0]-z[3][1]*tmp5[1]);
                  tab[nb_tmp3+2*nb3].imag(tmp30[1]+z[2][0]*tmp2[1]+z[2][1]*tmp2[0]+z[4][0]*tmp3[1]+z[4][1]*tmp3[0]+z[1][0]*tmp4[1]+z[1][1]*tmp4[0]+z[3][0]*tmp5[1]+z[3][1]*tmp5[0]);

                  tab[nb_tmp3+3*nb3].real(tmp30[0]+z[3][0]*tmp2[0]-z[3][1]*tmp2[1]+z[1][0]*tmp3[0]-z[1][1]*tmp3[1]+z[4][0]*tmp4[0]-z[4][1]*tmp4[1]+z[2][0]*tmp5[0]-z[2][1]*tmp5[1]);
                  tab[nb_tmp3+3*nb3].imag(tmp30[1]+z[3][0]*tmp2[1]+z[3][1]*tmp2[0]+z[1][0]*tmp3[1]+z[1][1]*tmp3[0]+z[4][0]*tmp4[1]+z[4][1]*tmp4[0]+z[2][0]*tmp5[1]+z[2][1]*tmp5[0]);

                  tab[nb_tmp3+4*nb3].real(tmp30[0]+z[4][0]*tmp2[0]-z[4][1]*tmp2[1]+z[3][0]*tmp3[0]-z[3][1]*tmp3[1]+z[2][0]*tmp4[0]-z[2][1]*tmp4[1]+z[1][0]*tmp5[0]-z[1][1]*tmp5[1]);
                  tab[nb_tmp3+4*nb3].imag(tmp30[1]+z[4][0]*tmp2[1]+z[4][1]*tmp2[0]+z[3][0]*tmp3[1]+z[3][1]*tmp3[0]+z[2][0]*tmp4[1]+z[2][1]*tmp4[0]+z[1][0]*tmp5[1]+z[1][1]*tmp5[0]);
                }
            }
        }
    }
///////////////////////////////////////////////////////////////////////

  void fun_fourier_transform_FFT_radix_7_stg_rest(std::complex<double> tab[],int nb1,int nb2,int nb3,int nb4,int nb5_stg_previous,double tmp99,double  (*z)[2])
    {
        double tmp1[2],tmp2[2],tmp3[2],tmp4[2],tmp5[2],tmp6[2],tmp7[2],tmp20[2],tmp30[2];
        double  w[7][2]={{1,0}};
        double tmp100=0.0;
        double tmp200=0.0;
        double tmp300=0.0;
        int nb_tmp1=0.0;
        int nb_tmp2=0.0;
        int nb_tmp3=0.0;

        for(int b=0;b<nb5_stg_previous;b=b+1)
        {
            tmp300=nb4*b*tmp99;
            tmp100=nb3*tmp300;
            for(int j=0;j<nb3;j=j+1)
            {
                tmp200=j*tmp300;
                w[0][0]=cos(0*tmp100+tmp200);
                w[0][1]=-sin(0*tmp100+tmp200);
                w[1][0]=cos(1*tmp100+tmp200);
                w[1][1]=-sin(1*tmp100+tmp200);
                w[2][0]=cos(2*tmp100+tmp200);
                w[2][1]=-sin(2*tmp100+tmp200);
                //w[2].real(cos(nb4*b*(2*nb3+j)*tmp99));
                w[3][0]=cos(3*tmp100+tmp200);
                w[3][1]=-sin(3*tmp100+tmp200);
                w[4][0]=cos(4*tmp100+tmp200);
                w[4][1]=-sin(4*tmp100+tmp200);
                w[5][0]=cos(5*tmp100+tmp200);
                w[5][1]=-sin(5*tmp100+tmp200);
                w[6][0]=cos(6*tmp100+tmp200);
                w[6][1]=-sin(6*tmp100+tmp200);

                nb_tmp1=b*nb2+j;
                for(int i=0;i<nb4;i=i+1)
                {
                    nb_tmp2=i*nb1;
                    nb_tmp3=nb_tmp1+nb_tmp2;

                  tmp1[0]=w[0][0]*tab[nb_tmp3+0*nb3].real()-w[0][1]*tab[nb_tmp3+0*nb3].imag();
                  tmp1[1]=w[0][0]*tab[nb_tmp3+0*nb3].imag()+w[0][1]*tab[nb_tmp3+0*nb3].real();
                  tmp2[0]=w[1][0]*tab[nb_tmp3+1*nb3].real()-w[1][1]*tab[nb_tmp3+1*nb3].imag();
                  tmp2[1]=w[1][0]*tab[nb_tmp3+1*nb3].imag()+w[1][1]*tab[nb_tmp3+1*nb3].real();
                  tmp3[0]=w[2][0]*tab[nb_tmp3+2*nb3].real()-w[2][1]*tab[nb_tmp3+2*nb3].imag();
                  tmp3[1]=w[2][0]*tab[nb_tmp3+2*nb3].imag()+w[2][1]*tab[nb_tmp3+2*nb3].real();
                  tmp4[0]=w[3][0]*tab[nb_tmp3+3*nb3].real()-w[3][1]*tab[nb_tmp3+3*nb3].imag();
                  tmp4[1]=w[3][0]*tab[nb_tmp3+3*nb3].imag()+w[3][1]*tab[nb_tmp3+3*nb3].real();
                  tmp5[0]=w[4][0]*tab[nb_tmp3+4*nb3].real()-w[4][1]*tab[nb_tmp3+4*nb3].imag();
                  tmp5[1]=w[4][0]*tab[nb_tmp3+4*nb3].imag()+w[4][1]*tab[nb_tmp3+4*nb3].real();
                  tmp6[0]=w[5][0]*tab[nb_tmp3+5*nb3].real()-w[5][1]*tab[nb_tmp3+5*nb3].imag();
                  tmp6[1]=w[5][0]*tab[nb_tmp3+5*nb3].imag()+w[5][1]*tab[nb_tmp3+5*nb3].real();
                  tmp7[0]=w[6][0]*tab[nb_tmp3+6*nb3].real()-w[6][1]*tab[nb_tmp3+6*nb3].imag();
                  tmp7[1]=w[6][0]*tab[nb_tmp3+6*nb3].imag()+w[6][1]*tab[nb_tmp3+6*nb3].real();

                  tmp30[0]=z[0][0]*tmp1[0]-z[0][1]*tmp1[1];
                  tmp30[1]=z[0][0]*tmp1[1]+z[0][1]*tmp1[0];
                  tmp20[0]=z[0][0]*(tmp1[0]+tmp2[0]+tmp3[0]+tmp4[0]+tmp5[0]+tmp6[0]+tmp7[0])-z[0][1]*(tmp1[1]+tmp2[1]+tmp3[1]+tmp4[1]+tmp5[1]+tmp6[1]+tmp7[1]);
                  tmp20[1]=z[0][0]*(tmp1[1]+tmp2[1]+tmp3[1]+tmp4[1]+tmp5[1]+tmp6[1]+tmp7[1])+z[0][1]*(tmp1[0]+tmp2[0]+tmp3[0]+tmp4[0]+tmp5[0]+tmp6[0]+tmp7[0]);

                 //radix-7
                  tab[nb_tmp3+0*nb3].real(tmp20[0]);
                  tab[nb_tmp3+0*nb3].imag(tmp20[1]);

                  tab[nb_tmp3+1*nb3].real(tmp30[0]+z[1][0]*tmp2[0]-z[1][1]*tmp2[1]+z[2][0]*tmp3[0]-z[2][1]*tmp3[1]+z[3][0]*tmp4[0]-z[3][1]*tmp4[1]+
                                          z[4][0]*tmp5[0]-z[4][1]*tmp5[1]+z[5][0]*tmp6[0]-z[5][1]*tmp6[1]+z[6][0]*tmp7[0]-z[6][1]*tmp7[1]);
                  tab[nb_tmp3+1*nb3].imag(tmp30[1]+z[1][0]*tmp2[1]+z[1][1]*tmp2[0]+z[2][0]*tmp3[1]+z[2][1]*tmp3[0]+z[3][0]*tmp4[1]+z[3][1]*tmp4[0]+
                                          z[4][0]*tmp5[1]+z[4][1]*tmp5[0]+z[5][0]*tmp6[1]+z[5][1]*tmp6[0]+z[6][0]*tmp7[1]+z[6][1]*tmp7[0]);

                  tab[nb_tmp3+2*nb3].real(tmp30[0]+z[2][0]*tmp2[0]-z[2][1]*tmp2[1]+z[4][0]*tmp3[0]-z[4][1]*tmp3[1]+z[6][0]*tmp4[0]-z[6][1]*tmp4[1]+
                                          z[1][0]*tmp5[0]-z[1][1]*tmp5[1]+z[3][0]*tmp6[0]-z[3][1]*tmp6[1]+z[5][0]*tmp7[0]-z[5][1]*tmp7[1]);
                  tab[nb_tmp3+2*nb3].imag(tmp30[1]+z[2][0]*tmp2[1]+z[2][1]*tmp2[0]+z[4][0]*tmp3[1]+z[4][1]*tmp3[0]+z[6][0]*tmp4[1]+z[6][1]*tmp4[0]+
                                          z[1][0]*tmp5[1]+z[1][1]*tmp5[0]+z[3][0]*tmp6[1]+z[3][1]*tmp6[0]+z[5][0]*tmp7[1]+z[5][1]*tmp7[0]);

                  tab[nb_tmp3+3*nb3].real(tmp30[0]+z[3][0]*tmp2[0]-z[3][1]*tmp2[1]+z[6][0]*tmp3[0]-z[6][1]*tmp3[1]+z[2][0]*tmp4[0]-z[2][1]*tmp4[1]+
                                          z[5][0]*tmp5[0]-z[5][1]*tmp5[1]+z[1][0]*tmp6[0]-z[1][1]*tmp6[1]+z[4][0]*tmp7[0]-z[4][1]*tmp7[1]);
                  tab[nb_tmp3+3*nb3].imag(tmp30[1]+z[3][0]*tmp2[1]+z[3][1]*tmp2[0]+z[6][0]*tmp3[1]+z[6][1]*tmp3[0]+z[2][0]*tmp4[1]+z[2][1]*tmp4[0]+
                                          z[5][0]*tmp5[1]+z[5][1]*tmp5[0]+z[1][0]*tmp6[1]+z[1][1]*tmp6[0]+z[4][0]*tmp7[1]+z[4][1]*tmp7[0]);

                  tab[nb_tmp3+4*nb3].real(tmp30[0]+z[4][0]*tmp2[0]-z[4][1]*tmp2[1]+z[1][0]*tmp3[0]-z[1][1]*tmp3[1]+z[5][0]*tmp4[0]-z[5][1]*tmp4[1]+
                                          z[2][0]*tmp5[0]-z[2][1]*tmp5[1]+z[6][0]*tmp6[0]-z[6][1]*tmp6[1]+z[3][0]*tmp7[0]-z[3][1]*tmp7[1]);
                  tab[nb_tmp3+4*nb3].imag(tmp30[1]+z[4][0]*tmp2[1]+z[4][1]*tmp2[0]+z[1][0]*tmp3[1]+z[1][1]*tmp3[0]+z[5][0]*tmp4[1]+z[5][1]*tmp4[0]+
                                          z[2][0]*tmp5[1]+z[2][1]*tmp5[0]+z[6][0]*tmp6[1]+z[6][1]*tmp6[0]+z[3][0]*tmp7[1]+z[3][1]*tmp7[0]);

                  tab[nb_tmp3+5*nb3].real(tmp30[0]+z[5][0]*tmp2[0]-z[5][1]*tmp2[1]+z[3][0]*tmp3[0]-z[3][1]*tmp3[1]+z[1][0]*tmp4[0]-z[1][1]*tmp4[1]+
                                          z[6][0]*tmp5[0]-z[6][1]*tmp5[1]+z[4][0]*tmp6[0]-z[4][1]*tmp6[1]+z[2][0]*tmp7[0]-z[2][1]*tmp7[1]);
                  tab[nb_tmp3+5*nb3].imag(tmp30[1]+z[5][0]*tmp2[1]+z[5][1]*tmp2[0]+z[3][0]*tmp3[1]+z[3][1]*tmp3[0]+z[1][0]*tmp4[1]+z[1][1]*tmp4[0]+
                                          z[6][0]*tmp5[1]+z[6][1]*tmp5[0]+z[4][0]*tmp6[1]+z[4][1]*tmp6[0]+z[2][0]*tmp7[1]+z[2][1]*tmp7[0]);

                  tab[nb_tmp3+6*nb3].real(tmp30[0]+z[6][0]*tmp2[0]-z[6][1]*tmp2[1]+z[5][0]*tmp3[0]-z[5][1]*tmp3[1]+z[4][0]*tmp4[0]-z[4][1]*tmp4[1]+
                                          z[3][0]*tmp5[0]-z[3][1]*tmp5[1]+z[2][0]*tmp6[0]-z[2][1]*tmp6[1]+z[1][0]*tmp7[0]-z[1][1]*tmp7[1]);
                  tab[nb_tmp3+6*nb3].imag(tmp30[1]+z[6][0]*tmp2[1]+z[6][1]*tmp2[0]+z[5][0]*tmp3[1]+z[5][1]*tmp3[0]+z[4][0]*tmp4[1]+z[4][1]*tmp4[0]+
                                          z[3][0]*tmp5[1]+z[3][1]*tmp5[0]+z[2][0]*tmp6[1]+z[2][1]*tmp6[0]+z[1][0]*tmp7[1]+z[1][1]*tmp7[0]);
                }
            }
        }
    }
///////////////////////////////////////////////////////////////////////

   void fun_fourier_transform_FFT_radix_11_stg_rest(std::complex<double> tab[],int nb1,int nb2,int nb3,int nb4,int nb5_stg_previous,double tmp99,double  (*z)[2])
    {
        double tmp1[2],tmp2[2],tmp3[2],tmp4[2],tmp5[2],tmp6[2],tmp7[2],tmp8[2],tmp9[2],tmp10[2],tmp11[2],tmp30[2],tmp20[2];
        double  w[11][2]={{1,0}};
        double tmp100=0.0;
        double tmp200=0.0;
        double tmp300=0.0;
        int nb_tmp1=0.0;
        int nb_tmp2=0.0;
        int nb_tmp3=0.0;

        for(int b=0;b<nb5_stg_previous;b=b+1)
        {
            tmp300=nb4*b*tmp99;
            tmp100=nb3*tmp300;
            for(int j=0;j<nb3;j=j+1)
            {
                tmp200=j*tmp300;
                w[0][0]=cos(0*tmp100+tmp200);
                w[0][1]=-sin(0*tmp100+tmp200);
                w[1][0]=cos(1*tmp100+tmp200);
                w[1][1]=-sin(1*tmp100+tmp200);
                w[2][0]=cos(2*tmp100+tmp200);
                w[2][1]=-sin(2*tmp100+tmp200);
                //w[2].real(cos(nb4*b*(2*nb3+j)*tmp99));
                w[3][0]=cos(3*tmp100+tmp200);
                w[3][1]=-sin(3*tmp100+tmp200);
                w[4][0]=cos(4*tmp100+tmp200);
                w[4][1]=-sin(4*tmp100+tmp200);
                w[5][0]=cos(5*tmp100+tmp200);
                w[5][1]=-sin(5*tmp100+tmp200);
                w[6][0]=cos(6*tmp100+tmp200);
                w[6][1]=-sin(6*tmp100+tmp200);
                w[7][0]=cos(7*tmp100+tmp200);
                w[7][1]=-sin(7*tmp100+tmp200);
                w[8][0]=cos(8*tmp100+tmp200);
                w[8][1]=-sin(8*tmp100+tmp200);
                w[9][0]=cos(9*tmp100+tmp200);
                w[9][1]=-sin(9*tmp100+tmp200);
                w[10][0]=cos(10*tmp100+tmp200);
                w[10][1]=-sin(10*tmp100+tmp200);

                nb_tmp1=b*nb2+j;
                for(int i=0;i<nb4;i=i+1)
                {
                    nb_tmp2=i*nb1;
                    nb_tmp3=nb_tmp1+nb_tmp2;

                  tmp1[0]=w[0][0]*tab[nb_tmp3+0*nb3].real()-w[0][1]*tab[nb_tmp3+0*nb3].imag();
                  tmp1[1]=w[0][0]*tab[nb_tmp3+0*nb3].imag()+w[0][1]*tab[nb_tmp3+0*nb3].real();
                  tmp2[0]=w[1][0]*tab[nb_tmp3+1*nb3].real()-w[1][1]*tab[nb_tmp3+1*nb3].imag();
                  tmp2[1]=w[1][0]*tab[nb_tmp3+1*nb3].imag()+w[1][1]*tab[nb_tmp3+1*nb3].real();
                  tmp3[0]=w[2][0]*tab[nb_tmp3+2*nb3].real()-w[2][1]*tab[nb_tmp3+2*nb3].imag();
                  tmp3[1]=w[2][0]*tab[nb_tmp3+2*nb3].imag()+w[2][1]*tab[nb_tmp3+2*nb3].real();
                  tmp4[0]=w[3][0]*tab[nb_tmp3+3*nb3].real()-w[3][1]*tab[nb_tmp3+3*nb3].imag();
                  tmp4[1]=w[3][0]*tab[nb_tmp3+3*nb3].imag()+w[3][1]*tab[nb_tmp3+3*nb3].real();
                  tmp5[0]=w[4][0]*tab[nb_tmp3+4*nb3].real()-w[4][1]*tab[nb_tmp3+4*nb3].imag();
                  tmp5[1]=w[4][0]*tab[nb_tmp3+4*nb3].imag()+w[4][1]*tab[nb_tmp3+4*nb3].real();
                  tmp6[0]=w[5][0]*tab[nb_tmp3+5*nb3].real()-w[5][1]*tab[nb_tmp3+5*nb3].imag();
                  tmp6[1]=w[5][0]*tab[nb_tmp3+5*nb3].imag()+w[5][1]*tab[nb_tmp3+5*nb3].real();
                  tmp7[0]=w[6][0]*tab[nb_tmp3+6*nb3].real()-w[6][1]*tab[nb_tmp3+6*nb3].imag();
                  tmp7[1]=w[6][0]*tab[nb_tmp3+6*nb3].imag()+w[6][1]*tab[nb_tmp3+6*nb3].real();
                  tmp8[0]=w[7][0]*tab[nb_tmp3+7*nb3].real()-w[7][1]*tab[nb_tmp3+7*nb3].imag();
                  tmp8[1]=w[7][0]*tab[nb_tmp3+7*nb3].imag()+w[7][1]*tab[nb_tmp3+7*nb3].real();
                  tmp9[0]=w[8][0]*tab[nb_tmp3+8*nb3].real()-w[8][1]*tab[nb_tmp3+8*nb3].imag();
                  tmp9[1]=w[8][0]*tab[nb_tmp3+8*nb3].imag()+w[8][1]*tab[nb_tmp3+8*nb3].real();
                  tmp10[0]=w[9][0]*tab[nb_tmp3+9*nb3].real()-w[9][1]*tab[nb_tmp3+9*nb3].imag();
                  tmp10[1]=w[9][0]*tab[nb_tmp3+9*nb3].imag()+w[9][1]*tab[nb_tmp3+9*nb3].real();
                  tmp11[0]=w[10][0]*tab[nb_tmp3+10*nb3].real()-w[10][1]*tab[nb_tmp3+10*nb3].imag();
                  tmp11[1]=w[10][0]*tab[nb_tmp3+10*nb3].imag()+w[10][1]*tab[nb_tmp3+10*nb3].real();

                  tmp30[0]=z[0][0]*tmp1[0]-z[0][1]*tmp1[1];
                  tmp30[1]=z[0][0]*tmp1[1]+z[0][1]*tmp1[0];
                  tmp20[0]=z[0][0]*(tmp1[0]+tmp2[0]+tmp3[0]+tmp4[0]+tmp5[0]+tmp6[0]+tmp7[0]+tmp8[0]+tmp9[0]+tmp10[0]+tmp11[0])-z[0][1]*(tmp1[1]+tmp2[1]+tmp3[1]+tmp4[1]+tmp5[1]+tmp6[1]+tmp7[1]+tmp8[1]+tmp9[1]+tmp10[1]+tmp11[1]);
                  tmp20[1]=z[0][0]*(tmp1[1]+tmp2[1]+tmp3[1]+tmp4[1]+tmp5[1]+tmp6[1]+tmp7[1]+tmp8[1]+tmp9[1]+tmp10[1]+tmp11[1])+z[0][1]*(tmp1[0]+tmp2[0]+tmp3[0]+tmp4[0]+tmp5[0]+tmp6[0]+tmp7[0]+tmp8[0]+tmp9[0]+tmp10[0]+tmp11[0]);

                 //radix-11
                  tab[nb_tmp3+0*nb3].real(tmp20[0]);
                  tab[nb_tmp3+0*nb3].imag(tmp20[1]);

                  tab[nb_tmp3+1*nb3].real(tmp30[0]+z[1][0]*tmp2[0]-z[1][1]*tmp2[1]+z[2][0]*tmp3[0]-z[2][1]*tmp3[1]+z[3][0]*tmp4[0]-z[3][1]*tmp4[1]+z[4][0]*tmp5[0]-z[4][1]*tmp5[1]+z[5][0]*tmp6[0]-z[5][1]*tmp6[1]+
                                            z[6][0]*tmp7[0]-z[6][1]*tmp7[1]+z[7][0]*tmp8[0]-z[7][1]*tmp8[1]+z[8][0]*tmp9[0]-z[8][1]*tmp9[1]+z[9][0]*tmp10[0]-z[9][1]*tmp10[1]+z[10][0]*tmp11[0]-z[10][1]*tmp11[1]);
                  tab[nb_tmp3+1*nb3].imag(tmp30[1]+z[1][0]*tmp2[1]+z[1][1]*tmp2[0]+z[2][0]*tmp3[1]+z[2][1]*tmp3[0]+z[3][0]*tmp4[1]+z[3][1]*tmp4[0]+z[4][0]*tmp5[1]+z[4][1]*tmp5[0]+z[5][0]*tmp6[1]+z[5][1]*tmp6[0]+
                                            z[6][0]*tmp7[1]+z[6][1]*tmp7[0]+z[7][0]*tmp8[1]+z[7][1]*tmp8[0]+z[8][0]*tmp9[1]+z[8][1]*tmp9[0]+z[9][0]*tmp10[1]+z[9][1]*tmp10[0]+z[10][0]*tmp11[1]+z[10][1]*tmp11[0]);

                  tab[nb_tmp3+2*nb3].real(tmp30[0]+z[2][0]*tmp2[0]-z[2][1]*tmp2[1]+z[4][0]*tmp3[0]-z[4][1]*tmp3[1]+z[6][0]*tmp4[0]-z[6][1]*tmp4[1]+z[8][0]*tmp5[0]-z[8][1]*tmp5[1]+z[10][0]*tmp6[0]-z[10][1]*tmp6[1]+
                                            z[1][0]*tmp7[0]-z[1][1]*tmp7[1]+z[3][0]*tmp8[0]-z[3][1]*tmp8[1]+z[5][0]*tmp9[0]-z[5][1]*tmp9[1]+z[7][0]*tmp10[0]-z[7][1]*tmp10[1]+z[9][0]*tmp11[0]-z[9][1]*tmp11[1]);
                  tab[nb_tmp3+2*nb3].imag(tmp30[1]+z[2][0]*tmp2[1]+z[2][1]*tmp2[0]+z[4][0]*tmp3[1]+z[4][1]*tmp3[0]+z[6][0]*tmp4[1]+z[6][1]*tmp4[0]+z[8][0]*tmp5[1]+z[8][1]*tmp5[0]+z[10][0]*tmp6[1]+z[10][1]*tmp6[0]+
                                            z[1][0]*tmp7[1]+z[1][1]*tmp7[0]+z[3][0]*tmp8[1]+z[3][1]*tmp8[0]+z[5][0]*tmp9[1]+z[5][1]*tmp9[0]+z[7][0]*tmp10[1]+z[7][1]*tmp10[0]+z[9][0]*tmp11[1]+z[9][1]*tmp11[0]);

                  tab[nb_tmp3+3*nb3].real(tmp30[0]+z[3][0]*tmp2[0]-z[3][1]*tmp2[1]+z[6][0]*tmp3[0]-z[6][1]*tmp3[1]+z[9][0]*tmp4[0]-z[9][1]*tmp4[1]+z[1][0]*tmp5[0]-z[1][1]*tmp5[1]+z[4][0]*tmp6[0]-z[4][1]*tmp6[1]+
                                            z[7][0]*tmp7[0]-z[7][1]*tmp7[1]+z[10][0]*tmp8[0]-z[10][1]*tmp8[1]+z[2][0]*tmp9[0]-z[2][1]*tmp9[1]+z[5][0]*tmp10[0]-z[5][1]*tmp10[1]+z[8][0]*tmp11[0]-z[8][1]*tmp11[1]);
                  tab[nb_tmp3+3*nb3].imag(tmp30[1]+z[3][0]*tmp2[1]+z[3][1]*tmp2[0]+z[6][0]*tmp3[1]+z[6][1]*tmp3[0]+z[9][0]*tmp4[1]+z[9][1]*tmp4[0]+z[1][0]*tmp5[1]+z[1][1]*tmp5[0]+z[4][0]*tmp6[1]+z[4][1]*tmp6[0]+
                                            z[7][0]*tmp7[1]+z[7][1]*tmp7[0]+z[10][0]*tmp8[1]+z[10][1]*tmp8[0]+z[2][0]*tmp9[1]+z[2][1]*tmp9[0]+z[5][0]*tmp10[1]+z[5][1]*tmp10[0]+z[8][0]*tmp11[1]+z[8][1]*tmp11[0]);

                  tab[nb_tmp3+4*nb3].real(tmp30[0]+z[4][0]*tmp2[0]-z[4][1]*tmp2[1]+z[8][0]*tmp3[0]-z[8][1]*tmp3[1]+z[1][0]*tmp4[0]-z[1][1]*tmp4[1]+z[5][0]*tmp5[0]-z[5][1]*tmp5[1]+z[9][0]*tmp6[0]-z[9][1]*tmp6[1]+
                                            z[2][0]*tmp7[0]-z[2][1]*tmp7[1]+z[6][0]*tmp8[0]-z[6][1]*tmp8[1]+z[10][0]*tmp9[0]-z[10][1]*tmp9[1]+z[3][0]*tmp10[0]-z[3][1]*tmp10[1]+z[7][0]*tmp11[0]-z[7][1]*tmp11[1]);
                  tab[nb_tmp3+4*nb3].imag(tmp30[1]+z[4][0]*tmp2[1]+z[4][1]*tmp2[0]+z[8][0]*tmp3[1]+z[8][1]*tmp3[0]+z[1][0]*tmp4[1]+z[1][1]*tmp4[0]+z[5][0]*tmp5[1]+z[5][1]*tmp5[0]+z[9][0]*tmp6[1]+z[9][1]*tmp6[0]+
                                            z[2][0]*tmp7[1]+z[2][1]*tmp7[0]+z[6][0]*tmp8[1]+z[6][1]*tmp8[0]+z[10][0]*tmp9[1]+z[10][1]*tmp9[0]+z[3][0]*tmp10[1]+z[3][1]*tmp10[0]+z[7][0]*tmp11[1]+z[7][1]*tmp11[0]);

                  tab[nb_tmp3+5*nb3].real(tmp30[0]+z[5][0]*tmp2[0]-z[5][1]*tmp2[1]+z[10][0]*tmp3[0]-z[10][1]*tmp3[1]+z[4][0]*tmp4[0]-z[4][1]*tmp4[1]+z[9][0]*tmp5[0]-z[9][1]*tmp5[1]+z[3][0]*tmp6[0]-z[3][1]*tmp6[1]+
                                            z[8][0]*tmp7[0]-z[8][1]*tmp7[1]+z[2][0]*tmp8[0]-z[2][1]*tmp8[1]+z[7][0]*tmp9[0]-z[7][1]*tmp9[1]+z[1][0]*tmp10[0]-z[1][1]*tmp10[1]+z[6][0]*tmp11[0]-z[6][1]*tmp11[1]);
                  tab[nb_tmp3+5*nb3].imag(tmp30[1]+z[5][0]*tmp2[1]+z[5][1]*tmp2[0]+z[10][0]*tmp3[1]+z[10][1]*tmp3[0]+z[4][0]*tmp4[1]+z[4][1]*tmp4[0]+z[9][0]*tmp5[1]+z[9][1]*tmp5[0]+z[3][0]*tmp6[1]+z[3][1]*tmp6[0]+
                                            z[8][0]*tmp7[1]+z[8][1]*tmp7[0]+z[2][0]*tmp8[1]+z[2][1]*tmp8[0]+z[7][0]*tmp9[1]+z[7][1]*tmp9[0]+z[1][0]*tmp10[1]+z[1][1]*tmp10[0]+z[6][0]*tmp11[1]+z[6][1]*tmp11[0]);

                  tab[nb_tmp3+6*nb3].real(tmp30[0]+z[6][0]*tmp2[0]-z[6][1]*tmp2[1]+z[1][0]*tmp3[0]-z[1][1]*tmp3[1]+z[7][0]*tmp4[0]-z[7][1]*tmp4[1]+z[2][0]*tmp5[0]-z[2][1]*tmp5[1]+z[8][0]*tmp6[0]-z[8][1]*tmp6[1]+
                                            z[3][0]*tmp7[0]-z[3][1]*tmp7[1]+z[9][0]*tmp8[0]-z[9][1]*tmp8[1]+z[4][0]*tmp9[0]-z[4][1]*tmp9[1]+z[10][0]*tmp10[0]-z[10][1]*tmp10[1]+z[5][0]*tmp11[0]-z[5][1]*tmp11[1]);
                  tab[nb_tmp3+6*nb3].imag(tmp30[1]+z[6][0]*tmp2[1]+z[6][1]*tmp2[0]+z[1][0]*tmp3[1]+z[1][1]*tmp3[0]+z[7][0]*tmp4[1]+z[7][1]*tmp4[0]+z[2][0]*tmp5[1]+z[2][1]*tmp5[0]+z[8][0]*tmp6[1]+z[8][1]*tmp6[0]+
                                            z[3][0]*tmp7[1]+z[3][1]*tmp7[0]+z[9][0]*tmp8[1]+z[9][1]*tmp8[0]+z[4][0]*tmp9[1]+z[4][1]*tmp9[0]+z[10][0]*tmp10[1]+z[10][1]*tmp10[0]+z[5][0]*tmp11[1]+z[5][1]*tmp11[0]);

                  tab[nb_tmp3+7*nb3].real(tmp30[0]+z[7][0]*tmp2[0]-z[7][1]*tmp2[1]+z[3][0]*tmp3[0]-z[3][1]*tmp3[1]+z[10][0]*tmp4[0]-z[10][1]*tmp4[1]+z[6][0]*tmp5[0]-z[6][1]*tmp5[1]+z[2][0]*tmp6[0]-z[2][1]*tmp6[1]+
                                            z[9][0]*tmp7[0]-z[9][1]*tmp7[1]+z[5][0]*tmp8[0]-z[5][1]*tmp8[1]+z[1][0]*tmp9[0]-z[1][1]*tmp9[1]+z[8][0]*tmp10[0]-z[8][1]*tmp10[1]+z[4][0]*tmp11[0]-z[4][1]*tmp11[1]);
                  tab[nb_tmp3+7*nb3].imag(tmp30[1]+z[7][0]*tmp2[1]+z[7][1]*tmp2[0]+z[3][0]*tmp3[1]+z[3][1]*tmp3[0]+z[10][0]*tmp4[1]+z[10][1]*tmp4[0]+z[6][0]*tmp5[1]+z[6][1]*tmp5[0]+z[2][0]*tmp6[1]+z[2][1]*tmp6[0]+
                                            z[9][0]*tmp7[1]+z[9][1]*tmp7[0]+z[5][0]*tmp8[1]+z[5][1]*tmp8[0]+z[1][0]*tmp9[1]+z[1][1]*tmp9[0]+z[8][0]*tmp10[1]+z[8][1]*tmp10[0]+z[4][0]*tmp11[1]+z[4][1]*tmp11[0]);

                  tab[nb_tmp3+8*nb3].real(tmp30[0]+z[8][0]*tmp2[0]-z[8][1]*tmp2[1]+z[5][0]*tmp3[0]-z[5][1]*tmp3[1]+z[2][0]*tmp4[0]-z[2][1]*tmp4[1]+z[10][0]*tmp5[0]-z[10][1]*tmp5[1]+z[7][0]*tmp6[0]-z[7][1]*tmp6[1]+
                                            z[4][0]*tmp7[0]-z[4][1]*tmp7[1]+z[1][0]*tmp8[0]-z[1][1]*tmp8[1]+z[9][0]*tmp9[0]-z[9][1]*tmp9[1]+z[6][0]*tmp10[0]-z[6][1]*tmp10[1]+z[3][0]*tmp11[0]-z[3][1]*tmp11[1]);
                  tab[nb_tmp3+8*nb3].imag(tmp30[1]+z[8][0]*tmp2[1]+z[8][1]*tmp2[0]+z[5][0]*tmp3[1]+z[5][1]*tmp3[0]+z[2][0]*tmp4[1]+z[2][1]*tmp4[0]+z[10][0]*tmp5[1]+z[10][1]*tmp5[0]+z[7][0]*tmp6[1]+z[7][1]*tmp6[0]+
                                            z[4][0]*tmp7[1]+z[4][1]*tmp7[0]+z[1][0]*tmp8[1]+z[1][1]*tmp8[0]+z[9][0]*tmp9[1]+z[9][1]*tmp9[0]+z[6][0]*tmp10[1]+z[6][1]*tmp10[0]+z[3][0]*tmp11[1]+z[3][1]*tmp11[0]);

                  tab[nb_tmp3+9*nb3].real(tmp30[0]+z[9][0]*tmp2[0]-z[9][1]*tmp2[1]+z[7][0]*tmp3[0]-z[7][1]*tmp3[1]+z[5][0]*tmp4[0]-z[5][1]*tmp4[1]+z[3][0]*tmp5[0]-z[3][1]*tmp5[1]+z[1][0]*tmp6[0]-z[1][1]*tmp6[1]+
                                            z[10][0]*tmp7[0]-z[10][1]*tmp7[1]+z[8][0]*tmp8[0]-z[8][1]*tmp8[1]+z[6][0]*tmp9[0]-z[6][1]*tmp9[1]+z[4][0]*tmp10[0]-z[4][1]*tmp10[1]+z[2][0]*tmp11[0]-z[2][1]*tmp11[1]);
                  tab[nb_tmp3+9*nb3].imag(tmp30[1]+z[9][0]*tmp2[1]+z[9][1]*tmp2[0]+z[7][0]*tmp3[1]+z[7][1]*tmp3[0]+z[5][0]*tmp4[1]+z[5][1]*tmp4[0]+z[3][0]*tmp5[1]+z[3][1]*tmp5[0]+z[1][0]*tmp6[1]+z[1][1]*tmp6[0]+
                                            z[10][0]*tmp7[1]+z[10][1]*tmp7[0]+z[8][0]*tmp8[1]+z[8][1]*tmp8[0]+z[6][0]*tmp9[1]+z[6][1]*tmp9[0]+z[4][0]*tmp10[1]+z[4][1]*tmp10[0]+z[2][0]*tmp11[1]+z[2][1]*tmp11[0]);

                  tab[nb_tmp3+10*nb3].real(tmp30[0]+z[10][0]*tmp2[0]-z[10][1]*tmp2[1]+z[9][0]*tmp3[0]-z[9][1]*tmp3[1]+z[8][0]*tmp4[0]-z[8][1]*tmp4[1]+z[7][0]*tmp5[0]-z[7][1]*tmp5[1]+z[6][0]*tmp6[0]-z[6][1]*tmp6[1]+
                                            z[5][0]*tmp7[0]-z[5][1]*tmp7[1]+z[4][0]*tmp8[0]-z[4][1]*tmp8[1]+z[3][0]*tmp9[0]-z[3][1]*tmp9[1]+z[2][0]*tmp10[0]-z[2][1]*tmp10[1]+z[1][0]*tmp11[0]-z[1][1]*tmp11[1]);
                  tab[nb_tmp3+10*nb3].imag(tmp30[1]+z[10][0]*tmp2[1]+z[10][1]*tmp2[0]+z[9][0]*tmp3[1]+z[9][1]*tmp3[0]+z[8][0]*tmp4[1]+z[8][1]*tmp4[0]+z[7][0]*tmp5[1]+z[7][1]*tmp5[0]+z[6][0]*tmp6[1]+z[6][1]*tmp6[0]+
                                            z[5][0]*tmp7[1]+z[5][1]*tmp7[0]+z[4][0]*tmp8[1]+z[4][1]*tmp8[0]+z[3][0]*tmp9[1]+z[3][1]*tmp9[0]+z[2][0]*tmp10[1]+z[2][1]*tmp10[0]+z[1][0]*tmp11[1]+z[1][1]*tmp11[0]);
                }

            }
        }
    }
///////////////////////////////////////////////////////////////////////


    void fun_inverse_table_FFT(int M,std::complex<double> tab[])
    {
        //author marcin matysek (r)ewertyn.PL
        //open-source
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
///////////////////////////////////////

     int radix_base(int N,int stg[])
        {
        int k=0;
        double M=(double)N;
        double epsilon1;
        stg[0]=1;
        //*flg2=0;
        //cout<<"M= "<<M<<endl;
        epsilon1=0.000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001;

        for(int j=0;j<200;j++)
        {
           if(fmod(M,11)<=epsilon1)
            {
                k++;
                M=M/11.0;
                stg[k]=11;
            }
            else if(fmod(M,7)<=epsilon1)
            {
                k++;
                M=M/7.0;
                stg[k]=7;
            }
            else if(fmod(M,5)<=epsilon1)
            {
                k++;
                M=M/5.0;
                stg[k]=5;
            }
            else if(fmod(M,4)<=epsilon1)
            {
                k++;
                M=M/4.0;
                stg[k]=4;
            }
            else if(fmod(M,3)<=epsilon1)
            {
                k++;
                M=M/3.0;
                stg[k]=3;
            }
            else if(fmod(M,2)<=epsilon1)
            {
                k++;
                M=M/2.0;
                stg[k]=2;
            }
            else if(M>=1.0-epsilon1&&M<=1.0+epsilon1)
            {
                //*flag=*flag+1;
                 //cout<<"*flag= "<<*flag<<" N= "<<N<<endl;
                break;
            }
            else
            {

                cout<<endl<< "Unsupported signal: N= "<<N<< endl;
                if(k>0)
                {
                    for(int m=1;m<=k;m++)
                    {
                     cout <<"stage:"<<m<<" = radix-"<<stg[m] << endl;
                    }
                }
                cout <<"stage:"<<k+1<<" = radix-??" << endl;

                k=0;
                break;
            }
            //*flg2=*flg2+1;
        }
       return k;
    }
/////////////////////////////////////////////////////////////


	void fun_fourier_transform_FFT_radix_2_stg_first(std::complex<double> tab[],int stg_first,double fi2,double (*z)[2])
    {
        double tmp1[2],tmp2[2],tmp20[2];
        double w[2][2]={{1,0}};

        w[0][0]=cos(0+fi2);
        w[0][1]=-sin(0+fi2);
        w[1][0]=cos(0+fi2);
        w[1][1]=-sin(0+fi2);

        for(int i=0;i<stg_first;i=i+1)
        {
            tmp1[0]=w[0][0]*tab[i+0*stg_first].real()-w[0][1]*tab[i+0*stg_first].imag();
            tmp1[1]=w[0][0]*tab[i+0*stg_first].imag()+w[0][1]*tab[i+0*stg_first].real();
            tmp2[0]=w[1][0]*tab[i+1*stg_first].real()-w[1][1]*tab[i+1*stg_first].imag();
            tmp2[1]=w[1][0]*tab[i+1*stg_first].imag()+w[1][1]*tab[i+1*stg_first].real();

          tmp20[0]=z[0][0]*(tmp1[0]+tmp2[0])-z[0][1]*(tmp1[1]+tmp2[1]);
          tmp20[1]=z[0][0]*(tmp1[1]+tmp2[1])+z[0][1]*(tmp1[0]+tmp2[0]);

          //radix-2
          tab[i+0*stg_first].real(tmp20[0]);
          tab[i+0*stg_first].imag(tmp20[1]);

          tab[i+1*stg_first].real(z[0][0]*tmp1[0]-z[0][1]*tmp1[1]+z[1][0]*tmp2[0]-z[1][1]*tmp2[1]);
          tab[i+1*stg_first].imag(z[0][0]*tmp1[1]+z[0][1]*tmp1[0]+z[1][0]*tmp2[1]+z[1][1]*tmp2[0]);
        }
    }
///////////////////////////////////////////////////////////////


   void fun_fourier_transform_FFT_radix_3_stg_first(std::complex<double> tab[],int stg_first,double fi2,double (*z)[2])
    {
        double tmp1[2],tmp2[2],tmp3[2],tmp30[2],tmp20[2];
        double  w[3][2]={{1,0}};

        w[0][0]=cos(0+fi2);
        w[0][1]=-sin(0+fi2);
        w[1][0]=cos(0+fi2);
        w[1][1]=-sin(0+fi2);
        w[2][0]=cos(0+fi2);
        w[2][1]=-sin(0+fi2);

        for(int i=0;i<stg_first;i=i+1)
        {
            tmp1[0]=w[0][0]*tab[i+0*stg_first].real()-w[0][1]*tab[i+0*stg_first].imag();
            tmp1[1]=w[0][0]*tab[i+0*stg_first].imag()+w[0][1]*tab[i+0*stg_first].real();
            tmp2[0]=w[1][0]*tab[i+1*stg_first].real()-w[1][1]*tab[i+1*stg_first].imag();
            tmp2[1]=w[1][0]*tab[i+1*stg_first].imag()+w[1][1]*tab[i+1*stg_first].real();
            tmp3[0]=w[2][0]*tab[i+2*stg_first].real()-w[2][1]*tab[i+2*stg_first].imag();
            tmp3[1]=w[2][0]*tab[i+2*stg_first].imag()+w[2][1]*tab[i+2*stg_first].real();

          tmp30[0]=z[0][0]*tmp1[0]-z[0][1]*tmp1[1];
          tmp30[1]=z[0][0]*tmp1[1]+z[0][1]*tmp1[0];
          tmp20[0]=z[0][0]*(tmp1[0]+tmp2[0]+tmp3[0])-z[0][1]*(tmp1[1]+tmp2[1]+tmp3[1]);
          tmp20[1]=z[0][0]*(tmp1[1]+tmp2[1]+tmp3[1])+z[0][1]*(tmp1[0]+tmp2[0]+tmp3[0]);


         //radix-3
          tab[i+0*stg_first].real(tmp20[0]);
          tab[i+0*stg_first].imag(tmp20[1]);

          tab[i+1*stg_first].real(tmp30[0]+z[1][0]*tmp2[0]-z[1][1]*tmp2[1]+z[2][0]*tmp3[0]-z[2][1]*tmp3[1]);
          tab[i+1*stg_first].imag(tmp30[1]+z[1][0]*tmp2[1]+z[1][1]*tmp2[0]+z[2][0]*tmp3[1]+z[2][1]*tmp3[0]);

          tab[i+2*stg_first].real(tmp30[0]+z[2][0]*tmp2[0]-z[2][1]*tmp2[1]+z[1][0]*tmp3[0]-z[1][1]*tmp3[1]);
          tab[i+2*stg_first].imag(tmp30[1]+z[2][0]*tmp2[1]+z[2][1]*tmp2[0]+z[1][0]*tmp3[1]+z[1][1]*tmp3[0]);
        }
    }
///////////////////////////////////////////////////////////////


   void fun_fourier_transform_FFT_radix_4_stg_first(std::complex<double> tab[],int stg_first,double fi2,double (*z)[2])
    {
        double tmp1[2],tmp2[2],tmp3[2],tmp4[2];
        double tmp101[2],tmp102[2],tmp103[2],tmp104[2];
        double tmp111[2],tmp112[2],tmp113[2],tmp114[2];
        double  w[4][2]={{1,0}};

        w[0][0]=cos(0+fi2);
        w[0][1]=-sin(0+fi2);
        w[1][0]=cos(0+fi2);
        w[1][1]=-sin(0+fi2);
        w[2][0]=cos(0+fi2);
        w[2][1]=-sin(0+fi2);
        w[3][0]=cos(0+fi2);
        w[3][1]=-sin(0+fi2);
        for(int i=0;i<stg_first;i=i+1)
        {
            tmp1[0]=w[0][0]*tab[i+0*stg_first].real()-w[0][1]*tab[i+0*stg_first].imag();
            tmp1[1]=w[0][0]*tab[i+0*stg_first].imag()+w[0][1]*tab[i+0*stg_first].real();
            tmp2[0]=w[1][0]*tab[i+1*stg_first].real()-w[1][1]*tab[i+1*stg_first].imag();
            tmp2[1]=w[1][0]*tab[i+1*stg_first].imag()+w[1][1]*tab[i+1*stg_first].real();
            tmp3[0]=w[2][0]*tab[i+2*stg_first].real()-w[2][1]*tab[i+2*stg_first].imag();
            tmp3[1]=w[2][0]*tab[i+2*stg_first].imag()+w[2][1]*tab[i+2*stg_first].real();
            tmp4[0]=w[3][0]*tab[i+3*stg_first].real()-w[3][1]*tab[i+3*stg_first].imag();
            tmp4[1]=w[3][0]*tab[i+3*stg_first].imag()+w[3][1]*tab[i+3*stg_first].real();


            tmp101[0]=tmp2[0]-tmp4[0];
            tmp101[1]=tmp2[1]-tmp4[1];
            tmp102[0]=tmp2[0]+tmp4[0];
            tmp102[1]=tmp2[1]+tmp4[1];
            tmp103[0]=tmp1[0]-tmp3[0];
            tmp103[1]=tmp1[1]-tmp3[1];
            tmp104[0]=tmp1[0]+tmp3[0];
            tmp104[1]=tmp1[1]+tmp3[1];

            tmp111[0]=z[0][0]*(tmp104[0]+tmp102[0])-z[0][1]*(tmp104[1]+tmp102[1]);
            tmp111[1]=z[0][0]*(tmp104[1]+tmp102[1])+z[0][1]*(tmp104[0]+tmp102[0]);
            tmp112[0]=z[0][0]*tmp103[0]-z[0][1]*tmp103[1];
            tmp112[1]=z[0][0]*tmp103[1]+z[0][1]*tmp103[0];
            tmp113[0]=z[0][0]*(tmp104[0]-tmp102[0])-z[0][1]*(tmp104[1]-tmp102[1]);
            tmp113[1]=z[0][0]*(tmp104[1]-tmp102[1])+z[0][1]*(tmp104[0]-tmp102[0]);
            tmp114[0]=z[0][0]*tmp103[0]-z[0][1]*tmp103[1];
            tmp114[1]=z[0][0]*tmp103[1]+z[0][1]*tmp103[0];

            //radix-4
            tab[i+0*stg_first].real(tmp111[0]);
            tab[i+0*stg_first].imag(tmp111[1]);
            tab[i+1*stg_first].real(tmp112[0]+z[1][0]*tmp101[0]-z[1][1]*tmp101[1]);
            tab[i+1*stg_first].imag(tmp112[1]+z[1][0]*tmp101[1]+z[1][1]*tmp101[0]);
            tab[i+2*stg_first].real(tmp113[0]);
            tab[i+2*stg_first].imag(tmp113[1]);
            tab[i+3*stg_first].real(tmp114[0]-(z[1][0]*tmp101[0]-z[1][1]*tmp101[1]));
            tab[i+3*stg_first].imag(tmp114[1]-(z[1][0]*tmp101[1]+z[1][1]*tmp101[0]));
        }
    }
/////////////////////////////////////////////////////////////////


void fun_fourier_transform_FFT_radix_5_stg_first(std::complex<double> tab[],int stg_first,double fi2,double  (*z)[2])
  {
        double tmp1[2],tmp2[2],tmp3[2],tmp4[2],tmp5[2],tmp30[2],tmp20[2];
        double  w[5][2]={{1,0}};

        w[0][0]=cos(0+fi2);
        w[0][1]=-sin(0+fi2);
        w[1][0]=cos(0+fi2);
        w[1][1]=-sin(0+fi2);
        w[2][0]=cos(0+fi2);
        w[2][1]=-sin(0+fi2);
        w[3][0]=cos(0+fi2);
        w[3][1]=-sin(0+fi2);
        w[4][0]=cos(0+fi2);
        w[4][1]=-sin(0+fi2);

        for(int i=0;i<stg_first;i++)
        {
          tmp1[0]=w[0][0]*tab[i+0*stg_first].real()-w[0][1]*tab[i+0*stg_first].imag();
          tmp1[1]=w[0][0]*tab[i+0*stg_first].imag()+w[0][1]*tab[i+0*stg_first].real();
          tmp2[0]=w[1][0]*tab[i+1*stg_first].real()-w[1][1]*tab[i+1*stg_first].imag();
          tmp2[1]=w[1][0]*tab[i+1*stg_first].imag()+w[1][1]*tab[i+1*stg_first].real();
          tmp3[0]=w[2][0]*tab[i+2*stg_first].real()-w[2][1]*tab[i+2*stg_first].imag();
          tmp3[1]=w[2][0]*tab[i+2*stg_first].imag()+w[2][1]*tab[i+2*stg_first].real();
          tmp4[0]=w[3][0]*tab[i+3*stg_first].real()-w[3][1]*tab[i+3*stg_first].imag();
          tmp4[1]=w[3][0]*tab[i+3*stg_first].imag()+w[3][1]*tab[i+3*stg_first].real();
          tmp5[0]=w[4][0]*tab[i+4*stg_first].real()-w[4][1]*tab[i+4*stg_first].imag();
          tmp5[1]=w[4][0]*tab[i+4*stg_first].imag()+w[4][1]*tab[i+4*stg_first].real();

          tmp30[0]=z[0][0]*tmp1[0]-z[0][1]*tmp1[1];
          tmp30[1]=z[0][0]*tmp1[1]+z[0][1]*tmp1[0];
          tmp20[0]=z[0][0]*(tmp1[0]+tmp2[0]+tmp3[0]+tmp4[0]+tmp5[0])-z[0][1]*(tmp1[1]+tmp2[1]+tmp3[1]+tmp4[1]+tmp5[1]);
          tmp20[1]=z[0][0]*(tmp1[1]+tmp2[1]+tmp3[1]+tmp4[1]+tmp5[1])+z[0][1]*(tmp1[0]+tmp2[0]+tmp3[0]+tmp4[0]+tmp5[0]);

         //radix-5
          tab[i+0*stg_first].real(tmp20[0]);
          tab[i+0*stg_first].imag(tmp20[1]);

          tab[i+1*stg_first].real(tmp30[0]+z[1][0]*tmp2[0]-z[1][1]*tmp2[1]+z[2][0]*tmp3[0]-z[2][1]*tmp3[1]+z[3][0]*tmp4[0]-z[3][1]*tmp4[1]+z[4][0]*tmp5[0]-z[4][1]*tmp5[1]);
          tab[i+1*stg_first].imag(tmp30[1]+z[1][0]*tmp2[1]+z[1][1]*tmp2[0]+z[2][0]*tmp3[1]+z[2][1]*tmp3[0]+z[3][0]*tmp4[1]+z[3][1]*tmp4[0]+z[4][0]*tmp5[1]+z[4][1]*tmp5[0]);

          tab[i+2*stg_first].real(tmp30[0]+z[2][0]*tmp2[0]-z[2][1]*tmp2[1]+z[4][0]*tmp3[0]-z[4][1]*tmp3[1]+z[1][0]*tmp4[0]-z[1][1]*tmp4[1]+z[3][0]*tmp5[0]-z[3][1]*tmp5[1]);
          tab[i+2*stg_first].imag(tmp30[1]+z[2][0]*tmp2[1]+z[2][1]*tmp2[0]+z[4][0]*tmp3[1]+z[4][1]*tmp3[0]+z[1][0]*tmp4[1]+z[1][1]*tmp4[0]+z[3][0]*tmp5[1]+z[3][1]*tmp5[0]);

          tab[i+3*stg_first].real(tmp30[0]+z[3][0]*tmp2[0]-z[3][1]*tmp2[1]+z[1][0]*tmp3[0]-z[1][1]*tmp3[1]+z[4][0]*tmp4[0]-z[4][1]*tmp4[1]+z[2][0]*tmp5[0]-z[2][1]*tmp5[1]);
          tab[i+3*stg_first].imag(tmp30[1]+z[3][0]*tmp2[1]+z[3][1]*tmp2[0]+z[1][0]*tmp3[1]+z[1][1]*tmp3[0]+z[4][0]*tmp4[1]+z[4][1]*tmp4[0]+z[2][0]*tmp5[1]+z[2][1]*tmp5[0]);

          tab[i+4*stg_first].real(tmp30[0]+z[4][0]*tmp2[0]-z[4][1]*tmp2[1]+z[3][0]*tmp3[0]-z[3][1]*tmp3[1]+z[2][0]*tmp4[0]-z[2][1]*tmp4[1]+z[1][0]*tmp5[0]-z[1][1]*tmp5[1]);
          tab[i+4*stg_first].imag(tmp30[1]+z[4][0]*tmp2[1]+z[4][1]*tmp2[0]+z[3][0]*tmp3[1]+z[3][1]*tmp3[0]+z[2][0]*tmp4[1]+z[2][1]*tmp4[0]+z[1][0]*tmp5[1]+z[1][1]*tmp5[0]);
        }
  }
/////////////////////////////////////////


    void fun_fourier_transform_FFT_radix_7_stg_first(std::complex<double> tab[],int stg_first,double fi2,double  (*z)[2])
  {
        double tmp1[2],tmp2[2],tmp3[2],tmp4[2],tmp5[2],tmp6[2],tmp7[2],tmp30[2],tmp20[2];
        double  w[7][2]={{1,0}};


        w[0][0]=cos(0+fi2);
        w[0][1]=-sin(0+fi2);
        w[1][0]=cos(0+fi2);
        w[1][1]=-sin(0+fi2);
        w[2][0]=cos(0+fi2);
        w[2][1]=-sin(0+fi2);
        w[3][0]=cos(0+fi2);
        w[3][1]=-sin(0+fi2);
        w[4][0]=cos(0+fi2);
        w[4][1]=-sin(0+fi2);
        w[5][0]=cos(0+fi2);
        w[5][1]=-sin(0+fi2);
        w[6][0]=cos(0+fi2);
        w[6][1]=-sin(0+fi2);

        for(int i=0;i<stg_first;i=i+1)
        {
          tmp1[0]=w[0][0]*tab[i+0*stg_first].real()-w[0][1]*tab[i+0*stg_first].imag();
          tmp1[1]=w[0][0]*tab[i+0*stg_first].imag()+w[0][1]*tab[i+0*stg_first].real();
          tmp2[0]=w[1][0]*tab[i+1*stg_first].real()-w[1][1]*tab[i+1*stg_first].imag();
          tmp2[1]=w[1][0]*tab[i+1*stg_first].imag()+w[1][1]*tab[i+1*stg_first].real();
          tmp3[0]=w[2][0]*tab[i+2*stg_first].real()-w[2][1]*tab[i+2*stg_first].imag();
          tmp3[1]=w[2][0]*tab[i+2*stg_first].imag()+w[2][1]*tab[i+2*stg_first].real();
          tmp4[0]=w[3][0]*tab[i+3*stg_first].real()-w[3][1]*tab[i+3*stg_first].imag();
          tmp4[1]=w[3][0]*tab[i+3*stg_first].imag()+w[3][1]*tab[i+3*stg_first].real();
          tmp5[0]=w[4][0]*tab[i+4*stg_first].real()-w[4][1]*tab[i+4*stg_first].imag();
          tmp5[1]=w[4][0]*tab[i+4*stg_first].imag()+w[4][1]*tab[i+4*stg_first].real();
          tmp6[0]=w[5][0]*tab[i+5*stg_first].real()-w[5][1]*tab[i+5*stg_first].imag();
          tmp6[1]=w[5][0]*tab[i+5*stg_first].imag()+w[5][1]*tab[i+5*stg_first].real();
          tmp7[0]=w[6][0]*tab[i+6*stg_first].real()-w[6][1]*tab[i+6*stg_first].imag();
          tmp7[1]=w[6][0]*tab[i+6*stg_first].imag()+w[6][1]*tab[i+6*stg_first].real();

          tmp30[0]=z[0][0]*tmp1[0]-z[0][1]*tmp1[1];
          tmp30[1]=z[0][0]*tmp1[1]+z[0][1]*tmp1[0];
          tmp20[0]=z[0][0]*(tmp1[0]+tmp2[0]+tmp3[0]+tmp4[0]+tmp5[0]+tmp6[0]+tmp7[0])-z[0][1]*(tmp1[1]+tmp2[1]+tmp3[1]+tmp4[1]+tmp5[1]+tmp6[1]+tmp7[1]);
          tmp20[1]=z[0][0]*(tmp1[1]+tmp2[1]+tmp3[1]+tmp4[1]+tmp5[1]+tmp6[1]+tmp7[1])+z[0][1]*(tmp1[0]+tmp2[0]+tmp3[0]+tmp4[0]+tmp5[0]+tmp6[0]+tmp7[0]);

         //radix-7
          tab[i+0*stg_first].real(tmp20[0]);
          tab[i+0*stg_first].imag(tmp20[1]);

          tab[i+1*stg_first].real(tmp30[0]+z[1][0]*tmp2[0]-z[1][1]*tmp2[1]+z[2][0]*tmp3[0]-z[2][1]*tmp3[1]+z[3][0]*tmp4[0]-z[3][1]*tmp4[1]+
                                  z[4][0]*tmp5[0]-z[4][1]*tmp5[1]+z[5][0]*tmp6[0]-z[5][1]*tmp6[1]+z[6][0]*tmp7[0]-z[6][1]*tmp7[1]);
          tab[i+1*stg_first].imag(tmp30[1]+z[1][0]*tmp2[1]+z[1][1]*tmp2[0]+z[2][0]*tmp3[1]+z[2][1]*tmp3[0]+z[3][0]*tmp4[1]+z[3][1]*tmp4[0]+
                                  z[4][0]*tmp5[1]+z[4][1]*tmp5[0]+z[5][0]*tmp6[1]+z[5][1]*tmp6[0]+z[6][0]*tmp7[1]+z[6][1]*tmp7[0]);

          tab[i+2*stg_first].real(tmp30[0]+z[2][0]*tmp2[0]-z[2][1]*tmp2[1]+z[4][0]*tmp3[0]-z[4][1]*tmp3[1]+z[6][0]*tmp4[0]-z[6][1]*tmp4[1]+
                                  z[1][0]*tmp5[0]-z[1][1]*tmp5[1]+z[3][0]*tmp6[0]-z[3][1]*tmp6[1]+z[5][0]*tmp7[0]-z[5][1]*tmp7[1]);
          tab[i+2*stg_first].imag(tmp30[1]+z[2][0]*tmp2[1]+z[2][1]*tmp2[0]+z[4][0]*tmp3[1]+z[4][1]*tmp3[0]+z[6][0]*tmp4[1]+z[6][1]*tmp4[0]+
                                  z[1][0]*tmp5[1]+z[1][1]*tmp5[0]+z[3][0]*tmp6[1]+z[3][1]*tmp6[0]+z[5][0]*tmp7[1]+z[5][1]*tmp7[0]);

          tab[i+3*stg_first].real(tmp30[0]+z[3][0]*tmp2[0]-z[3][1]*tmp2[1]+z[6][0]*tmp3[0]-z[6][1]*tmp3[1]+z[2][0]*tmp4[0]-z[2][1]*tmp4[1]+
                                  z[5][0]*tmp5[0]-z[5][1]*tmp5[1]+z[1][0]*tmp6[0]-z[1][1]*tmp6[1]+z[4][0]*tmp7[0]-z[4][1]*tmp7[1]);
          tab[i+3*stg_first].imag(tmp30[1]+z[3][0]*tmp2[1]+z[3][1]*tmp2[0]+z[6][0]*tmp3[1]+z[6][1]*tmp3[0]+z[2][0]*tmp4[1]+z[2][1]*tmp4[0]+
                                  z[5][0]*tmp5[1]+z[5][1]*tmp5[0]+z[1][0]*tmp6[1]+z[1][1]*tmp6[0]+z[4][0]*tmp7[1]+z[4][1]*tmp7[0]);

          tab[i+4*stg_first].real(tmp30[0]+z[4][0]*tmp2[0]-z[4][1]*tmp2[1]+z[1][0]*tmp3[0]-z[1][1]*tmp3[1]+z[5][0]*tmp4[0]-z[5][1]*tmp4[1]+
                                  z[2][0]*tmp5[0]-z[2][1]*tmp5[1]+z[6][0]*tmp6[0]-z[6][1]*tmp6[1]+z[3][0]*tmp7[0]-z[3][1]*tmp7[1]);
          tab[i+4*stg_first].imag(tmp30[1]+z[4][0]*tmp2[1]+z[4][1]*tmp2[0]+z[1][0]*tmp3[1]+z[1][1]*tmp3[0]+z[5][0]*tmp4[1]+z[5][1]*tmp4[0]+
                                  z[2][0]*tmp5[1]+z[2][1]*tmp5[0]+z[6][0]*tmp6[1]+z[6][1]*tmp6[0]+z[3][0]*tmp7[1]+z[3][1]*tmp7[0]);

          tab[i+5*stg_first].real(tmp30[0]+z[5][0]*tmp2[0]-z[5][1]*tmp2[1]+z[3][0]*tmp3[0]-z[3][1]*tmp3[1]+z[1][0]*tmp4[0]-z[1][1]*tmp4[1]+
                                  z[6][0]*tmp5[0]-z[6][1]*tmp5[1]+z[4][0]*tmp6[0]-z[4][1]*tmp6[1]+z[2][0]*tmp7[0]-z[2][1]*tmp7[1]);
          tab[i+5*stg_first].imag(tmp30[1]+z[5][0]*tmp2[1]+z[5][1]*tmp2[0]+z[3][0]*tmp3[1]+z[3][1]*tmp3[0]+z[1][0]*tmp4[1]+z[1][1]*tmp4[0]+
                                  z[6][0]*tmp5[1]+z[6][1]*tmp5[0]+z[4][0]*tmp6[1]+z[4][1]*tmp6[0]+z[2][0]*tmp7[1]+z[2][1]*tmp7[0]);

          tab[i+6*stg_first].real(tmp30[0]+z[6][0]*tmp2[0]-z[6][1]*tmp2[1]+z[5][0]*tmp3[0]-z[5][1]*tmp3[1]+z[4][0]*tmp4[0]-z[4][1]*tmp4[1]+
                                  z[3][0]*tmp5[0]-z[3][1]*tmp5[1]+z[2][0]*tmp6[0]-z[2][1]*tmp6[1]+z[1][0]*tmp7[0]-z[1][1]*tmp7[1]);
          tab[i+6*stg_first].imag(tmp30[1]+z[6][0]*tmp2[1]+z[6][1]*tmp2[0]+z[5][0]*tmp3[1]+z[5][1]*tmp3[0]+z[4][0]*tmp4[1]+z[4][1]*tmp4[0]+
                                  z[3][0]*tmp5[1]+z[3][1]*tmp5[0]+z[2][0]*tmp6[1]+z[2][1]*tmp6[0]+z[1][0]*tmp7[1]+z[1][1]*tmp7[0]);
        }

    }
///////////////////////////////////////////////////////////////////////

 void fun_fourier_transform_FFT_radix_11_stg_first(std::complex<double> tab[],int stg_first,double fi2,double (*z)[2])
  {
       double tmp1[2],tmp2[2],tmp3[2],tmp4[2],tmp5[2],tmp6[2],tmp7[2],tmp8[2],tmp9[2],tmp10[2],tmp11[2],tmp30[2],tmp20[2];
       double  w[11][2]={{1,0}};

        w[0][0]=cos(0+fi2);
        w[0][1]=-sin(0+fi2);
        w[1][0]=cos(0+fi2);
        w[1][1]=-sin(0+fi2);
        w[2][0]=cos(0+fi2);
        w[2][1]=-sin(0+fi2);
        w[3][0]=cos(0+fi2);
        w[3][1]=-sin(0+fi2);
        w[4][0]=cos(0+fi2);
        w[4][1]=-sin(0+fi2);
        w[5][0]=cos(0+fi2);
        w[5][1]=-sin(0+fi2);
        w[6][0]=cos(0+fi2);
        w[6][1]=-sin(0+fi2);
        w[7][0]=cos(0+fi2);
        w[7][1]=-sin(0+fi2);
        w[8][0]=cos(0+fi2);
        w[8][1]=-sin(0+fi2);
        w[9][0]=cos(0+fi2);
        w[9][1]=-sin(0+fi2);
        w[10][0]=cos(0+fi2);
        w[10][1]=-sin(0+fi2);


        for(int i=0;i<stg_first;i=i+1)
        {
          tmp1[0]=w[0][0]*tab[i+0*stg_first].real()-w[0][1]*tab[i+0*stg_first].imag();
          tmp1[1]=w[0][0]*tab[i+0*stg_first].imag()+w[0][1]*tab[i+0*stg_first].real();
          tmp2[0]=w[1][0]*tab[i+1*stg_first].real()-w[1][1]*tab[i+1*stg_first].imag();
          tmp2[1]=w[1][0]*tab[i+1*stg_first].imag()+w[1][1]*tab[i+1*stg_first].real();
          tmp3[0]=w[2][0]*tab[i+2*stg_first].real()-w[2][1]*tab[i+2*stg_first].imag();
          tmp3[1]=w[2][0]*tab[i+2*stg_first].imag()+w[2][1]*tab[i+2*stg_first].real();
          tmp4[0]=w[3][0]*tab[i+3*stg_first].real()-w[3][1]*tab[i+3*stg_first].imag();
          tmp4[1]=w[3][0]*tab[i+3*stg_first].imag()+w[3][1]*tab[i+3*stg_first].real();
          tmp5[0]=w[4][0]*tab[i+4*stg_first].real()-w[4][1]*tab[i+4*stg_first].imag();
          tmp5[1]=w[4][0]*tab[i+4*stg_first].imag()+w[4][1]*tab[i+4*stg_first].real();
          tmp6[0]=w[5][0]*tab[i+5*stg_first].real()-w[5][1]*tab[i+5*stg_first].imag();
          tmp6[1]=w[5][0]*tab[i+5*stg_first].imag()+w[5][1]*tab[i+5*stg_first].real();
          tmp7[0]=w[6][0]*tab[i+6*stg_first].real()-w[6][1]*tab[i+6*stg_first].imag();
          tmp7[1]=w[6][0]*tab[i+6*stg_first].imag()+w[6][1]*tab[i+6*stg_first].real();
          tmp8[0]=w[7][0]*tab[i+7*stg_first].real()-w[7][1]*tab[i+7*stg_first].imag();
          tmp8[1]=w[7][0]*tab[i+7*stg_first].imag()+w[7][1]*tab[i+7*stg_first].real();
          tmp9[0]=w[8][0]*tab[i+8*stg_first].real()-w[8][1]*tab[i+8*stg_first].imag();
          tmp9[1]=w[8][0]*tab[i+8*stg_first].imag()+w[8][1]*tab[i+8*stg_first].real();
          tmp10[0]=w[9][0]*tab[i+9*stg_first].real()-w[9][1]*tab[i+9*stg_first].imag();
          tmp10[1]=w[9][0]*tab[i+9*stg_first].imag()+w[9][1]*tab[i+9*stg_first].real();
          tmp11[0]=w[10][0]*tab[i+10*stg_first].real()-w[10][1]*tab[i+10*stg_first].imag();
          tmp11[1]=w[10][0]*tab[i+10*stg_first].imag()+w[10][1]*tab[i+10*stg_first].real();

          tmp30[0]=z[0][0]*tmp1[0]-z[0][1]*tmp1[1];
          tmp30[1]=z[0][0]*tmp1[1]+z[0][1]*tmp1[0];
          tmp20[0]=z[0][0]*(tmp1[0]+tmp2[0]+tmp3[0]+tmp4[0]+tmp5[0]+tmp6[0]+tmp7[0]+tmp8[0]+tmp9[0]+tmp10[0]+tmp11[0])-z[0][1]*(tmp1[1]+tmp2[1]+tmp3[1]+tmp4[1]+tmp5[1]+tmp6[1]+tmp7[1]+tmp8[1]+tmp9[1]+tmp10[1]+tmp11[1]);
          tmp20[1]=z[0][0]*(tmp1[1]+tmp2[1]+tmp3[1]+tmp4[1]+tmp5[1]+tmp6[1]+tmp7[1]+tmp8[1]+tmp9[1]+tmp10[1]+tmp11[1])+z[0][1]*(tmp1[0]+tmp2[0]+tmp3[0]+tmp4[0]+tmp5[0]+tmp6[0]+tmp7[0]+tmp8[0]+tmp9[0]+tmp10[0]+tmp11[0]);

         //radix-11
          tab[i+0*stg_first].real(tmp20[0]);
          tab[i+0*stg_first].imag(tmp20[1]);

          tab[i+1*stg_first].real(tmp30[0]+z[1][0]*tmp2[0]-z[1][1]*tmp2[1]+z[2][0]*tmp3[0]-z[2][1]*tmp3[1]+z[3][0]*tmp4[0]-z[3][1]*tmp4[1]+z[4][0]*tmp5[0]-z[4][1]*tmp5[1]+z[5][0]*tmp6[0]-z[5][1]*tmp6[1]+
                                    z[6][0]*tmp7[0]-z[6][1]*tmp7[1]+z[7][0]*tmp8[0]-z[7][1]*tmp8[1]+z[8][0]*tmp9[0]-z[8][1]*tmp9[1]+z[9][0]*tmp10[0]-z[9][1]*tmp10[1]+z[10][0]*tmp11[0]-z[10][1]*tmp11[1]);
          tab[i+1*stg_first].imag(tmp30[1]+z[1][0]*tmp2[1]+z[1][1]*tmp2[0]+z[2][0]*tmp3[1]+z[2][1]*tmp3[0]+z[3][0]*tmp4[1]+z[3][1]*tmp4[0]+z[4][0]*tmp5[1]+z[4][1]*tmp5[0]+z[5][0]*tmp6[1]+z[5][1]*tmp6[0]+
                                    z[6][0]*tmp7[1]+z[6][1]*tmp7[0]+z[7][0]*tmp8[1]+z[7][1]*tmp8[0]+z[8][0]*tmp9[1]+z[8][1]*tmp9[0]+z[9][0]*tmp10[1]+z[9][1]*tmp10[0]+z[10][0]*tmp11[1]+z[10][1]*tmp11[0]);

          tab[i+2*stg_first].real(tmp30[0]+z[2][0]*tmp2[0]-z[2][1]*tmp2[1]+z[4][0]*tmp3[0]-z[4][1]*tmp3[1]+z[6][0]*tmp4[0]-z[6][1]*tmp4[1]+z[8][0]*tmp5[0]-z[8][1]*tmp5[1]+z[10][0]*tmp6[0]-z[10][1]*tmp6[1]+
                                    z[1][0]*tmp7[0]-z[1][1]*tmp7[1]+z[3][0]*tmp8[0]-z[3][1]*tmp8[1]+z[5][0]*tmp9[0]-z[5][1]*tmp9[1]+z[7][0]*tmp10[0]-z[7][1]*tmp10[1]+z[9][0]*tmp11[0]-z[9][1]*tmp11[1]);
          tab[i+2*stg_first].imag(tmp30[1]+z[2][0]*tmp2[1]+z[2][1]*tmp2[0]+z[4][0]*tmp3[1]+z[4][1]*tmp3[0]+z[6][0]*tmp4[1]+z[6][1]*tmp4[0]+z[8][0]*tmp5[1]+z[8][1]*tmp5[0]+z[10][0]*tmp6[1]+z[10][1]*tmp6[0]+
                                    z[1][0]*tmp7[1]+z[1][1]*tmp7[0]+z[3][0]*tmp8[1]+z[3][1]*tmp8[0]+z[5][0]*tmp9[1]+z[5][1]*tmp9[0]+z[7][0]*tmp10[1]+z[7][1]*tmp10[0]+z[9][0]*tmp11[1]+z[9][1]*tmp11[0]);

          tab[i+3*stg_first].real(tmp30[0]+z[3][0]*tmp2[0]-z[3][1]*tmp2[1]+z[6][0]*tmp3[0]-z[6][1]*tmp3[1]+z[9][0]*tmp4[0]-z[9][1]*tmp4[1]+z[1][0]*tmp5[0]-z[1][1]*tmp5[1]+z[4][0]*tmp6[0]-z[4][1]*tmp6[1]+
                                    z[7][0]*tmp7[0]-z[7][1]*tmp7[1]+z[10][0]*tmp8[0]-z[10][1]*tmp8[1]+z[2][0]*tmp9[0]-z[2][1]*tmp9[1]+z[5][0]*tmp10[0]-z[5][1]*tmp10[1]+z[8][0]*tmp11[0]-z[8][1]*tmp11[1]);
          tab[i+3*stg_first].imag(tmp30[1]+z[3][0]*tmp2[1]+z[3][1]*tmp2[0]+z[6][0]*tmp3[1]+z[6][1]*tmp3[0]+z[9][0]*tmp4[1]+z[9][1]*tmp4[0]+z[1][0]*tmp5[1]+z[1][1]*tmp5[0]+z[4][0]*tmp6[1]+z[4][1]*tmp6[0]+
                                    z[7][0]*tmp7[1]+z[7][1]*tmp7[0]+z[10][0]*tmp8[1]+z[10][1]*tmp8[0]+z[2][0]*tmp9[1]+z[2][1]*tmp9[0]+z[5][0]*tmp10[1]+z[5][1]*tmp10[0]+z[8][0]*tmp11[1]+z[8][1]*tmp11[0]);

          tab[i+4*stg_first].real(tmp30[0]+z[4][0]*tmp2[0]-z[4][1]*tmp2[1]+z[8][0]*tmp3[0]-z[8][1]*tmp3[1]+z[1][0]*tmp4[0]-z[1][1]*tmp4[1]+z[5][0]*tmp5[0]-z[5][1]*tmp5[1]+z[9][0]*tmp6[0]-z[9][1]*tmp6[1]+
                                    z[2][0]*tmp7[0]-z[2][1]*tmp7[1]+z[6][0]*tmp8[0]-z[6][1]*tmp8[1]+z[10][0]*tmp9[0]-z[10][1]*tmp9[1]+z[3][0]*tmp10[0]-z[3][1]*tmp10[1]+z[7][0]*tmp11[0]-z[7][1]*tmp11[1]);
          tab[i+4*stg_first].imag(tmp30[1]+z[4][0]*tmp2[1]+z[4][1]*tmp2[0]+z[8][0]*tmp3[1]+z[8][1]*tmp3[0]+z[1][0]*tmp4[1]+z[1][1]*tmp4[0]+z[5][0]*tmp5[1]+z[5][1]*tmp5[0]+z[9][0]*tmp6[1]+z[9][1]*tmp6[0]+
                                    z[2][0]*tmp7[1]+z[2][1]*tmp7[0]+z[6][0]*tmp8[1]+z[6][1]*tmp8[0]+z[10][0]*tmp9[1]+z[10][1]*tmp9[0]+z[3][0]*tmp10[1]+z[3][1]*tmp10[0]+z[7][0]*tmp11[1]+z[7][1]*tmp11[0]);

          tab[i+5*stg_first].real(tmp30[0]+z[5][0]*tmp2[0]-z[5][1]*tmp2[1]+z[10][0]*tmp3[0]-z[10][1]*tmp3[1]+z[4][0]*tmp4[0]-z[4][1]*tmp4[1]+z[9][0]*tmp5[0]-z[9][1]*tmp5[1]+z[3][0]*tmp6[0]-z[3][1]*tmp6[1]+
                                    z[8][0]*tmp7[0]-z[8][1]*tmp7[1]+z[2][0]*tmp8[0]-z[2][1]*tmp8[1]+z[7][0]*tmp9[0]-z[7][1]*tmp9[1]+z[1][0]*tmp10[0]-z[1][1]*tmp10[1]+z[6][0]*tmp11[0]-z[6][1]*tmp11[1]);
          tab[i+5*stg_first].imag(tmp30[1]+z[5][0]*tmp2[1]+z[5][1]*tmp2[0]+z[10][0]*tmp3[1]+z[10][1]*tmp3[0]+z[4][0]*tmp4[1]+z[4][1]*tmp4[0]+z[9][0]*tmp5[1]+z[9][1]*tmp5[0]+z[3][0]*tmp6[1]+z[3][1]*tmp6[0]+
                                    z[8][0]*tmp7[1]+z[8][1]*tmp7[0]+z[2][0]*tmp8[1]+z[2][1]*tmp8[0]+z[7][0]*tmp9[1]+z[7][1]*tmp9[0]+z[1][0]*tmp10[1]+z[1][1]*tmp10[0]+z[6][0]*tmp11[1]+z[6][1]*tmp11[0]);

          tab[i+6*stg_first].real(tmp30[0]+z[6][0]*tmp2[0]-z[6][1]*tmp2[1]+z[1][0]*tmp3[0]-z[1][1]*tmp3[1]+z[7][0]*tmp4[0]-z[7][1]*tmp4[1]+z[2][0]*tmp5[0]-z[2][1]*tmp5[1]+z[8][0]*tmp6[0]-z[8][1]*tmp6[1]+
                                    z[3][0]*tmp7[0]-z[3][1]*tmp7[1]+z[9][0]*tmp8[0]-z[9][1]*tmp8[1]+z[4][0]*tmp9[0]-z[4][1]*tmp9[1]+z[10][0]*tmp10[0]-z[10][1]*tmp10[1]+z[5][0]*tmp11[0]-z[5][1]*tmp11[1]);
          tab[i+6*stg_first].imag(tmp30[1]+z[6][0]*tmp2[1]+z[6][1]*tmp2[0]+z[1][0]*tmp3[1]+z[1][1]*tmp3[0]+z[7][0]*tmp4[1]+z[7][1]*tmp4[0]+z[2][0]*tmp5[1]+z[2][1]*tmp5[0]+z[8][0]*tmp6[1]+z[8][1]*tmp6[0]+
                                    z[3][0]*tmp7[1]+z[3][1]*tmp7[0]+z[9][0]*tmp8[1]+z[9][1]*tmp8[0]+z[4][0]*tmp9[1]+z[4][1]*tmp9[0]+z[10][0]*tmp10[1]+z[10][1]*tmp10[0]+z[5][0]*tmp11[1]+z[5][1]*tmp11[0]);

          tab[i+7*stg_first].real(tmp30[0]+z[7][0]*tmp2[0]-z[7][1]*tmp2[1]+z[3][0]*tmp3[0]-z[3][1]*tmp3[1]+z[10][0]*tmp4[0]-z[10][1]*tmp4[1]+z[6][0]*tmp5[0]-z[6][1]*tmp5[1]+z[2][0]*tmp6[0]-z[2][1]*tmp6[1]+
                                    z[9][0]*tmp7[0]-z[9][1]*tmp7[1]+z[5][0]*tmp8[0]-z[5][1]*tmp8[1]+z[1][0]*tmp9[0]-z[1][1]*tmp9[1]+z[8][0]*tmp10[0]-z[8][1]*tmp10[1]+z[4][0]*tmp11[0]-z[4][1]*tmp11[1]);
          tab[i+7*stg_first].imag(tmp30[1]+z[7][0]*tmp2[1]+z[7][1]*tmp2[0]+z[3][0]*tmp3[1]+z[3][1]*tmp3[0]+z[10][0]*tmp4[1]+z[10][1]*tmp4[0]+z[6][0]*tmp5[1]+z[6][1]*tmp5[0]+z[2][0]*tmp6[1]+z[2][1]*tmp6[0]+
                                    z[9][0]*tmp7[1]+z[9][1]*tmp7[0]+z[5][0]*tmp8[1]+z[5][1]*tmp8[0]+z[1][0]*tmp9[1]+z[1][1]*tmp9[0]+z[8][0]*tmp10[1]+z[8][1]*tmp10[0]+z[4][0]*tmp11[1]+z[4][1]*tmp11[0]);

          tab[i+8*stg_first].real(tmp30[0]+z[8][0]*tmp2[0]-z[8][1]*tmp2[1]+z[5][0]*tmp3[0]-z[5][1]*tmp3[1]+z[2][0]*tmp4[0]-z[2][1]*tmp4[1]+z[10][0]*tmp5[0]-z[10][1]*tmp5[1]+z[7][0]*tmp6[0]-z[7][1]*tmp6[1]+
                                   z[4][0]*tmp7[0]-z[4][1]*tmp7[1]+z[1][0]*tmp8[0]-z[1][1]*tmp8[1]+z[9][0]*tmp9[0]-z[9][1]*tmp9[1]+z[6][0]*tmp10[0]-z[6][1]*tmp10[1]+z[3][0]*tmp11[0]-z[3][1]*tmp11[1]);
          tab[i+8*stg_first].imag(tmp30[1]+z[8][0]*tmp2[1]+z[8][1]*tmp2[0]+z[5][0]*tmp3[1]+z[5][1]*tmp3[0]+z[2][0]*tmp4[1]+z[2][1]*tmp4[0]+z[10][0]*tmp5[1]+z[10][1]*tmp5[0]+z[7][0]*tmp6[1]+z[7][1]*tmp6[0]+
                                   z[4][0]*tmp7[1]+z[4][1]*tmp7[0]+z[1][0]*tmp8[1]+z[1][1]*tmp8[0]+z[9][0]*tmp9[1]+z[9][1]*tmp9[0]+z[6][0]*tmp10[1]+z[6][1]*tmp10[0]+z[3][0]*tmp11[1]+z[3][1]*tmp11[0]);

          tab[i+9*stg_first].real(tmp30[0]+z[9][0]*tmp2[0]-z[9][1]*tmp2[1]+z[7][0]*tmp3[0]-z[7][1]*tmp3[1]+z[5][0]*tmp4[0]-z[5][1]*tmp4[1]+z[3][0]*tmp5[0]-z[3][1]*tmp5[1]+z[1][0]*tmp6[0]-z[1][1]*tmp6[1]+
                                    z[10][0]*tmp7[0]-z[10][1]*tmp7[1]+z[8][0]*tmp8[0]-z[8][1]*tmp8[1]+z[6][0]*tmp9[0]-z[6][1]*tmp9[1]+z[4][0]*tmp10[0]-z[4][1]*tmp10[1]+z[2][0]*tmp11[0]-z[2][1]*tmp11[1]);
          tab[i+9*stg_first].imag(tmp30[1]+z[9][0]*tmp2[1]+z[9][1]*tmp2[0]+z[7][0]*tmp3[1]+z[7][1]*tmp3[0]+z[5][0]*tmp4[1]+z[5][1]*tmp4[0]+z[3][0]*tmp5[1]+z[3][1]*tmp5[0]+z[1][0]*tmp6[1]+z[1][1]*tmp6[0]+
                                    z[10][0]*tmp7[1]+z[10][1]*tmp7[0]+z[8][0]*tmp8[1]+z[8][1]*tmp8[0]+z[6][0]*tmp9[1]+z[6][1]*tmp9[0]+z[4][0]*tmp10[1]+z[4][1]*tmp10[0]+z[2][0]*tmp11[1]+z[2][1]*tmp11[0]);

          tab[i+10*stg_first].real(tmp30[0]+z[10][0]*tmp2[0]-z[10][1]*tmp2[1]+z[9][0]*tmp3[0]-z[9][1]*tmp3[1]+z[8][0]*tmp4[0]-z[8][1]*tmp4[1]+z[7][0]*tmp5[0]-z[7][1]*tmp5[1]+z[6][0]*tmp6[0]-z[6][1]*tmp6[1]+
                                    z[5][0]*tmp7[0]-z[5][1]*tmp7[1]+z[4][0]*tmp8[0]-z[4][1]*tmp8[1]+z[3][0]*tmp9[0]-z[3][1]*tmp9[1]+z[2][0]*tmp10[0]-z[2][1]*tmp10[1]+z[1][0]*tmp11[0]-z[1][1]*tmp11[1]);
          tab[i+10*stg_first].imag(tmp30[1]+z[10][0]*tmp2[1]+z[10][1]*tmp2[0]+z[9][0]*tmp3[1]+z[9][1]*tmp3[0]+z[8][0]*tmp4[1]+z[8][1]*tmp4[0]+z[7][0]*tmp5[1]+z[7][1]*tmp5[0]+z[6][0]*tmp6[1]+z[6][1]*tmp6[0]+
                                    z[5][0]*tmp7[1]+z[5][1]*tmp7[0]+z[4][0]*tmp8[1]+z[4][1]*tmp8[0]+z[3][0]*tmp9[1]+z[3][1]*tmp9[0]+z[2][0]*tmp10[1]+z[2][1]*tmp10[0]+z[1][0]*tmp11[1]+z[1][1]*tmp11[0]);
        }
    }
///////////////////////////////////////////////////


//this is new in that method and fi:

//when you want to have equal results that are in false modificator in normal FFT then change this:
/*
 fun_fourier_transform_FFT_radix_4_N_256_official
{
    for(int j=0;j<N;j++)
    {
      tab[j].real() =tab[j].real()*2/N;
      tab[j].imag() =tab[j].imag()*2/N;
    }
}
//and:

fun_inverse_fourier_transform_FFT_radix_4_N_256_official
{
    for(int j=0;j<N;j++)
    {
      tab[j].real() =tab[j].real()*0.5;
      tab[j].imag() =tab[j].imag()*0.5;
    }
}

//for official modificator that is only in inverse FFT:

 fun_fourier_transform_FFT_radix_4_N_256_official
{

}
fun_inverse_fourier_transform_FFT_radix_4_N_256_official
{
    for(int i=0;i<N;i++)
    {
        tablica1[0][i]=tablica1[0][i]*1/(float)N;
        tablica1[1][i]=tablica1[1][i]*1/(float)N;
    }
}

//haven't try it with other function that cos(x)+jsin(x)=sin(x+pi/2)+jsin(x) that is mirror inverse
*/


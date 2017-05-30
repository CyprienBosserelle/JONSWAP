
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <cmath>
#include <math.h>
#include <fstream>
#include <algorithm>

#define pi 3.14159265

int nf,nt;
double Hs,Tp,Dir,Spread;
double alpha,beta,gamma;
double fp,fnyq,dfj,dang,mainang;
double *f,*x,*y,*theta,*Dd;
double maxy,sumy,sumDd;
char outfile[256];

template <class T> const T& max (const T& a, const T& b) {
  return (a<b)?b:a;     // or: return comp(a,b)?b:a; for version (2)
}

template <class T> const T& min (const T& a, const T& b) {
  return !(b<a)?a:b;     // or: return !comp(b,a)?a:b; for version (2)
}

template <class T> const T& round(const T& a)
  {
  return floor( a + 0.5 );
  }

float mod360(float x)
{
	float div=x/360.0f;
	float mod=(div-floor(div))*360;
	if (mod<0.0f)
	{
		mod=360+mod;
	}

	return mod;
}


int main(int argc, char **argv)
{
	
	//Read instruction file
	char opfile[]="jonswap.dat";
	FILE * fop;
	fop=fopen(opfile,"r");
		fscanf(fop,"%*s %s\t%*s",&outfile);
		fscanf(fop,"%lf\t%*s",&Hs);
		fscanf(fop,"%lf\t%*s",&Tp);
		fscanf(fop,"%lf\t%*s",&Dir);
		fscanf(fop,"%lf\t%*s",&Spread);
		fscanf(fop,"%lf\t%*s",&alpha);
		fscanf(fop,"%lf\t%*s",&beta);
		fscanf(fop,"%lf\t%*s",&gamma);
	fclose(fop);

	printf("Hs=%f\tTp=%f\tDir=%f\tSpread=%f\talpha=%f\tbeta=%f\tgamma=%f\n",Hs,Tp,Dir,Spread,alpha,beta,gamma);

	//Define frequency vector
	fp=1/Tp;
	dfj=1/400.0;
	nf=111;

	f= (double *)malloc(nf*sizeof(double ));
	x= (double *)malloc(nf*sizeof(double ));
	y= (double *)malloc(nf*sizeof(double ));


	//Define the frequency array
	for (int n=16; n<=126; n++)
	{
		f[n-16]=n*dfj;
		x[n-16]=n*dfj/fp;
	}



	//Calculate baseline JONSWAP spectra
	for (int n=0; n<nf; n++)
	{
		double sigma=0.09;
		if (x[n]<1.0f)
		{
			sigma=0.07;
		}

		double fac1=alpha*9.81*9.81*pow(x[n],-5);
		double fac2=exp(-1*beta*pow(x[n],-4));
		double fac3=pow(gamma,exp(-1*pow(x[n]-1,2)/(2*pow(sigma,2))));

		y[n]=fac1*fac2*fac3;

	}

	//scale the baseline spectra to one
	maxy=0.0;
	for (int n=0; n<nf; n++)
	{
		maxy=max(maxy,y[n]);
	}
	sumy=0.0;
	for (int n=0; n<nf; n++)
	{
		y[n]=y[n]/maxy;
		sumy=sumy+y[n];
	}

	//Scale the spectra to Hs
	for (int n=0; n<nf; n++)
	{
		y[n]=y[n]*pow(Hs/(4.0*sqrt(sumy*dfj)),2);
	}

	//Define the direction spectrum
	nt=360;
	dang=1;
	theta= (double *)malloc(nt*sizeof(double));
	Dd= (double *)malloc(nt*sizeof(double));
	for (int t=0; t<360; t++)
	{
		theta[t]=t*pi/180.0;
	}
	mainang=mod360(Dir)*pi/180;

	sumDd=0.0;
	for (int t=0; t<360; t++)
	{
		Dd[t]=pow(cos((theta[t]-mainang)/2.0),2.0*Spread);
		sumDd=sumDd+Dd[t];
	}
	//Scale to have a surface unity
	for (int t=0; t<360; t++)
	{
		Dd[t]=Dd[t]/(sumDd*dang);
	}

	FILE * ofile;
	ofile= fopen (outfile,"w");

	for (int i=0; i<nt; i++)
	{
		for (int ii=0; ii<nf; ii++)
		{
			fprintf(ofile,"%f\t%f\t%lf\n",f[ii],theta[i]*180/pi,y[ii]*Dd[i]);
		}
		
	}
	fclose (ofile);

	return 0;
}

#include <stdio.h> 
#include <math.h>
#include <stdlib.h> 

#define TINY 1.0e-20;
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

#define Dim 3              
#define La 200
#define Lb 200
#define XMAX 80.0
#define D 60
#include <iostream>
#include <limits> // 用于包含 std::numeric_limits 宏

void rk2(double h,double x[],double fvec[],double xmin[],double xmax[],long *point);
void force(double xp[],double fvec[]);
double ran1(long *idum);
double gasdev(long *idum); 
double findMin(double arr[], int size);

FILE *fp1,*fp2;
double a1;

main()
{  
    long r=1,*point;point=&r;
	char st1[20],st2[20];
	int i,m,k,index1,index2,A,B,filenum,ii,jj;
	int p[La][Lb]={0};
	double xc[2],xi[2],xf[2]; 
	double f_cyc[La][Lb]={0},f_cdk1[La][Lb]={0},fvec[Dim]={0};
	double iter=1.0,Tnf,tau,tau0;
	double x_min[Dim]={0,0,0},x_max[Dim]={150,150,150};;  
	


	char a_1[20]="S3";
	double Tn=1e10;//8
	double Tni=2e7;//2e6 
	double n_write=10000;
	double h=0.001,x[Dim]={1,1,1};;//53,34,25
	index1 = 1;
	index2 = 2;
		
	Tnf = Tn+Tni;
	 
//for(filenum=1;filenum<70;filenum++)
	{
		
	// a1=1+ (filenum-1.0)*1;
	 //a1 =11;
	   
	   
	  printf("filenum:%d	%f\n",filenum,a1);
    //撒初始值 
    
       /* for(ii=0;ii<30;ii++)
		{
		  for(jj=0;jj<30;jj++)
	     	{	 	
		      x[index1-1]=0+ii*0.1;
		      x[index2-1]=0+jj*0.1; 
            	printf("%d\n",k++);
            	printf("%f	%f\n",x[index1-1],x[index2-1]);
        */
        
   
       for(A=0;A<La;A++)                         
        {
	        for(B=0;B<Lb;B++)
	        {
		        p[A][B]= 0.0;
				f_cyc[A][B] = 0.0;
				f_cdk1[A][B] = 0.0;                  
			}
		}
   /*************** sprintf(st1,"bbistable%d%d_%s_%0.5fXm%0.1f.txt",index1,index2,a_1,a1,XMAX);            	                      
	if((fp1=fopen(st1,"w+"))==NULL)
	{
		printf("Cannot open file. \n");exit(0);
	}	
	tau0 = 0; 
***********/

	for(iter=1.0;iter<=Tnf;iter=iter+1)                                      
	{
		rk2(h,x,fvec,x_min,x_max,point); 
			if(fmod(iter,10000000)==0)  //function of fmod is remainder 
	           printf("%e	%f\n",iter,iter*h);  
       /************         
		if(iter>=Tni&iter<=Tni+n_write)    
		{
			fprintf(fp1,"%0.0f	%f",iter-Tni,tau0);  tau0=tau0+h;  
			for (i=0;i<Dim;i++)
			{
				fprintf(fp1,"   %f",x[i]);  
			}
			fprintf(fp1,"\n");                                         
		}
	       	************/
		xc[0]=x[index1-1]; xi[0]=x_min[index1-1]; xf[0]=x_max[index1-1];
		xc[1]=x[index2-1]; xi[1]=x_min[index2-1]; xf[1]=x_max[index2-1]; 
		
		
		

		if(iter>=Tni)               
		{
			A = (int)((xc[0]-xi[0])*La/(xf[0]-xi[0])); 
			B = (int)((xc[1]-xi[1])*Lb/(xf[1]-xi[1]));
			p[A][B] = p[A][B] + 1;
			f_cyc[A][B] = f_cyc[A][B] + fvec[index1-1];    
			f_cdk1[A][B] = f_cdk1[A][B] + fvec[index2-1];
		}	
	}
		sprintf(st2,"pp%d%d.txt",index1,index2);
		
		//	sprintf(st2,"pp%d%d.txt",index1,index2);
		        
		if ((fp2=fopen(st2,"w+"))==NULL)
		{
			printf("Cannot open file. \n");
			exit(0);
		}

		for(A=0;A<La;A++)
		{
			for(B=0;B<La;B++)
			{
					if(p[A][B]==0) 
					fprintf(fp2,"%d	%d	0	0.0	0.0\n",A,B);
					else 
					fprintf(fp2,"%d	%d	%d	%e	%e\n",A,B,p[A][B],f_cyc[A][B]/p[A][B],f_cdk1[A][B]/p[A][B]); 
				
			}         
		}
	fclose(fp1);
	fclose(fp2);
	

				  // 	}
			//	}//
 	}//
}
void rk2(double h,double x[],double q[],double xmin[],double xmax[],long *point)
{
	int i;
	double gauss[i],x1[Dim],q1[Dim];
		force(x,q);      //
		for(i=0;i<Dim;i++)
		{    
        	gauss[i]=gasdev(point);//
			x1[i]=x[i]+h*q[i]+sqrt(h)*gauss[i]*sqrt(2*D);//用heun方法的时候采用
		}
		force(x1,q1);   //获得下一个时刻力的表达式的值,用heun方法的时候采用 
		for(i=0;i<Dim;i++)
		{	
			//x[i]=x[i]+q[i]*h+sqrt(h)*sqrt(2*D)*gauss[i];//Euler method
			x[i]=x[i]+0.5*h*(q[i]+q1[i])+sqrt(h)*gauss[i]*sqrt(2*D);// Heun method
			if(x[i]<xmin[i]) x[i]=2*xmin[i]-x[i];//reflecting boundary condition
			else if(x[i]>xmax[i]) x[i]=2*xmax[i]-x[i]; 
			else x[i]=x[i];
		}	

}
void force(double xp[],double fvec[])        
{		  

double Kji[3][3] = {{1,0.6,0.3},
	                           {0.3,1,0.6},
	                           {0.6,0.3,1}}; //// //*****/ ///  Fig2A 

//double Cji[3][3] = {{0.07,0.04,0.04},
	          //                 {0.08,0.10,0.08},
	          //                 {0.10,0.10,0.14}};
	    double Cji[3][3] = {{0.04,0.044,0.07},
	                          {0.10,0.08,0.08},
	                         {0.10,0.14,0.10}};

  //double Cji[3][3] = {{a1,0.07,0.04},
	//                    {0.08,0.08,0.10},
	  //                     {0.14,0.10,0.10}};///Fig2B 
	                           
	      //  double Kji[3][3] = {{1,0.9,0.3},
	        //                   {0.3,1,0.9},
	       //                    {0.9,0.3,1}};  
	                           
	                           
        //	double Cji[3][3] = {{0.0585,0.07,0.04},
	      //                     {0.08,0.08,0.10},
	      //                     {0.14,0.10,0.10}};       ////Fig3A 
	        double Pji[3][3]={0};
	        double d=1;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
	 		double m1=0.25/d;                                                                                                                                                                                                                                                                                                                                                                       
	 		double m2=0.25/d;                                                                                                                                                                                                                                                                                                                                                                          
	 		double m3=0.25/d;                                                                                                                                                                                                                                                                                                                                                                       
		//	double D1=0.25/d;
 			double S1=6;                                                                                                                                                                                                                                                                                                                                                                               
 			double S2=10;                                                                                                                                                                                                                                                                                                                                                                              
 			double S3=13;
					 			
			double	r1=1/d;
			double	r2=1/d;
			double	r3=1/d;
			double R1,R2,R3;
		//	double P11,P12,P13,P21,P22,P23,P31,P32,P33;
			
			double rowMin[3]; // 用于存储每行的最小值
			int rows = 3;
    		int cols = 3;
			double minVals[cols]; // 用于存储每列的最小值
			
			R1= S1-Cji[0][0]*xp[0]-Cji[0][1]*xp[1]-Cji[0][2]*xp[2];
			R2= S2-Cji[1][0]*xp[0]-Cji[1][1]*xp[1]-Cji[1][2]*xp[2];
			R3= S3-Cji[2][0]*xp[0]-Cji[2][1]*xp[1]-Cji[2][2]*xp[2];	
			
		if (R1<=0)
            R1=0;
		if (R2<=0)
		    R2=0;
		if (R3<=0)
		    R3=0;

			
			Pji[0][0]=r1*R1/(Kji[0][0]+R1);
			Pji[0][1]=r2*R1/(Kji[0][1]+R1);
			Pji[0][2]=r3*R1/(Kji[0][2]+R1);
			
			Pji[1][0]=r1*R2/(Kji[1][0]+R2);
			Pji[1][1]=r2*R2/(Kji[1][1]+R2);
			Pji[1][2]=r3*R2/(Kji[1][2]+R2);
			
			Pji[2][0]=r1*R3/(Kji[2][0]+R3);
			Pji[2][1]=r2*R3/(Kji[2][1]+R3);
			Pji[2][2]=r3*R3/(Kji[2][2]+R3);
			
 // 寻找每列的最小值并输出
    for (int j = 0; j < cols; ++j) {
        double col[3];
        for (int i = 0; i < rows; ++i) {
            col[i] = Pji[i][j];
        }
        
        double minVal = findMin(col, rows);
        minVals[j] = minVal;
        //std::cout << "Minimum value in column " << j << ": " << minVal << std::endl;
    }
    
			fvec[0] = xp[0]*(minVals[0]-m1);                                         
			fvec[1] = xp[1]*(minVals[1]-m2);          
			fvec[2] = xp[2]*(minVals[2]-m3);
}

// 寻找数组中的最小值
double findMin(double arr[], int size) {
    double minVal = std::numeric_limits<double>::max();
    for (int i = 0; i < size; ++i) {
        if (arr[i] < minVal) {
            minVal = arr[i];
        }
    }
    return minVal;
}

double gasdev(long *idum) 
{
	 
	static int iset=0;
	static double gset;
	double fac,rsq,v1,v2;
	if (*idum < 0) iset=0;
	if  (iset == 0) 
	{
		do {
			v1=2.0*ran1(idum)-1.0; 
			v2=2.0*ran1(idum)-1.0;
			rsq=v1*v1+v2*v2;
		} 
		while (rsq >= 1.0 || rsq == 0.0);
		fac=sqrt(-2.0*log(rsq)/rsq);
		gset=v1*fac;
		iset=1;
		return v2*fac;
	} 
	else 
	{
		iset=0;
		return gset;
	}
}
double ran1(long *idum)   
{
	int j;
	long k;
	static long iy=0;
	static long iv[NTAB];
	double temp;
	if (*idum <= 0 || !iy) 
	{
		if (-(*idum) < 1) 
		*idum=1;
		else 
		*idum = -(*idum);
		for (j=NTAB+7;j>=0;j--) 
		{
			k=(*idum)/IQ;
			*idum=IA*(*idum-k*IQ)-IR*k;
			if (*idum < 0) 
			*idum += IM;
			if (j < NTAB)
			 iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ;
	*idum=IA*(*idum-k*IQ)-IR*k;
	if (*idum < 0) *idum += IM;
	j=iy/NDIV;
	iy=iv[j];
	iv[j] = *idum;
	if ((temp=AM*iy) > RNMX)
	 return RNMX;
	else 
	 return temp;
}

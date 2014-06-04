/*************************************************************************************************/
/*  complex_COSMOS  Software for steady-state analysis based on Biochemical systems theory	     */
/*	Calculating steady-state metabolite concentrtions,					                         */
/*		        rate constants and kinetic orders in S-system and GMA-system,Å@			         */
/*              logarithmic gains for metabolite concentrations and net fluxes in S-system,	     */
/*              logarithmic gains for metabolite concentrations and local fluxes in GMA-system,	 */
/*		        eigenvalues			                                                             */
/*  The complex variable method is used to differentiate a given function.                       */
/*	The calculated values are given within round-off errors.                                     */
/*	The output files can directly be used as data files for PLAS.                                */
/*	                                                                     (2014, March 14)        */
/*************************************************************************************************/
#include <windows.h>
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include <complex>
	using namespace std;	
#include "./InputFile/equation.h"
#define Step 1.0e-100
#define RADIX 2.0
#define SWAP(g,h) {y=(g);(g)=(h);(h)=y;}
#define SIGN(a,b) ((b)>=0.0 ? fabs(a) : -fabs(a))

void eigenvalue();
void balance(double **a, int n);
void elmhes(double **a, int n);
void hqr(double **a,int n,double *wr,double *wi);
void sort(int n,double *wr,double *wi);
void diff_depend(double *xp,double **df,double *fxp,double **vp,double **svp,int nd1,int ni,int cp,int cm,double sgm);
void prn_format();
void timecourse();
void timecourse2();

void steady();
void loggain();
void ssystem_finalwrite();
void gmasystem_finalwrite();
void formula_value(double *xp,double *fxp,double **vp,double **svp,int nd1,int ni,int cp,int cm);
void fo(int *nd,int *ni,int *cp,int *cm,double *q,int mode);
void fo1(int *nd,int *ni,int *cp,int *cm);
void fo2(double *t_start,double *t_end,double *t_samp,double *sub_div);
void fo3(int *iflagtime);
void spara_read_write(int nd,int na,double *as,double *bs,double **gs,double **hs,int mode);
void gmaparam_read_write(double **ag,double ***gg,int nd,int na,int cp,int cm,int ca,int mode);
void loggain_write(int nd,int na,int cp,int cm,int ca,double **lx,double **lv,double ***lsv,int mode);
void formula(complex<double> *x,complex<double> *fx,complex<double> **v,complex<double> **sv,int nd1,int ni,int cp,int cm);
void gj(double **a,int nd,int na);
void complex_step1(double *xp,double **df,double *fxp,double **vp,double **svp,int nd,int ni,int cp,int cm,double sgm);
void complex_step2(double *xp,double **dv0,double **dv1,double *fxp,double **vp,double **svp,int nd,int ni,int cp,int cm);
void complex_step3(double *xp,double ***dsv,double *fxp,double **vp,double **svp,int nd,int ni,int cp,int cm);
void coe_s(double *xp,double *fxp,double **vp,double **svp,double *as,double *bs,double **gs,double **hs,int nd,int ni,int cp,int cm);
void coe_g(double *xp,double *fxp,double **vp,double **svp,double **ag,double ***gg,int nd,int ni,int cp,int cm);
int* i_vector(int n);
double* d_vector(int n);
double** d_matrix(int m,int n);
double*** d_cube(int l,int m,int n);
complex<double>* C_d_vector(int n);
complex<double>** C_d_matrix(int m,int n);
void free_d_matrix(int m,double **a);
void free_i_matrix(int m,double **a);
void free_d_cube(int l,int m,double ***a);
void free_C_d_matrix(int m,complex<double> **a);
double** ope_matrix(int l,int m,int n,double **a,double **b);
double** rematrix(int n,double **a);
void gj2(double **a,int nd,int na);
void ffpowlaw();									

#define input_prn "./InputFile/print_format.dat"     
#define input_initial "./InputFile/initial.dat"
#define input_rungekutta "./InputFile/integration.dat"     
#define output1 "./OutputFile/ssystem.dat"
#define output2 "./OutputFile/gmasystem.dat"
#define output3 "./OutputFile/plasfile(ssystem).plc"
#define output4 "./OutputFile/plasfile(gmasystem).plc"
#define output5 "./OutputFile/eigenvalue.dat"

const complex<double> iu(0,1);			
                                        
int iflagtime;
int gmaprn,sprn,prns[4],prng[4];    

void main(){
	prn_format();
	fo3(&iflagtime);
	if(iflagtime==1) timecourse();
	if(iflagtime==2) timecourse2();
	steady();       
	loggain();      
	if(sprn==1) ssystem_finalwrite();     
	if(gmaprn==1) gmasystem_finalwrite(); 
	eigenvalue();
	ffpowlaw();					
}

void prn_format()
{
	int i=0;
	FILE *fr;

	if((fr=fopen(input_prn,"r"))==NULL){
		printf("The file 'print_format.dat' cannot be opened\n");
		puts("Program end.");
		puts("Pres any key.");
		getchar();
		exit(1);
	}
	while(fgetc(fr)!='=')
		;
	fscanf(fr,"%d",&sprn);
	while(fgetc(fr)!='=')
		;
	fscanf(fr,"%d",&gmaprn);

	for(i=0;i<4;i++){
		while(fgetc(fr)!='=')	           
			;
		fscanf(fr,"%d",&prns[i]);
	}
	for(i=0;i<4;i++){
		while(fgetc(fr)!='=')	           
			;
		fscanf(fr,"%d",&prng[i]);
	}
	fclose(fr);
}

// ***************************************************
void timecourse()
{
	double *xp,*fxp,**vp,**svp;
	double t_start,t_end,t_samp,sub_div;
	int i,j,iflag=0; 
	int nd,ni,cp,cm,nd1,na1,ca1;
	FILE *fpw;

	fo1(&nd,&ni,&cp,&cm);
	ca1=cp+cm+1;
	nd1=nd+1;
	na1=nd1+ni;
	xp=d_vector(na1);
	fo(&nd,&ni,&cp,&cm,xp,0);          
	fxp=d_vector(nd1);
	vp=d_matrix(nd1,2);
	svp=d_matrix(nd1,ca1);

	fo2(&t_start,&t_end,&t_samp,&sub_div);

	for(i=1;i<nd1;i++) printf("x[%d]=%lf\n",i,xp[i]);

	int  m;
	double s1,s2,s3,s4,s5,s6;
	double *k,*y,*q,*y0;
	double t,z,h;

	k=d_vector(nd1);
	y=d_vector(na1);
	y0=d_vector(nd1);
	q=d_vector(nd1);

	for(i=1;i<na1;i++) y[i]=xp[i];
	
	if((fpw=fopen("timecourse.dat","w"))==NULL){         
		printf("Output file cannot be opened.\n");
		exit(1);
	}

	s1=(2.0-sqrt(2.0))/2.0;       s2=(3.0*sqrt(2.0)-4.0)/2.0;
	s3=2.0-sqrt(2.0);             s4=(2.0+sqrt(2.0))/2.0;
	s5=-(3.0*sqrt(2.0)+4.0)/2.0;  s6=2.0+sqrt(2.0);

	h=t_samp/sub_div;
	z=t_start;                 

	printf("Time      \tRunge-Kutta\tNewton-Raphson\n");
	
	t=t_start;
	printf("%lf",t);
	fprintf(fpw,"%lf",t);
	for(i=1;i<nd1;i++){
		fprintf(fpw,"\t%12.5e",y[i]);
		printf("\t%23.15e",y[i]);
	}
	fprintf(fpw,"\n");
	printf("\n");

	while(fabs((z-t_end)/t_end)>1.0e-13){

		for(i=1;i<nd1;i++) y0[i]=y[i];

		for(m=0;m<sub_div;m++){	
			t=h*(double)m+z; 
			formula_value(y,fxp,vp,svp,nd1,ni,cp,cm); 
            for(j=1;j<nd1;j++){
				k[j]=h*fxp[j];
				y[j]=y[j]+0.5*k[j];
				q[j]=k[j];
			}

			t=0.5*h+h*(double)m+z;
			formula_value(y,fxp,vp,svp,nd1,ni,cp,cm);
			for(j=1;j<nd1;j++){
				k[j]=h*fxp[j];
				y[j]=y[j]+s1*(k[j]-q[j]);
				q[j]=s2*q[j]+s3*k[j];
			}
   
			t=0.5*h+h*(double)m+z;
			formula_value(y,fxp,vp,svp,nd1,ni,cp,cm);
			for(j=1;j<nd1;j++){
				k[j]=h*fxp[j];
				y[j]=y[j]+s4*(k[j]-q[j]);
				q[j]=s5*q[j]+s6*k[j];
			}
     
			t=h*(double)(m+1)+z;
			formula_value(y,fxp,vp,svp,nd1,ni,cp,cm);
			for(j=1;j<nd1;j++){
				k[j]=h*fxp[j];
				y[j]=y[j]+k[j]/6.0-q[j]/3.0;
			}
		}

		fprintf(fpw,"%lf",t);
		printf("%7.3lf",t);

		for(i=1;i<nd1;i++){
			fprintf(fpw,"\t%12.5e",y[i]);
			printf("\t%23.15e",y[i]);
		}
		fprintf(fpw,"\n");
		printf("\n");

		for(i=1;i<nd1;i++){
			if(fabs((y[i]-y0[i])/y0[i])<1.0e-2) //Use the value of less than 1.0e-2, if necessary.
				iflag=1;
			else {
				iflag=0; break;}
		}

		z+=t_samp;

		if(iflag==1) {
			printf("Initial guesses estimated.\nPress any key to start calcuration."); 
			getchar();
			break;
		}
	}
	
	if(iflag!=1){
			printf("The calculation was executed until the time end.\nPress any key to start calcuration."); 
			getchar();
	}

	fclose(fpw);

	if ((fpw=fopen("work4.dat","w"))==NULL){
		printf("Cannot open parameters file curent drive.\n");
		exit(1);
	}
	fprintf(fpw,"%d\n%d\n%d\n%d\n",nd,ni,cp,cm);
	for (i=1;i<na1;i++){
		fprintf(fpw,"%23.15e\n",y[i]);
	}
	fclose(fpw);

	free(k);
	free(q);
	free(y);
	free(y0);
	
	free(xp);
	free(fxp);
	free_d_matrix(nd1,svp);
	free_d_matrix(nd1,vp);
}

void timecourse2()
{
	double *xp,*fxp,**vp,**svp,**df;
	double t_start,t_end,t_samp,sub_div;
	int i,j,iflag=0; 
	int nd,ni,cp,cm,nd1,na1,ca1;
	FILE *fpw;

	fo1(&nd,&ni,&cp,&cm);
	ca1=cp+cm+1;
	nd1=nd+1;
	na1=nd1+ni;
	xp=d_vector(na1);
	fo(&nd,&ni,&cp,&cm,xp,0);          
	fxp=d_vector(nd1);
	vp=d_matrix(nd1,2);
	svp=d_matrix(nd1,ca1);

	fo2(&t_start,&t_end,&t_samp,&sub_div);

	int  k,m,trymax=20;
	double *yn,*y0,*y,*F,**rF;
	double t,z,h;
	double sgm=1,eps=10e-13;

	df=d_matrix(nd1,nd1+1);
	yn=d_vector(na1);
	y0=d_vector(nd1);
	y=d_vector(na1);
	F=d_vector(nd1);
	rF=d_matrix(nd1,nd1+1);


	for(i=1;i<na1;i++) {
		y[i]=xp[i];
		yn[i]=xp[i];
	}
	
	if((fpw=fopen("timecourse.dat","w"))==NULL){         
		printf("Output file cannot be opened.\n");
		exit(1);
	}

	h=t_samp/sub_div;
	z=t_start;                 

	printf("**** Solution by a Semi-implicit method ****\nTime      \tConcentration\n");
	
	t=t_start;
	printf("%lf",t);
	fprintf(fpw,"%lf",t);
	for(i=1;i<nd1;i++){
		fprintf(fpw,"\t%12.5e",y[i]);
		printf("\t%23.15e",y[i]);
	}
	fprintf(fpw,"\n");
	printf("\n");

	while(fabs((z-t_end)/t_end)>1.0e-13){
		for(m=1;m<=sub_div;m++){	
			for(i=1;i<nd1;i++) y0[i]=y[i];
			formula_value(y,fxp,vp,svp,nd1,ni,cp,cm);

			for(i=1;i<nd1;i++) yn[i]=y0[i]+fxp[i]*h;

			t=h*(double)m+z; 
			for (j=0;j<trymax;j++){                    //Start of the Newton-Raphson method
				complex_step1(yn,df,fxp,vp,svp,nd1,ni,cp,cm,sgm);						
				for(i=1;i<nd1;i++){
					for(k=1;k<nd1;k++){
						if(i==k)
							rF[i][i]=1.0-h*df[i][i];      
						else
							rF[i][k]=-h*df[i][k];     
					}
				}
				for(i=1;i<nd1;i++){
					F[i]=yn[i]-y0[i]-h*fxp[i];              
				}
				for(i=1;i<nd1;i++) rF[i][nd1]=F[i];   
								
				gj(rF,nd1,nd1+1);    
						
				for (i=1;i<nd1;i++)	yn[i]=yn[i]-rF[i][nd1];	

				for (i=1;i<nd1;i++){
					sgm=fabs(rF[i][nd1]/y0[i]);
					if (sgm>eps) break;
				}

				if(sgm<eps) break;
			}

			if (sgm>eps) printf("not convergence (trial=%d) \n",j);

			for(i=1;i<nd1;i++){
				y[i]=yn[i];
			}

			for(i=1;i<nd1;i++){
				if(fabs((yn[i]-y0[i])/y0[i])<1.0e-2) //Use the value of less than 1.0e-2, if necessary.
					iflag=1;
				else {
					iflag=0; break;}
			}
		}

		fprintf(fpw,"%lf",t);
		printf("%7.3lf",t);
		for(i=1;i<nd1;i++){
			fprintf(fpw,"\t%12.5e",y[i]);
			printf("\t%23.15e",y[i]);
		}
		fprintf(fpw,"\n");
		printf("\n");

		z+=t_samp;

		if(iflag==1) {
			printf("Initial guesses estimated.\nPress any key to start calcuration."); 
			getchar();
			break;
		}
	}

	if(iflag!=1){
			printf("The calculation was executed until the time end.\nPress any key to start calcuration."); 
			getchar();
	}

	fclose(fpw);

	if ((fpw=fopen("work4.dat","w"))==NULL){
		printf("Cannot open parameters file curent drive.\n");
		exit(1);
	}
	fprintf(fpw,"%d\n%d\n%d\n%d\n",nd,ni,cp,cm);
	for (i=1;i<na1;i++){
		fprintf(fpw,"%23.15e\n",y[i]);
	}
	fclose(fpw);

	free(xp);
	free(fxp);
	free_d_matrix(nd1,vp);
	free_d_matrix(nd1,svp);

	free_d_matrix(nd1,df);        
	free(yn);
	free(y0);
	free(y);
	free(F);
	free_d_matrix(nd1,rF); 
}

void formula_value(double *xp,double *fxp,double **vp,double **svp,int nd1,int ni,int cp,int cm)
{
	int i,j,na1,ca1;
	complex<double> *C_xp,*C_fxp,**C_vp,**C_svp;

	na1=nd1+ni;
	ca1=cp+cm+1;
	C_xp=C_d_vector(na1);
	C_fxp=C_d_vector(nd1);
	C_vp=C_d_matrix(nd1,2);
	C_svp=C_d_matrix(nd1,ca1);  

	for(j=1;j<na1;j++)
		C_xp[j]=xp[j];
		
	for(i=1;i<nd1;i++){
		formula(C_xp,C_fxp,C_vp,C_svp,nd1,ni,cp,cm);
		fxp[i]=real(C_fxp[i]);
		vp[i][0]=real(C_vp[i][0]);
		vp[i][1]=real(C_vp[i][1]);
	}
	
	free(C_xp);
	free(C_fxp);
	free_C_d_matrix(nd1,C_vp);
	free_C_d_matrix(nd1,C_svp);
}

void fo(int *nd,int *ni,int *cp,int *cm,double *x,int mode)
{
	double na1;
	int i;
	FILE *fp,*fpw;

	na1=*nd+*ni+1;
	if(mode==0){
		if ((fp=fopen(input_initial,"r"))==NULL){
			printf("The file 'initial.dat' cannot be opened.\n");
			exit(1);
		}
		fscanf(fp,"%d%d%d%d%d",&iflagtime,nd,ni,cp,cm);
		na1=*nd+*ni+1;
		for (i=1;i<na1;i++){
			fscanf(fp,"%lf",&x[i]);
		}
		fclose(fp);
	}
	if(mode==1){
		if ((fpw=fopen("work0.dat","w"))==NULL){
			printf("Cannot open parameters file curent drive.\n");
			exit(1);
		}
		fprintf(fpw,"%d\n%d\n%d\n%d\n",*nd,*ni,*cp,*cm);
		for (i=1;i<na1;i++){
			fprintf(fpw,"%23.15e\n",x[i]);
		}
		fclose(fpw);
	}
	if(mode==2){
		if ((fp=fopen("work4.dat","r"))==NULL){
			printf("Cannot open parameters file curent drive.\n");
			exit(1);
		}
		fscanf(fp,"%d%d%d%d",nd,ni,cp,cm);
		for (i=1;i<na1;i++){
			fscanf(fp,"%lf",&x[i]);
		}
		fclose(fp);
	}
}

void fo1(int *nd,int *ni,int *cp,int *cm)
{	
	FILE *fp;

	if ((fp=fopen(input_initial,"r"))==NULL){
		printf("The file 'initial.dat' cannot be opened.\n");
		exit(1);
	}
	fscanf(fp,"%d%d%d%d%d",&iflagtime,nd,ni,cp,cm);
	fclose(fp);
}

void fo2(double *t_start,double *t_end,double *t_samp,double *sub_div)
{	
	FILE *fp;
	
	if ((fp=fopen(input_rungekutta,"r"))==NULL){
		printf("The file 'rungekutta.dat' cannot be opened.\n");
		exit(1);
	}
	fscanf(fp,"%lf%lf%lf%lf%lf",t_start,t_end,t_samp,sub_div);  
	fclose(fp);
}

void fo3(int *iflagtime)    
{	
	FILE *fp;
	
	if ((fp=fopen(input_initial,"r"))==NULL){
		printf("The file 'initial.dat' cannot be opened.\n");
		exit(1);
	}
	fscanf(fp,"%d",iflagtime);
	fclose(fp);
}
// *************************************************

void steady()
{
	double *xp,*fxp,**df,**vp,**svp;
	double *as,*bs,**gs,**hs;
	double **ag,***gg;
	double sgm=1,eps=10e-15;
	int i,j; 
	int nd,nf,ni,na,cp,cm,nd1,na1,ca1,trymax=20;

	fo1(&nd,&ni,&cp,&cm);
	na=nd+ni;
	ca1=cp+cm+1;
	nd1=nd+1;
	na1=na+1;
	xp=d_vector(na1);
	
	if(iflagtime==1 || iflagtime==2) 
		fo(&nd,&ni,&cp,&cm,xp,2);
	else
		fo(&nd,&ni,&cp,&cm,xp,0);

	fxp=d_vector(nd1);
	vp=d_matrix(nd1,2);
	svp=d_matrix(nd1,ca1);
	df=d_matrix(nd1,nd1+1);
	nf=nd1;
	printf("nd=%d ni=%d cp=%d cm=%d\n",nd,ni,cp,cm);

	for (j=0;j<trymax;j++){
		complex_step1(xp,df,fxp,vp,svp,nd1,ni,cp,cm,sgm);
		for (i=1;i<nd1;i++){
			df[i][nf]=-fxp[i];
		}
		gj(df,nf,nf+1);				
		for (i=1;i<nf;i++){				
			xp[i]=xp[i]+df[i][nf];	
		}
		for (i=1;i<nf;i++){
			sgm=fabs(df[i][nf]/xp[i]);
			if (sgm>eps) break;
		}
		if(sgm<eps) break;
	}
	if (sgm>eps) printf("not convergence \n");
	printf("\n");
	printf("trial number %d\n\n",j);

	fo(&nd,&ni,&cp,&cm,xp,1);  
	free_d_matrix(nd1,df);
/*=====================================================*/
/*	s-system	*/
	as=d_vector(nd1);
	gs=d_matrix(nd1,na1);
	bs=d_vector(nd1);
	hs=d_matrix(nd1,na1);
	coe_s(xp,fxp,vp,svp,as,bs,gs,hs,nd1,ni,cp,cm);	
	spara_read_write(nd1,na1,as,bs,gs,hs,1);
	free(as);
	free(bs);
	free_d_matrix(nd1,gs);
	free_d_matrix(nd1,hs);

/*	gma	*/
	if(gmaprn==1){    
		ag=d_matrix(nd1,ca1);
		gg=d_cube(nd1,na1,ca1);
		coe_g(xp,fxp,vp,svp,ag,gg,nd1,ni,cp,cm);
		gmaparam_read_write(ag,gg,nd1,na1,cp,cm,ca1,1);
		free_d_matrix(nd1,ag);
		free_d_cube(nd1,na1,gg);
	}
/*=====================================================*/
	free(xp);
	free(fxp);
	free_d_matrix(nd1,svp);
	free_d_matrix(nd1,vp);
}

void loggain()
{
	double **a,**b,*as,*bs,**gs,**hs,**ag,***gg;
	double s,**lx,**lv,***lsv;
	int i,j,k,l;
	int nd,ni,cp,cm,nd1,na1,ca1;

	fo1(&nd,&ni,&cp,&cm);
	nd1=nd+1;
	na1=nd1+ni;
	ca1=cp+cm+1;
	as=d_vector(nd1);
	gs=d_matrix(nd1,na1);
	bs=d_vector(nd1);
	hs=d_matrix(nd1,na1);
	spara_read_write(nd1,na1,as,bs,gs,hs,0);
    a=d_matrix(nd1,na1);
    b=d_matrix(nd1,nd1);
    for (i=1;i<nd1;i++){
    	for (j=1;j<nd1;j++)
       		a[i][j]=gs[i][j]-hs[i][j];
	}
    b=rematrix(nd1,a);
    for (i=1;i<nd1;i++){
    	for (j=1;j<nd1;j++)
			a[i][j]=0.0;
    	for (j=nd1;j<na1;j++)
       		a[i][j]=-gs[i][j]+hs[i][j];
	}
	lx=d_matrix(nd1,na1);
	lx=ope_matrix(nd1,nd1,na1,b,a);
	lsv=d_cube(nd1,na1,ca1);
	lv=d_matrix(nd1,na1);

	loggain_write(nd1,na1,cp,cm,ca1,lx,lv,lsv,0);

	for (i=1;i<nd1;i++){
		for (j=1;j<na1;j++){
			s=0.0;
			for (k=1;k<nd1;k++)
				s=s+gs[i][k]*lx[k][j];
		lv[i][j]=gs[i][j]+s;
		}
	}
	loggain_write(nd1,na1,cp,cm,ca1,lx,lv,lsv,1);
	ag=d_matrix(nd1,ca1);
	gg=d_cube(nd1,na1,ca1);
	gmaparam_read_write(ag,gg,nd1,na1,cp,cm,ca1,0);

	for (i=1;i<nd1;i++){
		for (l=1;l<ca1;l++){
			for (k=1;k<na1;k++){
				s=0.0;
				for (j=1;j<nd1;j++)
					s=s+gg[i][j][l]*lx[j][k];
				lsv[i][k][l]=s+gg[i][k][l];
			}
		}
	}
	loggain_write(nd1,na1,cp,cm,ca1,lx,lv,lsv,2);
	free_d_cube(nd1,na1,lsv);
    free_d_matrix(nd1,lv);
	free_d_matrix(nd1,lx);
	free_d_matrix(nd1,ag);
	free_d_cube(nd1,na1,gg);
}

void coe_g(double *xp,double *fxp,double **vp,double **svp,double **ag,double ***gg,int nd1,int ni,int cp,int cm)
{
	int i,j,k,na1,ca1;
	double ***dsv,s;
	na1=nd1+ni;
	ca1=cp+cm+1;
	dsv=d_cube(nd1,na1,ca1);
	complex_step3(xp,dsv,fxp,vp,svp,nd1,ni,cp,cm);
	for (i=1;i<nd1;i++){
		for (k=1;k<ca1;k++){
			for (j=1;j<na1;j++){
			    if(svp[i][k]<10e-16) svp[i][k]=10e-16;
				gg[i][j][k]=dsv[i][j][k]*xp[j]/svp[i][k];
			}
		}	
	}
	for (i=1;i<nd1;i++){
		for (k=1;k<ca1;k++){
			s=1.0;
			for (j=1;j<na1;j++)
				s=s*pow(xp[j],gg[i][j][k]);
			ag[i][k]=svp[i][k]/s;
		}
	}
	free_d_cube(nd1,na1,dsv);
}

void coe_s(double *xp,double *fxp,double **vp,double **svp,double *as,double *bs,double **gs,double **hs,int nd1,int ni,int cp,int cm)
{
	int i,j,na1;    
	double **dv0,**dv1,s,t;
	na1=nd1+ni;
	dv0=d_matrix(nd1,na1);
	dv1=d_matrix(nd1,na1);
	complex_step2(xp,dv0,dv1,fxp,vp,svp,nd1,ni,cp,cm);
	for (i=1;i<nd1;i++){
		for (j=1;j<na1;j++){
			gs[i][j]=dv0[i][j]*xp[j]/vp[i][0];
			hs[i][j]=dv1[i][j]*xp[j]/vp[i][1];
		}
	}
	for (i=1;i<nd1;i++){
		s=1.0;
		t=1.0;
		for (j=1;j<na1;j++){
			s=s*pow(xp[j],gs[i][j]);
			t=t*pow(xp[j],hs[i][j]);
		}
		as[i]=vp[i][0]/s;
		bs[i]=vp[i][1]/t;
	}
	free_d_matrix(nd1,dv0);
	free_d_matrix(nd1,dv1);
}

void complex_step1(double *xp,double **df,double *fxp,double **vp,double **svp,int nd1,int ni,int cp,int cm,double sgm)
{
	int i,j,na1,ca1;
	double h=Step;
	complex<double> *C_xp,*C_fxp,**C_vp,**C_svp;

	na1=nd1+ni;
	ca1=cp+cm+1;
	C_xp=C_d_vector(na1);
	C_fxp=C_d_vector(nd1);
	C_vp=C_d_matrix(nd1,2);
	C_svp=C_d_matrix(nd1,ca1);  

	for(j=1;j<na1;j++)
		C_xp[j]=xp[j];
		
	for(i=1;i<nd1;i++){
		for(j=1;j<nd1;j++){	
			C_xp[j]=xp[j]+iu*h;
			formula(C_xp,C_fxp,C_vp,C_svp,nd1,ni,cp,cm);
			df[i][j]=imag(C_fxp[i])/h;
			C_xp[j]=xp[j];
		}
		formula(C_xp,C_fxp,C_vp,C_svp,nd1,ni,cp,cm);
		fxp[i]=real(C_fxp[i]);
		vp[i][0]=real(C_vp[i][0]);
		vp[i][1]=real(C_vp[i][1]);
	}
	
	free(C_xp);
	free(C_fxp);
	free_C_d_matrix(nd1,C_vp);
	free_C_d_matrix(nd1,C_svp);
}

void complex_step2(double *xp,double **dv0,double **dv1,double *fxp,double **vp,double **svp,int nd1,int ni,int cp,int cm)
{
	int i,j,na1,ca1;
	double h=Step;
	complex<double> *C_xp,*C_fxp,**C_vp,**C_svp;

	na1=nd1+ni;
	ca1=cp+cm+1;
	C_xp=C_d_vector(na1);
	C_fxp=C_d_vector(nd1);
	C_vp=C_d_matrix(nd1,2);
	C_svp=C_d_matrix(nd1,ca1);
	
	for(j=1;j<na1;j++)
		C_xp[j]=xp[j];
	for (i=1;i<nd1;i++){
		for (j=1;j<na1;j++){	
			C_xp[j]=xp[j]+iu*h;	
			formula(C_xp,C_fxp,C_vp,C_svp,nd1,ni,cp,cm);
			dv0[i][j]=imag(C_vp[i][0])/h;
			dv1[i][j]=imag(C_vp[i][1])/h;
			C_xp[j]=xp[j];
		}
		formula(C_xp,C_fxp,C_vp,C_svp,nd1,ni,cp,cm);
		fxp[i]=real(C_fxp[i]);
		vp[i][0]=real(C_vp[i][0]);
		vp[i][1]=real(C_vp[i][1]);
	}
	free(C_xp);
	free(C_fxp);
	free_C_d_matrix(nd1,C_vp);
	free_C_d_matrix(nd1,C_svp);
}

void complex_step3(double *xp,double ***dsv,double *fxp,double **vp,double **svp,int nd1,int ni,int cp,int cm)
{
	int i,j,kq,na1,ca1;
	double h=Step;
	complex<double> *C_xp,*C_fxp,**C_vp,**C_svp;

	na1=nd1+ni;
	ca1=cp+cm+1;
	C_xp=C_d_vector(na1);
	C_fxp=C_d_vector(nd1);
	C_vp=C_d_matrix(nd1,2);
	C_svp=C_d_matrix(nd1,ca1);  

	for(j=1;j<na1;j++)
		C_xp[j]=xp[j];

	for (i=1;i<nd1;i++){
		for (j=1;j<na1;j++){	
			C_xp[j]=xp[j]+iu*h;	
			formula(C_xp,C_fxp,C_vp,C_svp,nd1,ni,cp,cm);
			for (kq=1;kq<ca1;kq++){
				dsv[i][j][kq]=imag(C_svp[i][kq])/h;
				svp[i][kq]=real(C_svp[i][kq]);
			}
			C_xp[j]=xp[j];
		}
		formula(C_xp,C_fxp,C_vp,C_svp,nd1,ni,cp,cm);
		fxp[i]=real(C_fxp[i]);
		vp[i][0]=real(C_vp[i][0]);
		vp[i][1]=real(C_vp[i][1]);
	}
	free(C_xp);
	free(C_fxp);
	free_C_d_matrix(nd1,C_vp);
	free_C_d_matrix(nd1,C_svp);
}

void loggain_write(int nd1,int na1,int cp,int cm,int ca1,double **lx,double **lv,double ***lsv,int mode)
{	
	int i,j,k;
	FILE *fpw;

	if (mode==0){
		if ((fpw=fopen("work3.dat","w"))==NULL){
			printf("Cannot open parameters file curent drive.\n");
			exit(1);
		}
		for (i=1;i<nd1;i++){
			for (j=nd1;j<na1;j++)
				fprintf(fpw,"%23.15e\n",lx[i][j]);
		}
		fprintf(fpw,"\n");
		fclose(fpw);
	}
	if (mode==1){
		if ((fpw=fopen("work3.dat","a"))==NULL){
			printf("Cannot open parameters file curent drive.\n");
			exit(1);
		}
		for (i=1;i<nd1;i++){
			for (j=nd1;j<na1;j++)
				fprintf(fpw,"%23.15e\n",lv[i][j]);
		}
		fclose(fpw);
	}
	if(mode==2){
		if ((fpw=fopen("work3.dat","a"))==NULL){
			printf("Cannot open parameters file curent drive.\n");
			exit(1);
		}
		for (i=1;i<nd1;i++){
			for (k=1;k<ca1;k++){
				for (j=nd1;j<na1;j++)
					fprintf(fpw,"%23.15e\n",lsv[i][j][k]);
			}
		}
		fclose(fpw);
	}
}

void gmaparam_read_write(double **ag,double ***gg,int nd1,int na1,int cp,int cm,int ca1,int mode)
{	
	int i,j,k;
	FILE *fp,*fpw;

	if (mode ==0){
		if ((fp=fopen("work2.dat","r"))==NULL){
			printf("Cannot open parameters file curent drive.\n");
			exit(1);
		}
		for (i=1;i<nd1;i++){
			for (k=1;k<ca1;k++){
				for (j=1;j<na1;j++)
					fscanf(fp,"%lf",&gg[i][j][k]);
			}
		}
		for (i=1;i<nd1;i++){
			for (k=1;k<ca1;k++)
				fscanf(fp,"%lf",&ag[i][k]);
		}
		fclose(fp);
	}
	if(mode ==1){
		if ((fpw=fopen("work2.dat","w"))==NULL){
			printf("Cannot open parameters file curent drive.\n");
			exit(1);
		}
		for (i=1;i<nd1;i++){
			for (k=1;k<ca1;k++){
				for (j=1;j<na1;j++)
					fprintf(fpw,"%23.15e\n",gg[i][j][k]);
			}
		}
		for (i=1;i<nd1;i++){
			for (k=1;k<ca1;k++)
				fprintf(fpw,"%23.15e\n",ag[i][k]);
		}
		fclose(fpw);
	}
}

void spara_read_write(int nd1,int na1,double *as,double *bs,double **gs,double **hs,int mode)
{	
	int i,j;
	FILE *fp,*fpw;

	if (mode ==0){
		if ((fp=fopen("work1.dat","r"))==NULL){
			printf("Cannot open parameters file curent drive.\n");
			exit(1);
		}
		for (i=1;i<nd1;i++)
			fscanf(fp,"%lf",&as[i]);
		for (i=1;i<nd1;i++)
			fscanf(fp,"%lf",&bs[i]);
		for (i=1;i<nd1;i++){
			for (j=1;j<na1;j++)
				fscanf(fp,"%lf",&gs[i][j]);
		}
		for (i=1;i<nd1;i++){
			for (j=1;j<na1;j++)
				fscanf(fp,"%lf",&hs[i][j]);
		}
		fclose(fp);
	}
	if (mode==1){
		if ((fpw=fopen("work1.dat","w"))==NULL){
			printf("Cannot open parameters file curent drive.\n");
			exit(1);
		}
		for (i=1;i<nd1;i++)
			fprintf(fpw,"%23.15e\n",as[i]);
		for (i=1;i<nd1;i++)
			fprintf(fpw,"%23.15e\n",bs[i]);
		for (i=1;i<nd1;i++){
			for (j=1;j<na1;j++)
				fprintf(fpw,"%23.15e\n",gs[i][j]);
		}
		for (i=1;i<nd1;i++){
			for (j=1;j<na1;j++)
				fprintf(fpw,"%23.15e\n",hs[i][j]);
		}
		fclose(fpw);
	}
}

void ssystem_finalwrite()
{
	int nd,ni,cp,cm,ca1,nd1,na1;
	int i,j;    
	double *x;    
	double *as,*bs,**gs,**hs;
	double **lx,**lv;
	FILE *fp,*fpw;

	if ((fp=fopen("work0.dat","r"))==NULL){
		printf("Cannot open parameters file curent drive.\n");
		exit(1);
	}
	fscanf(fp,"%d%d%d%d",&nd,&ni,&cp,&cm);
	nd1=nd+1;
	na1=nd1+ni;
	ca1=cp+cm+1;
	x=d_vector(na1);
	for (i=1;i<na1;i++){
		fscanf(fp,"%lf",&x[i]);
	}
	fclose(fp);

	if ((fpw=fopen(output1,"w"))==NULL){
		printf("Cannot open parameters file curent drive.\n");
		exit(1);
	}
	fprintf(fpw,"nd=%2d ni=%2d \n",nd,ni);
	fprintf(fpw,"cp=%2d cm=%2d \n",cp,cm);
	fprintf(fpw,"\n");

	if(prns[0]==1){
		fprintf(fpw,"************ S-system ***********\n\n",nd);
		printf("\n************ S-system ***********\n\n",nd);
		printf("-------------------------------\n");
		printf(" i      X[i] \n");
		printf("-------------------------------\n");

		for(i=1;i<na1;i++){
			if(i==nd1) fprintf(fpw,"-------------------------------\n");
			if(i==nd1) printf("-------------------------------\n");
			printf("%2d %23.15e\n",i,x[i]);
			fprintf(fpw,"%2d %23.15e\n",i,x[i]);
		}
		fprintf(fpw,"-------------------------------\n");
		printf("-------------------------------\n");
		fprintf(fpw,"\n");
		printf("\n");
	}
	free(x);

	if(prns[1]==1){
		as=d_vector(nd1);
		bs=d_vector(nd1);
		gs=d_matrix(nd1,na1);
		hs=d_matrix(nd1,na1);

		if ((fp=fopen("work1.dat","r"))==NULL){
			printf("Cannot open parameters file curent drive.\n");
			exit(1);
		}
		for (i=1;i<nd1;i++)
			fscanf(fp,"%lf",&as[i]);
		for (i=1;i<nd1;i++)
			fscanf(fp,"%lf",&bs[i]);
		for (i=1;i<nd1;i++){
			for (j=1;j<na1;j++)
				fscanf(fp,"%lf",&gs[i][j]);
		}
		for (i=1;i<nd1;i++){
			for (j=1;j<na1;j++)
				fscanf(fp,"%lf",&hs[i][j]);
		}
		fclose(fp);

		fprintf(fpw,"alpha[%d]\n",nd);
		printf("alph[%d]\n",nd);
		for (i=1;i<nd1;i++){
			fprintf(fpw,"%23.15e ",as[i]);
			printf("%12.5e ",as[i]);
		}
		fprintf(fpw,"\n\n");
		printf("\n\n");
		fprintf(fpw,"beta[%d]\n",nd);
		printf("beta[%d]\n",nd);
		for (i=1;i<nd1;i++){
			fprintf(fpw,"%23.15e ",bs[i]);
			printf("%12.5e ",bs[i]);
		}
		fprintf(fpw,"\n\n");
		printf("\n\n");

		fprintf(fpw,"g[%d][%d]\n",nd,na1-1);
		printf("g[%d][%d]\n",nd,na1-1);
		for (i=1;i<nd1;i++){
			for (j=1;j<na1;j++){
				fprintf(fpw,"%23.15e ",gs[i][j]);
				printf("%12.5e ",gs[i][j]);
			}
			fprintf(fpw,"\n");
			printf("\n");
		}
		fprintf(fpw,"\n");
		printf("\n");

		fprintf(fpw,"h[%d][%d]\n",nd,na1-1);
		printf("h[%d][%d]\n",nd,na1-1);
		for (i=1;i<nd1;i++){
			for (j=1;j<na1;j++){
				fprintf(fpw,"%23.15e ",hs[i][j]);
			printf("%12.5e ",hs[i][j]);
			}
			fprintf(fpw,"\n");
			printf("\n");
		}
		fprintf(fpw,"\n");
		printf("\n");
	}
	free(as);
	free(bs);
	free_d_matrix(nd1,gs);
	free_d_matrix(nd1,hs);

	if(prns[2]==1){
		lx=d_matrix(nd1,na1);
		lv=d_matrix(nd1,na1);
		if ((fp=fopen("work3.dat","r"))==NULL){
			printf("Cannot open parameters file curent drive.\n");
			exit(1);
		}
		for (i=1;i<nd1;i++){
			for (j=nd1;j<na1;j++)
				fscanf(fp,"%lf",&lx[i][j]);
		}
		for (i=1;i<nd1;i++){
			for (j=nd1;j<na1;j++)
				fscanf(fp,"%lf",&lv[i][j]);
		}
		fclose(fp);

		fprintf(fpw,"Lx[%d][%d-%d]\n",nd,nd+1,na1-1);
		printf("Lx[%d][%d-%d]\n",nd,nd+1,na1-1);
		for (i=1;i<nd1;i++){
			for (j=nd1;j<na1;j++){
				fprintf(fpw,"%23.15e ",lx[i][j]);
				printf("%12.5e ",lx[i][j]);
			}
			fprintf(fpw,"\n");
			printf("\n");
		}
		fprintf(fpw,"\n");
		printf("\n");

		fprintf(fpw,"LV[%d][%d-%d]\n",nd,nd+1,na1-1);
		printf("LV[%d][%d-%d]\n",nd,nd+1,na1-1);
		for (i=1;i<nd1;i++){
			for (j=nd1;j<na1;j++){
				fprintf(fpw,"%23.15e ",lv[i][j]);
				printf("%12.5e ",lv[i][j]);
			}
			fprintf(fpw,"\n");
			printf("\n");
		}
		fprintf(fpw,"\n");
		printf("\n");
		free_d_matrix(nd1,lx);
		free_d_matrix(nd1,lv);
	}
	fclose(fpw);
}

void gmasystem_finalwrite()
{
	int nd,ni,cp,cm,nd1,na1,ca1;
	int i,j,k;    
	double *x;    
	double **ag,***gg;
	double **lx,**lv,***lsv;
	FILE *fp,*fpw;

	if ((fp=fopen("work0.dat","r"))==NULL){
		printf("Cannot open parameters file curent drive.\n");
		exit(1);
	}
	fscanf(fp,"%d%d%d%d",&nd,&ni,&cp,&cm);
	nd1=nd+1;
	na1=nd1+ni;
	ca1=cp+cm+1;
	x=d_vector(na1);
	for (i=1;i<na1;i++)
		fscanf(fp,"%lf",&x[i]);
	fclose(fp);

	if ((fpw=fopen(output2,"w"))==NULL){
		printf("Cannot open parameters file curent drive.\n");
		exit(1);
	}
	fprintf(fpw,"nd=%2d ni=%2d \n",nd,ni);
	fprintf(fpw,"cp=%2d cm=%2d \n",cp,cm);
	fprintf(fpw,"\n");

	if(prng[0]==1){
		fprintf(fpw,"\n************ GMA-system ***********\n\n",nd);
		printf("\n************ GMA-system ***********\n\n",nd);
		printf("-------------------------------\n");
		printf(" i      X[i] \n");
		printf("-------------------------------\n");

		for(i=1;i<na1;i++){
			if(i==nd1) fprintf(fp,"-------------------------------\n");
			if(i==nd1) printf("-------------------------------\n");
			printf("%2d %23.15e\n",i,x[i]);
			fprintf(fpw,"%2d %23.15e\n",i,x[i]);
		}
		fprintf(fpw,"-------------------------------\n");
		printf("-------------------------------\n");
		fprintf(fpw,"\n");
		printf("\n");
	}
	free(x);

	if(prng[1]==1){
		ag=d_matrix(nd1,ca1);
		gg=d_cube(nd1,na1,ca1);
		if ((fp=fopen("work2.dat","r"))==NULL){
			printf("Cannot open parameters file curent drive.\n");
			exit(1);
		}
		for (i=1;i<nd1;i++){
			for (k=1;k<ca1;k++){
				for (j=1;j<na1;j++)
					fscanf(fp,"%lf",&gg[i][j][k]);
			}
		}
		for (i=1;i<nd1;i++){
			for (k=1;k<ca1;k++)
				fscanf(fp,"%lf",&ag[i][k]);
		}
		fclose(fp);

		fprintf(fpw,"alpha[%d][%d]\n",nd,ca1-1);
		printf("alpha[%d][%d]\n",nd,ca1-1);
		for (i=1;i<nd1;i++){
			for (k=1;k<ca1;k++){
				fprintf(fpw,"%23.15e ",ag[i][k]);
				printf("%12.5e ",ag[i][k]);
			}
			fprintf(fpw,"\n");
			printf("\n");
		}
		fprintf(fpw,"\n");
		printf("\n");
	
		fprintf(fpw,"g[%d][%d][%d]\n",nd,na1-1,ca1-1);
		printf("g[%d][%d][%d]\n",nd,na1-1,ca1-1);
		for (i=1;i<nd1;i++){
			for (k=1;k<ca1;k++){
				for (j=1;j<na1;j++){
					fprintf(fpw,"%23.15e ",gg[i][j][k]);
					printf("%12.5e ",gg[i][j][k]);
				}
				fprintf(fpw,"\n");
				printf("\n");
			}
			fprintf(fpw,"\n");
			printf("\n");
		}

		free_d_matrix(nd1,ag);
		free_d_cube(nd1,na1,gg);
	}

	if(prng[2]==1){
		lx=d_matrix(nd1,na1);
		lv=d_matrix(nd1,na1);
		lsv=d_cube(nd1,na1,ca1);

		if ((fp=fopen("work3.dat","r"))==NULL){
			printf("Cannot open parameters file curent drive.\n");
			exit(1);
		}
		for (i=1;i<nd1;i++){
			for (j=nd1;j<na1;j++)
				fscanf(fp,"%lf",&lx[i][j]);
		}
		for (i=1;i<nd1;i++){
			for (j=nd1;j<na1;j++)
				fscanf(fp,"%lf",&lv[i][j]);
		}
		for (i=1;i<nd1;i++){
			for (k=1;k<ca1;k++){
				for (j=nd1;j<na1;j++)
					fscanf(fp,"%lf",&lsv[i][j][k]);
			}
		}
		fclose(fp);

		fprintf(fpw,"Lx[%d][%d-%d]\n",nd,nd1,na1-1);
		printf("Lx[%d][%d-%d]\n",nd,nd1,na1-1);
		for (i=1;i<nd1;i++){
			for (j=nd1;j<na1;j++){
				fprintf(fpw,"%23.15e ",lx[i][j]);
				printf("%12.5e ",lx[i][j]);
			}
			fprintf(fpw,"\n");
			printf("\n");
		}
		fprintf(fpw,"\n");
		printf("\n");

		fprintf(fpw,"Lv[%d][%d-%d][%d]\n",nd,nd1,na1-1,ca1-1);
		printf("Lv[%d][%d-%d][%d]\n",nd,nd1,na1-1,ca1-1);
		for (i=1;i<nd1;i++){
			for (k=1;k<ca1;k++){
				for (j=nd1;j<na1;j++){
					fprintf(fpw,"%23.15e ",lsv[i][j][k]);
					printf("%12.5e ",lsv[i][j][k]);
				}
				fprintf(fpw,"\n");
				printf("\n");
			}
			fprintf(fpw,"\n");
			printf("\n");
		}
		free_d_matrix(nd1,lx);
		free_d_matrix(nd1,lv);
		free_d_cube(nd1,na1,lsv);
	}
	fclose(fpw);
}

void ffpowlaw()	
{	
	int nd,ni,cp,cm,ca1,nd1,na1;
	int i,j,k;    
	double *x;    
	double *as,*bs,**gs,**hs;
	double **ag,***gg;
	FILE *fp,*fpw;

	if ((fp=fopen("work0.dat","r"))==NULL){
		printf("Cannot open parameters file curent drive.\n");
		exit(1);
	}
	fscanf(fp,"%d%d%d%d\n",&nd,&ni,&cp,&cm);
	nd1=nd+1;
	na1=nd1+ni;
	ca1=cp+cm+1;
	x=d_vector(na1);
	for (i=1;i<na1;i++){
		fscanf(fp,"%lf\n",&x[i]);
	}
	fclose(fp);	

	if(prns[3]==1){
		as=d_vector(nd1);
		bs=d_vector(nd1);
		gs=d_matrix(nd1,na1);
		hs=d_matrix(nd1,na1);

		if ((fp=fopen("work1.dat","r"))==NULL){
			printf("Cannot open parameters file curent drive.\n");
			exit(1);
		}
		for (i=1;i<nd1;i++)
			fscanf(fp,"%lf",&as[i]);
		for (i=1;i<nd1;i++)
			fscanf(fp,"%lf",&bs[i]);
		for (i=1;i<nd1;i++){
			for (j=1;j<na1;j++)
				fscanf(fp,"%lf",&gs[i][j]);
		}
		for (i=1;i<nd1;i++){
			for (j=1;j<na1;j++)
				fscanf(fp,"%lf",&hs[i][j]);
		}
		fclose(fp);

		if ((fpw=fopen(output3,"w"))==NULL){
			printf("Cannot open parameters file curent drive.\n");
			exit(1);
		}
		for (i=1;i<nd1;i++){
			fprintf(fpw,"X%d' = %lf ",i,as[i]);					
			for (j=1;j<na1;j++){
				if(gs[i][j] != 0)
					fprintf(fpw,"X%d^%lf ",j,gs[i][j]);		
			}

			fprintf(fpw,"- %lf ",bs[i]);						
			for (j=1;j<na1;j++){
				if(hs[i][j] != 0)
					fprintf(fpw,"X%d^%lf ",j,hs[i][j]);	         
			}
			fprintf(fpw,"\n");
		}

		fprintf(fpw,"\n\nSet values of dependent variables\n\n");
		for(i=1;i<nd1;i++)                                   
			fprintf(fpw,"X%d=%lf\n",i,x[i]);
		fprintf(fpw,"\n\n");
		fprintf(fpw,"Set values of independent variables\n\n");
		for(i=nd1;i<na1;i++)                                     
			fprintf(fpw,"X%d=%lf\n",i,x[i]);
		fprintf(fpw,"\n\n");

		fprintf(fpw,"&& ");										
		for(i=1;i<ni+1;i++)
			fprintf(fpw,"X%d ",i+nd);
		fprintf(fpw,"\n\nSet solution parameters\n\n");
		fprintf(fpw,"t0 = 0.0\nhr = 0.01\ntf = 10.0\n");				
		fclose(fpw);

		free(as);
		free(bs);
		free_d_matrix(nd1,gs);
		free_d_matrix(nd1,hs);
	}
	
	if(prng[3]==1){
		ag=d_matrix(nd1,ca1);
		gg=d_cube(nd1,na1,ca1);
		if ((fp=fopen("work2.dat","r"))==NULL){
			printf("Cannot open parameters file curent drive.\n");
			exit(1);
		}
		for (i=1;i<nd1;i++){
			for (k=1;k<ca1;k++){
				for (j=1;j<na1;j++)
					fscanf(fp,"%lf",&gg[i][j][k]);
			}
		}
		for (i=1;i<nd1;i++){
			for (k=1;k<ca1;k++)
				fscanf(fp,"%lf",&ag[i][k]);
		}
		fclose(fp);

		if ((fpw=fopen(output4,"w"))==NULL){
			printf("Cannot open parameters file curent drive.\n");
			exit(1);
		}

		for (i=1;i<nd1;i++){
			fprintf(fpw,"X%d' = ",i);
			for (k=1;k<ca1;k++){
				if(ag[i][k] > 1.0e-14){                             
					if(k==1) 
						fprintf(fpw," %lf ",ag[i][k]);				
					else{
						if(k<=cp) fprintf(fpw," + %lf ",ag[i][k]);		
						if(k>cp) fprintf(fpw," - %lf ",ag[i][k]);      
					}													
					
					for (j=1;j<na1;j++){
							if(gg[i][j][k] != 0.0)
								fprintf(fpw,"X%d^%lf ",j,gg[i][j][k]);	
					}
				}
			}
			fprintf(fpw,"\n");
		}
		
		fprintf(fpw,"\n\nSet values of dependent variables\n\n");
		for(i=1;i<nd1;i++)                                    
			fprintf(fpw,"X%d=%lf\n",i,x[i]);
		fprintf(fpw,"\n\n");
		fprintf(fpw,"Set values of independent variables\n\n");
		for(i=nd1;i<na1;i++)                                    
			fprintf(fpw,"X%d=%lf\n",i,x[i]);
		fprintf(fpw,"\n\n");

		fprintf(fpw,"&& ");								
		for(i=1;i<ni+1;i++)
			fprintf(fpw,"X%d ",i+nd);
		fprintf(fpw,"\n\nSet solution parameters\n\n");
		fprintf(fpw,"t0 = 0.0\nhr = 0.01\ntf = 10.0\n");			
		fclose(fpw);
	
		free_d_matrix(nd1,ag);
		free_d_cube(nd1,na1,gg);
		free(x);
	}
}

double** rematrix(int n1,double **a)
{
	int i,j;  
	double **b,**c;

	b=d_matrix(n1,n1+n1);
	for(i=1;i<n1;i++){
		for(j=1;j<n1;j++)
			b[i][j]=a[i][j];
		for(j=n1;j<n1+n1;j++){
			if (i==j-n1) b[i][j]=1.0;
			else b[i][j]=0.0;
		}
	}
	gj2(b,n1,n1+n1);
	for(i=1;i<n1;i++){
		for(j=1;j<n1;j++)
			b[i][j]=b[i][j+n1];
	}
	c=d_matrix(n1,n1);
	c=ope_matrix(n1,n1,n1,a,b);

	return(b);
	free_d_matrix(n1,c);
}

double** ope_matrix(int l1,int m1,int n1,double **a,double **b)
{
	int i,j,k;
	double **c,s;

	c=d_matrix(l1,n1);
	for (i=1;i<l1;i++){
		for (k=1;k<n1;k++){
			s=0.0;
			for (j=1;j<m1;j++)
				s=s+a[i][j]*b[j][k];
			c[i][k]=s;
		}
	}
	return(c);
}

void gj2(double **a,int nd1,int na1)
{
int i,j,k,p,q,k1,*l;   
double a2,a1,a0,d,eps=10e-4;
	l=i_vector(nd1);
	for (i=1;i<nd1;i++)
		l[i]=i;
	d=1.0;
	for (k=1;k<nd1;k++){
		a1=fabs(a[k][k]);
		p=k;
		q=k;
		for (j=k;j<nd1;j++){
			for (i=k;i<nd1;i++){
				if (a1<fabs(a[i][j])){
					a1=fabs(a[i][j]);
					p=i;
					q=j;
				}
			}
		}
		if (a1<eps){
			printf("singular matrix\n");
			exit(1);
		}
		if (k!=p){
			d=-d;
			for (j=k;j<na1;j++){
				a0=a[k][j];
				a[k][j]=a[p][j];
				a[p][j]=a0;
			}
		}
		if (k!=q){
			d=-d;
			for (i=1;i<nd1;i++){
				a0=a[i][k];
				a[i][k]=a[i][q];
				a[i][q]=a0;
			}
			j=l[k];l[k]=l[q];l[q]=j;
		}
		a1=a[k][k];
		d=d*a1;
		k1=k+1;
		for (j=k1;j<na1;j++)
			a[k][j]=a[k][j]/a1;
		for (i=1;i<nd1;i++){
			if (i!=k){
			a2=a[i][k];
			for (j=k1;j<na1;j++)
				a[i][j]=a[i][j]-a2*a[k][j];
			}
		}
	}

	for (j=nd1;j<na1;j++){
		for (i=1;i<nd1;i++){
			p=l[i];
			a[p][nd1-1]=a[i][j];
		}
		for (i=1;i<nd1;i++)
			a[i][j]=a[i][nd1-1];
	}
	free(l);
}

void gj(double **a,int nd1,int nd2)
{
	int i,j,k,p,q,k1,*l,n1;   
	double a2,a1,a0,d,eps=10e-5;         //change if necessary.
	l=i_vector(nd1);
	for (i=1;i<nd1;i++)
		l[i]=i;
	d=1.0;
	for (k=1;k<nd1;k++){
		a1=fabs(a[k][k]);
		p=k;
		q=k;
		for (j=k;j<nd1;j++){
			for (i=k;i<nd1;i++){
				if (a1<fabs(a[i][j])){
					a1=fabs(a[i][j]);
					p=i;
					q=j;
				}
			}
		}
		if (a1<eps){
			printf("singular matrix\n");
			exit(1);
		}
		if (k!=p){
			d=-d;
			for (j=k;j<nd2;j++){
				a0=a[k][j];
				a[k][j]=a[p][j];
				a[p][j]=a0;
			}
		}
		if (k!=q){
			d=-d;
			for (i=1;i<nd1;i++){
				a0=a[i][k];
				a[i][k]=a[i][q];
				a[i][q]=a0;
			}
			j=l[k];l[k]=l[q];l[q]=j;
		}
		a1=a[k][k];
		d=d*a1;
		k1=k+1;
		for (j=k1;j<nd2;j++)
			a[k][j]=a[k][j]/a1;
		for (i=1;i<nd1;i++){
			if (i!=k){
			a2=a[i][k];
			for (j=k1;j<nd2;j++)
				a[i][j]=a[i][j]-a2*a[k][j];
			}
		}
	}
	n1=nd2;
	for (j=n1-1;j<nd2;j++){				
		for (i=1;i<nd1;i++){
			p=l[i];
			a[p][nd1-1]=a[i][j];
		}
		for (i=1;i<nd1;i++)
			a[i][j]=a[i][nd1-1];
	}
	free(l);
}

double* d_vector(int n)
{
    double *a;
    if((a=(double*)malloc(n*sizeof(double)))==NULL){
	printf("memory unavailable\n");
	exit(1);
    }
    return(a);
}

double** d_matrix(int m,int n)
{
double **a;
int i;
	if((a=(double**)malloc(m*sizeof(double*)))==NULL){
		printf("memory unavailable\n");
		exit(1);
	}
	for (i=0;i<m;i++){
		if((a[i]=(double*)malloc(n*sizeof(double)))==NULL){
			printf("memory unavailable\n");
			exit(1);
		}
	}
	return(a);
}

double*** d_cube(int l,int m,int n)
{
double ***a;
int i,j;
	if((a=(double***)malloc(l*sizeof(double**)))==NULL){
		printf("memory unavailable\n");
		exit(1);
	}
	for (i=0;i<l;i++){
		if((a[i]=(double**)malloc(m*sizeof(double*)))==NULL){
			printf("memory unavailable\n");
			exit(1);
		}
	}
	for (i=0;i<l;i++){
		for (j=0;j<m;j++){
			if((a[i][j]=(double*)malloc(n*sizeof(double)))==NULL){
				printf("memory unavailable\n");
				exit(1);
			}
		}
	}
	return(a);
}

int* i_vector(int n)
{
    int *a;
    if((a=(int*)malloc(n*sizeof(int)))==NULL){
		printf("memory unavailable\n");
		exit(1);
    }
    return(a);
}

void free_d_matrix(int m,double **a)
{
	int i;
	for (i=0;i<m;i++)
		free(a[i]);
	free(a);
}

void free_i_matrix(int m,double **a)
{
	int i;
	for (i=0;i<m;i++)
		free(a[i]);
	free(a);
}

void free_d_cube(int l,int m,double ***a)
{
	int i,j;
	for (i=0;i<l;i++){
		for (j=0;j<m;j++)
			free(a[i][j]);
	}
	for (i=0;i<l;i++)
		free(a[i]);
	free(a);
}

void free_C_d_matrix(int m,complex<double> **a)
{
	int i;
	for (i=0;i<m;i++)
		free(a[i]);
	free(a);
}

void eigenvalue()
{
	double **a,*wr,*wi;
	double *xp,*fxp,**df,**vp,**svp;
	double sgm=1,eps=10e-15;
	int i; 
	int nd,nf,ni,na,cp,cm,nd1,na1,ca1,dummy;
	FILE *fp,*fpw;

	fo1(&nd,&ni,&cp,&cm);    
	na=nd+ni;
	ca1=cp+cm+1;
	nd1=nd+1;
	na1=na+1;
	xp=d_vector(na1);
	fo(&nd,&ni,&cp,&cm,xp,0);
	fxp=d_vector(nd1);
	vp=d_matrix(nd1,2);
	svp=d_matrix(nd1,ca1);
	df=d_matrix(nd1,nd1+1);
	nf=nd1;
	printf("nd=%d ni=%d cp=%d cm=%d\n",nd,ni,cp,cm);
	
	wr=d_vector(nd1);
	wi=d_vector(nd1);
	a=d_matrix(nd1,nd1);

	if ((fp=fopen("work0.dat","r"))==NULL){
		printf("Cannot open parameters file curent drive.\n");
		exit(1);
	}
	fscanf(fp,"%d%d%d%d",&dummy,&dummy,&dummy,&dummy);
	na1=nd+ni+1;
	for (i=1;i<na1;i++){
		fscanf(fp,"%lf",&xp[i]);
	}
	fclose(fp);

	complex_step1(xp,df,fxp,vp,svp,nd1,ni,cp,cm,sgm);   

	balance(df,nd);
	elmhes(df,nd);
	hqr(df,nd,wr,wi);
	sort(nd,wr,wi);

	if ((fpw=fopen(output5,"w"))==NULL){
			printf("Cannot open parameters file curent drive.\n");
			exit(1);
	}

	printf("Eigenvalues----------\n");
	for(i=1;i<=nd;i++){
		printf("%23.15e+%23.15ei\n",wr[i],wi[i]);
		fprintf(fpw,"%23.15e+%23.15ei\n",wr[i],wi[i]);
	}

	fclose(fpw);
	free_d_matrix(nd1,a);
    free(wr);    
	free(wi);
	free(xp);
	free(fxp);
	free_d_matrix(nd1,svp);
	free_d_matrix(nd1,vp);
	free_d_matrix(nd1,df);
}

void balance(double **a, int n)
{
	int last,j,i;
	double s,r,g,f,c,sqrdx;

	sqrdx=RADIX*RADIX;
	last=0;
	while (last==0) {
		last=1;
		for (i=1;i<=n;i++) {
			r=0.0,c=0.0;
			for (j=1;j<=n;j++)
				if (j != i) {
					c += fabs(a[j][i]);
					r += fabs(a[i][j]);
				}
			if (c && r) {
				g=r/RADIX;
				f=1.0;
				s=c+r;
				while (c<g) {
					f *= RADIX;
					c *= sqrdx;
				}
				g=r*RADIX;
				while (c>g) {
					f /= RADIX;
					c /= sqrdx;
				}
				if ((c+r)/f < 0.95*s) {
					last=0;
					g=1.0/f;
					for (j=1;j<=n;j++) a[i][j] *= g;
					for (j=1;j<=n;j++) a[j][i] *= f;
				}
			}
		}
	}
}


void elmhes(double **a, int n)
{
	int m,j,i;
	double y,x;

	for (m=2;m<n;m++) {
		x=0.0;
		i=m;
		for (j=m;j<=n;j++) {
			if (fabs(a[j][m-1]) > fabs(x)) {
				x=a[j][m-1];
				i=j;
			}
		}
		if (i != m) {
			for (j=m-1;j<=n;j++) SWAP(a[i][j],a[m][j]);
			for (j=1;j<=n;j++) SWAP(a[j][i],a[j][m]);
		}
		if (x) {
			for (i=m+1;i<=n;i++) {
				if ((y=a[i][m-1]) != 0.0) {
					y /= x;
					a[i][m-1]=y;
					for (j=m;j<=n;j++) 
						a[i][j] -= y*a[m][j];
					for (j=1;j<=n;j++) 
						a[j][m] += y*a[j][i];
				}
			}
		}
	}
}


void hqr(double **a,int n,double *wr,double *wi)
{
	int nn,m,l,k,j,its,i,mmin;
	double z,y,x,w,v,u,t,s,r,q,p,anorm;

	anorm=fabs(a[1][1]);
	for (i=2;i<=n;i++)
		for (j=(i-1);j<=n;j++)
			anorm += fabs(a[i][j]);
	nn=n;
	t=0.0;
	while (nn >= 1) {
		its=0;
		do {
			for (l=nn;l>=2;l--) {
				s=fabs(a[l-1][l-1])+fabs(a[l][l]);
				if (s == 0.0) s=anorm;
				if ((fabs(a[l][l-1]) +s) == s) break;
			}
			x=a[nn][nn];
			if (l == nn) {
				wr[nn]=x+t;
				wi[nn--]=0.0;
			} else {
				y=a[nn-1][nn-1];
				w=a[nn][nn-1]*a[nn-1][nn];
				if (l == (nn-1)) {
					p=0.5*(y-x);
					q=p*p+w;
					z=sqrt(fabs(q));
					x += t;
					if (q >= 0.0) {
						z=p+SIGN(z,p);
						wr[nn-1]=wr[nn]=x+z;
						if (z) wr[nn]=x-w/z;
						wi[nn-1]=wi[nn]=0.0;
					} else {
						wr[nn-1]=wr[nn]=x+p;
						wi[nn-1]=-(wi[nn]=z);
					}
					nn -= 2;
				} else {
					if (its == 30) {printf("Too many iterations in hqr");exit(0);}
					if (its == 10 || its == 20) {
						t += x;
						for (i=1;i<=nn;i++) a[i][i] -= x;
						s=fabs(a[nn][nn-1])+fabs(a[nn-1][nn-2]);
						y=x=0.75*s;
						w = -0.4375*s*s;
					}
					++its;
					for (m=(nn-2);m>=l;m--) {
						z=a[m][m];
						r=x-z;
						s=y-z;
						p=(r*s-w)/a[m+1][m]+a[m][m+1];
						q=a[m+1][m+1]-z-r-s;
						r=a[m+2][m+1];
						s=fabs(p)+fabs(q)+fabs(r);
						p /= s;
						q /= s;
						r /= s;
						if (m == l) break;
						u=fabs(a[m][m-1])*(fabs(q)+fabs(r));
						v=fabs(p)*(fabs(a[m-1][m-1])+fabs(z)+fabs(a[m+1][m+1]));
						if ((double)(u+v) <= v) break;
					}
					for (i=m+2;i<=nn;i++) {
						a[i][i-2]=0.0;
						if (i != (m+2)) a[i][i-3]=0.0;
					}
					for (k=m;k<=nn-1;k++) {
						if (k != m) {
							p=a[k][k-1];
							q=a[k+1][k-1];
							r=0.0;
							if (k != (nn-1)) r=a[k+2][k-1];
							if ((x=fabs(p)+fabs(q)+fabs(r)) != 0.0) {
								p /= x;
								q /= x;
								r /= x;
							}
						}
						if ((s=SIGN(sqrt(p*p+q*q+r*r),p)) != 0.0) {
							if (k == m) {
								if (l != m)
								a[k][k-1] = -a[k][k-1];
							} else
								a[k][k-1] = -s*x;
							p += s;
							x=p/s;
							y=q/s;
							z=r/s;
							q /= p;
							r /= p;
							for (j=k;j<=nn;j++) {
								p=a[k][j]+q*a[k+1][j];
								if (k != (nn-1)) {
									p += r*a[k+2][j];
									a[k+2][j] -= p*z;
								}
								a[k+1][j] -= p*y;
								a[k][j] -= p*x;
							}
							mmin = nn < k+3 ? nn : k+3;
							for (i=l;i<=mmin;i++) {
								p=x*a[i][k]+y*a[i][k+1];
								if (k != (nn-1)) {
									p += z*a[i][k+2];
									a[i][k+2] -= p*r;
								}
								a[i][k+1] -= p*q;
								a[i][k] -= p;
							}
						}
					}
				}
			}
		} while (l < (nn-1));
	}
}

void sort(int n,double *wr,double *wi)
{
	int i,j,s;
	double x,y,t;

	for (i=1;i<n;i++) {
		x=wr[i];
		y=wi[i];
		s=i;
		for (j=i+1;j<=n;j++) {
			if (fabs(wr[j]) > fabs(x)) {
				x=wr[j];
				y=wi[j];
				s=j;
			}
		}
		t=wr[i];wr[i]=wr[s];wr[s]=t;
		t=wi[i];wi[i]=wi[s];wi[s]=t;
	}
}

complex<double>* C_d_vector(int n)
{
    complex<double> *a;
    if((a=(complex<double>*)malloc(n*sizeof(complex<double>)))==NULL){
	printf("memory unavailable\n");
	exit(1);
    }
    return(a);
}

complex<double>** C_d_matrix(int m,int n)
{
complex<double> **a;
int i;
	if((a=(complex<double>**)malloc(m*sizeof(complex<double>*)))==NULL){
		printf("memory unavailable\n");
		exit(1);
	}
	for (i=0;i<m;i++){
		if((a[i]=(complex<double>*)malloc(n*sizeof(complex<double>)))==NULL){
			printf("memory unavailable\n");
			exit(1);
		}
	}
	return(a);
}

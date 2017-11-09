#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<gsl/gsl_eigen.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_permutation.h>
#include<gsl/gsl_linalg.h>

#include<gsl/gsl_odeiv.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include<gsl/gsl_odeiv.h>
#include <gsl/gsl_errno.h>
#include<gsl/gsl_math.h>
#include<gsl/gsl_eigen.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_sort.h>
#include<gsl/gsl_sort_vector.h>

#define id 7
#define idm 7
#define idm1 idm*idm
#define M 1
#define N id
#define niter 1000000
#define npart 2000000
#define period 365


class sys_para{
public:
    double ep,a,a0,b,d,tpi,w,NN,R0,nu,b1,bet1,bet2,b1bet1,b1bet2,b1bet3,muM0,muM,muB,dB,PIM,PIB,NB,NM,Rb,muMav,PIM0;
double R01,PIMav,R02,lam,lamMB,lamMH,PIM01,PIM1,peakPIM,fracPIM,lamMHav,lamMBav;
double PIH,muH,bet3,b2,alpha,del,dH,tau,c,q,NBsq,R0d,Corrfactr,Omeg,tauB,sigB,sigM,dM,etaM,etaB,q1,q2,q3,q4,b0,pc1,pc2,R0D,RcD,Corrfact2;
int con,zit;
};

int
monod (double t, const double ym0[], double dym[],
void *params);
int
jac (double t, const double y0[], double *dfdy,
double dfdt[], void *params);
int
jacm (double t, const double ym0[], double *dfmdy,
double dfmdt[], void *params);

double BRN(sys_para &sp,double lam);
double PIM0fn(sys_para &sp,double PIM0);
double b10fn(sys_para &sp,double b10);
void root(sys_para &sp);
void rootPIM0(sys_para &sp);
void write(sys_para &sp,FILE *fp[4],double zz[10001][id+1]);


int
main (void)
{
int i,j,l,s;
double nrm,maxnorm,max,min,zz[10001][id+1];
double ea1,ea2,ea3,ea4,ea5,es1,es2,es3,es4,es5,es6,es7,er=10e-7;

FILE *fp[4];
fp[0]=fopen("bifSEIRSnSEIstableforcedCorrfactr_SInSI_tauB_p0_alph_p0_d_p3.dat","w");
fp[1]=fopen("bifSEIRSnSEIunstableforcedCorrfactr_SInSI_tauB_p0_alph_p0_d_p3.dat","w");

sys_para sp;


sp.NN=1;
sp.NB=200000;
sp.NM=2*pow(10,5);
sp.nu=10;

//----------------------------PARAMETERS FOR SEIR and SEI model near WNV model-----------------------------------------------//
sp.b1=0.1;//biting rate of mosquitoes on birds
sp.bet1=0.16;//transmission probability from bird to mosquito
sp.bet2=0.88;//transmission probability from mosquito to bird
sp.sigB=1000;//1000;//a large value of sigB means very less or zero amount of time is spent in latent period of birds
sp.sigM=1000;//a large value of sigM means very less or zero amount of time spent in latent period fo mosquitoes
sp.PIB=1000.0;// recruitment rate of  suceptible birds per days
sp.PIM=30000;//recruitment rate of uninfected mosquitoes (per day)-- variable parameter
sp.muB=0.001;//1/1000.0;// 1/sp.muB is average life span of a bird (in days)
sp.dB=0.005; //WNV-induced death rate of birds (fraction per day)
sp.dM=0;//WNV-induced death rate of mosquitoes
sp.etaB=0.0;
sp.etaM=0.0;
sp.muM=1/16.0;// 1/sp.muM is average life span of a mosquito (in days)--variable parameter
sp.tauB=0.0000001;//001428;//recovery rate for birds
sp.w=0.00;//Rate of loss of immunity
//-----------------------------------------------------------------------------------------------------------------------------------//

sp.q1=sp.sigB+sp.muB;
sp.q2=sp.tauB+sp.muB+sp.dB;
sp.q3=sp.sigM+sp.muM;
sp.q4=sp.muM+sp.dM;
sp.tpi=8*atan(1.0); 
sp.Omeg=sp.tpi/period;  
    
    
sp.ep=0.0;


double y0[N+1],y1[N+1];

double m11,m12,m13,m14,m15,m16,m17,m21,m22,m23,m24,m25,m26,m27,m31,m32,m33,m34,m35,m36,m37,m41,m42,m43,m44,m45,m46,m47;
double m130,m140,m230,m240,m330,m340,m430,m440,m51,m52,m53,m54,m55,m56,m57,m61,m62,m63,m64,m65,m66,m67,m71,m72,m73,m74,m75,m76,m77;
double gg1,gg2,fe,t;

sp.con=1;//controlling vector control
sp.fracPIM=1;//1-fracPIM is fraction of mosquitoes killed during RMS of the peak;
sp.d=0.3;
sp.PIM=sp.PIM0;



gsl_matrix *m1=gsl_matrix_alloc(sqrt(idm1),sqrt(idm1));
gsl_matrix *Minv=gsl_matrix_alloc(sqrt(idm1),sqrt(idm1));



gsl_vector *c1=gsl_vector_alloc(sqrt(idm1)),*c0=gsl_vector_alloc(sqrt(idm1)),*f0=gsl_vector_alloc(sqrt(idm1));

gg1=(sp.dB+sp.muB+sp.tauB)*sp.muM;
 gg2=(sp.Omeg*sp.Omeg)+((sp.dB+sp.muB+sp.tauB+sp.muM)*(sp.dB+sp.muB+sp.tauB+sp.muM)); 
fe=(1-((gg1*sp.d*sp.d*0.5)/gg2));
sp.Corrfactr=fe;
sp.R0d=0.65;
for(sp.R0d=1.1;sp.R0d>=0.6;sp.R0d-=0.0005)
{
sp.R02=sp.R0d;
sp.b0=sp.muB*sp.PIB*sp.q1*sp.q2*sp.q3*sp.q4*(sp.b1*sp.bet1*sp.muB*(sp.q2*sp.etaB+sp.sigB)+sp.muM*(2*(sp.muB*sp.q2+sp.muB*sp.sigB+sp.sigB*sp.tauB)-sp.q1*sp.q2*sp.R0d));

sp.a0=sp.PIB*sp.q3*sp.q4*(sp.q2*sp.muB+sp.sigB*(sp.muB+sp.tauB))*(sp.muB*sp.b1*sp.bet1*(sp.etaB*sp.q2+sp.sigB)+sp.muB*sp.muM*(sp.q2+sp.sigB)+sp.muM*sp.sigB*sp.tauB);

sp.RcD=sqrt(1-(sp.b0*sp.b0)/(4*sp.a0*sp.q1*sp.q1*sp.q2*sp.q2*sp.q3*sp.q4*sp.muB*sp.muB*sp.muM*sp.PIB));

sp.pc1=(-sp.Omeg*sp.Omeg+sp.q1*sp.q2+sp.q1*sp.q3+sp.q1*sp.q4+sp.q2*sp.q3+sp.q2*sp.q4+sp.q3*sp.q4);
sp.pc2=sp.Omeg*sp.Omeg*(sp.q1+sp.q2+sp.q3+sp.q4)-sp.q1*sp.q2*(sp.q3+sp.q4)-sp.q3*sp.q4*(sp.q1+sp.q2);

sp.Corrfact2=(1-(sp.d*sp.d*0.5*sp.pc1*sp.q1*sp.q2*sp.q3*sp.q4)/((sp.Omeg*sp.Omeg*sp.pc1*sp.pc1)+sp.pc2*sp.pc2));

//sp.R02=sqrt((sp.b1*sp.b1*sp.bet2*sp.bet1*sp.muB*sp.PIM0)/(sp.muM*sp.muM*(sp.muB+sp.dB)*sp.PIB));

//rootPIM0(sp);
sp.PIM0=(sp.R0d*sp.Corrfact2*sp.q1*sp.q2*sp.q3*sp.q4*sp.muM*sp.PIB)/(sp.b1*sp.b1*sp.bet1*sp.bet2*sp.muB*(sp.etaB*sp.q2+sp.sigB)*(sp.etaM*sp.q4+sp.sigM));
//sp.PIM0=(sp.R02*sp.R02*sp.muM*sp.muM*(sp.muB+sp.dB)*sp.PIB)/(sp.b1*sp.b1*sp.bet2*sp.bet1*sp.muB);

printf("%lf %lf %lf\n",sp.PIM0,sp.R02,sp.Corrfact2);
//for(double il=1;il<=100;il++)
//for(double il=0.1;il<=1;il+=0.01)
double il=1;
while(il<=50)
{
y0[1]=sp.PIM0/sp.muM*drand48()+0.1;
y0[2]=sp.PIM0/sp.muM*drand48()+0.1;
y0[3]=sp.PIM0/sp.muM*drand48()+0.1;
y0[4]=sp.PIB/sp.muB*drand48()+0.1;
y0[5]=sp.PIB/sp.muB*drand48()+0.1;
y0[6]=sp.PIB/sp.muB*drand48()+0.1;
y0[7]=sp.PIB/sp.muB*drand48()+0.1;
//for(t=0;t<=period;t+=0.1)
{

//printf("%lf %lf %lf %lf %lf %lf\n",sp.PIM0,sp.R02,y0[1],y0[2],y0[3],y0[4]);
es1=10;es2=10;es3=10;es4=10;es5=10;es6=10;es7=10;
while(es1>er && es2>er && es3>er && es4>er && es5>er && es6>er && es7>er)
{
  gsl_permutation *p=gsl_permutation_alloc(sqrt(idm1));

sp.NB=y0[4]+y0[5]+y0[6]+y0[7];
sp.NBsq=sp.NB*sp.NB;

m11=-sp.b1*sp.bet1*(sp.etaB*y0[5]+y0[6])/sp.NB-sp.muM;
m12=0;
m130=0;
m140=sp.b1*sp.bet1*y0[1]*(sp.etaB*y0[5]+y0[6])/sp.NBsq;
m15=sp.b1*sp.bet1*y0[1]*(sp.etaB*y0[5]+y0[6])/sp.NBsq - (sp.b1*sp.bet1*y0[1]*sp.etaB)/sp.NB;
m16=-sp.b1*sp.bet1*y0[1]/sp.NB + sp.b1*sp.bet1*y0[1]*(sp.etaB*y0[5]+y0[6])/sp.NBsq;
m17=sp.b1*sp.bet1*y0[1]*(sp.etaB*y0[5]+y0[6])/sp.NBsq;

gsl_matrix_set(m1,0,0,m11);gsl_matrix_set(m1,0,1,m12);
gsl_matrix_set(m1,0,2,m130);gsl_matrix_set(m1,0,3,m140);
gsl_matrix_set(m1,0,4,m15);gsl_matrix_set(m1,0,5,m16);
gsl_matrix_set(m1,0,6,m17);

m21=sp.b1*sp.bet1*(sp.etaB*y0[5]+y0[6])/sp.NB;
m22=-sp.muM-sp.sigM;
m230=0;
m240=-sp.b1*sp.bet1*y0[1]*(sp.etaB*y0[5]+y0[6])/sp.NBsq;
m25=-sp.b1*sp.bet1*y0[1]*(sp.etaB*y0[5]+y0[6])/sp.NBsq + (sp.b1*sp.bet1*y0[1]*sp.etaB)/sp.NB;
m26=sp.b1*sp.bet1*y0[1]/sp.NB-sp.b1*sp.bet1*y0[1]*(sp.etaB*y0[5]+y0[6])/sp.NBsq;
m27=-sp.b1*sp.bet1*y0[1]*(sp.etaB*y0[5]+y0[6])/sp.NBsq;

gsl_matrix_set(m1,1,0,m21);gsl_matrix_set(m1,1,1,m22);
gsl_matrix_set(m1,1,2,m230);gsl_matrix_set(m1,1,3,m240);
gsl_matrix_set(m1,1,4,m25);gsl_matrix_set(m1,1,5,m26);
gsl_matrix_set(m1,1,6,m27);

m31=0;
m32=sp.sigM;
m33=-(sp.muM+sp.dM);

m330=m33;
m340=0;
m35=0;
m36=0;
m37=0;

gsl_matrix_set(m1,2,0,m31);gsl_matrix_set(m1,2,1,m32);
gsl_matrix_set(m1,2,2,m330);gsl_matrix_set(m1,2,3,m340);
gsl_matrix_set(m1,2,4,m35);gsl_matrix_set(m1,2,5,m36);
gsl_matrix_set(m1,2,6,m37);

m41=0;
m42=-sp.b1*sp.bet2*sp.etaM*y0[4]/sp.NB;
m43=-sp.b1*sp.bet2*y0[4]/sp.NB;
m44=-sp.b1*sp.bet2*(sp.etaM*y0[2]+y0[3])/sp.NB+(sp.b1*sp.bet2*(sp.etaM*y0[2]+y0[3])*y0[4]/sp.NBsq)-sp.muB;
m45=sp.b1*sp.bet2*(sp.etaM*y0[2]+y0[3])*y0[4]/sp.NBsq;
m46=m45;
m47=m45+sp.w;
m430=m43;
m440=m44;

gsl_matrix_set(m1,3,0,m41);gsl_matrix_set(m1,3,1,m42);
gsl_matrix_set(m1,3,2,m430);gsl_matrix_set(m1,3,3,m440);
gsl_matrix_set(m1,3,4,m45);gsl_matrix_set(m1,3,5,m46);
gsl_matrix_set(m1,3,6,m47);

m51=0;
m52=sp.b1*sp.bet2*sp.etaM*y0[4]/sp.NB;
m53=sp.b1*sp.bet2*y0[4]/sp.NB;
m54=sp.b1*sp.bet2*(sp.etaM*y0[2]+y0[3])/sp.NB-sp.b1*sp.bet2*(sp.etaM*y0[2]+y0[3])*y0[4]/sp.NBsq;
m55=-sp.b1*sp.bet2*(sp.etaM*y0[2]+y0[3])*y0[4]/sp.NBsq-sp.muB-sp.sigB;
m56=-sp.b1*sp.bet2*(sp.etaM*y0[2]+y0[3])*y0[4]/sp.NBsq;
m57=-sp.b1*sp.bet2*(sp.etaM*y0[2]+y0[3])*y0[4]/sp.NBsq;

gsl_matrix_set(m1,4,0,m51);gsl_matrix_set(m1,4,1,m52);
gsl_matrix_set(m1,4,2,m53);gsl_matrix_set(m1,4,3,m54);
gsl_matrix_set(m1,4,4,m55);gsl_matrix_set(m1,4,5,m56);
gsl_matrix_set(m1,4,6,m57);

m61=0;
m62=0;
m63=0;
m64=0;
m65=sp.sigB;
m66=-(sp.muB+sp.dB+sp.tauB);
m67=0;


gsl_matrix_set(m1,5,0,m61);gsl_matrix_set(m1,5,1,m62);
gsl_matrix_set(m1,5,2,m63);gsl_matrix_set(m1,5,3,m64);
gsl_matrix_set(m1,5,4,m65);gsl_matrix_set(m1,5,5,m66);
gsl_matrix_set(m1,5,6,m67);

m71=0;
m72=0;
m73=0;
m74=0;
m75=0;
m76=sp.tauB;
m77=-sp.muB;


gsl_matrix_set(m1,6,0,m71);gsl_matrix_set(m1,6,1,m72);
gsl_matrix_set(m1,6,2,m73);gsl_matrix_set(m1,6,3,m74);
gsl_matrix_set(m1,6,4,m75);gsl_matrix_set(m1,6,5,m76);
gsl_matrix_set(m1,6,6,m77);

  gsl_linalg_LU_decomp(m1,p,&s);
  gsl_linalg_LU_invert(m1,p,Minv);
//sp.PIM=sp.PIM0*(1+sp.d*sin(sp.tpi*t/period+0.1));
gsl_vector_set(f0,0,sp.PIM0-sp.b1*sp.bet1*y0[1]*(sp.etaB*y0[5]+y0[6])/sp.NB-sp.muM*y0[1]);
gsl_vector_set(f0,1,sp.b1*sp.bet1*y0[1]*(sp.etaB*y0[5]+y0[6])/sp.NB-(sp.sigM+sp.muM)*y0[2]);
gsl_vector_set(f0,2,sp.sigM*y0[2]-sp.muM*y0[3]-sp.dM*y0[3]);
gsl_vector_set(f0,3,sp.PIB-sp.b1*sp.bet2*(sp.etaM*y0[2]+y0[3])*y0[4]/sp.NB-sp.muB*y0[4]+sp.w*y0[7]);
gsl_vector_set(f0,4,sp.b1*sp.bet2*(sp.etaM*y0[2]+y0[3])*y0[4]/sp.NB-sp.muB*y0[5]-sp.sigB*y0[5]);
gsl_vector_set(f0,5,sp.sigB*y0[5]-sp.muB*y0[6]-sp.dB*y0[6]-sp.tauB*y0[6]);
gsl_vector_set(f0,6,sp.tauB*y0[6]-sp.muB*y0[7]-sp.w*y0[7]);
       gsl_blas_dgemv(CblasNoTrans,1.0,Minv,f0,0.0,c0);


ea1=fabs((gsl_vector_get(c0,0)-y0[1]))/gsl_vector_get(c0,0);
ea2=fabs((gsl_vector_get(c0,1)-y0[2]))/gsl_vector_get(c0,1);


y1[1]=y0[1]-gsl_vector_get(c0,0);
y1[2]=y0[2]-gsl_vector_get(c0,1);
y1[3]=y0[3]-gsl_vector_get(c0,2);
y1[4]=y0[4]-gsl_vector_get(c0,3);
y1[5]=y0[5]-gsl_vector_get(c0,4);
y1[6]=y0[6]-gsl_vector_get(c0,5);
y1[7]=y0[7]-gsl_vector_get(c0,6);
//printf("%lf %lf %lf %lf\n",es1,es2,es3,es4);

es1=fabs((y1[1]-y0[1]));///y1[1]);
es2=fabs((y1[2]-y0[2]));///y1[2]);
es3=fabs((y1[3]-y0[3]));///y1[3]);
es4=fabs((y1[4]-y0[4]));///y1[4]);
es5=fabs((y1[5]-y0[5]));
es6=fabs((y1[6]-y0[6]));
es7=fabs((y1[7]-y0[7]));
//printf("%lf\n",il);
//printf("%lf %lf %lf %lf\n",es1,es2,es3,es4);
y0[1]=y1[1];
y0[2]=y1[2];
y0[3]=y1[3];
y0[4]=y1[4];
y0[5]=y1[5];
y0[6]=y1[6];
y0[7]=y1[7];
sp.NB=y0[4]+y0[5]+y0[6]+y0[7];
sp.NBsq=sp.NB*sp.NB;
}
sp.lamMH=sp.b1*sp.bet2*(sp.etaM*y0[2]+y0[3])/sp.NB;
if(sp.lamMH>=0 && y0[1]>=0 && y0[2]>=0 && y0[3]>=0 && y0[4]>=0 && y0[5]>=0 && y0[6]>=0 && y0[7]>=0 && y0[2]>0.851 || y0[2]<0.000001 && sp.R02<1)
{
fprintf(fp[0],"%lf %lf %lf %lf %lf %lf %lf\n",sp.PIM0,sp.R02,sp.lamMH,y0[1],y0[2],y0[3],y0[4]);
il++;
}
if(sp.lamMH>=0 && y0[1]>=0 && y0[2]>=0 && y0[3]>=0 && y0[4]>=0 && y0[5]>=0 && y0[6]>=0 && y0[7]>=0 && y0[2]>0.851 && sp.R02>1)
{
fprintf(fp[0],"%lf %lf %lf %lf %lf %lf %lf\n",sp.PIM0,sp.R02,sp.lamMH,y0[1],y0[2],y0[3],y0[4]);
il++;
}
if(sp.lamMH>=0 && y0[1]>=0 && y0[2]>=0 && y0[3]>=0 && y0[4]>=0 && y0[5]>=0 && y0[6]>=0 && y0[7]>=0 && y0[2]<0.851  && y0[2]>0.000001 && sp.R02 <1)
{
il++;
fprintf(fp[1],"%lf %lf %lf %lf %lf %lf %lf\n",sp.PIM0,sp.R02,sp.lamMH,y0[1],y0[2],y0[3],y0[4]);
//printf("%lf\n",il);
}
if(sp.lamMH>=0 && y0[1]>=0 && y0[2]>=0 && y0[3]>=0 && y0[4]>=0 && y0[5]>=0 && y0[6]>=0 && y0[7]>=0 && y0[2]<0.851 && sp.R02 >1)
{
il++;
fprintf(fp[1],"%lf %lf %lf %lf %lf %lf %lf\n",sp.PIM0,sp.R02,sp.lamMH,y0[1],y0[2],y0[3],y0[4]);
//printf("%lf\n",il);
}

//printf("%lf %lf %lf\n",il,t,sp.Corrfact2);
}


}
}
fclose(*fp);
}


int
monod (double t, const double ym0[], double dym[],
void *params)
{
//double mu = *(double *)params;
int i,j;
  sys_para sp=*(sys_para*)params;

//monodromy matrix
  dym[1]=-sp.muM*ym0[1]+((sp.b1bet1*sp.muB*sp.PIM)/(sp.muM*sp.PIB*sp.lam))*ym0[3];
  dym[2]=-sp.muM*ym0[2]+((sp.b1bet1*sp.muB*sp.PIM)/(sp.muM*sp.PIB*sp.lam))*ym0[4];
  dym[3]=((sp.b1bet2)/sp.lam)*ym0[1] - (sp.dB+sp.muB)*ym0[3];
  dym[4]=((sp.b1bet2)/sp.lam)*ym0[2] - (sp.dB+sp.muB)*ym0[4];

  

 
return GSL_SUCCESS;
}


void write(sys_para &sp,FILE *fp[4],double zz[10001][id+1])
{
int jj=12,yes1=0,yes2=0,yes3=0,yes4=0;
//fprintf(fp,"%lf %lf %lf\n",zit*1.0,y0[1],y0[2]);

while(jj<=sp.zit-11)
{
if(jj%50==0)
{
  //  fprintf(fp[0],"%lf  %lf   %lf  %lf %lf  %lf %lf\n",sp.b1,sp.R01,sp.PIM0,sp.zz[jj][4],sp.zz[jj][1],sp.lamMHav,sp.lamMBav);
    //fprintf(fp[1],"%lf  %lf   %lf  %lf %lf  %lf %lf\n",sp.b1,sp.R01,sp.PIM0,sp.zz[jj][4],sp.zz[jj][2],sp.lamMHav,sp.lamMBav);

 if(zz[jj][1]>=zz[jj-1][1] && zz[jj][1]>=zz[jj+1][1])
   {
    fprintf(fp[2],"%lf  %lf   %lf  %lf %lf  %lf %lf %lf\n",sp.b1,sp.R01,sp.PIM0,sp.PIMav,zz[jj][4],zz[jj][1],sp.lamMHav,sp.lamMBav);
   }

 if(zz[jj][2]>=zz[jj-1][2] && zz[jj][2]>=zz[jj+1][2])
   fprintf(fp[3],"%lf  %lf  %lf  %lf  %lf  %lf %lf %lf\n",sp.b1,sp.R01,sp.PIM0,sp.PIMav,zz[jj][4],zz[jj][2],sp.lamMHav,sp.lamMBav);
}
jj++;
}

for(jj=1;jj<=100011;jj++)
   {
    zz[jj][1]=0;
    zz[jj][2]=0;
   }

}


int
jac (double t, const double y0[], double *dfdy,
double dfdt[], void *params)
{

//double mu = *(double *)params;
//sys_para sp=*(sys_para*)params;
//gsl_matrix_view dfdy_mat
//= gsl_matrix_view_array (dfdy, 2, 2);
//gsl_matrix * m = &dfdy_mat.matrix;

//gsl_matrix_set (m, 0, 0, -sp.a*y[1]-sp.b);
//gsl_matrix_set (m, 0, 1, -sp.a*y[0]+sp.c);

//gsl_matrix_set (m, 1, 0, sp.a*y[1]);
//gsl_matrix_set (m, 1, 1, sp.a*y[0]-sp.b-sp.c);

//dfdt[0] = 0.0;
//dfdt[1] = 0.0;
return NULL;//GSL_SUCCESS;
}
int
jacm (double t, const double ym0[], double *dfmdy,
double dfmdt[], void *params)
{

//double mu = *(double *)params;
//sys_para sp=*(sys_para*)params;
//gsl_matrix_view dfdy_mat
//= gsl_matrix_view_array (dfdy, 2, 2);
//gsl_matrix * m = &dfdy_mat.matrix;

//gsl_matrix_set (m, 0, 0, -sp.a*y[1]-sp.b);
//gsl_matrix_set (m, 0, 1, -sp.a*y[0]+sp.c);

//gsl_matrix_set (m, 1, 0, sp.a*y[1]);
//gsl_matrix_set (m, 1, 1, sp.a*y[0]-sp.b-sp.c);

//dfdt[0] = 0.0;
//dfdt[1] = 0.0;
return NULL;//GSL_SUCCESS;
}


double BRN(sys_para &sp,double lam)
{
int i,j,l;
double nrm,maxnorm,max,min;
const gsl_odeiv_step_type * T
= gsl_odeiv_step_rkck;

gsl_matrix *m1=gsl_matrix_alloc(2,2);

double t = 0.0, tt = period;
double h = 0.0031,hmin=1e-3;
double ym0[idm+1]= {0,0.1, 0.2,0.3,0.4};
gsl_odeiv_system sysm = {monod, jacm, idm, &sp};

h = 0.0031;hmin=1e-3;
gsl_odeiv_step * sm
= gsl_odeiv_step_alloc (T, idm+1);
gsl_odeiv_control * cm
= gsl_odeiv_control_y_new (3e-14, 0.0);
gsl_odeiv_evolve * em
= gsl_odeiv_evolve_alloc (idm+1);
t=0;

sp.lam=lam;
sp.PIM01=sp.PIM0;
ym0[1]=1;ym0[2]=0;
ym0[3]=0;ym0[4]=1;
min=0;max=1;
while (t <tt)
{
sp.PIM=sp.PIM0*(1+sp.d*sin(sp.tpi*t/365));
if(sp.d>0)
{
sp.PIM1=sp.PIM01*(1+sp.d*sin(sp.tpi*t/365));
if(sp.PIM1>=sp.PIM01*sqrt(1+(sp.d*sp.d)/2.0) && sp.con==1)//
sp.PIM=sp.PIM0*sp.fracPIM*(1+sp.d*sin(sp.tpi*t/365));
}

if(sp.d==0)
{
sp.PIM1=sp.PIM01*(1+sin(sp.tpi*t/365));
if(sp.PIM1>=sp.PIM01*sqrt(1+0.5) && sp.con==1)//
sp.PIM=sp.PIM0*sp.fracPIM*(1+sp.d*sin(sp.tpi*t/365));
}

max=sp.PIM;
int status = gsl_odeiv_evolve_apply (em, cm, sm,
&sysm,
&t, tt,
&h, ym0);
}
gsl_vector_complex *eval=gsl_vector_complex_alloc(2);
gsl_matrix_complex *evec=gsl_matrix_complex_alloc(2,2);
gsl_eigen_nonsymmv_workspace *ww=gsl_eigen_nonsymmv_alloc(2);


gsl_matrix_set(m1,0,0,ym0[1]);gsl_matrix_set(m1,0,1,ym0[2]);
gsl_matrix_set(m1,1,0,ym0[3]);gsl_matrix_set(m1,1,1,ym0[4]);

gsl_eigen_nonsymmv(m1,eval,evec,ww);

gsl_vector *norm=gsl_vector_alloc(idm/2.0);
gsl_permutation *p=gsl_permutation_alloc(idm/2.0);
//printf("%lf ",sp.lam);
for(l=0;l<idm/2.0;l++)
{
gsl_complex eval_0=gsl_vector_complex_get(eval,l);
nrm=sqrt((GSL_REAL(eval_0)*GSL_REAL(eval_0))+(GSL_IMAG(eval_0)*GSL_IMAG(eval_0)));
  gsl_vector_set(norm,l,nrm);
//printf("%lf ",nrm);
}
//printf("\n");
gsl_sort_vector_index(p,norm);

  l=gsl_permutation_get(p,(idm/2.0)-1);
maxnorm=gsl_vector_get(norm,l);

gsl_odeiv_evolve_free (em);
gsl_odeiv_control_free (cm);
gsl_odeiv_step_free (sm);

gsl_eigen_nonsymmv_free(ww);
return (maxnorm-1);
}


void root(sys_para &sp)
{
double yl,yr,ym;
double acc=1e-9;
double Llam,Rlam,red;
red=0.00019804;
Llam=sp.R02+red;
Rlam=sp.R02-pow(10,-1)*.5;
{
while(BRN(sp,Llam)<0)
{
Llam/=1.0001;
}
Rlam=Llam*1.0001;


yl=Llam;
yr=Rlam;
do
{

  ym=(yl+yr)/2.;
if(BRN(sp,ym)>0)
{
yl=ym;
yr=yr;
}
if(BRN(sp,ym)<0)
{
yr=ym;
yl=yl;
}
 }
 while(fabs((yr-yl)/(yr+yl))>=acc);
sp.lam=ym;
}

//if(sp.d<=0.0000000000)
//sp.lam=sp.R02;
}

double PIM0fn(sys_para &sp,double PIM0)
{
double fx;
sp.PIM0=PIM0;
sp.b1bet1=sp.b1*sp.bet1;
sp.b1bet2=sp.b1*sp.bet2;
sp.R02=sqrt((sp.b1bet1*sp.b1bet2*sp.muB*sp.PIM0)/(sp.muM*sp.muM*(sp.muB+sp.dB)*sp.PIB));
root(sp);

sp.R01=sp.lam;
fx=sp.R01-sp.R0d;

return (fx);
}

double b10fn(sys_para &sp,double b10)
{
double fx;
sp.b1=b10;
sp.b1bet1=sp.b1*sp.bet1;
sp.b1bet2=sp.b1*sp.bet2;
sp.R02=sqrt((sp.b1*sp.bet1*sp.b1*sp.bet2*sp.muB*sp.PIM0)/(sp.muM*sp.muM*(sp.muB+sp.dB)*sp.PIB));
root(sp);

sp.R01=sp.lam;
fx=sp.R01-sp.R0d;

return (fx);
}



void rootPIM0(sys_para &sp)
{
double yl,yr,ym;
double acc=1e-9;
double LPIM,RPIM;
LPIM=5;
RPIM=80000000;
if(PIM0fn(sp,LPIM)*PIM0fn(sp,RPIM)<0)
{

if(PIM0fn(sp,LPIM)<0)
{
yl=LPIM;
yr=RPIM;
}
if(PIM0fn(sp,LPIM)>0)
{
yl=RPIM;
yr=LPIM;
}
do
{
  ym=(yl+yr)/2.;
if(PIM0fn(sp,ym)<0)
{
yl=ym;
yr=yr;
}
if(PIM0fn(sp,ym)>0)
{
yr=ym;
yl=yl;
}
 }
 while(fabs((yr-yl)/(yr+yl))>=acc);
}
sp.PIM0=ym;
}

void rootb10(sys_para &sp)
{
double yl,yr,ym;
double acc=1e-9;
double Lb1,Rb1;
Lb1=0.001;
Rb1=1;
if(b10fn(sp,Lb1)*b10fn(sp,Rb1)<0)
{

if(b10fn(sp,Lb1)<0)
{
yl=Lb1;
yr=Rb1;
}
if(b10fn(sp,Lb1)>0)
{
yl=Rb1;
yr=Lb1;
}
do
{
  ym=(yl+yr)/2.;
if(b10fn(sp,ym)<0)
{
yl=ym;
yr=yr;
}
if(b10fn(sp,ym)>0)
{
yr=ym;
yl=yl;
}
 }
 while(fabs((yr-yl)/(yr+yl))>=acc);
}
sp.b1=ym;
}





void gs(double *y,double *sum)
{
  double a[M+1];
  double bnorm;
  int ii,itt,it,iii;

  for(ii=1;ii<=M;ii++){
    a[ii]=0.0;
  }

  for(itt=1;itt<=M;itt++){
    bnorm=0.0;
    if(itt==1){
      for(it=1;it<=id;it++){
    bnorm=bnorm+pow(y[it+id],2.0);
      }
      for(it=1;it<=id;it++){
    y[it+id]=y[it+id]/sqrt(bnorm);
      }
      sum[itt]=sum[itt]+0.5*log(bnorm);
      bnorm=0.0;
    }
    else{
      for(iii=1;iii<=(itt-1);iii++){
    a[iii]=0.0;
    for(it=1;it<=id;it++){
      a[iii]=a[iii]+y[it+iii*id]*y[it+itt*id];
    }
    for(it=1;it<=id;it++){
      y[it+itt*id]=y[it+itt*id]-a[iii]*y[it+iii*id];
        }
      }
      for(it=1;it<=id;it++){
    bnorm=bnorm+pow(y[it+itt*id],2.0);
      }
      for(it=1;it<=id;it++){
    y[it+itt*id]=y[it+itt*id]/sqrt(bnorm);
      }
      sum[itt]=sum[itt]+0.5*log(bnorm);
      bnorm=0.0;
    }}}



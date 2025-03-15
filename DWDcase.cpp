#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <math.h>
#include <fstream>
#include <iostream>
#include <sys/time.h>
#include <sys/timeb.h>
#include <cmath>
#include "VBBinaryLensingLibrary.h"
using std::cout;
using std::endl;
using std::cin;
///=============================================================================
const double RA=180.0/M_PI;
const double KP=3.08568025*pow(10.,19); // in meter.
const double G= 6.67384*pow(10.,-11.0);// in [m^3/s^2*kg].
const double velocity=299792458.0;//velosity of light  m/s
const double Msun=1.98892*pow(10.,30.0); //in [kg].
const double Mjupiter=1.898*pow(10,27.0); 
const double Mearth=  5.9722*pow(10.0,24.0);
const double AU=1.495978707*pow(10.0,11.0);
const double Rsun=6.9634*pow(10.0,8.0); ///solar radius [meter]
const double Teffsun= 5772.0;//kelvin
const double Avks=double(8.20922);
const double cadence=double(2.0/60.0/24.0);  
const double KBol=1.380649*pow(10.0,-23.0);
const double Hplank=6.62607015*pow(10.0,-34.0);
const double const1= sqrt(4.0*G*Msun)/velocity; 
const double const2= Hplank*velocity*pow(10.0,9.0)/KBol; 


///=============================================================================
const int    Nli=int(798);  
const int    NB=int(1000); 
const int    nbh=int(30);  
const int    Nl=110; /// wavelength TESS throughput 
const int    nx=1051; 
const int    ny=1051;  
const double thre=0.001; 
const double wave[4]= {0.673,0.7865,0.532,0.797};//https://en.wikipedia.org/wiki/Photometric_system  G, T, BP, RP
const double AlAv[4]= {0.791986524645539, 0.617155245862836  ,1.0386670640894,0.601810722049874};
const double sigma[4]={0.017, 0.02,0.02, 0.02};// G, T, BP, RP  Table (2) Cardeli
///=============================================================================
struct source{
    int    num;  
    double Ds;
    double nsbl, blend, magb, Ai, Mab, Map;
    double Teff, Rstar, mass, Logg;
    double ros, limb, grav, lon, lat;
    double magG, magBP, magRP;
    double cdp, Lumi;  
};
struct limbD{
    double Logg[Nli],  Tef[Nli],  Limb[Nli],  grav[Nli];  
};

struct lens{
    int    num;  
    double ecen, inc, Lumi, tet, tp, period, a; 
    double phi, RE, limb, grav;   
    double magG, magBP, magRP;
    double ratio, q, MBH, RBH, Map, Mab, Teff, Logg;
    double dx,dy,xi,yi,xsi, ysi, xsc, ysc, num0, num1, Dls, Dl;
    int    flag;  
};
struct doppler{
  double wave[Nl], throu[Nl]; 
  double waven, Fpl0, Fpl1;  
};
///=============================================================================
double ErrorTESS(double maga); 
double RandN(double sigma,double);
double RandR(double down, double up);
double Fluxlimb(double limb, double rstar);  
double Kepler(double phi, double ecen); 
double Bessel(int n,double x); 
double CDPPM(source & s, double); 
double DopplerF(doppler & dp, double , double); 
double Fplank(double , double);  
double Ellipsoid(double, double, double, double, double, double);
void   LensEq2(source & s, lens & l, double, double );
void   FiniteLens(source & s, lens & l, double , double , double, double);
int LimbF(limbD & li, double, double); 

time_t  _timeNow;
unsigned int _randVal;
unsigned int _dummyVal;
FILE * _randStream;
///===========================================================================//
///                                                                           //
///                  Main program                                             //
///                                                                           //
///===========================================================================//	
int main()
{


   time(&_timeNow);
   _randStream = fopen("/dev/urandom", "r");
   _dummyVal = fread(&_randVal, sizeof(_randVal), 1, _randStream);
   srand(_randVal);
   time( &_timeNow);
   printf("START time:   %s",ctime(&_timeNow));
      
   
   VBBinaryLensing vbb;
   vbb.Tol=1.e-5;
   vbb.a1 =0.0;  
   vbb.LoadESPLTable("./ESPL.tbl");
  
  
   source s;
   lens l;
   doppler dp;
   limbD li;  
   FILE* film;  
   FILE* distr;
   FILE* dopp;
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
   char   filenam0[40], filenam2[40];
   double ksi, x0, y0, x1, y1, z1, dis, b1,b2; 
   double phase, u, us, Astar, As, disp, dt, dt0;   
   double Rho, RE, proj, tim, finl, SourceA, Astar2, Roche;
   double cosi2, li3, ampl, gc, L0, Delf2;
   double Dist1, Av, Atot;  
   int    self, nstep, Flag, Nl, Ns;    
   int    dwd=6;  
   nstep=0;
   
   
  
   dopp=fopen("./files/TESS_Throught.txt","r");
   for(int i=0;  i<Nl;  ++i){
   fscanf(dopp,"%lf  %lf\n",&dp.wave[i], &dp.throu[i]);}
   fclose(dopp); 
   cout<<"********** File TESS_Throught.txt was read ************"<<endl;    
   
   
   dopp=fopen("./files/WDCoefC2020.dat","r");
   for(int i=0; i<Nli;  ++i){
   fscanf(dopp,"%lf  %lf  %lf   %lf\n",&li.Logg[i], &li.Tef[i], &li.Limb[i], &li.grav[i]);}
   fclose(dopp);     
   cout<<"********** File TESS_Limb_gravity_coefficient.txt was read ************"<<endl;    
   
 
   sprintf(filenam0,"./files/light/lcF0/files/%c%c%d.dat",'C','_',1);
   distr=fopen(filenam0,"a+"); 
   fclose(distr);              
   sprintf(filenam2,"./files/light/lcF0/files/%c%c%d.dat",'M','_',dwd);
   film=fopen(filenam2,"a+");  


///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH  
if(dwd==1){//J1717+6757_TESS_ID:219868627  
l.MBH=0.9;//0.185;
l.RBH=0.0085;//0.1;
l.Logg=8.22;//5.67;
l.Teff=15500.0;//14900.0;
s.mass=0.185;///0.9;
s.Rstar=0.1;//0.0085; 
s.Logg=5.67;// 8.22;
s.Teff=14900.0;///15500.0; 
s.lon=98.475950;
s.lat=33.796074;
l.period=0.246137; 
l.inc=fabs(90.0-86.75)*M_PI/180.0;    
Dist1= 178.61022*0.001;}  
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
if(dwd==2){//J1557+2823  TESS_ID: 11015992282
l.MBH=0.52;//0.461;
l.RBH=0.01125*sqrt(pow(l.MBH/1.454,-2.0/3.0)-pow(l.MBH/1.454,2.0/3.0));///#[Rsun];//0.0147;
l.Logg=7.9;//7.762;
l.Teff=11000.0;//12560.0;
s.mass=0.461;//0.52; 
s.Rstar=0.0147;// 0.01125*sqrt(pow(s.mass/1.454,-2.0/3.0)-pow(s.mass/1.454,2.0/3.0));///#[Rsun]
s.Logg=7.762;//7.9;
s.Teff=12560.0;//11000.0; 
s.lon=45.8954759;
s.lat=49.1562009;
l.period=0.40741; 
l.inc=fabs(90.0-60.0)*M_PI/180.0;    
Dist1= 247.013151719299*0.001;}
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
if(dwd==3){//LP400-22;  TESS_ID: 2002564035  
l.MBH=0.48;//0.19;
l.RBH=0.01125*sqrt(pow(l.MBH/1.454,-2.0/3.0)-pow(l.MBH/1.454,2.0/3.0));
l.Logg=7.8;// 6.42;
l.Teff=10000.0;//11140;
s.mass=0.19;//0.48;
s.Rstar=0.01125*sqrt(pow(s.mass/1.454,-2.0/3.0)-pow(s.mass/1.454,2.0/3.0));
s.Logg=6.42;//7.8;
s.Teff=11140.0;//10000.0;
s.lon=86.3614383;
s.lat=-30.5842598;
l.period=1.01016;
l.inc=fabs(90.0-83.0)*M_PI/180.0;    
Dist1= 365.796395381606*0.001;}
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
if(dwd==4){//J1449+1717; TESS_ID: 1101113865
l.MBH=0.81;//0.168;
l.RBH=0.01125*sqrt(pow(l.MBH/1.454,-2.0/3.0)-pow(l.MBH/1.454,2.0/3.0));
l.Logg=7.9;// 6.08;
l.Teff=10000.0;//9700.0;
s.mass=0.168;//0.81;
s.Rstar=0.01125*sqrt(pow(s.mass/1.454,-2.0/3.0)-pow(s.mass/1.454,2.0/3.0));
s.Logg=6.08;//7.9;
s.Teff=9700.0;//10000.0;
s.lon=19.3641738;
s.lat= 60.9492427;
l.period=0.29075;//days 
l.inc=fabs(90.0-60.0)*M_PI/180.0;    
Dist1= 613.386335188751*0.001;}
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
if(dwd==5){//J2132+0754;  TESS_ID: 2000073295
l.MBH=1.33;//0.17;
l.RBH=0.01125*sqrt(pow(l.MBH/1.454,-2.0/3.0)-pow(l.MBH/1.454,2.0/3.0));
l.Logg=9.2;// 5.995;
l.Teff=15000.0;// 13700.0;
s.mass=0.17;// 1.33; //That should be more massive WD
s.Rstar=0.01125*sqrt(pow(s.mass/1.454,-2.0/3.0)-pow(s.mass/1.454,2.0/3.0));
s.Logg=5.995;// 9.2;
s.Teff=13700.0;//15000.0;
s.lon=61.6602807;
s.lat=-30.4603525;
l.period=0.25056;//days 
l.inc=fabs(90.0-60.0)*M_PI/180.0;    
Dist1= 1221.31745*0.001;}
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
if(dwd==6){//J2151+1614;  TESS_ID: 2000525780
l.MBH=0.66;//0.176;
l.RBH=0.01125*sqrt(pow(l.MBH/1.454,-2.0/3.0)-pow(l.MBH/1.454,2.0/3.0));///#[Rsun];
l.Logg=8.0;// 6.24;
l.Teff=11000.0;// 12300.0;
s.mass=0.176;// 0.66;
s.Rstar=0.01125*sqrt(pow(s.mass/1.454,-2.0/3.0)-pow(s.mass/1.454,2.0/3.0));///#[Rsun]
s.Logg=6.24;// 8.0;
s.Teff=12300.0;//11000.0;
s.lon= 72.4831190;
s.lat=-28.5585599;
l.period=0.59152;//day
l.inc=fabs(90.0-60.0)*M_PI/180.0;    
Dist1=391.21190*0.001;}
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
   
   l.tet=0.0;    
   l.tp=-0.2*l.period;    
   l.ecen=0.0; 
   s.blend=1.0; 
   Ns=LimbF(li, s.Teff, s.Logg);
   s.limb=li.Limb[Ns]; 
   s.grav=li.grav[Ns];
   Nl=LimbF(li, l.Teff, l.Logg);
   l.limb=li.Limb[Nl]; 
   l.grav=li.grav[Nl];
   dt=l.period/1500.0;
   l.dx=0.0915563545; 
   l.dy=0.0915563545;
   dt0=dt; 
   Flag=-1;
   l.q=double(l.MBH/s.mass);  
   l.Lumi=1.0;
   s.Lumi=pow(s.Rstar/l.RBH,2.0)*pow(s.Teff/l.Teff,4.0); 
   l.ratio=l.Lumi/s.Lumi; 
   b1=l.ratio/(1.0+l.ratio);//for lens  
   b2=    1.0/(1.0+l.ratio);//for source 
   l.a= pow(pow(l.period*24.0*3600.0,2.0)*G*Msun*(s.mass+l.MBH)/(4.0*M_PI*M_PI),1.0/3.0 );//meter
   Roche=l.a-l.a*0.49*pow(l.q,2.0/3.0)/(0.6*pow(l.q,2.0/3.0)+log(1.0+pow(l.q,1.0/3.0)));
   
   
   cout<<dwd<<"\t\t *********************************************"<<endl;   
   cout<<"Mass_Lens:   "<<l.MBH<<"\t Mass_source:   "<<s.mass<<endl;
   cout<<"Radius_Lens: "<<l.RBH<<"\t Radius_source: "<<s.Rstar<<endl;
   cout<<"Teff_lens:   "<<l.Teff<<"\tTeff_source:   "<<s.Teff<<endl;
   cout<<"logg_lens:   "<<l.Logg<<"\tLogg_source:   "<<s.Logg<<endl;
   cout<<"limb_lens:   "<<l.limb<<"\tlimb_source:   "<<s.limb<<endl;
   cout<<"grav_lens:   "<<l.grav<<"\t grav_source:  "<<s.grav<<endl;
   cout<<"\nRoche/Radius_low_mass_WD:    "<<Roche/l.RBH<<endl;
   cout<<"low_mass_Radius/Rsun: "<<l.RBH<<"\t Roche/Rsun:  "<<Roche/Rsun<<endl;
   cout<<"Roche/RBH:   "<<Roche/(l.RBH*Rsun)<<endl;
   cout<<"Semi/Rsun:  "<<l.a/Rsun<<"\t ratio:  "<<l.ratio<<endl;
   cout<<"ratio: "<<l.ratio<<"\t b_lens: "<<b1<<"\t b_source: "<<b2<<endl;
   int uue; cin>>uue;   
   cout<<"*******************************************************"<<endl;    
   
   distr=fopen(filenam0,"a+"); 
   fprintf(distr,"%d %.4lf %.4lf  %.4lf  %.4lf %.4lf %.4lf %.1lf  %.3lf %.4lf %.4lf %.4lf %.1lf %.3lf %.6lf  %.4lf  %.4lf %.4lf %.2lf\n",
   dwd,s.lat,s.lon,Dist1,l.MBH,l.RBH,l.Logg,l.Teff,l.limb,
   s.mass,s.Rstar,s.Logg,s.Teff,s.limb,l.period,l.a/Rsun,l.ratio,Roche/Rsun/l.RBH,double(90.0-l.inc*180.0/M_PI));//19
   fclose(distr);   
   
   
   
   for(tim=0.0; tim<=l.period; tim+=dt0){
   l.phi=double((tim-l.tp)*2.0*M_PI/l.period); 
   if(l.ecen<0.01) ksi=l.phi;
   else            ksi=Kepler(l.phi,l.ecen);
   x0=l.a*(cos(ksi)-l.ecen);//major axis[m] source with respect to the lens object
   y0=l.a*sin(ksi)*sqrt(1.0-l.ecen*l.ecen);//minor axis [m]
   y1=              y0*cos(l.tet)+x0*sin(l.tet);//[m]
   x1= cos(l.inc)*(-y0*sin(l.tet)+x0*cos(l.tet));//[m]
   z1=-sin(l.inc)*(-y0*sin(l.tet)+x0*cos(l.tet));//[m]
   dis= sqrt(x1*x1 + y1*y1 + z1*z1)+1.0e-50;//[m] 
   disp=sqrt(y1*y1 + z1*z1)+1.0e-50;///meter
   phase=acos(-x1/dis);//radian
   Astar=Astar2=Atot=1.0;
   l.num0=l.num1=finl=0.0;
   self=0; 

   
   
   
   //Ellipsiodal_Variation_Source_star 
   gc=s.grav;//gravity-darkening 
   cosi2=cos(l.inc)*cos(l.inc);
   li3= double(3.0-s.limb);
   ampl=double(s.Rstar*Rsun/dis);
   L0=fabs(1.0+(15.0+s.limb)*(1.0+gc)*pow(ampl,3.0)*(2.0+5.0*l.q)*(2.0-3.0*cosi2)/(60.0*(3.0-s.limb)) +
       9.0*(1.0-s.limb)*(3.0+gc)*pow(ampl,5.0)*l.q*(8.0-40.0*cosi2+35.0*cosi2*cosi2)/(256.0*(3.0-s.limb)));
   Delf2=double((15*s.limb*(2.0+gc)*pow(ampl,4.0)*l.q*(4.0*cos(l.inc)-5.0*pow(cos(l.inc),3))/(32.0*li3))*cos(phase) + 
   (-3.0*(15.0+s.limb)*(1.0+gc)*pow(ampl,3.0)*l.q*cosi2/(20.0*li3)-15.0*(1.0-s.limb)*(3.0+gc)*pow(ampl,5.0)*l.q*(6.0*cosi2-7.0*cosi2*cosi2)/(64.0*li3))*cos(2.0*phase)+ 
   (-25.0*s.limb*(2.0+gc)*pow(ampl,4.0)*l.q*pow(cos(l.inc),3.0)/(32.0*li3))*cos(3.0*phase) + 
   (105.0*(1.0-s.limb)*(3.0+gc)*pow(ampl,5.0)*l.q*pow(cos(l.inc),4.0)/(256.0*li3) )*cos(4.0*phase))/L0;//with zero average
   cout<<"L0: "<<L0<<"\t Delf2:  "<<Delf2<<endl;



   /*Doppler effect 
    DF0=0.0;  
    DF1=0.0;  
    for(int i=0; i<int(Nl-1); ++i){
    dw=double(dp.wave[i+1]-dp.wave[i])*pow(10.0,-9);
    wave1= dp.wave[i]*pow(10.0,-9.0);              
    wave2=wave1*(1.0+vx/velocity); 
    //wave1=wave1*pow(10.0,-9.0);//#[m]  
    con=Hplank*velocity/(KBol*Tstar*wave1); 
    Fplank1=2.0*Hplank*velocity*velocity*pow(wave1,-5.0)/(exp(con)-1.0);
    con=Hplank*velocity/(KBol*Tstar*wave2); 
    Fplank2= double(2.0*Hplank*velocity*velocity*pow(wave2,-5.0)/(exp(con)-1.0)  ) ;
    Delf1[m]+= double(Fplank2 - Fplank1)*throu[i,1]*dw ;      
    DF0+=     Fplank1 * throu[i,1] * dw; } 
    Delf1[m]=float(Delf1[m]/DF0) ;*/    


   
        
   if(x1<-0.0001){//Lensing of source star by masive lens object
   l.Dl=Dist1*KP;//meter
   s.Ds=Dist1*KP+fabs(x1);//meter
   l.Dls=fabs(x1);//[m]
   proj=double(l.Dl/s.Ds);
   l.RE=sqrt(4.0*G*Msun*l.MBH)*sqrt(l.Dls*proj)/velocity+1.0e-50;
   s.ros=fabs(s.Rstar*Rsun*proj/l.RE); 
   SourceA=double(M_PI*s.ros*s.ros*(1.0-s.limb/3.0)/(l.dx*l.dy));
   u=fabs(disp/l.RE);
   l.xsc=double(y1/l.RE);  
   l.ysc=double(z1/l.RE);
   if(u<float(15.0*s.ros)){
   self=1;
   if(Flag!=1){dt0=dt/35.0;  Flag=1;   if(dt0>cadence) dt0=double(cadence);} 
   if(s.ros>100.0){//Valerio_Magnification
   if(u<s.ros){Astar=double(1.0+2.0/s.ros/s.ros);        }
   else       {Astar=double(u*u+2.0)/sqrt(u*u*(u*u+4.0));}}
   else{vbb.a1=s.limb;  Astar=vbb.ESPLMag2(u,s.ros);}
   FiniteLens(s,l, double(l.MBH), double(l.RBH), double(s.limb), double(u) );
   finl=   double(l.num1/SourceA);
   Astar2= double(l.num0/SourceA);}
   //us=double(disp-l.RBH*Rsun-s.Rstar*Rsun*proj);
   //Occul=1.0; 
   //if(us<=0.0){
   //if(disp<=fabs(s.Rstar*Rsun*s.Ds*KP/(s.Ds*KP-x1)-l.RBH*Rsun)){frac=1.0;  frac0=1.0;}
   //else{   
   //frac=0.0;  frac0=0.0;    
   //for(int i=0; i<nbh; ++i){
   //for(int j=0; j<nbh; ++j){
   //yb=double(-l.RBH*Rsun + i*l.stepb);   
   //zb=double(-l.RBH*Rsun + j*l.stepb); 
   //zlim= sqrt(l.RBH*l.RBH*Rsun*Rsun-yb*yb); 
   //if(fabs(zb)<=zlim){
   //yc=y1-yb;    
   //zc=z1-zb;  
   //rstar=sqrt(yc*yc+zc*zc)/(s.Rstar*Rsun*s.Ds*KP/(s.Ds*KP-x1)); 
   //frac0+=1.0; 
   //if(rstar<=1.0) frac+=1.0;
   //}}}}   
   //Occul=double(1.0-frac/frac0);}
   Atot=(Astar-finl+Delf2)*b2 + b1;
   fprintf(film,"%.11lf   %.9lf   %.9lf  %.9lf  %.7lf  %.8lf  %.8lf  %.8lf %.8lf   %d\n",//10
   tim, Astar*b2+b1, Astar2*b2+b1, (1.0-finl)*b2+b1, l.RE/Rsun, s.ros, u/s.ros, Atot, (1.0+Delf2)*b2+b1, self);}
   
   
   
  
   else{//Lensing of Lens by source star_x1>0.0
   l.Dl=Dist1*KP-fabs(x1);//[m]
   s.Ds=Dist1*KP;//[m]
   l.Dls=fabs(x1);//[m]
   proj=double(l.Dl/s.Ds);
   l.RE=sqrt(4.0*G*Msun*s.mass)*sqrt(l.Dls*proj)/velocity+1.0e-50;
   s.ros=fabs(l.RBH*Rsun*proj/l.RE); 
   SourceA=double(M_PI*s.ros*s.ros*(1.0-l.limb/3.0)/(l.dx*l.dy));
   u=fabs(disp/l.RE);
   l.xsc=double(-y1/l.RE);  
   l.ysc=double(-z1/l.RE);
   if(u<float(15.0*s.ros)){
   self=2;
   if(Flag!=2){dt0=dt/35.0; Flag=2;   if(dt0>cadence) dt0=double(cadence);}  
   if(s.ros>100.0){//Valerio_Magnification
   if(u<s.ros){Astar=double(1.0+2.0/s.ros/s.ros);        }
   else       {Astar=double(u*u+2.0)/sqrt(u*u*(u*u+4.0));}}
   else{vbb.a1=l.limb; Astar=vbb.ESPLMag2(u,s.ros);}
   FiniteLens(s,l,double(s.mass),double(s.Rstar),double(l.limb),double(u));
   finl=   double(l.num1/SourceA);
   Astar2= double(l.num0/SourceA);}
   //us=double(disp-l.RBH*Rsun-s.Rstar*Rsun*proj);
   //Occul=1.0; 
   //if(us<=0.0){
   //if(disp<=fabs(s.Rstar*Rsun*s.Ds*KP/(s.Ds*KP-x1)-l.RBH*Rsun)){frac=1.0;  frac0=1.0;}
   //else{   
   //frac=0.0;  frac0=0.0;    
   //for(int i=0; i<nbh; ++i){
   //for(int j=0; j<nbh; ++j){
   //yb=double(-l.RBH*Rsun + i*l.stepb);   
   //zb=double(-l.RBH*Rsun + j*l.stepb); 
   //zlim= sqrt(l.RBH*l.RBH*Rsun*Rsun-yb*yb); 
   //if(fabs(zb)<=zlim){
   //yc=y1-yb;    
   //zc=z1-zb;  
   //rstar=sqrt(yc*yc+zc*zc)/(s.Rstar*Rsun*s.Ds*KP/(s.Ds*KP-x1)); 
   //frac0+=1.0; 
   //if(rstar<=1.0) frac+=1.0;
   //}}}}   
   //Occul=double(1.0-frac/frac0);}
   Atot=(Astar-finl)*b1+(1.0+Delf2)*b2;
   //Atot= double((Astar-finl+Delf2)*l.Lumi+s.Lumi)/(s.Lumi+l.Lumi);
   fprintf(film,"%.11lf   %.9lf   %.9lf  %.9lf  %.7lf  %.8lf  %.8lf  %.8lf %.8lf   %d\n",//10
   tim, Astar*b1+b2, Astar2*b1+b2, (1.0-finl)*b1+b2, l.RE/Rsun, s.ros, u/s.ros, Atot, (1.0+Delf2)*b2+b1 ,self);}
   if(self==0) dt0=dt; 
   nstep+=1;
   
   
   
   
   if(nstep%1==0){  
   cout<<"nstep:  "<<nstep<<"\t dt0: "<<dt0/dt<<"\t period  "<<l.period<<endl;
   cout<<"self: "<<self<<"\t ro_star: "<<s.ros<<"\t u/ros: "<<u/s.ros<<endl;
   cout<<"A_VBB:  "<<Astar<<"\t A_IRS: "<<Astar2<<"\t FinL: "<<finl<<endl;
   cout<<"num0:  "<<l.num0<<"\t num1: "<<l.num1<<"\t sourceA: "<<SourceA<<endl;
   cout<<"u:  "<<u<<"\t RE:   "<<l.RE/(s.Rstar*Rsun)<<"\t inc(deg):  "<<l.inc*RA<<endl;
   cout<<"semi/Rsun: "<<l.a/Rsun<<"\t tim/period: "<<tim/l.period<<"\t x1: "<<x1/l.a<<endl;
   cout<<"*********************************************************"<<endl;
   if((x1<-0.01 and fabs(phase)>M_PI/2.0) or l.RE<-0.01 or s.ros<-0.01 or dis<-0.01 or x1>dis or disp<-0.01 or 
   phase<-0.00001 or phase>float(M_PI*1.01) or s.limb>1.01 or Astar<0.9 or l.num0<0.0 or s.limb<0.0 or 
   l.limb<0.0 or l.limb>1.01){ 
   cout<<"ERROx1:  "<<x1<<"\t phase: "<<phase<<"\t ros:  "<<s.ros<<endl;
   cout<<"RE:    "<<RE<<"\t dis:  "<<dis<<"disp:  "<<disp<<endl;
   cout<<"u:  "<<u<<"\t Astar:  "<<Astar<<"\t limb:  "<<s.limb<<"\t l.limb: "<<l.limb<<endl;}}
   }//time loop
   fclose(film);    
  
   cout<<"=============================================================="<<endl;
   cout<<"latit:  "<<s.lat<<"\t longt: "<<s.lon<<"\t nstep:  "<<nstep<<endl;
   cout<<"l.dx:   "<<l.dx<<"\t l.dy:  "<<l.dy<<endl;
   cout<<"s.Lumi: "<<s.Lumi<<"\t l.Lumi:  "<<l.Lumi<<"\t Ds:  "<<Dist1<<endl;
   cout<<"inc(deg):    "<<l.inc*RA<<"\t  teta(deg):  "<<l.tet*RA<<"\t dt:  "<<dt<<endl;
   cout<<"l.Mab:  "<<l.Mab<<"\t l.Map: "<<l.Map<<"\t ratio:  "<<l.ratio<<endl;
   cout<<"s.Mab:  "<<s.Mab<<"\t s.Map: "<<s.Map<<"\t period: "<<l.period<<endl;
   cout<<"l.MBH:  "<<l.MBH<<"\t l.RBH: "<<l.RBH<<"\t semi/Rsun: "<<l.a/Rsun<<endl;
   cout<<"mass:   "<<s.mass<<"\t s.Rstar:  "<<s.Rstar<<"\t blend:  "<<s.blend<<endl;
   cout<<"l.a/AU: "<<l.a/AU<<"\t Ratio:  "<<l.ratio<<"\t ecen:  "<<l.ecen<<endl;
   cout<<"=============================================================="<<endl;
   fclose(_randStream);
   return(0);
}
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
int LimbF(limbD & li, double tstar, double logg){
    int num=-1,num2=-1,i;double distm, dist;
    
    if(logg>li.Logg[Nli-1])   num=int(Nli-1); 
    else if(logg<li.Logg[0])  num=0;  
    else{
    for(int j=1; j<Nli; ++j){
    if(double((logg-li.Logg[j])*(logg-li.Logg[j-1]))<=0.0){num=int(j-1); break;}}}

    distm=100000000.0; 
    for(int j=-40; j<40; ++j){
    i=int(num+j);
    if(i>=0 and i<Nli){
    dist=sqrt(pow(tstar-li.Tef[i],2.0)+pow(logg-li.Logg[i],2.0)); 
    if(dist<distm){distm=dist;  num2=i; }}}
    
    
    if(num<0 or num2<0 or num>=Nli or num2>=Nli or li.Limb[num2]<0.0 or 
    li.Limb[num2]>1.0 or fabs(tstar-li.Tef[num2])>2500){
    cout<<"Error num: "<<num<<"\t  num2:  "<<num2<<endl;  
    cout<<"TeffC:    "<<li.Tef[num2]<<"\tLoggC:  "<<li.Logg[num2]<<endl;
    cout<<"Tstar:   "<<tstar<<"\t logg:  "<<logg<<endl;  
    int yye;  cin>>yye;}
    
    cout<<"Tstar:   "<<tstar<<"\t       logg:  "<<logg<<endl;
    cout<<"TeffC:   "<<li.Tef[num2]<<"\tLoggC: "<<li.Logg[num2]<<endl;
    cout<<"limb_Drakening:   "<<li.Limb[num2]<<"\tnum2: "<<num2<<endl;
    cout<<"****************************************************"<<endl;
    return(int(num2));  
}
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
double Kepler(double phi, double ecen){  
    double ksi=0;   
    double term, term0;  
    phi=double(phi*RA); 
    while(phi>360.0) phi=phi-360.0; 
    while(phi<0.0)   phi=phi+360.0;      
    if(phi>180)      phi=double(phi-360.0);
    if(phi<-181.0 or phi>181.0){ 
    cout<<"Error :  Phi:  "<<phi<<"\t ecent:  "<<ecen<<endl;   int yye;  cin>>yye;}
    phi=double(phi/RA);
    ksi=phi; 
    for(int i=1; i<NB; ++i){
    term= Bessel(i,i*ecen)*sin(i*phi)*2.0/i;  
    ksi+=term; 
    if(i==1) term0=fabs(term); 
    if(fabs(term)<double(thre*term0) and i>5)  break;}        
    return(ksi); 
}    
///#############################################################################
double ErrorTESS(double maga){
   double emt=-1.0, m;     
   
   if(maga<7.5)        emt=double(0.22*maga-5.850); 
   else if(maga<12.5)  emt=double(0.27*maga-6.225);  
   else                emt=double(0.31*maga-6.725);    
   emt=emt+RandN(0.1,3.0);
   if(emt<-5.0) emt=-5.0;  
   emt=pow(10.0,emt);
   if(emt<0.00001 or emt>0.5 or maga<0.0){
   cout<<"Error emt:  "<<emt<<"\t maga:  "<<maga<<endl;}
   return(emt); 
}
///#############################################################################
double RandN(double sigma, double nn){
   double rr,f,frand;
   do{
   rr=RandR(-sigma*nn , sigma*nn); ///[-N sigma:N sigma]
   f= exp(-0.5*rr*rr/(sigma*sigma));
   frand=RandR(0.0 , 1.0);
   }while(frand>f);
   return(rr);
}
///#############################################################################
double RandR(double down, double up){
   double p =(double)rand()/((double)(RAND_MAX)+(double)(1.0));
   return(p*(up-down)+down);
}
///#############################################################################
double Fluxlimb(double limb, double rstar){
    return ( double(1.0-limb*(1.0-sqrt(fabs(1.0-rstar*rstar)))) );
}
///#############################################################################
double Bessel(int n,double x){
    double j1=0.00000001,tet;
    int kmax=10000;
    for(int k=0; k<kmax; ++k){
    tet=double(k*M_PI/kmax);
    j1+=double(M_PI/kmax/1.0)*cos(n*tet-x*sin(tet)); }
    return(j1/M_PI);
}    
///#############################################################################
void FiniteLens(source & s, lens & l, double lmass, double lradius, double limbs, double u){  
    double diss, mu, Fsour, lim1, lim2, x, y;  
    l.num0=0.0; l.num1=0.0;
    lim1=double(u-s.ros-sqrt((u-s.ros)*(u-s.ros)+4.0))/2.0-s.ros;///negative[RE]
    lim2=double(u+s.ros+sqrt((u+s.ros)*(u+s.ros)+4.0))/2.0+s.ros;//positive[RE]
    l.xsc=u;
    l.ysc=0.0;
    for(x=lim1; x<=lim2;  x=x+l.dx){
    for(y=0.0;  y<=double(1.0+2.5*s.ros); y=y+l.dy){
    l.xi=double(x*l.RE);//meter 
    l.yi=double(y*l.RE);//meter    
    LensEq2(s,l, lmass, lradius);
    diss=sqrt((l.xsi-l.xsc)*(l.xsi-l.xsc)+(l.ysi-l.ysc)*(l.ysi-l.ysc));//##[RE]
    if(diss<s.ros or diss==s.ros){
    mu=sqrt(fabs(1.0-diss*diss/(s.ros*s.ros)));   
    Fsour=1.0-limbs*(1.0-mu);//-s.limb2*(1.0-mu)*(1.0-mu); 
    l.num0+=Fsour*2.0;//area_of_images
    if(l.flag<1) l.num1+=Fsour*2.0;//part of images' area occultated
    if(Fsour<0.0 or mu>1.0 or s.ros<=0.0 or diss<0.0 or l.flag<0){ 
    cout<<"diss:  "<<diss<<"\t ros:  "<<s.ros<<"\t mu:  "<<mu<<"\t Fosour:  "<<Fsour<<endl;
    cout<<"num0:  "<<l.num0<<"\t num1:  "<<l.num1<<"\t flag: "<<l.flag<<endl; 
    cout<<"xsi:  "<<l.xsi<<"\t xsc:  "<<l.xsc<<endl;
    cout<<"ysi:  "<<l.ysi<<"\t ysc:  "<<l.ysc<<endl;  int uue;    cin>>uue;}}}}  
} 
//////##########################################################################
void LensEq2(source & s , lens & l, double lmass, double lradius){
    l.flag=-1; 
    double b, tanb, beta, xs0, ys0, d2, angle, tant, ttm;  
    
    b=    sqrt(l.xi*l.xi+l.yi*l.yi)+1.0e-50;//Impact parameter[m]
    angle=fabs(4.0*G*lmass*Msun/(velocity*velocity*b));//#radian  
    tant= double(b/l.Dl);
    ttm=  double(tan(angle)-tant)/(1.0+tan(angle)*tant );//tan(alfa-teta)
    tanb=tant-(tant + ttm)*l.Dls/s.Ds;
    beta= double(tanb*s.Ds);//#[m] in lens plane
    l.xsi=double(beta*l.xi/b)/l.RE;//[RE]
    l.ysi=double(beta*l.yi/b)/l.RE;//[RE]
    d2=   double(b*b/(l.RE*l.RE)); 
    xs0=  double(l.xi/l.RE-l.xi/l.RE/d2);
    ys0=  double(l.yi/l.RE-l.yi/l.RE/d2);
    /*
    if(double(angle*180.0/M_PI)>10.0 or fabs(l.xsi-xs0)>2.1 or fabs(l.ysi-ys0)>2.1){
    cout<<"new xs:  "<<l.xsi<<"\t ys:   "<<l.ysi<<endl;
    cout<<"old xs0: "<<xs0<<"\t   ys0: "<<ys0<<endl;
    cout<<"alpha: "<<angle*180.0/M_PI<<"\t b: "<<b<<endl;
    cout<<"xi/RE:  "<<l.xi/l.RE<<"\t yi/RE:  "<<l.yi/l.RE<<endl;
    cout<<"Dls:  "<<l.Dls<<"\t Ds(KP):  "<<s.Ds/KP<<"\t Dls/Ds "<<l.Dls/s.Ds<<endl;
    cout<<"lmass:  "<<lmass<<"\t lradius:  "<<lradius<<endl;
    cout<<"RE/Rstar: "<<l.RE/(lradius*Rsun)<<"\t Delta:  "<<fabs(tan(angle)-angle)<<endl;
    cout<<"tanb:  "<<tanb<<"\t tant:   "<<tant<<endl;} */
    //int yye; cin>>yye;}
    if(b<double(lradius*Rsun) or b==(lradius*Rsun)) l.flag=0;///##Occultation
    else                                             l.flag=1;
    //cout<<"LENSEQ xsi:  "<<l.xsi<<"\t xsc:  "<<l.xsc<<endl;
    //cout<<"ysi:  "<<l.ysi<<"\t ysc:  "<<l.ysc<<endl;
}
///#############################################################################

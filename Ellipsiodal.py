import numpy as np 
from numpy import conj
import matplotlib.pyplot as plt
import matplotlib
import pylab as py 
from matplotlib import rcParams
from matplotlib.ticker import FormatStrFormatter
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
rcParams["font.size"] = 11.5
rcParams["font.family"] = "sans-serif"
rcParams["font.sans-serif"] = ["Computer Modern Sans"]
rcParams["text.usetex"] = True
rcParams["text.latex.preamble"] = r"\usepackage{cmbright}"
rcParams['text.usetex'] = True
matplotlib.rc('text', usetex=True)
rcParams["text.latex.preamble"].join([r"\usepackage{dashbox}",r"\setmainfont{xcolor}",])
cmap=plt.get_cmap('viridis')
import scipy.special as ss
import warnings
warnings.filterwarnings("ignore")
cmap=plt.get_cmap('viridis')
from matplotlib import colors

RA=float(np.pi/180.0)
G= 6.67430*pow(10.0,-11.0)
AU=1.495978707*pow(10.0,11)
Msun=1.989*pow(10.0,30.0)
Rsun =6.9634*pow(10.0,8.0)
KPC= 3.0857*pow(10.0,19.0)## meter 
velocity=299792458.0##m/s
ni=90
nt=180; 
#dphase=float(180.0/nt)##in degree
###############################################################################
'''
dwd=  3#LP400-22; TESS_ID: 2002564035
limb= 0.3045 
MBH=  0.19;
#RBH=  0.01125*np.sqrt(np.power(MBH/1.454,-2.0/3.0)-np.power(MBH/1.454,2.0/3.0));
period=1.01016;#days
gc=   0.291 
li3= float(3.0-limb);

RBH=np.array([0.02, 0.1, 0.2, 0.15])

#mass>0.37

###############################################################################

dwd=5  #J2132+0754;  TESS_ID: 2000073295
limb=0.2712
MBH =0.17;
#RBH =0.01125*np.sqrt(np.power(MBH/1.454,-2.0/3.0)-np.power(MBH/1.454,2.0/3.0));
period=0.25056;##days
gc= 0.3408 
li3= float(3.0-limb);
RBH=np.array([0.02, 0.1, 0.2, 0.15])
#M2>0.95
'''
###############################################################################

dwd=6#J2151+1614;  TESS_ID: 2000525780
limb=0.2836
MBH =0.176;
#RBH =0.01125*np.sqrt(np.power(MBH/1.454,-2.0/3.0)-np.power(MBH/1.454,2.0/3.0));
period=0.59152;##days
gc= 0.3419 
li3= float(3.0-limb);
RBH=np.array([0.02, 0.1, 0.2, 0.15])
#M2>0.49
###############################################################################
fr=np.zeros((ni,5,5,3))
fij=open("./Ellipoid{0:d}.dat".format(dwd),"w")
fij.close(); 
Mass=np.zeros((10))

for l in range(3):
    for k in range(5): 
        Mass[k]=float(0.49+0.4*k)
        q=Mass[k]/MBH
        dis=np.power(np.power(period*24.0*3600.0,2.0)*G*Msun*(Mass[k]+MBH)/(4.0*np.pi*np.pi),1.0/3.0)
        ampl=float(RBH[l]*Rsun/dis);
        for i in range(ni): 
            inc=float(i*RA)#[radian]
            cosi2=np.sin(inc)*np.sin(inc);
            Delf2=np.zeros((nt,2))
            for j in range(nt):
                phase=float(j*1.0)*RA;#[radian]  
                L0=abs(1.0+(15.0+limb)*(1.0+gc)*pow(ampl,3.0)*(2.0+5.0*q)*(2.0-3.0*cosi2)/(60.0*(3.0-limb))+9.0*(1.0-limb)*(3.0+gc)*np.power(ampl,5.0)*q*(8.0-40.0*cosi2+35.0*cosi2*cosi2)/(256.0*(3.0-limb)));
                Delf2[j,1]=float((15*limb*(2.0+gc)*pow(ampl,4.0)*q*(4.0*np.sin(inc)-5.0*pow(np.sin(inc),3))/(32.0*li3))*np.cos(phase)+(-3.0*(15.0+limb)*(1.0+gc)*pow(ampl,3.0)*q*cosi2/(20.0*li3)-15.0*(1.0-limb)*(3.0+gc)*pow(ampl,5.0)*q*(6.0*cosi2-7.0*cosi2*cosi2)/(64.0*li3))*np.cos(2.0*phase)+(-25.0*limb*(2.0+gc)*pow(ampl,4.0)*q*pow(np.sin(inc),3.0)/(32.0*li3))*np.cos(3.0*phase)+(105.0*(1.0-limb)*(3.0+gc)*pow(ampl,5.0)*q*pow(np.sin(inc),4.0)/(256.0*li3) )*np.cos(4.0*phase))/L0;
                Delf2[j,0]=float(phase/RA)
        ###########################################################################      
            Delf=float(np.max(np.abs(Delf2[:,1])))
            DD=q*pow(np.sin(inc),2.0)*pow(ampl,3.0)*0.15*(15.0+limb)*(1.0+gc)/(3.0-limb)
            
            Dave=3.0*np.pi*np.pi*(15.0+limb)*(1.0+gc)*Mass[k]*pow(RBH[l]*Rsun,3.0)*np.sin(inc)*np.sin(inc)
            Dave=Dave/(5.0*pow(period*24.0*3600.0,2.0)*(3.0-limb)*G*MBH*(MBH+Mass[k])*Msun )
            
            test=np.array([ i, inc/RA, Delf, DD, Dave ])
            fr[i,0,k,l], fr[i,1,k,l], fr[i,2,k,l], fr[i,3,k,l], fr[i,4,k,l]= i, inc/RA, Delf, DD, Dave

            fij=open("./Ellipoid{0:d}.dat".format(dwd),"a+")
            np.savetxt(fij,test.reshape((1,5)),fmt="%d    %.5f     %.10f  %.10f   %.10f") 
            fij.close() 
    ###########################################################################
            if(k==4 and l==2 and int(i)%10==0):  
                plt.clf()
                plt.cla()
                fig=plt.figure(figsize=(8,6))
                ax1=fig.add_subplot(111)
                plt.plot(Delf2[:,0], Delf2[:,1], "k--", lw=1.9, label=r"$\rm{inclination}(\rm{deg})=$"+str(round(inc/RA,1)) )
                plt.xlabel(r"$\rm{Phase}(\rm{deg})$",  fontsize=18)
                plt.ylabel(r"$\rm{Ellipsiodal}~\rm{variation}$", fontsize=18)
                plt.xticks(fontsize=19, rotation=0)
                plt.yticks(fontsize=19, rotation=0)
                ax1.legend(prop={"size":17})
                ax1.grid("True")
                ax1.grid(linestyle='dashed')
                fig=plt.gcf()
                plt.subplots_adjust(hspace=.0)
                fig.savefig("./ellip{0:d}_{1:d}.jpg".format(dwd,i),dpi=200)
        ###########################################################################

cm=colors.ListedColormap(['k', 'y', 'red','blue', 'c', 'lightblue', 'm', 'g', 'orange', 'gray'])
plt.clf()
plt.cla()
fig=plt.figure(figsize=(8,6))
ax1=fig.add_subplot(111)

for k in range(5):
    plt.plot(fr[:,1,k,0], fr[:,2,k,0], "-", color=cm(k), lw=2.1, label=r"$M_{2}(M_{\odot})=~$"+str(round(Mass[k],2)) )

for k in range(5):
    plt.plot(fr[:,1,k,1], fr[:,2,k,1], "--", color=cm(k), lw=2.1)##, label=r"$R_{1}(R_{\odot})=~$"+str(round(RBH[1],2)) )
    
for k in range(5):
    plt.plot(fr[:,1,k,2], fr[:,2,k,2], ".-", color=cm(k), lw=2.1)##, label=r"$R_{1}(M_{\odot})=~$"+str(round(RBH[2],2)) )
        
    
plt.xlabel(r"$\rm{Inclination}(\rm{deg})$",  fontsize=18)
plt.ylabel(r"$\rm{Normalized}~\rm{ellipsiodal}~\rm{amplitude}$", fontsize=17)
plt.xticks(fontsize=19, rotation=0)
plt.yticks(fontsize=19, rotation=0)
plt.title(
r"$M_{1}(M_{\odot})=$"+'{0:.2f}'.format(MBH)+
#r"$;~R_{1}(R_{\odot})=$"+'{0:.2f}, {1:.2f}, {2:.2f}'.format(RBH[0], RBH[1], RBH[2])+
r"$;~\beta_{1}=$"+'{0:.2f}'.format(gc)+
r"$;~\Gamma_{1}=$"+'{0:.2f}'.format(limb)+r"$;~T(\rm{days})=$"+'{0:.2f}'.format(period), fontsize=17.0,color='k')
plt.xlim(0.0,90.0)
plt.yscale('log')
plt.ylim(1e-8, 2e-3)
plt.text(70.0, 8e-4, r"$R_{1}(R_{\odot})=$"+str(round(RBH[2],2)),fontsize =16)
plt.text(70.0, 1e-4, r"$R_{1}(R_{\odot})=$"+str(round(RBH[1],2)),fontsize =16)
plt.text(70.0, 2e-6, r"$R_{1}(R_{\odot})=$"+str(round(RBH[0],2)),fontsize =16)
#ax1.legend(prop={"size":14})
#ax1.grid("True")
#ax1.grid(linestyle='dashed')
legend=ax1.legend(prop={"size":22},loc='best',frameon=True, fancybox = True,shadow=True,framealpha=1.0)
legend.get_frame().set_facecolor('gray')
plt.legend(loc='best',fancybox=True, shadow=True)
if(dwd==3): dd=ax1.legend(title=r"$\rm{\bf LP400-22}$", prop={"size":15})
if(dwd==5): dd=ax1.legend(title=r"$\rm{\bf J2132+0754}$",prop={"size":15})
if(dwd==6): dd=ax1.legend(title=r"$\rm{\bf J2151+1614}$", prop={"size":15})
#ax1.legend(prop={"size":15})
#plt.setp(legend.get_title(),color="k")
fig=plt.gcf()
fig.tight_layout()
#plt.subplots_adjust(hspace=.0)
fig.savefig("./AmpliEllip{0:d}.jpg".format(dwd),dpi=200)   
###############################################################################



















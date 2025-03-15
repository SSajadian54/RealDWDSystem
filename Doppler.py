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
KPC= 3.0857*pow(10.0,19.0)## meter 
velocity=299792458.0##m/s
KBol= 1.380649*pow(10.0,-23.0)
Hplank=6.62607015*pow(10.0,-34.0)
const= np.sqrt(4.0*G*Msun)/velocity
const2=float(Hplank*velocity*np.power(10.0,9.0)/KBol )
ni=30
nt=180;#time 

dinc= float(90.0/ni)

tmod=  np.zeros((nt));
phi=   np.zeros((nt));  
x1=    np.zeros((nt));  
x0=    np.zeros((nt));
y1=    np.zeros((nt));  
y0=    np.zeros((nt));
z1=    np.zeros((nt));  
Vx=    np.zeros((nt));
Delf1= np.zeros((nt));
tp=0.0

for i in range(nt):  
    tmod[i]=float(-0.5+i/nt/1.0)
    phi[i] =(tmod[i]-tp)*2.0*np.pi 


###############################################################################
dwd=  3#LP400-22; TESS_ID: 2002564035
limb= 0.3045 
MBH=  0.19;
RBH=  0.01125*np.sqrt(np.power(MBH/1.454,-2.0/3.0)-np.power(MBH/1.454,2.0/3.0));
period=1.01016;#days
Tstar= 11140.0
gc=   0.291 
li3= float(3.0-limb);
#mass>0.37

###############################################################################
'''
dwd=5  #J2132+0754;  TESS_ID: 2000073295
limb=0.2712
MBH =0.17;
#RBH =0.01125*np.sqrt(np.power(MBH/1.454,-2.0/3.0)-np.power(MBH/1.454,2.0/3.0));
period=0.25056;##days
gc= 0.3408 
li3= float(3.0-limb);
RBH=np.array([0.02, 0.1, 0.2, 0.15])
#M2>0.95
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
'''
###############################################################################

Nw=110
throu=np.zeros((Nw , 2))
throu=np.loadtxt("./TESS_Throught.txt")


def Fplank(wave):
    wave=wave*pow(10.0,-9.0)#[m]  
    con= Hplank*velocity/(KBol*Tstar*wave)
    return(2.0*Hplank*velocity*velocity*pow(wave,-5.0)/(np.exp(con)-1.0)  )  

def Doppler(vx):
    for m in range(nt): 
        DF0=0.0
        Delf1[m]=0.0
        for s in range(Nw-1):
            wave= throu[s,0]             
            dw=float(throu[s+1,0]-throu[s,0])*pow(10.0,-9)
            waven= wave+wave*vx[m]/velocity
            Delf1[m]+= float(Fplank(waven) - Fplank(wave))*throu[s,1]*dw       
            DF0+=Fplank(wave) * throu[s,1] * dw
        Delf1[m]=float(Delf1[m]/DF0)     
    return(Delf1)



###############################################################################
fr=np.zeros((ni,4,10))
fij=open("./Doppler{0:d}.dat".format(dwd),"w")
fij.close(); 
Mass=np.zeros((10))
ecen=0.0
teta=0.0
for l in range(1):## primary radius
    for k in range(10):## secondary mass 
        Mass[k]=float(0.37+0.63*k)
        dis=np.power(np.power(period*24.0*3600.0,2.0)*G*Msun*(Mass[k]+MBH)/(4.0*np.pi*np.pi),1.0/3.0)#meter
        for i in range(ni):
            print ("**************, i", i, k) 
            inc=float(i*dinc)*RA#[radian
            x0=dis*(np.cos(phi)-ecen)
            y0=dis*np.sin(phi)*np.sqrt(1.0-ecen**2.0)
            y1=               y0*np.cos(teta)+x0*np.sin(teta)
            x1= np.cos(inc)*(-y0*np.sin(teta)+x0*np.cos(teta))
            z1=-np.sin(inc)*(-y0*np.sin(teta)+x0*np.cos(teta)) 
            ###################################################################
            for s1 in range(nt-1):   
                 Vx[s1]=float(x1[s1+1]-x1[s1])/(tmod[s1+1]-tmod[s1])/(period*24.0*3600.0)##[m/s]
            Vx[nt-1]=Vx[nt-2]+(Vx[nt-2]-Vx[nt-3])*(tmod[nt-1]-tmod[nt-2])/(tmod[nt-2]-tmod[nt-3])
            Delf1=Doppler(Vx)
            Delf=np.max(np.abs(Delf1))
            vmax=np.max(np.abs(Vx))
            test=np.array([ i, inc/RA, Delf, vmax ])
            fr[i,0,k], fr[i,1,k], fr[i,2,k], fr[i,3,k]= i, inc/RA, Delf, vmax
            fij=open("./Doppler{0:d}.dat".format(dwd),"a+")
            np.savetxt(fij,test.reshape((1,4)),fmt="%d    %.5f     %.10f      %.10f ") 
            fij.close() 
            ###########################################################################
            if(int(i)%10==0):  
                plt.clf()
                plt.cla()
                fig=plt.figure(figsize=(8,6))
                ax1=fig.add_subplot(111)
                plt.plot(tmod, Delf1, "k--", lw=1.9, label=r"$\rm{inclination}(\rm{deg})=$"+str(round(inc/RA,1)) )
                plt.xlabel(r"$\rm{Phase}(\rm{deg})$",  fontsize=18)
                plt.ylabel(r"$\rm{Doppler}~\rm{Boosting}$", fontsize=18)
                plt.xticks(fontsize=19, rotation=0)
                plt.yticks(fontsize=19, rotation=0)
                ax1.legend(prop={"size":17})
                ax1.grid("True")
                ax1.grid(linestyle='dashed')
                fig=plt.gcf()
                plt.subplots_adjust(hspace=.0)
                fig.savefig("./Doppler{0:d}_{1:d}_{2:d}.jpg".format(dwd,i, k),dpi=200)
            ###########################################################################

cm=colors.ListedColormap(['k', 'y', 'red','blue', 'c', 'lightblue', 'm', 'g', 'orange', 'gray'])
plt.clf()
plt.cla()
fig=plt.figure(figsize=(8,6))
ax1=fig.add_subplot(111)
for k in range(8):
    plt.plot(90.0-fr[:,1,k], fr[:,2,k], "-", color=cm(k), lw=2.1, label=r"$M_{2}(M_{\odot})=~$"+str(round(Mass[k],2)) )
plt.xlabel(r"$\rm{Inclination}(\rm{deg})$",  fontsize=18)
plt.ylabel(r"$\rm{Normalized}~\rm{Doppler}~\rm{amplitude}$", fontsize=17)
plt.xticks(fontsize=19, rotation=0)
plt.yticks(fontsize=19, rotation=0)
plt.title(
r"$M_{1}(M_{\odot})=$"+'{0:.2f}'.format(MBH)+
r"$;~R_{1}(R_{\odot})=$"+'{0:.2f}'.format(RBH)+
#r"$;~a (R_{\odot})=$"+'{0:.2f}'.format(dis/Rsun)+
r"$;~T(\rm{days})=$"+'{0:.2f}'.format(period), fontsize=17.0,color='k')
plt.xlim(0.0,90.0)
#plt.yscale('log')
plt.ylim(0.0,0.0036)
#plt.text(70.0, 8e-4, r"$R_{1}(R_{\odot})=$"+str(round(RBH[2],2)),fontsize =16)
#plt.text(70.0, 1e-4, r"$R_{1}(R_{\odot})=$"+str(round(RBH[1],2)),fontsize =16)
#plt.text(70.0, 2e-6, r"$R_{1}(R_{\odot})=$"+str(round(RBH[0],2)),fontsize =16)
#ax1.legend(prop={"size":14})
#ax1.grid("True")
#ax1.grid(linestyle='dashed')
legend=ax1.legend(prop={"size":20},loc='best',frameon=True, fancybox = True,shadow=True,framealpha=1.0)
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
fig.savefig("./AmpliDoppler{0:d}.jpg".format(dwd),dpi=200)   
###############################################################################



















from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib import gridspec
import numpy as np
import pandas as pd
import matplotlib as mpl
import pylab
rcParams["font.size"] = 18
rcParams["font.family"] = "sans-serif"
rcParams["font.sans-serif"] = ["Computer Modern Sans"]
rcParams["text.usetex"] = True
rcParams["text.latex.preamble"] = r"\usepackage{cmbright}"
from matplotlib.ticker import FormatStrFormatter
from matplotlib.gridspec import GridSpec
import pandas as pd
###########################################################################################    
#source_id     ra     dec     parallax     parallax_error     pmra     pmra_error     pmdec     pmdec_error     phot_g_mean_mag     bp_rp     bp_g     g_rp     radial_velocity     radial_velocity_error

##source_id,ra,dec,parallax,parallax_error,pmra,pmra_error,pmdec,pmdec_error,phot_g_mean_mag,bp_rp,bp_g,g_rp,radial_velocity,radial_velocity_error,phot_variable_flag,target_id,target_ra,target_dec,target_separation (deg) 

df = pd.read_csv("./LP40022.csv")
nrow= len(df.source_id)
print ("n_row:  ",  nrow)
Ra=  np.zeros((nrow))
Dec= np.zeros((nrow))
Mag= np.zeros((nrow))
col= np.zeros((nrow))
GG = np.zeros((nrow))

ra0 = 15.0*(22.0+36.0/60.0+30.10/3600.0) 
dec0=  1.0*(22.0+32.0/60.0+24.00/3600.0) 
print ("Source star position:  ",   ra0, dec0)
mind= 100000.0
for i in range(nrow): 
    Ra[i]=   df.ra[i]
    Dec[i]= df.dec[i]
    Mag[i]= (21.1- df.phot_g_mean_mag[i])*26.5##10.0**(-0.4*df.phot_g_mean_mag[i]) *1000000000.0
    col[i]= df.bp_rp[i]
    GG[i]=  df.phot_g_mean_mag[i]
    #print "counter:  ",  Ra[i],     Dec[i],        Mag[i],   col[i],  GG[i]
    dis=np.sqrt( (Ra[i]-ra0)**2.0 + (Dec[i]-dec0)**2.0 )
    if(dis<mind):  
        mind=dis
        id0= i
id1=0 
id2=0
for i in range(nrow): 
    if(df.source_id[i]==1874523804732334464): ## target itself 
        id1=i
    if(df.source_id[i]==1874523770372593664): ## variable star  1874523770372593664
        id2=i  
    dis=np.sqrt((Ra[i]-ra0)**2.0 + (Dec[i]-dec0)**2.0)*3600.0##  arcsec
    if(dis<21.0 and i!=id0):  
        print( " *** RA:  ",  Ra[i],   "Dec:  ", Dec[i],  "  distance:  ", dis, "    counter:  ",  i,  id0)
        
print(np.max(df.phot_g_mean_mag[:]))          
        
print( "***********************************************" )
print( "Min Dis:   ",  mind,  " iD:  ",    df.source_id[id0],  "    counter:  ",  id0)
print( "parallax_erro: ",  df.parallax_error[id0] ) 
print( "parallax (marcs):   ", df.parallax[id0], " distance(kpc):   ",    float(1.0/df.parallax[id0]) )
print(df.phot_g_mean_mag[id0],  df.phot_g_mean_mag[id1],  np.sqrt((Ra[id1]-ra0)**2.0 + (Dec[id1]-dec0)**2.0)*3600.0 )
print(df.phot_g_mean_mag[id0],  df.phot_g_mean_mag[id2],  np.sqrt((Ra[id2]-ra0)**2.0 + (Dec[id2]-dec0)**2.0)*3600.0 )
print( "***********************************************")        
########################################################       

plt.clf()
fig, ax = plt.subplots()
plt.figure(figsize=(8,6))
plt.scatter( (Ra-Ra[id0])*60.0 ,     (Dec-Dec[id0])*60.0 ,     s=Mag*3, marker="o", facecolors='pink', edgecolors='m')##,  "ro")
plt.scatter( (Ra[id0]-Ra[id0])*60.0, (Dec[id0]-Dec[id0])*60.0, s=Mag[id0]*3, marker="*", color="g")##,  "ro")
plt.scatter( (Ra[id2]-Ra[id0])*60.0, (Dec[id2]-Dec[id0])*60.0, s=Mag[id2]*3, marker="^", color="b")##,  "ro")
plt.xticks(fontsize=18, rotation=0)
plt.yticks(fontsize=18, rotation=0)
plt.xlabel(r"$\rm{RA(arcm)}-339.1^{\circ}$", fontsize=18.0)
plt.ylabel(r"$\rm{Dec(arcm)}-22.5^{\circ}$", fontsize=18.0)
#plt.text(18.4,7.11,r"$-28^{\circ}$", fontsize=18.0)
plt.text(-0.2,-0.15, r"$\rm{LP400}-22$",fontsize=18.0)
plt.text(-0.8,-0.4,r"$\rm{Variable}~\rm{star}$",fontsize=18.0)
plt.xlim([-1.5, 1.0])
plt.ylim([-1.5, 1.0])
plt.grid("True")
plt.grid(linestyle='dashed')
fig=plt.gcf()
fig.tight_layout()
fig.savefig("gaiaLP40022.jpg",dpi=200)
##########################################################


plt.clf()
fig, ax = plt.subplots()
plt.figure(figsize=(8,6))
plt.scatter(col, GG , marker="o", color="r")##,  "ro")
plt.scatter(col[id0], GG[id0] , marker="*", color="b")##,  "ro")
plt.xticks(fontsize=16, rotation=0)
plt.yticks(fontsize=16, rotation=0)
plt.xlabel(r"$\rm{bp-rp} (mag)$",  fontsize=18.0)
plt.ylabel(r"$\rm{G}-\rm{mag}$", fontsize=18.0)
#plt.xlim([0.31*60.0,0.285*60.0])
#plt.ylim([(-27.898+28.0)*60.0, (28.0-27.877)*60.0 ])
#plt.xlim([18.4,  17.15])
#plt.ylim([6.2, 7.1])
#ax.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
plt.grid("True")
plt.grid(linestyle='dashed')
fig=plt.gcf()
fig.savefig("LP40022B.jpg",dpi=200)

##########################################################

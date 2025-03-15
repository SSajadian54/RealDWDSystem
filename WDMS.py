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

Rsun=6.9634*np.power(10.0,8.0)###; ///solar radius [meter]
G=6.67384*np.power(10.,-11.0);##// in [m^3/s^2*kg].
Msun=1.98892*np.power(10.,30.0);## //in [kg].
velocity=299792458.0;##//velosity of light  m/s
AU=1.495978707*np.power(10.0,11.0);
year=float(365.2422)
parcs=float(30856775814913673.0)
mm= AU*0.001/(year*24.0*3600.0) 

Mlens=0.6 
Msource=0.6
period=30.0
Msum=(Mlens+Msource)*Msun
semi=np.power(G*Msum*period*period*24.0*24.0*3600.0*3600.0/(4.0*np.pi*np.pi),1.0/3.0)## meter

RE=np.sqrt(4.0*G*Mlens*Msun)/velocity *np.sqrt(semi)+1.0e-50;#meter 
print(RE/Rsun/0.01,  semi/AU)

input("Enter a number ")

#################################################
#SELECT  
#tw.source_id,   tw.period,  tw.mass_ratio,  tw.eccentricity,  tw.inclination,  
#tw.a_thiele_innes, tw.b_thiele_innes,  tw.f_thiele_innes,  tw.g_thiele_innes
#FROM   gaiadr3.nss_two_body_orbit AS tw
#where 

fil2=open("IDS2.txt","w")
fil2.close()

#################################################
#RA, DEC, dis, source_id, simbadid, sed_fit, parallax, e_parallax, dist, pmra, e_pmra, pmdec, e_pmdec,  g_mag, gbp_mag,  grp_mag, ruwe, excess_noise, e_excess_noise, teff_wd, e_teff_wd, logg_wd, e_logg_wd, lbol_wd, e_lbol_wd, r_wd, e_r_wd, m_wd, e_m_wd, teff_ms, e_teff_ms, lbol_ms, e_lbol_ms, r_ms, e_r_ms, flag
df = pd.read_csv("./wdw4_1725261998.csv")##http://svocats.cab.inta-csic.es/wdw4/index.php?action=search
No= len(df.source_id)
print ("n_row_WDMS_cataloge:  ",  No)
#################################################

#Cataloge TwoBodyOrbit table of Gaia Data Release3(gaiadr3.nss_two_body_orbit) 
df2 = pd.read_csv("./TwoBodyOrbit.csv")#source_id,period,mass_ratio,eccentricity,inclination,a_thiele_innes,b_thiele_innes,f_thiele_innes,g_thiele_innes
Nob= len(df2.source_id)
print ("n_row:  ",  Nob)
#################################################
##https://ui.adsabs.harvard.edu/abs/2024MNRAS.527.6100N/abstract
#Index,Gaia-ID, -ID, RA,DEC,Teff_MS,e_Teff_MS,Rad1_MS,eRad1_MS,Lbol_MS,Lberr_MS,Teff_WD, e_Teff_WD,Rad1_WD,eRad1_WD,Lbol_WD,Lberr_WD,wd_mass,wd_mass_error,RedChi2,Vgfb,frac_excess_FUV,frac_excess_NUV
df3 = pd.read_csv("./well_fitted_WDMS_binaries.csv")
Noc= len(df3.GaiaID)

print ("n_row:  ",  Noc)

#################################################

for i in range(Noc): 
    ids=df3.GaiaID[i]
    #fil2= open("IDS2.txt","a+")
    #np.savetxt(fil2,ids.reshape(-1,1), fmt="gs.source_id='%d'  OR ")
    #fil2.close() 
    for j in range(Nob):
        idw=df2.source_id[j] 
        if(ids==idw): 
            print ("one case is sound, Source_IDs:  ",    ids,     idw )















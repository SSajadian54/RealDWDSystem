import numpy as np

################################################################################
N1=3192; 
#http://cdsarc.u-strasbg.fr/viz-bin/cat/J/A+A/634/A93#/browse
#logg   Teff      ZR      u(Ke)    u(Te)   u(Gbp)    u(Gg)    u(Grp)
#logg   Teff      ZR      (Kesig)  (Tesig) (Gbpsig)  (Ggsig)  (Grpsig)
#logg   Teff      ZR      I(Ke)    I(Te)   I(Gbp)    I(Gg)    I(Grp)
#logg   Teff      ZR      be(Ke)   be(Te)  be(Gbp)   be(Gg)   be(Grp)
dat=np.zeros((N1,8))
dat=np.loadtxt('./OriginalTab/TABLE99C')
file1=open("WDLimbC2020.dat","w")
file1.close()
nn=np.argsort(dat[:,1])
for i in range(N1): 
    if(i%4==0):
        par=np.array([ dat[i,0], dat[i,1], dat[i,4]])## logg,   Teff,   u(Te)
        file1=open("WDLimbC2020.dat","a+")
        np.savetxt(file1,par.reshape(1,3),fmt='%.2f   %.1f   %.5f') 
        file1.close()
###############################################################################
#http://cdsarc.u-strasbg.fr/viz-bin/cat/J/A+A/634/A93#/browse
#logg    Teff     ZR      y1(Ke)   y1(Te)   y1(Gbp)  y1(Gg)   y1(Grp)
#logg    Teff     ZR      y2(Ke)   y2(Te)   y2(Gbp)  y2(Gg)   y2(Grp)
N2=1596
dat=np.zeros((N2,8))
dat=np.loadtxt('./OriginalTab/TABLE105C')## gravity-daekening 
file2=open("WDGravC2020.dat","w")
file2.close()
nn=np.argsort(dat[:,3])
for i in range(N2): 
    if(i%2==0):
        par=np.array([ dat[i,0], dat[i,1], dat[i,4]])## logg, Teff,  gravityD
        file2=open("WDGravC2020.dat","a+")
        np.savetxt(file2,par.reshape(1,3),fmt='%.2f   %.1f   %.5f') 
        file2.close()
        
###############################################################################     
N3=798           
dat=np.zeros((N3,6))
dat[:,:3]=np.loadtxt('./WDLimbC2020.dat')
dat[:,3:]=np.loadtxt('./WDGravC2020.dat')        
file3=open("WDCoefC2020.dat","w")       
file3.close()  
        
for i in range(N3): 
    if(dat[i,0]!= dat[i,3] or dat[i,1]!=dat[i,4]): 
        print ("Error i : ", i, dat[i,: ])
        input("Enter a number ")
    else:
        file3=open("WDCoefC2020.dat","a+")   
        par=np.array([dat[i,0],  dat[i,1],  dat[i,2],  dat[i,5] ])
        np.savetxt(file3,par.reshape(1,4),fmt='%.2f   %.1f   %.5f   %.5f')       
        file3.close()        
        
###############################################################################                
        
        
        
        
        


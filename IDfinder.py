import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.timeseries import BoxLeastSquares
from matplotlib import rcParams
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib import rcParams 
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
import warnings
warnings.filterwarnings("ignore")
rcParams["font.size"] = 13
rcParams["font.family"] = "sans-serif"
rcParams["font.sans-serif"] = ["Computer Modern Sans"]
rcParams["text.usetex"] = True
rcParams["text.latex.preamble"] = r"\usepackage{cmbright}"
################################################################################
nm=[
'0001', '0002', '0003', '0004', '0005', '0006', '0007', '0008', '0009', '0010', '0011', '0012', '0013', '0014', '0015', '0016', '0017', '0018', '0019', '0020', '0021', '0022', '0023', '0024', '0025', '0026', '0027', '0028', '0029', '0030', '0031', '0032', '0033', '0034', '0035', '0036', '0037', '0038', '0039', '0040', '0041', '0042', '0043', '0044', '0045', '0046', '0047', '0048', '0049', '0050', '0051', '0052', '0053', '0054', '0055', '0056', '0057', '0058', '0059', '0060', '0061', '0062', '0063', '0064', '0065', '0066', '0067', '0068', '0069', '0070', '0071', '0072', '0073', '0074', '0075', '0076', '0077', '0078', '0079', '0080', '0081', '0082', '0083', '0084', '0085', '0086', '0087', '0088', '0089', '0090',
'0091', '0092', '0093', '0094', '0095', '0096', '0097', '0098', '0099']



RA0=np.zeros((100));  DEC0=np.zeros((100))
RA0[0], DEC0[0]= float(92.48600123120592),    float(-61.35807261978761)  ## target_1 variable star
RA0[1], DEC0[1]= float(86.46431482068473),    float(-60.546206368911164) ## target_2  
RA0[2], DEC0[2]= float(98.62039453444818),    float(-74.28581516547928) ## target_3
RA0[3], DEC0[3]= float(92.649519  ),          float(-70.871977)         ## target_4  
RA0[4], DEC0[4]= float(96.28613363517778),    float(-70.17620296135674) ## target_5
RA0[5], DEC0[5]= float(97.68466962238169),    float(-62.29245205993212) ## target_6
RA0[6], DEC0[6]= float(85.82503109011775),    float(-73.63105439941128) ## target_7
RA0[7], DEC0[7]= float(86.74908360771985),    float(-71.76723742488748) ## target_8
RA0[8], DEC0[8]= float(83.13470215435092),    float(-67.2331821107543) ## target_9
RA0[9], DEC0[9]= float(80.25340141517492),    float(-64.33199021560245) ## target_10
RA0[10], DEC0[10]= float(92.22775836932337),  float(-70.74242005500108) ## target_11
RA0[11], DEC0[11]= float(91.97987838272857),  float(-69.69265402800174)  ##12   variable
RA0[12], DEC0[12]= float(99.28053337552106),  float(-68.47080880934016) ##13
RA0[13], DEC0[13]= float(84.58558420365227),  float(-69.3515603598234)  ##14
RA0[14], DEC0[14]= float(89.90246546092416),  float(-66.58763373340597) ##15
RA0[15], DEC0[15]= float(87.47292034744957),  float(-64.39040231695888)##16
RA0[16], DEC0[16]= float(84.4025349566182),   float(-64.42439227225557)   ##17
RA0[17], DEC0[17]= float(84.30658267065951),  float(-64.22107106059525) ##18
RA0[18], DEC0[18]= float(82.41372411719799),  float(-64.0367758881953)  ##19
RA0[19], DEC0[19]= float(82.19339537284809),  float(-62.0890234645653)  ##20
RA0[20], DEC0[20]= float(82.04680346478071),  float(-61.26959446183596) ##21
RA0[21], DEC0[21]= float(80.26299241734958),  float(-62.67877108080558)##22
RA0[22], DEC0[22]= float(92.83338395288197),  float(-69.97479703223877) ##23
RA0[23], DEC0[23]= float(93.9897147742992),   float(-69.90880264005166)  ##24
RA0[24], DEC0[24]= float(98.416667910905),    float(-68.74134662621121)   ##25
RA0[25], DEC0[25]= float(97.38074539079668),  float(-64.72966889784003)  ##26
RA0[26], DEC0[26]= 100.4774007995452,  -67.1070466544051  ##27 SB1  
#5280389169171864448 100.4774007995452 -67.1070466544051 1.3988924 95.35859558816462 NOT_AVAILABLE SB1
RA0[27], DEC0[27]= 105.63980855791702, -56.62827883161682
RA0[28], DEC0[28]=101.8269703026828, -65.5728905729796
RA0[29], DEC0[29]=98.79424171766726, -70.68182432045731
RA0[30], DEC0[30]= 86.59601638101688, -55.936853544056774

RA0[31], DEC0[31]=327.9967083  ,  16.2468667  



co=np.zeros((4))

for i in range(1):
    i=i+31
    fil=open("./target{0:d}.dat".format(i+1),"w")
    fil.close()
    print(" HHHHHHHHHHHHHHHHHHHH Target:  ",  i+1)
    for k in range(70): 
        sector=k+1
        name= "./sec"+str(int(k+1))+"/s"+str(nm[k])+".csv"
        df=pd.read_csv(name , usecols=[0, 1, 2] )
        nd=int(len(df['RA']))
        ra=np.zeros((nd))
        dec=np.zeros((nd))
        Tic=['']*nd ## an string array
        ra, dec= df['RA'], df['DEC']
        Tic= df['#TIC_ID']
        dis= np.sqrt((ra-RA0[i])**2.0+(dec-DEC0[i])**2.0)
        nn=int(np.argmin(dis))
        print("file_name: ",  name)
        print("Sector, dis, No_row,   TIC_number:  ", sector,   dis[nn], nn, Tic[nn]) 
        print("RA, DEC", RA0[i],   ra[nn],      DEC0[i],    dec[nn])    
        print("*******************")
        co=np.array([int(sector), int(nn), int(Tic[nn])  ,  dis[nn] ])
        fil=open("./target{0:d}.dat".format(i+1),"a+")
        np.savetxt(fil,co.reshape((1,4)),fmt ="%d    %d        %d        %.5f")
        fil.close()




























import numpy as np
from math import *
from numpy import *
#from scipy import *
import pandas as pd
import matplotlib.pyplot as plt
     
a=pd.read_csv('HF1.csv',header=0)
dt=pd.DataFrame(a)
print('\nESTE ES EL DATA FRAME\n',a)
tiempo=dt.loc[0:len(a),'time']
ax=dt.loc[0:len(a),'Ax']    
ay=dt.loc[0:len(a),'Ay']     
az=dt.loc[0:len(a),'Az']    
Aceb=np.array([ax,ay,az])   #matriz de aceleración

mx=dt.loc[0:len(a),'Mx'] 
my=dt.loc[0:len(a),'My'] 
mz=dt.loc[0:len(a),'Mz'] 
Mb=np.array([mx,my,mz])     #matriz campo magnético
Acei=np.array([0,0,-9.81])
print(Acei)

print('\nMb\n',Mb,np.shape(Mb))
print('Aceb\n',Aceb)
print('\nAceb[0]\n',Aceb[0])
print('\nMb\n',Mb)

M_b=[]
A_b=[]
Time=[]
Time.append(tiempo)
RR=[]
for i in range(2):
    Mb=np.array([mx[i],my[i],mz[i]])
    Aceb=np.array([ax[i],ay[i],az[i]])
    M_b.append(Mb)
    A_b.append(Aceb)

    print('\nAceb[i]\n',Aceb[i],mx[i],my[i],mz[i])
    print('\nMb[i]\n',Mb[i])
    Acei=np.array([0,0,-9.81])
    Mi=np.array([0.377,-0.499,-1.124])
    #print('Mi\n',Mi)

    AceB=Aceb/linalg.norm(Aceb)     #normaliza vectores
    MB=Mb/linalg.norm(Mb)
    AceI=Acei/linalg.norm(Acei)
    MI=Mi/linalg.norm(Mi)

    print('AceB\n',AceB,'\nMB\n',MB)
    t1b=AceB
    t1i=AceI
    t2B=np.cross(AceB,MB)       
    t2I=np.cross(AceI,MI)
    t2b=t2B/np.linalg.norm(t2B)
    t2i=t2I/np.linalg.norm(t2I)
    #print('\nEste es t2b\n',t2b)
    #print('\nEste es t2i\n',t2i)
    t3b=np.cross(t1b,t2b)
    t3i=np.cross(t1i,t2i)
    #print('\nEste es t3b\n',t3b)
    #print('\nEste es t3i\n',t3i)
    Tb=np.array([t1b,t2b,t3b])
    Ti=np.array([t1i,t2i,t3i])
    Tb_trans=np.transpose(Tb)      #checar transposicion
    #print('\nTb\n',Tb,'\nTb transpuesta\n',Tb_trans,'\nTi\n',Ti)
    R=Tb_trans@Ti
    RR.append(R)
print(R)   
print('\nMatrices de rotación R\n',RR)
print(np.shape(RR))
print('\nRR[0]\n',RR[0])
print('\n',np.transpose(RR[0]))
#print('\nAceb iterado\n',A_b,'\nMb iterado\n',M_b) #,'\nTiempo\n',Time)
#print('\nM_b[0]\n',M_b[0])
np.savetxt('text.txt',RR,fmt='%.2f')

df=pd.DataFrame(RR)
#df.columns=['
ruta="E:\EMIDSS\DATOSTRIAD1.csv"
df.to_csv(ruta)

np.savetxt('text.csv',RR,fmt,'%.2f')
               







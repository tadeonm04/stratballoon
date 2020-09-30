from math import *
from numpy import *
from pandas import *
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
#from pyquaternion import Quaternion
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import style

a=pd.read_csv('HF1.csv',header=0)
df=pd.DataFrame(a)
ax=df.loc[0:len(a),'Ax']
ay=df.loc[0:len(a),'Ay']
az=df.loc[0:len(a),'Az']
V_1b=np.array([ax,ay,az])

mx=df.loc[0:len(a),'Mx']
my=df.loc[0:len(a),'My']
mz=df.loc[0:len(a),'Mz']
V_2b=np.array([mx,my,mz])

V_1i=np.array([0,0,-9.81])  #se supone que ya están transpuestos (vector fila originalmemte columna)
V_2i=np.array([0.377,-4.88,-1.124])

Time=[]
Q=[]
V1b=[]
V2b=[]
x=[]
y=[]
z=[]
n=[]
for i in range(10):
    V1B=[ax[i],ay[i],az[i]]
    V2B=[mx[i],my[i],mz[i]]
    V_1b=V1B/np.linalg.norm(V1B)
    V_1b=np.transpose(V_1b).reshape(3,1)
    V_2b=V2B/np.linalg.norm(V2B)
    V_2b=np.transpose(V_2b).reshape(3,1)
    V1b.append(V_1b)
    V2b.append(V_2b)
    V_1i=V_1i/np.linalg.norm(V_1i)
    V_1itrans=np.transpose(V_1i).reshape(1,3)
    V_2i=V_2i/np.linalg.norm(V_2i)
    V_2itrans=np.transpose(V_2i).reshape(1,3) #originalmente 3x1 se transpone
    
    rho_k=2
    lambda_opt=rho_k+rho_k

    B=rho_k*(V_1b@V_1itrans)+rho_k*(V_2b@V_2itrans)
    #print(B)
    BT=np.transpose(B).reshape(3,3)
    BB=np.transpose(BT) #matriz igual a "B" pero ya manipulable por entradas
    #print(BB[0][0])
    #print(BB)
    S=B+BT
    #print(S)
    
    Z=np.array([BB[1][2]-BB[2][1],BB[2][0]-BB[0][2],BB[0][1]-BB[1][0]]) #se resta "1" a pos original, porque el arreglo comienza en 0 no en 1
    ZZ=np.vstack(Z)
    sigma=np.trace(BB)
    tao=(np.linalg.inv((lambda_opt+sigma)*np.identity(3)-S))@ZZ 
    #print('es tao\n',tao)
    TAO=np.transpose(tao).reshape(1,3) #vector 1x3
    #print(np.shape(TAO))
    taoo=np.transpose(TAO) #vector 3x1
    #print(np.shape(taoo))
    T=np.append(1,taoo) #vectror (,4)
    T=np.expand_dims(T,axis=0) #vector 1x4
    #print(np.shape(T))
    #print(taoo, T)
    q=(1/sqrt(1+TAO@taoo)@T)     
    x.append(q[0][1])   #primera coordenada imaginaria de cuaternión
    y.append(q[0][2])
    z.append(q[0][3]) 
    Q.append(q)
    N=q[0][0]**2+q[0][1]**2+q[0][2]**2+q[0][3]**2 #suma de cuadrados, cuaternión
    n.append(N)

print('\ncuaterniones\n',Q)
print('\nPRIMER CUATERNION',Q[0])
print(np.shape(Q[0]))
print('\nX\n',x,'\nY\n',y,'\nZ\n',z)
print('\nSUMA DE CUADRADOS\n',n)
print(min(x),max(x))
print(min(y),max(y))
print(min(z),max(z))

fig=plt.figure()
ax=fig.add_subplot(111,projection='3d')
#ax.plot_wireframe(x,y,z,restride=1,cestride=1, cmap='Blues')
ax.plot(x,y,z)
#plt.xlim()
#plt.ylim()
#plt.zlim()
plt.show()





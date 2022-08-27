import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from pathlib import Path
import glob

#from mpl_toolkits.mplot3d import Axes3D
xres = 2000
x = np.linspace(0,1,196)



omega=np.loadtxt("output/omega.txt",delimiter=" ")
plt.plot(omega[:,0],omega[:,1])
#%%    
#fig = plt.figure(figsize=(6,6))
for count in range(6):
    file=glob.glob("output/field_"+str(count+5)+".csv",recursive=True)
    w=np.loadtxt(file[0],delimiter=",",dtype=float)
    print(file)
    
    #plt.clf()
    x=np.sin(w[:,0])*(1+0.01*(w[:,1]+w[:,2]))
    y=np.cos(w[:,0])*(1+0.01*(w[:,1]+w[:,2]))
    plt.axis('equal')
    plt.plot(x,y)
    plt.title("lanslide number="+str(count))
    plt.pause(0.5)
 #%%   
fig = plt.figure(figsize=(6,6))   
for count in range(len(omega)):

    file=glob.glob("output/field_"+str(count)+".csv",recursive=True)
    w=np.loadtxt(file[0],delimiter=",",dtype=float)
    print(file)
    
    plt.clf()
    x=np.sin(w[:,0])*(1+0.01*(w[:,1]+w[:,2]))
    y=np.cos(w[:,0])*(1+0.01*(w[:,1]+w[:,2]))
    plt.plot(w[:,0],w[:,1])
    plt.pause(0.5)  
    
    
# textFile = open("uout2.txt")
# lineu1 = textFile.readlines()
# for line in lineu1:
#     datu.append(line.split(" "))
    
# textFile = open("vout2.txt")
# linev1 = textFile.readlines()
# for line in linev1:
#     datv.append(line.split(" "))
      
# textFile = open("tout2.txt")
# lino = textFile.readlines()
# for line in lino:
#     dato.append(line.split(" "))
    
# textFile = open("sout2.txt")
# linu = textFile.readlines()
# for line in linu:
#     dats.append(line.split(" "))
        
# textFile = open("mout2.txt")
# linm = textFile.readlines()
# for line in linm:
#     datm.append(line.split(" "))
    
# textFile = open("kout2.txt")
# link = textFile.readlines()
# for line in link:
#     datk.append(line.split(" "))
# 	


# for i in range (len(dato)):#
#     datao[i]=float(dato[i][0])
#     datas[i]=float(dats[i][0])
#     datak[i]=float(datk[i][0])
#     datam[i]=float(datm[i][0])


# #plt.plot(datas,datak,'-')

# #axh.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
# #axh.plot(x,datahi)
# #len(dath)/no of avalanches/xres
    
# fig = plt.figure(figsize=(8,4))

# c=0
# for i in range(len(dath)):#len(dath1)-xres*10,
#     datah[i%xres]=float(dath[i][0])
#     datau[i%xres]=float(datu[i][0])
#     datav[i%xres]=float(datv[i][0])
    
#     if (i%xres==(xres-1) and i!=0):
#         c=c+1
#         if (c%4040==0):
#             '''
#             axh = fig.add_axes([0.15,0.8,0.3,0.15],xlim=(0,1))
# #            plt.title("$\Omega$=2,  $\delta$=25$^\circ$,  timesteps="+str(c*20))
#             plt.ylabel('h')
#             axu = fig.add_axes([0.15,0.5,0.3,0.15],xlim=(0,1))
#             plt.ylabel('u')
#             axv = fig.add_axes([0.15,0.2,0.3,0.15],xlim=(0,1))
#             plt.xlabel('s')
#             plt.ylabel('v')
#             axh.plot(x,datah)
#             axu.plot(x,datau)
#             axv.plot(x,datav)
#             axomega = fig.add_axes([0.6,0.65,0.3,0.15])
#             plt.ylabel('$\Omega$')
#             axmass = fig.add_axes([0.6,0.3,0.3,0.15])
#             plt.ylabel('mass')
#             plt.xlabel('t')
#             axmass.plot(datas[0:c-1],datam[0:c-1],'-')
#             axomega.plot(datas[0:c-1],datao[0:c-1])
#             plt.draw()
#             plt.pause(0.01)
# #            plt.clf()
            
#             '''
#             axh = fig.add_axes([0.1,0.2,0.3,0.75],xlim=(0,0.99))
#             #axh = fig.add_axes([0.1,0.2,0.3,0.75],xlim=(0.925,1.025),ylim=(-0.01,0.01))
#             #axh.grid();
#             plt.ylabel('h')
#             plt.xlabel('s')
#             axh.plot(x,datah)
#             # axomega = fig.add_axes([0.5,0.2,0.45,0.35])
#             # plt.ylabel('$\Omega$')
#             # plt.xlabel('t')
#             # axomega.plot(datas[0:c-1],datao[0:c-1])
#             # axmass = fig.add_axes([0.5,0.6,0.45,0.35])
#             # plt.ylabel('mass')
#             # axmass.plot(datas[0:c-1],datam[0:c-1],'-')
#             plt.draw()
#             plt.pause(0.25)
#             #plt.clf()
#             # if (c>2550):
#             #     break

# #fig = plt.figure(figsize=(5,5))
# #axh = fig.add_axes([0.1,0.1,0.7,0.7],xlim=(0,1))
# #plt.ylabel('h')
# #plt.xlabel('s')
# #axh.plot(x,datah)
# #plt.axis('equal')

# axomega = fig.add_axes([0.5,0.2,0.45,0.35])
# plt.ylabel('$\Omega$')
# plt.xlabel('t')
# axomega.plot(datas[0:c-1],datao[0:c-1])
# axmass = fig.add_axes([0.5,0.6,0.45,0.35])
# plt.ylabel('mass')
# axmass.plot(datas[0:c-1],datam[0:c-1],'-')
# plt.draw()

# #fig = plt.figure(figsize=(2.5,1))
# #axp = fig.add_axes([0.20,0.25,0.7,0.65])
# #axp.plot(datas[0:c-1],datam[0:c-1],'-')

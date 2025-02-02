import os
import sys
import pickle
import argparse
import numpy as np
import pandas as pd
import csv
from Bio import SeqIO
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.cm as cm
import matplotlib.patches as patches
from matplotlib.colors import LinearSegmentedColormap
#import tensorflow as tf
from tensorflow.python.keras.layers import Input,Dense,Flatten,Dropout,Conv1D,Conv1DTranspose,Conv2D,Conv2DTranspose,Conv3D,Conv3DTranspose,Dense,MaxPooling1D,MaxPooling2D,Activation,ReLU,MaxPooling3D,Reshape
from tensorflow.python.keras import initializers, models 
from tensorflow.python.keras.models import Model, load_model

import idr_parm
from idr_parm import pearson_r, pearson_r_loss

aa1=["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","U","X"]
cc1=["-","1","2","3","4","5","6","7","8","9","a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","t","u","v","w","x","y","z","*"]

def generate_cmap(colors,vmax):
    values = range(vmax+1)
    color_list = []
    for v, c in zip(values, colors):
        color_list.append( ( v/ vmax, c) )
    return LinearSegmentedColormap.from_list('custom_cmap', color_list)

parser = argparse.ArgumentParser()
parser.add_argument('infile', help='input fasta file')
parser.add_argument('-ist', default=idr_parm.threshold, help='threshold for interacting site (0.0-1.0)')
parser.add_argument('-pgt', default=idr_parm.threshold_pg, help='threshold for protogroup (0.0-1.0)')
args = parser.parse_args()
if args.ist:
  is_threshold = float(args.ist)
#else:
#  is_threshold = idr_parm.threshold 
if args.pgt:
  pg_threshold = float(args.pgt) 
#else:
#  pg_threshold = idr_parm.threshold_pg

la= []
k=c=0
csize=15
para2=idr_parm.para2

model_1=load_model("idr_aetrain_pg_240420.hd5", custom_objects={"pearson_r": pearson_r, "pearson_r_loss": pearson_r_loss})
model_2 = Model(inputs=model_1.input, outputs=model_1.get_layer('dense_1').output)
model_3 = Model(inputs=model_1.input, outputs=model_1.get_layer('dropout').output)

with open('idr_pca_model.pkl', 'rb') as f0:
  pca = pd.read_pickle(f0)
with open('idr_pca_norm.pkl', 'rb') as f0:
  norm = pd.read_pickle(f0)
v3=np.zeros(idr_parm.numpg+1,dtype=float)
with open('idr_map_landscape.pkl', 'rb') as f2:
  data = pickle.load(f2)
cm0=generate_cmap(['white','dodgerblue'],1)
fig, ax = plt.subplots(dpi=512)
xmin,xmax,ymin,ymax,xbin,ybin,zbin=data[0:7,0,12]
im = ax.imshow(data[:,::-1,0].transpose(1,0),cmap=cm0,interpolation='none',extent=(xmin,xmax,ymin,ymax))
if os.path.isfile('idr_map_pg_scorefreq_name.pkl'):
  with open('idr_map_pg_scorefreq_name.pkl', 'rb') as f0:
    pgl = pickle.load(f0)
  for i in range(0,idr_parm.numpg):
    v3[i]=float(pgl[i][0])
  mv3=sum(v3)
  for i in range(0,idr_parm.numpg):
    v3[i]=v3[i]/max(mv3,1.0)
for record in SeqIO.parse(args.infile,"fasta"):
  name=record.id
  line=record.seq
  c+=1
  if c==1: 
    c=n=o=0
    q=1
    l= np.zeros((20,20,para2), dtype=float)
    lp= np.zeros(len(line),dtype=float)
    lp0=np.zeros(len(line),dtype=float)
    lp1=np.zeros(len(line),dtype=float)
    for i in range(0,len(line)):
      for j in range(i,len(line)):
        if j-i<para2 and aa1.index(line[i])<20 and aa1.index(line[j])<20 :
          l[aa1.index(line[i])][aa1.index(line[j])][j-i]+=1
          n+=1
          for s in range(max(1,j-i-idr_parm.diff_range),min(para2,j-i+idr_parm.diff_range)):
            if s!=(j-i):
              l[aa1.index(line[i])][aa1.index(line[j])][s]+=idr_parm.diff_rate**abs(s-j+i)
          if l[aa1.index(line[i])][aa1.index(line[j])][j-i] > q and i!=j :
            q = l[aa1.index(line[i])][aa1.index(line[j])][j-i]
    if n>0:
      if idr_parm.norm1 == 1:
        n=1
      if idr_parm.norm1 == 2:
        n=q
      if idr_parm.norm3 == 1:
        o=1
      for i in range(0,20):
        for j in range(0,20):
          for m in range(0,para2-1):
            l[i][j][m]=l[i][j][m]/n
            if m==0 and i==j:
              l[i][j][m]=l[i][j][m]*idr_parm.norm4
    l=np.array(l,dtype=float)
    l=l.reshape(1,20*20*para2)
    la0=model_2.predict(l).reshape(1,-1)
    for i in range(norm.shape[1]):
      if norm['std'][i]>0:
        la0[0][i]=(la0[0][i]-norm['mean'][i])/norm['std'][i]
    feature = pca.transform(la0)
    print('seq_',k+1,'\t',"name",'\t','>',name,sep='') 
    print('seq_',k+1,'\t',"seq",'\t',line.rstrip('\n'),sep='')
    print('seq_',k+1,'\t',"ve",sep='',end='\t')
    for i in range(len(la0[0])):
      print("{:.3f}".format(la0[0][i]),end='\t')
    print('\n',end='')
    print('seq_',k+1,'\t',"vepca",sep='',end='\t')
    for i in range(feature.shape[1]):
      print("{:.3f}".format(feature[0][i]),end='\t')
    print('\n',end='')
    ax.add_patch(patches.Circle(xy=(feature[0][(idr_parm.pc1-1)],feature[0][(idr_parm.pc2-1)]),radius=0.005,color='darkred'))
    if idr_parm.showno==1:
       ax.text(feature[0][(idr_parm.pc1-1)],feature[0][(idr_parm.pc2-1)],k+1, ha='left',va='center',fontsize=6,color='darkblue')
    if idr_parm.showno==2:
       ax.text(feature[0][(idr_parm.pc1-1)],feature[0][(idr_parm.pc2-1)],name,ha='left',va='center',fontsize=6,color='darkblue')
    la1=model_3.predict(l)
    for i in range(0,len(line)):
      for j in range(i,len(line)):
        if j-i<para2:
          lp[i]+=la1[0][aa1.index(line[i])*20*para2+aa1.index(line[j])*para2+(j-i)]
          lp[j]+=la1[0][aa1.index(line[i])*20*para2+aa1.index(line[j])*para2+(j-i)]
    for i in range(0,len(line)-1):
      lp0[i]=lp[i]
    for i in range(0,len(line)-1):
      for s in range(max(0,i-idr_parm.diff_range),min(len(line)-1,i+idr_parm.diff_range)):
        lp[s]+=lp0[i]*(idr_parm.diff_rate3**abs(s-i))
    ssm=0.
    for i in range(0,len(line)):
      if (lp[i]*idr_parm.weight)>ssm:
        ssm=lp[i]*idr_parm.weight
    if ssm==0.:
      ssm=1.
    print('seq_',k+1,'\t',"issum",sep='',end='\t')
    for i in range(0,len(line)):
      ss=0.
      if max(lp)-min(lp)!=0:
        ss=(((lp[i]*idr_parm.weight)/(ssm))*idr_parm.scorescale) 
        if lp[i]==0:
          ss=0
      if ss<=0:
        print('-',end="")
      elif ss>=is_threshold:  #36
        print('*',end="")
      else:
        print(cc1[int((ss/is_threshold)*36)],end="")    
    print('\n',end='')
    m=1
    lps=np.argsort(-lp)
    print('seq_',k+1,'\t',"iscol",'\t','rnk','\t','res','\t','score',sep='',end='\n')
    for i in range(0,len(line)):
      if ((lp[lps[i]]*idr_parm.weight)/(ssm))*idr_parm.scorescale < is_threshold:
        break
      print('seq_',k+1,'\t',"islst",'\t',m,'\t',lps[i]+1,line[int(lps[i])],'\t',"{:.3f}".format(((lp[lps[i]]*idr_parm.weight)/(ssm))*idr_parm.scorescale),sep='',end='\n')
      m+=1
    la2=model_1.predict(l)
    mv1=max(sum(la2[0]),1)
    ssm=0.0
    for i in range(0,min(len(la2[0])-1,idr_parm.numpg)):
      if (la2[0][i+1]*idr_parm.weightpg) > ssm:
        ssm=(la2[0][i+1]*idr_parm.weightpg)
    if ssm==0.:
      ssm=1.
    print('seq_',k+1,'\t',"pgsum",sep='',end='\t')
    for i in range(0,min(len(la2[0])-1,idr_parm.numpg)):
      ss=0
      if max(la2[0])-min(la2[0])!=0:
        ss=(((la2[0][i+1]*idr_parm.weightpg)/(ssm))*idr_parm.scorescale)
        if la2[0][i+1]==0:
          ss=0
      if ss<=0:
        print('-',end="")
      elif ss>=pg_threshold: #36
        print('*',end="")
      else:
        print(cc1[int((ss/pg_threshold)*36)],end="")
    m=1
    print('\n',end='')
    la2[0][0]=0.
    la2s=np.argsort(-la2[0]) #@@@@@
    print('seq_',k+1,'\t',"pgcol",'\t','rnk','\t','pg','\t','score','\t','name','\t','smiles',sep="",end="\n")
    for i in range(0,min(len(la2[0])-1,idr_parm.numpg)):
      if ((la2[0][la2s[i]]*idr_parm.weightpg)/(ssm))*idr_parm.scorescale < pg_threshold:
        break
      print('seq_',k+1,'\t',"pglst",'\t',m,'\t',pgl[la2s[i]-1][2],'_',la2s[i],'\t',"{:.3f}".format(((la2[0][la2s[i]]*idr_parm.weightpg)/(ssm))*idr_parm.scorescale),'\t',pgl[la2s[i]-1][3],'\t',pgl[la2s[i]-1][4],sep="",end="\n")
      m+=1
    if 0:
      print("@\t",name.split(' ')[0],"\t",end='')
      for i in range(feature.shape[1]):
        print("{:.5f}".format(feature[0][i]),end='\t')
      print(":\t",end='')
      for i in range(len(la0[0])):
        print("{:.5f}".format(la0[0][i]),end='\t')
      print(":\t",end='')
      print(line.rstrip('\n'))
    k+=1

fig.savefig("idr_map_landscape_01.png")

exit(0)

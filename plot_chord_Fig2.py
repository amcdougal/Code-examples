"""
Create chord plot for Figure 2 in Truebenbach & Darling, 2018, arXiv, 180807130.

Created by: Alex Truebenbach
Last Editted: Sept 2017
"""
from pylab import *
import pyfits
import plotly.plotly as py
#from plotly.offline import download_plotlyjs,init_notebook_mode,plot,iplot
from plotly.graph_objs import *
import numpy as np

#read in objects
ivs,ra,dec,pmRA,pmDE,d,glon,glat,sglon,sglat=loadtxt('100Mpcbin_objects_dHubble1.csv',usecols=[0,1,2,3,4,8,10,11,12,13],dtype={'names':('ivs','ra','dec','pmRA','pmDE','dist','glon','glat','sglon','sglat'),'formats':('S10','f4','f4','f4','f4','f4','f4','f4','f4','f4')},skiprows=1,unpack=True,delimiter=',')

#read in full pairs list
f=pyfits.open('100Mpcbin_PM_dHubble1.fits')
data=f[1].data
ivs_1=data.field('IVS_1')
ivs_2=data.field('IVS_2')
ra_1=data.field('RA_1')
dec_1=data.field('DE_1')
pmRA_1=data.field('pmRA_1')
pmDE_1=data.field('pmDE_1')
dist_1=data.field('dist_1')
ra_2=data.field('RA_2')
dec_2=data.field('DE_2')
pmRA_2=data.field('pmRA_2')
pmDE_2=data.field('pmDE_2')
dist_2=data.field('dist_2')
dtheta=data.field('dtheta')
theta=data.field('theta')

#create list of tuples where the elements are indices of objects in the pair.
E=[]
edge_colors=[]
N=np.zeros(len(ivs)) #number of pairs the object is connected with

for ind in range(len(ivs_1)):
    a=np.where(ivs==ivs_1[ind])
    b=np.where(ivs==ivs_2[ind])
    
    E.append([a[0][0],b[0][0]])
    
    #ribbon colors, blue for converging, red for diverging
    if (dtheta[ind]/sin(theta[ind]*pi/180)) > 0:
        edge_colors.append('rgba(255,0,0,0.75)')
    else:
        edge_colors.append('rgba(0,0,205,0.75)')
    N[a[0][0]]+=1
    N[b[0][0]]+=1

#Label ends as object name
labels=ivs

#make evenly spaced locations around the unit circle
layt=[]
layt_2=[] #for the number of connection in a point. 3 deg lower than the main label at normal layt location
phi=2*pi/len(labels)
theta=phi
for i in range(len(labels)):
    layt.append([cos(theta),sin(theta)])
    print  theta, labels[i]
    if theta<pi/2. or theta > 3*pi/2:
        layt_2.append([cos(theta-0.05),sin(theta-0.05)])
    else:
        layt_2.append([cos(theta+0.05),sin(theta+0.05)])
    theta+=phi

L=len(layt)



def dist (A,B):
    return np.linalg.norm(np.array(A)-np.array(B))

Dist=[0, dist([1,0], 2*[np.sqrt(2)/2]), np.sqrt(2),
      dist([1,0],  [-np.sqrt(2)/2, np.sqrt(2)/2]), 2.0]
params=[1.2, 1.5, 1.8, 2.1]

def get_idx_interv(d, D):
    k=0
    while(d>D[k]): 
        k+=1
    return  k-1


class InvalidInputError(Exception):
    pass

def deCasteljau(b,t): 
    N=len(b) 
    if(N<2):
        raise InvalidInputError("The  control polygon must have at least two points")
    a=np.copy(b) #shallow copy of the list of control points 
    for r in range(1,N): 
        a[:N-r,:]=(1-t)*a[:N-r,:]+t*a[1:N-r+1,:]                             
    return a[0,:]

def BezierCv(b, nr=5):
    t=np.linspace(0, 1, nr)
    return np.array([deCasteljau(b, t[k]) for k in range(nr)])


#Calculate the text angle and position for labels
angles=[]
pos_text=[]
pos_text2=[]
for k in range(L):
    v=layt[k]
    ang=-(180*np.arctan([v[1]/v[0]])/pi)
    angles.append(ang[0])
    pos_text.append([1.2*v[0],1.2*v[1]])
    v=layt_2[k]
    pos_text2.append([1.2*v[0],1.2*v[1]])


import plotly.plotly as py
from plotly.graph_objs import *

line_color=['rgba(0,51.181.0.85)' for v in labels]
node_color=['black' for v in labels]

Xn=[layt[k][0] for k in range(L)]
Yn=[layt[k][1] for k in range(L)]

lines=[]# the list of dicts defining   edge  Plotly attributes
edge_info=[]# the list of points on edges where  the information is placed


for j, e in enumerate(E):
    A=np.array(layt[e[0]])
    B=np.array(layt[e[1]])
    d=dist(A, B)
    K=get_idx_interv(d, Dist)
    b=[A, A/params[K], B/params[K], B]
    color=edge_colors[j]
    pts=BezierCv(b, nr=5)
    mark=deCasteljau(b,0.9)
    edge_info.append(Scatter(x=mark[0], 
                             y=mark[1], 
                             mode='markers', 
                             marker=Marker( size=0.5,  color=edge_colors)
                             )
                    )
    lines.append(Scatter(x=pts[:,0],
                         y=pts[:,1],
                         mode='lines',
                         line=Line(color=color, 
                                  shape='spline',
                                  width=1
                                 ), 
                        hoverinfo='none' 
                       )
                )
    
trace2=Scatter(x=Xn,
           y=Yn,
           mode='markers',
           name='',
           marker=Marker(symbol='dot',
                         size=15, 
                         color=node_color, 
                         line=Line(color=line_color, width=0.5)
                         ),
           text=labels,
           hoverinfo='text',
           )

axis=dict(showline=False, # hide axis line, grid, ticklabels and  title
          zeroline=False,
          showgrid=False,
          showticklabels=False,
          title='' 
          )

width=800
height=850
layout=Layout(font= Font(size=12),
              showlegend=False,
              autosize=False,
              width=width,
              height=height,
              xaxis=XAxis(axis),
              yaxis=YAxis(axis),          
              margin=Margin(l=40,
                            r=40,
                            b=85,
                            t=100,
                          ),
              hovermode='closest'
              )

for k in range(L):
    layout['annotations']+=[Annotation(x=pos_text[k][0],
                      y=pos_text[k][1],
                      text=labels[k],
                      textangle=angles[k],
                      font=Font(size=15,color='rgb(10,10,10)'),
                      showarrow=False)]
    layout['annotations']+=[Annotation(x=pos_text2[k][0],
                      y=pos_text2[k][1],
                      text='N='+str(int(N[k])),
                      textangle=angles[k],
                      font=Font(size=12,color='rgb(10,10,10)'),
                      showarrow=False)]
                            
data=Data(lines+edge_info+[trace2])
fig=Figure(data=data, layout=layout)
py.iplot(fig, filename='pairs_chord_dHubble1') 
py.image.save_as(fig,filename='pairs_chord_dHubble1.pdf')

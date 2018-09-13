"""
Plots the final AIPS CLEANed images. Include object positions and contour maps over images.

Created by: Alex Truebenbach
Last Editted: April 2017
"""
import pyfits,pywcs, csv
from pylab import *
import numpy as np
from scipy import ndimage
from matplotlib.patches import Ellipse
import matplotlib.patches as patches


image1=['0216+011','2305+198','2208-373'] #BT131A
image2=['1304-318','1232+366','1441+522','1625+582'] #BT131B
image3=['1716+686','1721-146','1726+769','1755+055','1839+548','1843+356','1915+657'] #BT131C
image4=['0651+410','0836+182','M81','NGC3862','1145+268','1306+360','M84','1013+127','MRK180','NGC5675','1144+352','1424+240','1623+578','NGC4261','1456+044','1727+502','1754+159','3C390.3'] #BT134A
image5=['2319+444','2324+151','NGC0315','0111+021','0241+622','NGC1218','0415+379'] #BT134B

image_dir='images_temp/'
output_dir='/duanestorage/data/student/altr4657/real-time_cosmology/PMs/papergraphs/'

scale=9.0e-5 #.09 milliarcsec/pix
#Open file of beam sizes
reader=csv.reader(open('beam_sizes.txt','rb'),delimiter=',')
data=[]
data.extend(reader)


i=0
j=0
k=0
for obj in image1:
    #Get beam data
    for row in data:
        if '#BT131A1' in row:
            j=1
        if j==1 and row[0]==obj:
            maj1=float(row[4])
            min1=float(row[5])
            pa1=float(row[6])
            j=0
        if '#BT131A2' in row:
            k=1
        if k==1 and row[0]==obj:
            maj2=float(row[4])
            min2=float(row[5])
            pa2=float(row[6])
            k=0
    
    f=pyfits.open('BT131A1/'+obj+'_image.fits')
    img1=f[0].data
    h1=f[0].header
    img1=img1[0][0]
    f=pyfits.open('BT131A2/'+obj+'_image.fits')
    img2=f[0].data
    h2=f[0].header
    img2=img2[0][0]
    
    
    #Crop image -> 250x250 pixels, 22.5 milliarcsec per side
    img1=img1[131:381,131:381]
    img2=img2[131:381,131:381]
    img1_max=img1.max(axis=None)
    img2_max=img2.max(axis=None)
    img1_sigma=img1.std(axis=None)
    img2_sigma=img2.std(axis=None)
    r1=.8*img1_max/(3*img1_sigma)
    r2=.8*img2_max/(3*img2_sigma)
    cont_levels1=[3*img1_sigma,(r1**.25)*3*img1_sigma,(r1**.5)*3*img1_sigma,(r1**.75)*3*img1_sigma,r1*3*img1_sigma]
    cont_levels2=[3*img2_sigma,(r2**.25)*3*img2_sigma,(r2**.5)*3*img2_sigma,(r2**.75)*3*img2_sigma,r2*3*img2_sigma]

    #Find image center
    ra_c1=h1['crval1']
    de_c1=h1['crval2']
    #convert to sexigesimal
    ra=ra_c1/15 #hours
    min,sec=divmod(ra*3600,60)
    hr,min=divmod(min,60)
    if hr<10:
        ra1='0'+str(int(hr))
    else:
        ra1=str(int(hr))
    if min<10:
        ra1=ra1+':0'+str(int(min))
    else:
        ra1=ra1+':'+str(int(min))
    if sec<10:
        ra1=ra1+':0'+str(round(sec,6)).ljust(8,'0')
    else:
        ra1=ra1+':'+str(round(sec,6)).ljust(9,'0')
    #convert to sexigesimal
    if de_c1<0:
        dec1='-'
    else:
        dec1='+'
    min,sec=divmod(fabs(de_c1)*3600,60)
    hr,min=divmod(min,60)
    if hr<10:
        dec1=dec1+'0'+str(int(hr))
    else:
        dec1=dec1+str(int(hr))
    if min<10:
        dec1=dec1+':0'+str(int(min))
    else:
        dec1=dec1+':'+str(int(min))
    if sec<10:
        dec1=dec1+':0'+str(round(sec,6)).ljust(8,'0')
    else:
        dec1=dec1+':'+str(round(sec,6)).ljust(9,'0')

        
    #Set up plot
    f=plt.figure(1)
    ax1=plt.subplot(121,aspect='equal')
    ax2=plt.subplot(122,aspect='equal')

    #Plot image
    ax1.imshow(img1,cmap=cm.gray_r)
    ax2.imshow(img2,cmap=cm.gray_r)

    cs=ax1.contour(range(250),range(250),img1,levels=cont_levels1,colors='black')
    cs=ax2.contour(range(250),range(250),img2,levels=cont_levels2,colors='black')
    
    #Add arcsec labels
    ax1.set_xticks([13.9,69.4,125,180.6,236.11])
    ax1.set_xticklabels([10,5,0,-5,-10])
    ax1.set_yticks([13.9,69.4,125,180.6,236.11])
    ax1.set_yticklabels([10,5,0,-5,-10])
    ax2.set_xticks([13.9,69.4,125,180.6,236.11])
    ax2.set_xticklabels([10,5,0,-5,-10])
    ax2.set_yticks([13.9,69.4,125,180.6,236.11])
    ax2.set_yticklabels([10,5,0,-5,-10])
    ax1.set_xlabel('Milliarcseconds')
    ax2.set_xlabel('Milliarcseconds')
    ax1.set_ylabel('Milliarcseconds')
    
    ax1.text(10,25,obj,color='blue',fontsize='large')
    ax1.text(10,-5,'Center: RA '+ra1+'   Dec '+dec1,fontsize='small')
    if i==0:
        ax1.text(220,25,'A',color='blue',fontsize='large')
        ax2.text(220,25,'B',color='blue',fontsize='large')

    #Add ellipse of beam size
    if pa1 >= 90:
        y=pa1-90
    else:
        y=pa1
    x=np.max([sin(y*pi/180)*maj1,cos(y*pi/180)*maj1])
    ax1.add_patch(
        patches.Rectangle(
            (0,250-1.5*x),
            1.5*x,
            1.5*x,
            fill=False
            )
        )
    ell=Ellipse(xy=(.75*x,250-.75*x),width=maj1,height=min1,angle=-1*(90-pa1))
    ax1.add_artist(ell)
    ell.set_facecolor('white')

    #Add ellipse of beam size
    if pa2 >= 90:
        y=pa2-90
    else:
        y=pa2
    x=np.max([sin(y*pi/180)*maj2,cos(y*pi/180)*maj2])
    ax2.add_patch(
        patches.Rectangle(
            (0,250-1.5*x),
            1.5*x,
            1.5*x,
            fill=False
            )
        )
    ell2=Ellipse(xy=(.75*x,250-.75*x),width=maj2,height=min2,angle=-1*(90-pa2))
    ax2.add_artist(ell2)
    ell2.set_facecolor('white')
    
    savefig(output_dir+obj+'.eps')
    show()
    clf()
    i+=1

    

for obj in image2:
    #Get beam data
    for row in data:
        if '#BT131B1' in row:
            j=1
        if j==1 and row[0]==obj:
            maj1=float(row[4])
            min1=float(row[5])
            pa1=float(row[6])
            j=0
        if '#BT131B2' in row:
            k=1
        if k==1 and row[0]==obj:
            maj2=float(row[4])
            min2=float(row[5])
            pa2=float(row[6])
            k=0
    f=pyfits.open('BT131B1/'+obj+'_image.fits')
    img1=f[0].data
    h1=f[0].header
    img1=img1[0][0]
    f=pyfits.open('BT131B2/'+obj+'_image.fits')
    img2=f[0].data
    h2=f[0].header
    img2=img2[0][0]
    
    
    #Crop image -> 250x250 pixels, 22.5 milliarcsec per side
    img1=img1[131:381,131:381]
    img2=img2[131:381,131:381]
    img1_max=img1.max(axis=None)
    img2_max=img2.max(axis=None)
    img1_sigma=img1.std(axis=None)
    img2_sigma=img2.std(axis=None)
    r1=.8*img1_max/(3*img1_sigma)
    r2=.8*img2_max/(3*img2_sigma)
    cont_levels1=[3*img1_sigma,(r1**.25)*3*img1_sigma,(r1**.5)*3*img1_sigma,(r1**.75)*3*img1_sigma,r1*3*img1_sigma]
    cont_levels2=[3*img2_sigma,(r2**.25)*3*img2_sigma,(r2**.5)*3*img2_sigma,(r2**.75)*3*img2_sigma,r2*3*img2_sigma]

    #Find image center
    ra_c1=h1['crval1']
    de_c1=h1['crval2']
    #convert to sexigesimal
    ra=ra_c1/15 #hours
    min,sec=divmod(ra*3600,60)
    hr,min=divmod(min,60)
    if hr<10:
        ra1='0'+str(int(hr))
    else:
        ra1=str(int(hr))
    if min<10:
        ra1=ra1+':0'+str(int(min))
    else:
        ra1=ra1+':'+str(int(min))
    if sec<10:
        ra1=ra1+':0'+str(round(sec,6)).ljust(8,'0')
    else:
        ra1=ra1+':'+str(round(sec,6)).ljust(9,'0')
    #convert to sexigesimal
    if de_c1<0:
        dec1='-'
    else:
        dec1='+'
    min,sec=divmod(fabs(de_c1)*3600,60)
    hr,min=divmod(min,60)
    if hr<10:
        dec1=dec1+'0'+str(int(hr))
    else:
        dec1=dec1+str(int(hr))
    if min<10:
        dec1=dec1+':0'+str(int(min))
    else:
        dec1=dec1+':'+str(int(min))
    if sec<10:
        dec1=dec1+':0'+str(round(sec,6)).ljust(8,'0')
    else:
        dec1=dec1+':'+str(round(sec,6)).ljust(9,'0')

        
    #Set up plot
    f=plt.figure(1)
    ax1=plt.subplot(121,aspect='equal')
    ax2=plt.subplot(122,aspect='equal')

    #Plot image
    ax1.imshow(img1,cmap=cm.gray_r)
    ax2.imshow(img2,cmap=cm.gray_r)

    cs=ax1.contour(range(250),range(250),img1,levels=cont_levels1,colors='black')
    cs=ax2.contour(range(250),range(250),img2,levels=cont_levels2,colors='black')
    
    #Add arcsec labels
    ax1.set_xticks([13.9,69.4,125,180.6,236.11])
    ax1.set_xticklabels([10,5,0,-5,-10])
    ax1.set_yticks([13.9,69.4,125,180.6,236.11])
    ax1.set_yticklabels([10,5,0,-5,-10])
    ax2.set_xticks([13.9,69.4,125,180.6,236.11])
    ax2.set_xticklabels([10,5,0,-5,-10])
    ax2.set_yticks([13.9,69.4,125,180.6,236.11])
    ax2.set_yticklabels([10,5,0,-5,-10])
    ax1.set_xlabel('Milliarcseconds')
    ax2.set_xlabel('Milliarcseconds')
    ax1.set_ylabel('Milliarcseconds')
    
    ax1.text(10,25,obj,color='blue',fontsize='large')
    ax1.text(10,-5,'Center: RA '+ra1+'   Dec '+dec1,fontsize='small')

    #Add ellipse of beam size
    if pa1 >= 90:
        y=pa1-90
    else:
        y=pa1
    x=np.max([sin(y*pi/180)*maj1,cos(y*pi/180)*maj1])
    ax1.add_patch(
        patches.Rectangle(
            (0,250-1.5*x),
            1.5*x,
            1.5*x,
            fill=False
            )
        )
    ell=Ellipse(xy=(.75*x,250-.75*x),width=maj1,height=min1,angle=-1*(90-pa1))
    ax1.add_artist(ell)
    ell.set_facecolor('white')

    #Add ellipse of beam size
    if pa2 >= 90:
        y=pa2-90
    else:
        y=pa2
    x=np.max([sin(y*pi/180)*maj2,cos(y*pi/180)*maj2])
    ax2.add_patch(
        patches.Rectangle(
            (0,250-1.5*x),
            1.5*x,
            1.5*x,
            fill=False
            )
        )
    ell2=Ellipse(xy=(.75*x,250-.75*x),width=maj2,height=min2,angle=-1*(90-pa2))
    ax2.add_artist(ell2)
    ell2.set_facecolor('white')
    
    savefig(output_dir+obj+'.eps')
    #show()
    clf()
    i+=1
    
for obj in image3:
    #Get beam data
    for row in data:
        if '#BT131C1' in row:
            j=1
        if j==1 and row[0]==obj:
            maj1=float(row[4])
            min1=float(row[5])
            pa1=float(row[6])
            j=0
        if '#BT131C2' in row:
            k=1
        if k==1 and row[0]==obj:
            maj2=float(row[4])
            min2=float(row[5])
            pa2=float(row[6])
            k=0
    f=pyfits.open('BT131C1/'+obj+'_image.fits')
    img1=f[0].data
    h1=f[0].header
    img1=img1[0][0]
    f=pyfits.open('BT131C2/'+obj+'_image.fits')
    img2=f[0].data
    h2=f[0].header
    img2=img2[0][0]
    
    
    #Crop image -> 250x250 pixels, 22.5 milliarcsec per side
    img1=img1[131:381,131:381]
    img2=img2[131:381,131:381]
    img1_max=img1.max(axis=None)
    img2_max=img2.max(axis=None)
    img1_sigma=img1.std(axis=None)
    img2_sigma=img2.std(axis=None)
    r1=.8*img1_max/(3*img1_sigma)
    r2=.8*img2_max/(3*img2_sigma)
    cont_levels1=[3*img1_sigma,(r1**.25)*3*img1_sigma,(r1**.5)*3*img1_sigma,(r1**.75)*3*img1_sigma,r1*3*img1_sigma]
    cont_levels2=[3*img2_sigma,(r2**.25)*3*img2_sigma,(r2**.5)*3*img2_sigma,(r2**.75)*3*img2_sigma,r2*3*img2_sigma]

    #Find image center
    ra_c1=h1['crval1']
    de_c1=h1['crval2']
    #convert to sexigesimal
    ra=ra_c1/15 #hours
    min,sec=divmod(ra*3600,60)
    hr,min=divmod(min,60)
    if hr<10:
        ra1='0'+str(int(hr))
    else:
        ra1=str(int(hr))
    if min<10:
        ra1=ra1+':0'+str(int(min))
    else:
        ra1=ra1+':'+str(int(min))
    if sec<10:
        ra1=ra1+':0'+str(round(sec,6)).ljust(8,'0')
    else:
        ra1=ra1+':'+str(round(sec,6)).ljust(9,'0')
    #convert to sexigesimal
    if de_c1<0:
        dec1='-'
    else:
        dec1='+'
    min,sec=divmod(fabs(de_c1)*3600,60)
    hr,min=divmod(min,60)
    if hr<10:
        dec1=dec1+'0'+str(int(hr))
    else:
        dec1=dec1+str(int(hr))
    if min<10:
        dec1=dec1+':0'+str(int(min))
    else:
        dec1=dec1+':'+str(int(min))
    if sec<10:
        dec1=dec1+':0'+str(round(sec,6)).ljust(8,'0')
    else:
        dec1=dec1+':'+str(round(sec,6)).ljust(9,'0')

        
    #Set up plot
    f=plt.figure(1)
    ax1=plt.subplot(121,aspect='equal')
    ax2=plt.subplot(122,aspect='equal')

    #Plot image
    ax1.imshow(img1,cmap=cm.gray_r)
    ax2.imshow(img2,cmap=cm.gray_r)

    cs=ax1.contour(range(250),range(250),img1,levels=cont_levels1,colors='black')
    cs=ax2.contour(range(250),range(250),img2,levels=cont_levels2,colors='black')
    
    #Add arcsec labels
    ax1.set_xticks([13.9,69.4,125,180.6,236.11])
    ax1.set_xticklabels([10,5,0,-5,-10])
    ax1.set_yticks([13.9,69.4,125,180.6,236.11])
    ax1.set_yticklabels([10,5,0,-5,-10])
    ax2.set_xticks([13.9,69.4,125,180.6,236.11])
    ax2.set_xticklabels([10,5,0,-5,-10])
    ax2.set_yticks([13.9,69.4,125,180.6,236.11])
    ax2.set_yticklabels([10,5,0,-5,-10])
    ax1.set_xlabel('Milliarcseconds')
    ax2.set_xlabel('Milliarcseconds')
    ax1.set_ylabel('Milliarcseconds')
    
    ax1.text(10,25,obj,color='blue',fontsize='large')
    ax1.text(10,-5,'Center: RA '+ra1+'   Dec '+dec1,fontsize='small')
    

    #Add ellipse of beam size
    if pa1 >= 90:
        y=pa1-90
    else:
        y=pa1
    x=np.max([sin(y*pi/180)*maj1,cos(y*pi/180)*maj1])
    ax1.add_patch(
        patches.Rectangle(
            (0,250-1.5*x),
            1.5*x,
            1.5*x,
            fill=False
            )
        )
    ell=Ellipse(xy=(.75*x,250-.75*x),width=maj1,height=min1,angle=-1*(90-pa1))
    ax1.add_artist(ell)
    ell.set_facecolor('white')

    #Add ellipse of beam size
    if pa2 >= 90:
        y=pa2-90
    else:
        y=pa2
    x=np.max([sin(y*pi/180)*maj2,cos(y*pi/180)*maj2])
    ax2.add_patch(
        patches.Rectangle(
            (0,250-1.5*x),
            1.5*x,
            1.5*x,
            fill=False
            )
        )
    ell2=Ellipse(xy=(.75*x,250-.75*x),width=maj2,height=min2,angle=-1*(90-pa2))
    ax2.add_artist(ell2)
    ell2.set_facecolor('white')
    
    savefig(output_dir+obj+'.eps')
    #show()
    clf()
    i+=1
    

for obj in image4:
    #Get beam data
    for row in data:
        if '#BT134A2' in row:
            j=1
        if j==1 and row[0]==obj:
            maj1=float(row[4])
            min1=float(row[5])
            pa1=float(row[6])
            j=0
        if '#BT134A3' in row:
            k=1
        if k==1 and row[0]==obj:
            maj2=float(row[4])
            min2=float(row[5])
            pa2=float(row[6])
            k=0
    
    f=pyfits.open('BT134A2/'+obj+'_image.fits')
    img1=f[0].data
    h1=f[0].header
    img1=img1[0][0]
    f=pyfits.open('BT134A3/'+obj+'_image.fits')
    img2=f[0].data
    h2=f[0].header
    img2=img2[0][0]
    
    
    #Crop image -> 250x250 pixels, 22.5 milliarcsec per side
    img1=img1[131:381,131:381]
    img2=img2[131:381,131:381]
    img1_max=img1.max(axis=None)
    img2_max=img2.max(axis=None)
    img1_sigma=img1.std(axis=None)
    img2_sigma=img2.std(axis=None)
    r1=.8*img1_max/(3*img1_sigma)
    r2=.8*img2_max/(3*img2_sigma)
    cont_levels1=[3*img1_sigma,(r1**.25)*3*img1_sigma,(r1**.5)*3*img1_sigma,(r1**.75)*3*img1_sigma,r1*3*img1_sigma]
    cont_levels2=[3*img2_sigma,(r2**.25)*3*img2_sigma,(r2**.5)*3*img2_sigma,(r2**.75)*3*img2_sigma,r2*3*img2_sigma]

    #Find image center
    ra_c1=h1['crval1']
    de_c1=h1['crval2']
    #convert to sexigesimal
    ra=ra_c1/15 #hours
    min,sec=divmod(ra*3600,60)
    hr,min=divmod(min,60)
    if hr<10:
        ra1='0'+str(int(hr))
    else:
        ra1=str(int(hr))
    if min<10:
        ra1=ra1+':0'+str(int(min))
    else:
        ra1=ra1+':'+str(int(min))
    if sec<10:
        ra1=ra1+':0'+str(round(sec,6)).ljust(8,'0')
    else:
        ra1=ra1+':'+str(round(sec,6)).ljust(9,'0')
    #convert to sexigesimal
    if de_c1<0:
        dec1='-'
    else:
        dec1='+'
    min,sec=divmod(fabs(de_c1)*3600,60)
    hr,min=divmod(min,60)
    if hr<10:
        dec1=dec1+'0'+str(int(hr))
    else:
        dec1=dec1+str(int(hr))
    if min<10:
        dec1=dec1+':0'+str(int(min))
    else:
        dec1=dec1+':'+str(int(min))
    if sec<10:
        dec1=dec1+':0'+str(round(sec,6)).ljust(8,'0')
    else:
        dec1=dec1+':'+str(round(sec,6)).ljust(9,'0')

        
    #Set up plot
    f=plt.figure(1)
    ax1=plt.subplot(121,aspect='equal')
    ax2=plt.subplot(122,aspect='equal')

    #Plot image
    ax1.imshow(img1,cmap=cm.gray_r)
    ax2.imshow(img2,cmap=cm.gray_r)

    cs=ax1.contour(range(250),range(250),img1,levels=cont_levels1,colors='black')
    cs=ax2.contour(range(250),range(250),img2,levels=cont_levels2,colors='black')
    
    #Add arcsec labels
    ax1.set_xticks([13.9,69.4,125,180.6,236.11])
    ax1.set_xticklabels([10,5,0,-5,-10])
    ax1.set_yticks([13.9,69.4,125,180.6,236.11])
    ax1.set_yticklabels([10,5,0,-5,-10])
    ax2.set_xticks([13.9,69.4,125,180.6,236.11])
    ax2.set_xticklabels([10,5,0,-5,-10])
    ax2.set_yticks([13.9,69.4,125,180.6,236.11])
    ax2.set_yticklabels([10,5,0,-5,-10])
    ax1.set_xlabel('Milliarcseconds')
    ax2.set_xlabel('Milliarcseconds')
    ax1.set_ylabel('Milliarcseconds')
    
    ax1.text(10,25,obj,color='blue',fontsize='large')
    ax1.text(10,-5,'Center: RA '+ra1+'   Dec '+dec1,fontsize='small')
    

    #Add ellipse of beam size
    if pa1 >= 90:
        y=pa1-90
    else:
        y=pa1
    x=np.max([sin(y*pi/180)*maj1,cos(y*pi/180)*maj1])
    ax1.add_patch(
        patches.Rectangle(
            (0,250-1.5*x),
            1.5*x,
            1.5*x,
            fill=False
            )
        )
    ell=Ellipse(xy=(.75*x,250-.75*x),width=maj1,height=min1,angle=-1*(90-pa1))
    ax1.add_artist(ell)
    ell.set_facecolor('white')

    #Add ellipse of beam size
    if pa2 >= 90:
        y=pa2-90
    else:
        y=pa2
    x=np.max([sin(y*pi/180)*maj2,cos(y*pi/180)*maj2])
    ax2.add_patch(
        patches.Rectangle(
            (0,250-1.5*x),
            1.5*x,
            1.5*x,
            fill=False
            )
        )
    ell2=Ellipse(xy=(.75*x,250-.75*x),width=maj2,height=min2,angle=-1*(90-pa2))
    ax2.add_artist(ell2)
    ell2.set_facecolor('white')
    
    savefig(output_dir+obj+'.eps')
    #show()
    clf()
    i+=1


for obj in image5:
    #Get beam data
    for row in data:
        if '#BT134B1' in row:
            j=1
        if j==1 and row[0]==obj:
            maj1=float(row[4])
            min1=float(row[5])
            pa1=float(row[6])
            j=0
        if '#BT134B2' in row:
            k=1
        if k==1 and row[0]==obj:
            maj2=float(row[4])
            min2=float(row[5])
            pa2=float(row[6])
            k=0
    
    f=pyfits.open('BT134B1/'+obj+'_image.fits')
    img1=f[0].data
    h1=f[0].header
    img1=img1[0][0]
    """
    f=pyfits.open('BT134B2/'+obj+'_image.fits')
    img2=f[0].data
    h2=f[0].header
    img2=img2[0][0]
    """
    
    #Crop image -> 250x250 pixels, 22.5 milliarcsec per side
    img1=img1[131:381,131:381]
    #img2=img2[131:381,131:381]
    img1_max=img1.max(axis=None)
    #img2_max=img2.max(axis=None)
    img1_sigma=img1.std(axis=None)
    #img2_sigma=img2.std(axis=None)
    r1=.8*img1_max/(3*img1_sigma)
    #r2=.8*img2_max/(3*img2_sigma)
    cont_levels1=[3*img1_sigma,(r1**.25)*3*img1_sigma,(r1**.5)*3*img1_sigma,(r1**.75)*3*img1_sigma,r1*3*img1_sigma]
    #cont_levels2=[3*img2_sigma,(r2**.25)*3*img2_sigma,(r2**.5)*3*img2_sigma,(r2**.75)*3*img2_sigma,r2*3*img2_sigma]

    #Find image center
    ra_c1=h1['crval1']
    de_c1=h1['crval2']
    #convert to sexigesimal
    ra=ra_c1/15 #hours
    min,sec=divmod(ra*3600,60)
    hr,min=divmod(min,60)
    if hr<10:
        ra1='0'+str(int(hr))
    else:
        ra1=str(int(hr))
    if min<10:
        ra1=ra1+':0'+str(int(min))
    else:
        ra1=ra1+':'+str(int(min))
    if sec<10:
        ra1=ra1+':0'+str(round(sec,6)).ljust(8,'0')
    else:
        ra1=ra1+':'+str(round(sec,6)).ljust(9,'0')
    #convert to sexigesimal
    if de_c1<0:
        dec1='-'
    else:
        dec1='+'
    min,sec=divmod(fabs(de_c1)*3600,60)
    hr,min=divmod(min,60)
    if hr<10:
        dec1=dec1+'0'+str(int(hr))
    else:
        dec1=dec1+str(int(hr))
    if min<10:
        dec1=dec1+':0'+str(int(min))
    else:
        dec1=dec1+':'+str(int(min))
    if sec<10:
        dec1=dec1+':0'+str(round(sec,6)).ljust(8,'0')
    else:
        dec1=dec1+':'+str(round(sec,6)).ljust(9,'0')

        
    #Set up plot
    f=plt.figure(1)
    ax1=plt.subplot(121,aspect='equal')
    #ax2=plt.subplot(122,aspect='equal')

    #Plot image
    ax1.imshow(img1,cmap=cm.gray_r)
    #ax2.imshow(img2,cmap=cm.gray_r)

    cs=ax1.contour(range(250),range(250),img1,levels=cont_levels1,colors='black')
    #cs=ax2.contour(range(250),range(250),img2,levels=cont_levels2,colors='black')
    
    #Add arcsec labels
    ax1.set_xticks([13.9,69.4,125,180.6,236.11])
    ax1.set_xticklabels([10,5,0,-5,-10])
    ax1.set_yticks([13.9,69.4,125,180.6,236.11])
    ax1.set_yticklabels([10,5,0,-5,-10])
    #ax2.set_xticks([13.9,69.4,125,180.6,236.11])
    #ax2.set_xticklabels([10,5,0,-5,-10])
    #ax2.set_yticks([13.9,69.4,125,180.6,236.11])
    #ax2.set_yticklabels([10,5,0,-5,-10])
    ax1.set_xlabel('Milliarcseconds')
    #ax2.set_xlabel('Milliarcseconds')
    ax1.set_ylabel('Milliarcseconds')
    
    ax1.text(10,25,obj,color='blue',fontsize='large')
    ax1.text(10,-5,'Center: RA '+ra1+'   Dec '+dec1,fontsize='small')
    

    #Add ellipse of beam size
    if pa1 >= 90:
        y=pa1-90
    else:
        y=pa1
    x=np.max([sin(y*pi/180)*maj1,cos(y*pi/180)*maj1])
    ax1.add_patch(
        patches.Rectangle(
            (0,250-1.5*x),
            1.5*x,
            1.5*x,
            fill=False
            )
        )
    ell=Ellipse(xy=(.75*x,250-.75*x),width=maj1,height=min1,angle=-1*(90-pa1))
    ax1.add_artist(ell)
    ell.set_facecolor('white')
    """
    #Add ellipse of beam size
    if pa2 >= 90:
        y=pa2-90
    else:
        y=pa2
    x=np.max([sin(y*pi/180)*maj2,cos(y*pi/180)*maj2])
    ax2.add_patch(
        patches.Rectangle(
            (0,250-1.5*x),
            1.5*x,
            1.5*x,
            fill=False
            )
        )
    ell2=Ellipse(xy=(.75*x,250-.75*x),width=maj2,height=min2,angle=-1*(90-pa2))
    ax2.add_artist(ell2)
    ell2.set_facecolor('white')
    """
    savefig(output_dir+obj+'.eps')
    #show()
    clf()
    i+=1

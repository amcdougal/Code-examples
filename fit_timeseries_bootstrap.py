"""
Calculate proper motions with linear regression from time series data. Eliminate outliers through bootstrapping and 3-sigma clipping. Creates Table 6 and Figure 3 in Truebenbach & Darling, 2017, ApJS, 233, 3.

This program reads in the time series solution output by Calc/Solve and fits a proper motion vector to each source.  Excludes the 39 sources with nonlinear PMs listed in the file exclude.tbl.  Also excludes all sessions before 1990, all sessions with fewer than 3 group delays and sources with <10 sessions and/or  <10 yrs of observations. Fits a line by analytically calculating the least squares estimators.  Then do 500 rounds of bootstrapping to assess influence of individual points. Output RA pms are multiplied by cos(dec).

Created by: Alex Truebenbach
Last Editted: Jan 2017
"""
from scipy import *
from pylab import *
import csv,os
import numpy as np
from scipy import optimize
from scipy.stats import chi2,gaussian_kde
from numpy import random
import matplotlib.pyplot as plt
from astropy.time import Time
import matplotlib.ticker as ticker

file='gsf2017a.csv'
Ts=5

def convert_coord(ra,dec):
    #RA and Dec are each a list of 3 numbers, [h,m,s] for ra, [d,m,s] for dec
    r=(float(ra[0])+(float(ra[1])+(float(ra[2])/60))/60)*15
    if float(dec[0])>0:
        d=float(dec[0])+(float(dec[1])+(float(dec[2])/60))/60
    elif float(dec[0])==0.0:
        if '-' in str(dec[0]):
            d=-1*(fabs(float(dec[0]))+(float(dec[1])+(float(dec[2])/60.))/60.)
        else:
            d=float(dec[0])+(float(dec[1])+(float(dec[2])/60.))/60.
    else:
        d=-1*(fabs(float(dec[0]))+(float(dec[1])+(float(dec[2])/60))/60)
    return r,d

def comp_chisq(obs,model,w):
    """
    Calculate chi squared for a model, given a set of data and errors.
    Inputs
    ------
    model=array of model values
    obs=array of observed valurd
    w=array of errors for each observed value

    Output
    -----
    cs=float, the chi squared value for this model.
    Notes
    ------
    model,obs, and w must all be the same length.
    """
    if len(obs)!=len(model) or len(model)!=len(w) or len(obs)!=len(w):
        print 'ERROR! Arrays input in comp_chisq are different lengths!'
    cs=0
    for findex in range(len(model)):
        cs+=((obs[findex]-model[findex])**2/w[findex]**2)

    return cs

#Bunch of functions for finding the MLE for a line.
s= lambda y_e: np.sum(1/y_e**2)
sx= lambda y_e, x: np.sum(x/y_e**2)
sy= lambda y_e, y: np.sum(y/y_e**2)
sx2= lambda y_e, x: np.sum(x**2/y_e**2)
sxy= lambda y_e, x,y: np.sum((x*y)/y_e**2)
delta= lambda y_e, x: s(y_e)*sx2(y_e,x)-sx(y_e,x)**2
fit=lambda y_e,x,y: [
    (s(y_e)*sxy(y_e,x,y)-sx(y_e,x)*sy(y_e,y))/delta(y_e,x),
    (sx2(y_e,x)*sy(y_e,y)-sx(y_e,x)*sxy(y_e,x,y))/delta(y_e,x)] #slope,yinter in units of arcsec/day and arcsec
e_fit=lambda y_e,x,y: [
    [s(y_e)/delta(y_e,x), -1*sx(y_e,x)/delta(y_e,x)],
    [-1*sx(y_e,x)/delta(y_e,x), sx2(y_e,x)/delta(y_e,x)]] #covariance matrix

func=lambda m,yint, x: m*x+yint
avg=lambda x: np.sum(x)/len(x)
cov=lambda x,y: np.sum((x-avg(x))*(y-avg(y)))/(len(x)-1)



#Read in the list of objects to be excluded
reader=csv.reader(open('exclude.tbl','rb'))
data=[]
data.extend(reader)
exclude=[]
data=data[2:]
for row in data:
    exclude.append(row[0])

#Read in objects I've observed with the VLBA
reader=csv.reader(open('VLBA_PMs.csv','rb'))
data=[]
data.extend(reader)
new_id=[]
new_epoch=[]
new_ra=[]
new_era=[]
new_de=[]
new_ede=[]
data.pop(0)
for row in data:
    if len(row)==0:
        continue
    if '#' in row[0]:
        continue
    line=row[0].split(' ')
    line=filter(None,line)
    new_id.append(line[0])
    new_epoch.append(float(line[1]))
    new_ra.append(float(line[2]))
    new_era.append(float(line[3]))
    new_de.append(float(line[4]))
    new_ede.append(float(line[5]))
    
#Open time series file.    
reader=csv.reader(open(file,'rb'))
data=[]
data.extend(reader)
data.pop(0)
ivs=[]
pm_ra=[] #muas/yr
pm_dec=[]
epm_ra=[]
epm_dec=[]

total_objects = 0

#Open output file
f=open('2017a_PMs.csv','wb')
writer=csv.writer(f,delimiter=',')
writer.writerow(['IVS','RA','DE','e_RA','e_DE','PM_RA','e_PM_RA','RA_cs','PTE_RA','N_RA','PM_DE','e_PM_DE','DE_cs','PTE_DE','N_DE','RA0','DE0','Length','lastobs','cov_ra', 'cov_de','Source']) #RA0 and DE0 are the y intercepts for the RA and Dec.  Aka the objects location at MJD=0.  RA and DE and the most recent coordinates of the object. Lastobs is julian date of last observation. cov is the entire covariance matrix for the fit, they are in the original fit units.

#Open file to record rejected objects
g=open('2017a_rejects.csv','wb')
w2=csv.writer(g,delimiter=',')
w2.writerow(['IVS','RA','DE','N','lastobs','Length'])

ivs_i=''
time=[]
ra=[]
dec=[]
e_ra=[]
e_dec=[]
corr=[]
n=0

for row in data:
    line=row
    if line[0]==ivs_i: #while still on the same object
        #Exclude dates before 1990 and sessions with <3 group delays and objects in the exclude list
        if float(line[11])<47892.5 or float(line[12])<3 or line[0] in exclude:
            continue
        time.append(float(line[11])) #days
        r=[float(line[2]),float(line[3]),float(line[4])]
        d=[float(line[5]),float(line[6]),float(line[7])]
        r,d=convert_coord(r,d)
        ra.append(r*3600) #arcsec
        dec.append(d*3600)
        e_ra.append(float(line[8])*15) #arcssec
        e_dec.append(float(line[9]))
        corr.append(float(line[10]))
        n+=1
    #Fit a line to the object's position and date
    else:
        total_objects +=1
        ra=np.array(ra)
        dec=np.array(dec)
        e_ra=np.array(e_ra)
        e_dec=np.array(e_dec)
        time=np.array(time)
        if ivs_i in new_id: #object in new observations
            ra_temp=0
            de_temp=0
            t_time=0
            era_temp=0
            ede_temp=0
            n_temp=0
            for index in range(len(new_id)):
                if new_id[index] == ivs_i:
                    #average together the two observations
                    ra_temp+=new_ra[index]*3600
                    de_temp+=new_de[index]*3600
                    t_time+=new_epoch[index]
                    era_temp+=new_era[index]**2
                    ede_temp+=new_ede[index]**2
                    n_temp+=1
            ra=np.append(ra,[ra_temp/n_temp])
            dec=np.append(dec,[de_temp/n_temp])
            e_ra=np.append(e_ra,[sqrt(era_temp)/n_temp])
            e_dec=np.append(e_dec,[sqrt(ede_temp)/n_temp])
            time=np.append(time,[t_time/n_temp])

            #Fit line to data. No clipping or bootstrapping
            fit_ra=fit(e_ra,time,ra)
            e_fit_ra=e_fit(e_ra,time,ra)
            fit_de=fit(e_dec,time,dec)
            e_fit_de=e_fit(e_dec,time,dec)

            
            #Plot fits
            dec_temp=(dec[dec.size-1]/3600)*pi/180
            
            
            #Plot timeseries with bestfit line
            clf()
            rcParams['axes.linewidth']=2
            #convert mjd to years
            time_yr=Time(time,format='mjd').decimalyear
            #find new y intercept with these new units.
            def find_b_ra(x,b):
                return fit_ra[0]*10**6*365.25*x + b
            def find_b_de(x,b):
                return fit_de[0]*10**6*365.25*x + b
            b_ra,c=optimize.curve_fit(find_b_ra,time_yr,(ra-ra[-1])*10**6)
            b_de,c=optimize.curve_fit(find_b_de,time_yr,(dec-dec[-1])*10**6)
            
            fig, (ax1,ax2)=plt.subplots(2,sharex=True)
            
            ax1.errorbar(time_yr[0:-1],(ra[0:-1]-ra[-1])*10**6,yerr=e_ra[0:-1]*10**6,fmt='.',mfc='blue',ecolor="black",elinewidth=2,ms=15)
            
            x=np.array(range(int(time_yr.min()),int(time_yr.max())+2))
            ax1.plot(x,func(fit_ra[0]*10**6*365.25,b_ra,x),color="red",linewidth=2)
            ax1.plot([time_yr[-1]],[0],'o',markersize=15, mfc='none',mec='green',mew=2)
            ax1.errorbar([time_yr[-1]],[0],yerr=[e_ra[-1]*10**6],fmt='',ecolor="green",elinewidth=2,capsize=10,capthick=2)
            ax1.xaxis.set_major_locator(ticker.MaxNLocator(integer=True))
            xticks(fontsize='medium')
            yticks(fontsize='medium')
            ax1.set_title(ivs_i,fontsize=20)
            ax1.set_ylabel('$\Delta$ RA ($\mu$as)',fontsize=20)
            plt.setp(ax1.get_xticklabels(),visible=False)
            
            
            #locator_params(axis='x',nbins=4)
            ax2.errorbar(time_yr[0:-1],(dec[0:-1]-dec[-1])*10**6,yerr=e_dec[0:-1]*10**6,fmt='.',mfc='blue',ecolor="black",elinewidth=2,ms=15)
            x=np.array(range(int(time_yr.min()),int(time_yr.max())+2))
            ax2.plot(x,func(fit_de[0]*10**6*365.25,b_de,x),color="red",linewidth=2)
            ax2.plot([time_yr[-1]],[0],'o',markersize=15,mfc='none',mec='green',mew=2)
            ax2.errorbar([time_yr[-1]],[0],yerr=[e_dec[-1]*10**6],fmt='',ecolor="green",elinewidth=3,capsize=10,capthick=2)
            xlim([int(time_yr.min())-1,int(time_yr.max())+2])
            xticks(fontsize='medium')
            yticks(fontsize='medium')
           
            ax2.xaxis.set_major_locator(ticker.MaxNLocator(integer=True))
            ax2.set_xlabel('Year',fontsize=20)
            ax2.set_ylabel('$\Delta$ Dec ($\mu$as)',fontsize=20)
            ax2.set_yticks(ax2.get_yticks()[:-1])
            fig.subplots_adjust(hspace=0)
            savefig('papergraphs/'+ivs_i+'_timeseries.eps')
            #show()
            
            #calc chisq
            model_ra=func(fit_ra[0],fit_ra[1],time)
            model_de=func(fit_de[0],fit_de[1],time)
            cs_ra=comp_chisq(ra,model_ra,e_ra)/(len(ra)-2) #reduced chisq
            cs_de=comp_chisq(dec,model_de,e_dec)/(len(dec)-2)
            #calc number of standard deviations chisq is from expected (1)
            if len(ra)>2:
                pte_ra=abs(cs_ra-1)/(2./(float(len(ra))-2))
            else:
                pte_ra=999
            if len(dec)>2:
                pte_de=abs(cs_de-1)/(2./(float(len(dec))-2))
            else:
                pte_de=999

            dec_temp=(dec[dec.size-1]/3600)*pi/180
            #last time index
            tindex=np.where(np.array(time)==max(time))
            tindex=tindex[0][0]
            
            writer.writerow([ivs_i,ra[tindex]/3600,dec[tindex]/3600,e_ra[tindex],e_dec[tindex],fit_ra[0]*10**6*365.25*cos(dec_temp),sqrt(e_fit_ra[0][0])*10**6*365.25*cos(dec_temp),cs_ra,pte_ra,ra.size,fit_de[0]*10**6*365.25,sqrt(e_fit_de[0][0])*10**6*365.25,cs_de,pte_de,dec.size,fit_ra[1]/3600,fit_de[1]/3600,(max(time)-min(time))/365.25,max(time),cov_ra,cov_de,'Truebenbach&Darling2017'])

            
        elif len(ra)>0 and n>=10  and fabs(max(time)-min(time))>=3652.5: #Only include quasars w/ >=10 sessions over a span of more than 10 years
            #do bootstrap to find mean proper motion and estimate errors
            its=500
            fit_ra_slope=[]
            fit_ra_int=[]
            fit_de_slope=[]
            fit_de_int=[]
            e_fit_ra_slope=[]
            e_fit_ra_int=[]
            e_fit_ra_cov=[]
            e_fit_de_slope=[]
            e_fit_de_int=[]
            e_fit_de_cov=[]
            for iter in range(its):
                selrandom=random.random_integers(0,len(ra)-1,len(ra))
                selrandom_de=random.random_integers(0,len(dec)-1,len(dec))
                #fit pm   
                fit_ra=fit(e_ra[selrandom],time[selrandom],ra[selrandom])
                e_fit_ra=e_fit(e_ra[selrandom],time[selrandom],ra[selrandom])
                fit_de=fit(e_dec[selrandom_de],time[selrandom_de],dec[selrandom_de])
                e_fit_de=e_fit(e_dec[selrandom_de],time[selrandom_de],dec[selrandom_de])
                fit_ra_slope.append(fit_ra[0])
                fit_ra_int.append(fit_ra[1])
                fit_de_slope.append(fit_de[0])
                fit_de_int.append(fit_de[1])
                e_fit_ra_slope.append(e_fit_ra[0][0])
                e_fit_ra_int.append(e_fit_ra[1][1])
                e_fit_ra_cov.append(e_fit_ra[1][0])
                e_fit_de_slope.append(e_fit_de[0][0])
                e_fit_de_int.append(e_fit_de[1][1])
                e_fit_de_cov.append(e_fit_de[1][0])

            
    
            #pick average PM and y-intercept
            avgslope_ra=median(fit_ra_slope)
            yint_ra=median(fit_ra_int)
            avgslope_de=median(fit_de_slope)
            yint_de=median(fit_de_int)

            #estimate covariance
            e_avgslope_ra=sqrt(cov(fit_ra_slope,fit_ra_slope)) #arcsec/day
            cov_ra=cov(fit_ra_slope,fit_ra_int)
            e_yint_ra=sqrt(cov(fit_ra_int,fit_ra_int))
            cov_ra=[[e_avgslope_ra**2,cov_ra],[cov_ra,e_yint_ra**2]] #full covariance matrix
            e_avgslope_de=sqrt(cov(fit_de_slope,fit_de_slope)) #arcsec/day
            cov_de=cov(fit_de_slope,fit_de_int)
            e_yint_de=sqrt(cov(fit_de_int,fit_de_int))
            cov_de=[[e_avgslope_de**2,cov_de],[cov_de,e_yint_de**2]]
            
            dec_temp=(dec[dec.size-1]/3600)*pi/180
            
            
            #calc chisq
            model_ra=func(avgslope_ra,yint_ra,time)
            model_de=func(avgslope_de,yint_de,time)
            cs_ra=comp_chisq(ra,model_ra,e_ra)/(len(ra)-2) #reduced chisq
            cs_de=comp_chisq(dec,model_de,e_dec)/(len(dec)-2)
            #calc number of standard deviations chisq is from expected (1)
            if len(ra)>2:
                pte_ra=abs(cs_ra-1)/(2./(float(len(ra))-2))
            else:
                pte_ra=999
            if len(dec)>2:
                pte_de=abs(cs_de-1)/(2./(float(len(dec))-2))
            else:
                pte_de=999
            
            dec_temp=(dec[dec.size-1]/3600)*pi/180
            
            tindex=np.where(np.array(time)==max(time))
            tindex=tindex[0][0]
            
            writer.writerow([ivs_i,ra[tindex]/3600,dec[tindex]/3600,e_ra[tindex],e_dec[tindex],avgslope_ra*10**6*365.25*cos(dec_temp),e_avgslope_ra*10**6*365.25*cos(dec_temp),cs_ra,pte_ra,n,avgslope_de*10**6*365.25,e_avgslope_de*10**6*365.25,cs_de,pte_de,n,yint_ra/3600,yint_de/3600,(max(time)-min(time))/365.25,max(time),e_fit_ra,e_fit_de,'GSFC'])

        elif len(ra)>0 and ivs_i not in new_id: #add to reject list
            tindex=np.where(np.array(time)==max(time))
            tindex=tindex[0][0]
            
            w2.writerow([ivs_i,ra[tindex]/3600,dec[tindex]/3600,n,max(time),fabs(max(time)-min(time))])
        
            
        #Get ready for next object
        ivs_i=line[0]
        if float(line[11])<47892.5 or float(line[12])<3:
            time=[]
            ra=[]
            dec=[]
            e_ra=[]
            e_dec=[]
            corr=[]
            n=0
        else:
            time=[float(line[11])] #MJD
            r=[float(line[2]),float(line[3]),float(line[4])]
            d=[float(line[5]),float(line[6]),float(line[7])]
            r,d=convert_coord(r,d)
            ra=[r*3600] #arcsec
            dec=[d*3600]
            e_ra=[float(line[8])*15] #arcsec
            e_dec=[float(line[9])]
            corr=[float(line[10])] #unitless
            n=1 #number of sessions

f.close()
g.close()

print total_objects

#Make a version of the PM catalog that has redshifts
os.system('java -jar $HOME/bin/stilts/stilts.jar tmatch2 \
            in1=2017a_PMs.csv ifmt1=csv in2=../pairs/redshift_list.csv \
            ifmt2=csv ofmt=csv out=2017a_PMs_withz.csv \
            matcher=sky \
            values1="RA DE" values2="RA_deg DE_deg" \
            params=1.0 join=all1 find=best1 \
            ocmd="delcols Name" ocmd="delcols RA_deg" \
            ocmd="delcols DE_deg"  \
            ocmd="delcols Separation"')
            

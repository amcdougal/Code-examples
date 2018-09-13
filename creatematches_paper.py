"""
This program creates a large catalog by cross-matching and filtering several online catalogs. Uses STILTS, a catalog comparison tool. Used to create catalog published in Truebenbach & Darling, 2017, MNRAS, 468, 196.

Compares a local portion of the FIRST catalog to the online SDSS DR9 catalog and the online AllWISE catalog.  The program ultimately produces a list of FIRST sources with WISE counterparts but without SDSS DR9 optical counterparts.  It removes sources that are outside of the SDSS DR9 footprint, are near bright stars (r mag < 10), and sources that aren't detected in several WISE bands with SNR>5 or in one WISE band with SNR>7.  WISE and FIRST sources must be 6arcsec apart to be considered for a match, and a WISE-FIRST source must not have a SDSS source within the largest positional uncertainty of WISE/FIRST to be considered 'without optical counterpart.' The largest positional uncertainty is also used to determine the initial WISE-FIRST match.  The AllWISE coordinates are used in all catalog comparisons if possible because these are the most accurate.  The program also calculates the probability of association between the WISE and FIRST sources using their positional uncertainties. The positional uncertainties are the 1sigma uncertainties for WISE and the 90% uncertainties for FIRST.

Steps to run:
1)Download FIRST catalog from http://sundog.stsci.edu/cgi-bin/searchfirst
   Note: Entire FIRST catalog available for download elsewhere.
2)Edit FIRST catalog to be ascii readable.  There must be only one header row, with no space btw first column and #.  One header text block per column.  Label RA and Dec columns as RAh RAm RAs Decd Decm Decs.
3)Modify file names immediately below.
4)Run this program >>python creatematches_variable.py


By: Alex Truebenbach
Last Editted: Feb 2015
"""

import csv,os,urllib
from pylab import *

#Make the program nice
os.nice(30)

#File names, modify these for a new matching.

direc='allsky_PAPERversion'
first_data='catalog_14mar04.bin'
brightstar_catalog='brightstars_sdssdr9.csv' #Catalog of bright stars that we'll steer clear of.
brightstar_catalog_new='brightstars_wradii.csv'
matchwise_combined='matchir_allskyV.tbl'
matchwise_conserv='matchir_allskyV_conserv.tbl'
footprint_check='matchir_allskyV_sdsscheck.tbl'
nosdss_dr9='noopDR9_allskyV.tbl'
nosdss_atall='noopt_allskyV.tbl'
final_list_almost='final_match_allskyV_almost.tbl'
final_list='final_match_allskyV.tbl'
final_list_2mass='final_match_allskyV.no2mass.tbl'
final_real='final_match_allskyV.no2mass_FINAL.tbl'
final_short='final_match_allskyV.no2mass_FINAL.short.tbl'
######

###Functions###
def convert_coords(file,output):
    #Converts RA and Dec from sexigesimal to degrees
    reader=csv.reader(open(file,'rb'))
    data=[]
    data.extend(reader)
    header=data[0]
    header=header[0]
    header=header.split()
    rah_loc=header.index('#RAh')
    ram_loc=header.index('RAm')
    ras_loc=header.index('RAs')
    decd_loc=header.index('Decd')
    decm_loc=header.index('Decm')
    decs_loc=header.index('Decs')
    data.pop(0)
    f=open(output,'wb')
    writer=csv.writer(f,delimiter=',')
    writer.writerow(header+['RA_deg','Dec_deg'])
    for row in data:
        line=row[0]
        line=line.split()
        rah=float(line[rah_loc])
        ram=float(line[ram_loc])
        ras=float(line[ras_loc])
        decm=float(line[decm_loc])
        decd=line[decd_loc]
        decs=float(line[decs_loc])
        ra=(rah+(ram+(ras/60))/60)*15
        x=str(format((decm+(decs/60))/60,'.10f'))
        x=x.split('.')
        dec=float(decd+'.'+x[1])
        writer.writerow(line+[ra,dec])
    f.close()
    
    
def trim_first():
    #Removes sections of FIRST that aren't covered by SLOAN
    sdss_cover=csv.reader(open('allrunsdr7db.par','rb'))
    data=[]
    data.extend(sdss_cover)
    data=data[36:]
    mu_dict={}
    

    dec=[]
    ra=[]
    for row in data:
        line=row[0]
        line=line.split()
        mu_up=float(line[6])
        mu_lower=float(line[5])
        stripe=float(line[3])
        #Create dictionary of mu upper and lower ranges for the main group of stripes
        if stripe > 9 and stripe < 38:
            if stripe not in mu_dict.keys():
                mu_dict[stripe]=[mu_lower,mu_up]
            elif mu_lower < mu_dict[stripe][0]:
                mu_dict[stripe][0]=mu_lower
            elif mu_up > mu_dict[stripe][1]:
                mu_dict[stripe][1]=mu_up
                
        #Add in RA and DECs of the partial stripes around the main group    
        elif stripe == 9 or (stripe > 37 and stripe <= 44):
            mu_range=arange(mu_lower+1.25,mu_up-1.25,1.25)
            i=((float(line[3])-10)*2.5)*(pi/180)            
            if line[4] == 'N':
                nu=0.625*(pi/180)
            else:
                nu=-0.625*(pi/180)
            for j in mu_range:
                dec_new=math.asin(sin(nu)*cos(i)+sin((j-95)*(pi/180))*cos(nu)*sin(i))*(180/pi)
                ra_new=math.acos((cos((j-95)*(pi/180))*cos(nu))/cos(dec_new*(pi/180)))*(180/pi)+95
                if dec_new not in dec:
                    dec.append(dec_new)
                    ra.append(ra_new)

    #Add in RA and Decs of stripes in dictionary.  Assume all points between the highest and lowest mu's for that stripe are covered by sloan.  Implied by the image on the website but not actually shown in the data file.
    for stripe in mu_dict.keys():
        mu_lower=mu_dict[stripe][0]
        mu_up=mu_dict[stripe][1]
        mu_range=arange(mu_lower+1.25,mu_up-1.25,1.25)
        i=((stripe-10)*2.5)*(pi/180)            
        nu1=0.625*(pi/180)
        nu2=-0.625*(pi/180)
        for j in mu_range:
            dec_new1=math.asin(sin(nu1)*cos(i)+sin((j-95)*(pi/180))*cos(nu1)*sin(i))*(180/pi)
            dec_new2=math.asin(sin(nu2)*cos(i)+sin((j-95)*(pi/180))*cos(nu2)*sin(i))*(180/pi)
            ra_new1=math.acos((cos((j-95)*(pi/180))*cos(nu1))/cos(dec_new1*(pi/180)))*(180/pi)+95
            ra_new2=math.acos((cos((j-95)*(pi/180))*cos(nu2))/cos(dec_new2*(pi/180)))*(180/pi)+95
            if dec_new1 not in dec:
                dec.append(dec_new1)
                ra.append(ra_new1)
            if dec_new2 not in dec:
                dec.append(dec_new2)
                ra.append(ra_new2)

    #Output as a list of RA and Dec
    f=open('sdss_coverage.csv','w')
    writer=csv.writer(f,delimiter=',')
    writer.writerow(['RA','Dec'])
    for i,d in enumerate(dec):
        writer.writerow([ra[i],dec[i]])
    f.close()

def split_first(input):
    #Splits the FIRST catalog into 10 subcatalogs that are each matched with WISE individually.  The main purpose of this is to be able to recover part of the matching process if there is an error during the stilts matching.
    #Returns a list of FIRST files that the catalog has been divided up into.
    output_root=input.split('.')
    output_root=output_root[0]

    reader=csv.reader(open(input,'rb'))
    data=[]
    data.extend(reader)
    header=data[0]
    data.pop(0)
    
    i=1 #File number counter
    modnum=int(len(data)/10) #Number of sources in each file

    filelist=[]
    
    f=open(output_root+'_0.tbl','wb')
    writer=csv.writer(f,delimiter=',')
    writer.writerow(header)
    filelist.append(output_root+'_0.tbl')
    for index,row in enumerate(data):
        writer.writerow(row)
        if fmod(index,modnum)==0 and index!=0 and index<len(data)-5:
            f.close()
            f=open(output_root+'_'+str(i)+'.tbl','wb')
            writer=csv.writer(f,delimiter=',')
            writer.writerow(header)
            filelist.append(output_root+'_'+str(i)+'.tbl')
            i+=1
            
    f.close()
    return filelist
    
def w4detect(file,output):
    #Select those WISE sources with ONLY W4 detections   

    reader = csv.reader(open(file,'rb'))   
    data = []
    data.extend(reader)
    header=data[0]
    w1sn_loc=header.index('snr1')
    w2sn_loc=header.index('snr2')
    w3sn_loc=header.index('snr3')
    w4sn_loc=header.index('snr4')
    data.pop(0)

    f = open(output,'w')
    writer=csv.writer(f,delimiter=',')
    writer.writerow(header)
    for index,row in enumerate(data):
        if row[w1sn_loc]=='':
            w1sn=0
        else:
            w1sn=float(row[w1sn_loc])
        if row[w2sn_loc]=='':
            w2sn=0
        else:
            w2sn=float(row[w2sn_loc])
        if row[w3sn_loc]=='':
            w3sn=0
        else:
            w3sn=float(row[w3sn_loc])
        if row[w4sn_loc]=='':
            w4sn=0
        else:
            w4sn=float(row[w4sn_loc])
        if w1sn<=5 and w2sn<=5 and w3sn<=5 and w4sn>=7:
            writer.writerow(row)

def rmdup(file,output):
    """
    This function removes duplicates under the assumption that if two objects have the same GroupID, then the one with the smallest Separation is the best object to keep.  
    """
    reader=csv.reader(open(file,'rb'))
    data=[]
    data.extend(reader)
    header=data[0]
    sep=header.index('Separation')
    id=header.index('GroupID')
    dup_dict={}
    for index,row in enumerate(data):
        if index != 0 and row[id]!='':
            id_val=float(row[id])
            sep_val=float(row[sep])
            if id_val not in dup_dict.keys():
                dup_dict[id_val]=[sep_val]
            else:
                if sep_val not in dup_dict[id_val]:
                    dup_dict[id_val].append(sep_val)
                    
    f=open(output,'w')
    writer=csv.writer(f,delimiter=',')
    for index,row in enumerate(data):
        if index !=0 and row[id]!='':
            id_val=float(row[id])
            sep_val=float(row[sep])
            if sep_val == min(dup_dict[id_val]):
                writer.writerow(row)
        else:
            writer.writerow(row)


def create_sdssmatches(file,nomatches,stars):
    #From a list of all matches with SDSS DR7, including nonmatches, create a) a list of FIRST sources with no optical counterparts (these will have no Separation parameter listed) and b) a list of FIRSt sources whose optical counterpart is listed as a star. Only include stars that are >0.4" away from the FIRST source; these are the most likely to be misidentifications.
    reader=csv.reader(open(file,'rb'))
    data=[]
    data.extend(reader)
    header=data[0]
    separation_loc=header.index('Separation')
    cl_loc=header.index('cl_cone')
    data.pop(0)

    f=open(nomatches,'w')
    writer=csv.writer(f,delimiter=',')
    g=open(stars,'w')
    writer_star=csv.writer(g,delimiter=',')
    writer.writerow(header)
    writer_star.writerow(header)
    for row in data:
        if row[separation_loc]=='':
            writer.writerow(row)
        if row[cl_loc] == '6' and float(row[separation_loc])>=1.11e-4:
            writer_star.writerow(row)
    f.close()
    g.close()
    
def rm_2mass(file,output):
    #Remove those wise sources that have 2mass detections.  These represent galaxies with optical counterparts that slipped through the filtering process.
    reader=csv.reader(open(file,'rb'))
    data=[]
    data.extend(reader)
    header=data[0]
    j_loc=header.index('Jmag')
    h_loc=header.index('Hmag')
    k_loc=header.index('Kmag')
    jsigma_loc=header.index('e_Jmag')
    hsigma_loc=header.index('e_Hmag')
    ksigma_loc=header.index('e_Kmag')
    data.pop(0)
    
    #Create table of sources with 2MASS detections
    writer=csv.writer(open(direc+'/tmp_2mass.tbl','wb'),delimiter=',')
    writer.writerow(header)
    for row in data:
        j=row[j_loc]
        h=row[h_loc]
        k=row[k_loc]
        if (j!='' and row[jsigma_loc]!='') or (h!='' and row[hsigma_loc]!='') or (k!='' and row[ksigma_loc]!=''):
            writer.writerow(row)

def snr_limit(input,output):
    """
    This program only keeps those sources that are either detected in multiple bands at SNR > 5.0 or in one band at SNR > 7.0.
    """
    reader=csv.reader(open(input,'rb'),delimiter=',')
    data=[]
    data.extend(reader)
    header=data[0]
    w1sn_loc=header.index('snr1')
    w2sn_loc=header.index('snr2')
    w3sn_loc=header.index('snr3')
    w4sn_loc=header.index('snr4')
    data.pop(0)
    
    file=open(output,'w')
    writer=csv.writer(file,delimiter=',')
    writer.writerow(header)
    
    for row in data:
        if row[w1sn_loc]!='':
            w1_sn=float(row[w1sn_loc])
        else:
            w1_sn=0
        if row[w2sn_loc]!='':
            w2_sn=float(row[w2sn_loc])
        else:
            w2_sn=0
        if row[w3sn_loc]!='':
            w3_sn=float(row[w3sn_loc])
        else:
            w3_sn=0
        if row[w4sn_loc]!='':
            w4_sn=float(row[w4sn_loc])
        else:
            w4_sn=0
        if w1_sn>7.0 or w2_sn>7.0 or w3_sn>7.0 or w4_sn>7.0 or (w1_sn>5.0 and w2_sn>5.0) or (w1_sn>5.0 and w3_sn>5.0) or (w1_sn>5.0 and w4_sn>5.0) or (w2_sn>5.0 and w3_sn>5.0) or (w2_sn>5.0 and w4_sn>5.0)or (w3_sn>5.0 and w4_sn>5.0):
            writer.writerow(row)
    file.close()

def sep_uncertainty(input):
    #Calculates the uncertainty in the separation parameter based on the WISE and FIRST positional uncertainties.  See notes on 7/10/13 for derivation.  
    #Also outputs the positional uncertainities of the WISE and FIRST coordinates.

    reader=csv.reader(open(input,'rb'))
    data=[]
    data.extend(reader)
    header=data[0]
    wra_loc=header.index('RAJ2000')
    wdec_loc=header.index('DEJ2000')
    wmaj_loc=header.index('eeMaj')
    wmin_loc=header.index('eeMin')
    wpa_loc=header.index('eePA')
    fra_loc=header.index('RA_deg')
    fdec_loc=header.index('Dec_deg')
    fmaj_loc=header.index('fMaj')
    fmin_loc=header.index('fMin')
    fpa_loc=header.index('fPA')
    fpeak_loc=header.index('Fpeak')
    rms_loc=header.index('RMS')
    data.pop(0)

    f=open(input,'wb')
    writer=csv.writer(f,delimiter=',')
    writer.writerow(header+['sigma_sep','WISE_RAsig','WISE_Decsig','FIRST_RAsig','FIRST_Decsig'])
    for row in data:
        wra=float(row[wra_loc])*(pi/180) #positions all in degrees, convert to radians for below
        wdec=float(row[wdec_loc])*(pi/180)
        wmaj=float(row[wmaj_loc]) #arcsec
        wpa=float(row[wpa_loc])*(pi/180)
        fra=float(row[fra_loc])*(pi/180)
        fdec=float(row[fdec_loc])*(pi/180)
        fmaj=float(row[fmaj_loc]) #arcsec
        fmin=float(row[fmin_loc])
        fpa=float(row[fpa_loc])
        fpeak=float(row[fpeak_loc])
        rms=float(row[rms_loc])

        #Calculate WISE ra and dec 1sigma uncertainites and covariance
        wsig_ra=wmaj**2*math.sin(wpa)**2+wmin**2*math.cos(wpa)**2 #arcsec^2
        wsig_dec=wmaj**2*math.cos(wpa)**2+wmin**2*math.sin(wpa)**2 #arcsec^2

    	#Calculate the derivatives of r wrt the coordinates.
    	x=sin(wdec)*sin(fdec)+cos(fdec)*cos(wdec)*cos(wra-fra)
    	d_wdec=-1*(cos(wdec)*sin(fdec)-sin(wdec)*cos(fdec)*cos(wra-fra))/(math.sqrt(1-x**2))
    	d_fdec=-1*(sin(wdec)*cos(fdec)-cos(wdec)*sin(fdec)*cos(wra-fra))/(math.sqrt(1-x**2))
    	d_wra=(cos(wdec)*cos(fdec)*sin(wra-fra))/math.sqrt(1-x**2)
    	d_fra=-1*(cos(wdec)*cos(fdec)*sin(wra-fra))/math.sqrt(1-x**2)

        #Calculate wise covariance	
        if wpa==0 or wpa==90 or wpa==180:
            #RA and Dec aren't correlated.
            wsigradec=0
        else:
           wsigradec=(wmaj**2-wsig_ra)*math.tan(wpa) #arcsec^2

        #Calculate FIRST ra and dec 90% confidence uncertainites and covariance.
        snr=(fpeak-.25)/rms
        err_fmaj=fmaj*(1./snr+1.0/20) #arcsec
        err_fmin=fmin*(1./snr+1.0/20) #arcsec
        fsig_ra=(err_fmaj**2*math.sin(fpa*(pi/180))**2+err_fmin**2*math.cos(fpa*(pi/180))**2) #arcsec^2
        fsig_dec=(err_fmaj**2*math.cos(fpa*(pi/180))**2+err_fmin**2*math.sin(fpa*(pi/180))**2) #arcsec^2

        #Calculate FIRST covariance
        if fpa==0 or fpa==90 or fpa==180:
            #RA and dec aren't correlated
            fsigradec=0
        else:
            fsigradec=(err_fmaj**2-fsig_ra)*math.tan(fpa*(pi/180)) #arcsec^2


        #Calculate sigma_r
        sigr=math.sqrt(d_wdec**2*wsig_dec+2*d_wdec*d_wra*wsigradec+d_wra**2*wsig_ra+d_fdec**2*fsig_dec+d_fra**2*fsig_ra+2*d_fra*d_fdec*fsigradec) #arcsec
        #Print out row
        writer.writerow(row+[sigr,sqrt(wsig_ra),sqrt(wsig_dec),sqrt(fsig_ra),sqrt(fsig_dec)])
    f.close()

        
def filter_conservative(file,matchwise_conserv):
    #Filter the FIRST-WISE matches to only include matches with a match radius that's less than the greates    #IMPORTANT: The separation uncertainty calculated here is wrong.  See next function for correct calculation.    #IMPORTANT: The separation uncertainty calculated here is wrong.  See next function for correct calculation.t positional uncertainity.  So, if r=5" and the FIRST positional uncertainty (in either RA or DEC) is 4", reject it!  Similarly, if r=5" and the WISE positional uncertaintiy is 3", reject it!
    #Filtered output is matchwise_conserv
     #Also Choose which search radius to use for matching the FIRST-WISE catalog with SDSS.  If FIRST has a larger positional error, use that error as the search radius and use the FIRST coordinates to compare with SDSS.  Otherwise, if WISE has a larger positional error, use that error as the search radius and the WISE coordinates to compare with SDSS.
    #Outputs same file, but with 3 additional columns: search_ra, search_dec, and search_r.  The first two are the object coordinates used to compare with SDSS in degrees.  The last is the search radius to be used, in degrees.
    reader=csv.reader(open(file,'rb'))
    data=[]
    data.extend(reader)
    header=data[0]
    wra_loc=header.index('RAJ2000')
    wdec_loc=header.index('DEJ2000')
    fra_loc=header.index('RA_deg')
    fdec_loc=header.index('Dec_deg')
    w_raerr=header.index('WISE_RAsig')
    w_decerr=header.index('WISE_Decsig')
    f_raerr=header.index('FIRST_RAsig')
    f_decerr=header.index('FIRST_Decsig')
    r_loc=header.index('Separation')
    data.pop(0)

    f=open(matchwise_conserv,'wb')
    writer=csv.writer(f,delimiter=',')
    writer.writerow(header+['search_ra','search_dec','search_r'])
    for row in data:
        wra=float(row[wra_loc])
        wdec=float(row[wdec_loc])
        fra=float(row[fra_loc])
        fdec=float(row[fdec_loc])
        first_raerr=float(row[f_raerr])
        first_decerr=float(row[f_decerr])
        wise_raerr=float(row[w_raerr])
        wise_decerr=float(row[w_decerr])
        r=float(row[r_loc])*3600
        errors=[first_raerr,first_decerr,wise_raerr,wise_decerr]
        maxerr=max(errors)
        errtype=errors.index(maxerr)
        #Convert maxerr from arcsec to degrees
        
        if r<maxerr:
            maxerr=maxerr/3600
            if errtype==0 or errtype==1:
                writer.writerow(row+[fra,fdec,maxerr])
            else:
                writer.writerow(row+[wra,wdec,maxerr])
    f.close()
    
    
    
        
def getsdssimages(file):
    #From a list of FIRST sources whose optical counterparts are stars, download SDSS DR7 images of these sources centered on their FIRST coordinates.  Images will be placed in images_sdss.
    
    ra = []
    dec = []
    reader = csv.reader(open(file,'rb'))
    data=[]
    data.extend(reader)
    header=data[0]
    ra_loc=header.index('RA_deg')
    dec_loc=header.index('Dec_deg')
    data.pop(0)

    #Create list of FIRST coordinates of the matches
    for row in data:
        ra.append(row[ra_loc])
        dec.append(row[dec_loc])


    #Examine images_sdss folder for existing files.  Create a list of existing images.  Only download images that aren't already in images_sdss.
    curr_images=os.listdir('images_sdss')

    for i,r in enumerate(ra):
        filename=str(r)+'_'+str(dec[i])+'.jpeg'
        if filename not in curr_images:
            f = urllib.urlretrieve('http://casjobs.sdss.org/ImgCutoutDR7/getjpeg.aspx?ra='+str(r)+'&dec='+str(dec[i])+'&scale=0.09903&width=400&height=400&opt=GLP','images_sdss/'+filename) 
    

###Main###
first_data_mod=first_data.split('.bin')
first_data_mod=first_data_mod[0]+'_mod.tbl'

convert_coords(first_data,first_data_mod)
#Remove any FIRST sources that are in regions of the sky not covered by SDSS DR9.  The region is not covered by SDSS if there are no SDSS sources within 1'.

curr_files=os.listdir('.')
if 'sdss_coverage.tbl' not in curr_files:
    print 'Error, no coverage table.  Run coverage_map_paper.py'
os.system('java -jar $HOME/bin/stilts/stilts.jar tmatch2 \
     ifmt1=csv ifmt2=csv in1=sdss_coverage.tbl in2='+first_data_mod+' \
     ofmt=csv out='+first_data_mod+' \
     ocmd="delcols RA" ocmd="delcols Dec" ocmd="delcols GroupID" \
     ocmd="delcols GroupSize" ocmd="delcols Separation" \
     matcher=sky values1="RA Dec" values2="RA_deg Dec_deg" params=600 \
     join=1and2 find=best2')

#Break the first list into 10 sublists.
filelist=split_first(first_data_mod)

#Create list of WISE sources with FIRST counterparts. r=12 arcsec.
#Do this for each sublist.
for index,file in enumerate(filelist):
    index=index
    os.system('java -jar $HOME/bin/stilts/stilts.jar coneskymatch \
          icmd=progress ofmt=fits parallel=5 \
          ifmt=csv find=all ra=RA_deg dec=Dec_deg sr=0.003333 \
          erract=retry10 in='+file+' \
          ocmd="delcols "*_1"" \
          serviceurl="http://vizier.u-strasbg.fr/viz-bin/votable/-A?-source=II/328&-out.all&" \
          out='+direc+'/firsttmp_'+str(index)+'.tbl')

#ocmd="delcols "ccf"" \
# ocmd="delcols "wy"" ocmd="delcols "Jmag"" ocmd="delcols "e_Jmag"" \
 #         ocmd="delcols "Hmag"" ocmd="delcols "e_Hmag"" ocmd="delcols "Kmag"" \
  #        ocmd="delcols "e_Kmag"" ocmd="delcols "2Mkey"" ocmd="delcols "d2M""
#Combine the lists
if len(filelist)!=10:
    print 'Warning, combining the lists wont work / may not join all of the lists!'
   
os.system('java -jar $HOME/bin/stilts/stilts.jar tcatn \
          nin=10 in1='+direc+'/firsttmp_0.tbl in2='+direc+'/firsttmp_1.tbl \
          in3='+direc+'/firsttmp_2.tbl in4='+direc+'/firsttmp_3.tbl \
          in5='+direc+'/firsttmp_4.tbl in6='+direc+'/firsttmp_5.tbl \
          in7='+direc+'/firsttmp_6.tbl in8='+direc+'/firsttmp_7.tbl \
          in9='+direc+'/firsttmp_8.tbl in10='+direc+'/firsttmp_9.tbl \
          out='+direc+'/tmp_combined.tbl \
          ifmt1=fits ifmt2=fits ifmt3=fits ifmt4=fits ifmt5=fits ifmt6=fits \
          ifmt7=fits ifmt8=fits ifmt9=fits ifmt10=fits ofmt=csv' )

#Remove duplicates
os.system('java -jar $HOME/bin/stilts/stilts.jar tmatch1 \
           action=keep1 matcher=exact values="AllWISE" \
           ifmt=csv in='+direc+'/tmp_combined.tbl \
           out='+matchwise_combined+' ofmt=csv')

#Calculate the uncertainty in the separation between FIRST and WISE
sep_uncertainty(matchwise_combined)

#Filter FIRST+WISE catalog to create a conservative list of matches.  Only keeps those with r < largest positional uncertainty.  Also put in a search radius for SDSS vs. FIRST+WISE comparison.  Use largest positional uncertainty for the radius.  Use corresponding FIRST or WISE coordinates for matching.  (columns are search_ra,search_dec, and search_r).
filter_conservative(matchwise_combined,matchwise_conserv)

#Only keep sources that are detected in multiple bands at SNR>5 or are detected in one band at SNR>7 and all other bands at SNR<5.
snr_limit(matchwise_conserv,matchwise_conserv)

#First, do another check that all of these objects are in SDSS footprint.  Remove objects that dont have an SDSS object within 1'. Also convert Separation_wf and its error to arcsec.  
os.system('java -jar $HOME/bin/stilts/stilts.jar coneskymatch \
     ifmt=csv in='+matchwise_conserv+' \
     icmd=progress icmd="colmeta -name "Separation_wf" "Separation"" \
     icmd="replacecol Separation_wf Separation_wf*3600" \
     icmd="replacecol sigma_sep sigma_sep*3600" \
     out='+footprint_check+' ofmt=csv \
     ra=search_ra dec=search_dec sr=0.01667 \
     fixcols=all suffix0="" suffix1="_cone" \
     find=each servicetype=cone parallel=5 \
     ocmd="delcols "*_cone*"" \
     ocmd="delcols Separation" erract=retry10 \
     serviceurl="http://vizier.u-strasbg.fr/viz-bin/votable/-A?-source=V/139&"')


#Compare FIRST+WISE catalog to online SDSS DR9 catalog (vizier).  Search radius is 2xlargest FIRST/WISE positional uncertainty.

os.system('java -jar $HOME/bin/stilts/stilts.jar coneskymatch \
     ifmt=csv in='+footprint_check+' \
     icmd=progress \
     out='+nosdss_dr9+' ofmt=csv \
     icmd="addcol sdss_searchr search_r*2" \
     ra=search_ra dec=search_dec sr=sdss_searchr \
     fixcols=all suffix0="" suffix1="_cone" \
     find=each servicetype=cone parallel=5 \
     ocmd="delcols "*_cone*"" ocmd="select NULL_Separation"\
     ocmd="delcols Separation" erract=retry10 \
     serviceurl="http://vizier.u-strasbg.fr/viz-bin/votable/-A?-source=V/139&"')
#Redo with DR7 catalog
os.system('java -jar $HOME/bin/stilts/stilts.jar coneskymatch \
     ifmt=csv in='+nosdss_dr9+' \
     icmd=progress  \
     out='+nosdss_atall+' ofmt=csv \
     ra=search_ra dec=search_dec sr=sdss_searchr \
     fixcols=all suffix0="" suffix1="_cone" \
     find=each servicetype=cone parallel=5 \
     ocmd="delcols "*_cone*"" ocmd="select NULL_Separation"\
     ocmd="delcols Separation" erract=retry10 \
     serviceurl="http://vizier.u-strasbg.fr/viz-bin/votable/-A?-source=II/294&"')

###Remove all sources that are near a bright star.
radius=[40,60,120] #radius of circle to be removed for r mags between 10-9, 9-8, <8, respectively
#Some stars have mag = -9999, an error. These stars are almost always very very faint or w radius<5arcsec.  Dont include these in the bright star catalog since it's unlikely they'll align with a WISE-FIRST object by chance (bc of the small radius)

#Read in bright star catalog.  Add in a new column for search radius
reader=csv.reader(open(brightstar_catalog,'rb'))
data=[]
data.extend(reader)
header=data[0]
data.pop(0)

f=open(brightstar_catalog_new,'wb')
writer=csv.writer(f,delimiter=',')
writer.writerow(header+['radius'])

for row in data:
    if float(row[2])>9:
        writer.writerow(row+[radius[0]])
    elif float(row[2])>8:
        writer.writerow(row+[radius[1]])
    elif float(row[2])<-999:
        a=0 #i.e. do nothing
    else:
        writer.writerow(row+[radius[2]])

f.close()

#Match this new catalog with WISE-FIRST catalog. Remove all WISE-FIRST catalog entries that are within these new radii of a bright star.
#Do this in 2 steps.  First, find sep between catalog and nearest bright star. Use 120 arcsec because thats the largest area that will potentially have to be removed
os.system('java -jar $HOME/bin/stilts/stilts.jar tmatch2 \
     ifmt1=csv ifmt2=csv in1='+brightstar_catalog_new+' in2='+nosdss_atall+' \
     ofmt=csv out='+final_list_almost+' \
     matcher=sky values1="ra dec" values2="search_ra search_dec" params=120 \
     join=all2 find=best2')


#Then remove all objects that have a bright star match whose separation is < the bright star's removal radius
reader=csv.reader(open(final_list_almost,'rb'))
data=[]
data.extend(reader)
header=data[0]
sep_loc=header.index('Separation')
star_r_loc=header.index('radius')
data.pop(0)

f=open(final_list,'wb')
writer=csv.writer(f,delimiter=',')
writer.writerow(header[4:len(header)-3])

for row in data:
    sep=row[sep_loc]
    if sep=='': #No bright star nearby. Can keep object
        writer.writerow(row[4:len(header)-3])
    else:
        star_r=float(row[star_r_loc])
        if float(sep)>star_r: #object is far enough away from the bright star to be kept
            writer.writerow(row[4:len(header)-3])
        

f.close()

#Remove objects that have a match in 2MASS.  Use the sdss search radius (2x objects largest positional uncertainty.  Basically a 3sigma cut since largest uncertainty is almost always from FIRST, which has 90% errors.
os.system('java -jar $HOME/bin/stilts/stilts.jar coneskymatch \
     ifmt=csv in='+final_list+' \
     icmd=progress \
     out='+final_list_2mass+' ofmt=csv \
     ra=search_ra dec=search_dec sr=sdss_searchr \
     fixcols=all suffix0="" suffix1="_cone" \
     find=each servicetype=cone parallel=5 \
     ocmd="delcols "*_cone*"" ocmd="select NULL_Separation"\
     ocmd="delcols Separation" erract=retry10 \
     serviceurl="http://vizier.u-strasbg.fr/viz-bin/votable/-A?-source=II/246&"')

#Remove objects that I know shouldnt be in the catalog.  
os.system('java -jar $HOME/bin/stilts/stilts.jar tmatch2 \
           ifmt1=csv ifmt2=csv in1='+final_list_2mass+' \
           in2=exclude.tbl matcher=sky values1="search_ra search_dec" \
           values2="ra dec" find=best1 join=1not2 ofmt=csv \
           out='+final_real+' params=.1')

reader=csv.reader(open(final_real,'rb'))
data=[]
data.extend(reader)
f=open(final_short,'wb')
writer=csv.writer(f,delimiter=',')

for header in data:
    writer.writerow([header[0],header[1],header[2],header[3],header[4],header[5],header[6],header[7],header[8],header[27],header[28],header[30],header[70],header[31],header[32]]+header[37:45]+[header[54],header[56],header[58],header[60]]+header[106:116])

f.close()

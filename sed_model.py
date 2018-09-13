"""
This program uses parallel processing with python to create a list of best fitting models and redshifts for a set of observations, based on the flux each models produces at each redshift through the same filters as the observations. The models' flux information is obtained from files in the subdirectory 'modelflux', and are created by the program 'calc_modelflux.py'.  The observations information file is user created and has a very specific format (see below).  In this file, the user may specify which filters have definite observations and which (if any) set upper limits on the flux of the object.  The program calculates which models are the best fit by computing a chi squared statistic using all filters in which the object was detected and rejecting those models that have fluxes above the limiting filters.  Program assumes H0=70, omega_m=0.3, omega_l=0.7
Note: This program requires the accompanying program sed_model_functions.py to run.  This file must be in the same directory as the current file.  It also requires that luminosity distances for all redshifts are stored in z.save (created by make_redshifts.py)
Note: In order for the program to run, all redshifts examined by this program must have a corresponding model flux file in the subdirectory 'modelflux', and each of those files must have model information at all of the filters specified in the input file.
Note: If the redshift grid is finer than dz=0.0001, then the lookup for luminosity distances will need to be changed in the code.

This code has been parallelized.  To execute, type >>mpirun -np 8(number of processors to run on) python sed_model.py.

Input File Format:
-----------------
Two lines per object.  First line: object_name flux1 err_flux1 flux2 err_flux2 etc.
Second line: one number for each filter. 0=do not use this filter, 1=this filter has a definite flux measurement.  Use it in calculating chi squared, 2=this filter has an upper limit flux measurement.  Reject models with fluxes above this flux at this filter location.
Notes: -If 1 is set for a filter in the second line, then an error must be given for that flux.  If 0 or 2 is set for the filter, than the error is not used in the program and may be set to 0.
-The fluxes in line one must follow the order u,g,r,i,z,J,H,Ks,W1,W2,W3,W4.  This is the order of ascending wavelengths.  Other filters may be added to the code but it would need to be modified somewhat and an RSR file for that filter obtained.
-It is assumed that all fluxes and their errors are in Jy and that they have been corrected for galactic attenuation.  All errors are 1 sigma errors.

Output File Format:
------------------
Files are output to the subdirectory 'results', with one file per object. All models with chi squared less than 50 are output.
Format: one line for each model: chi squared for that model, degrees of freedom, index of the optical model within its binary file, index of the IR model within its binary file (-99 if no IR model was used), a (the scaling factor used to scale the models to the data. Effectively the stellar mass of the model (in Msolar)), the redshift of the model, log of the model's stellar mass [M*], log of the SSFR averaged over the last 0.1 Gyrs [yr^-1], Ldust [Lsolar/Msolar], and log of the dust mass [M*].

Note:the output model parameters (dust mass, ssfr haven't been scaled by a yet).  Actually dust mass = a*ldust*mdust, actual ssfr=ssfr*a

Description of Program:
----------------------
1)Read in the user created input file with information about the observations.
2)For each object, create a list of fluxes to be used in calculating chi squared and fluxes to be used as upper limits to reject models.
3)At each redshift, read in the model fluxes and other model information in the corresponding redshift file. Also convert the user detected and limiting fluxes from Jy to Lsolar/Hz in the emitted rest frame. 
4)For each optical model, if the model flux != -99 (this implies the model contained stars that were older than the age of the universe at the current redshift), find each IR model that can go with this optical model (ie., the optical and IR models must have Fmu within 0.15 of each other).  This means that the two models have roughly the same amount of dust luminosity from the ISM and stellar birth clouds.
5)For each of the optical - IR model pairs, compute the flux at each filter.  If the filter's central wavelength at the current redshift is <2.5 microns, the filter is all optical so just use the optical flux from the file.  If the filter is at >10 microns, the filter is all IR so use the IR flux from the file, scaled by the current optical model's total dust luminosity.  This preserved energy balance between the two models.  If the filter is between 2.5 and 10 microns, add the optical and IR flux (scaled by the total dust luminosity) from the file.
6)Compute a, the scale factor for this optical-IR model.
7)Check if the scaled model's fluxes are above the user-specified limiting fluxes in any of the filters.  Reject models that are.
8)For all models that pass this test, compute the chi squared value for this model.  
9)Print out model information for all models with a chi squared > 50.  Models with chi squareds greater have a <1e-5 chance of being the true model.  Note that this may change if more degrees of freedom (>12) are added to the program.
10)Finally, the program reads in the optical and IR binary files and looks up additional parameters for each acceptable model (e.g., dust mass, ssfr).  It converts the files to .fits for easier reading.


Created by: Alex Truebenbach
Last Editted: June 2014
"""
from sed_model_functions import *
from pylab import *
import csv,os,sys,array,struct,pickle,time
import numpy as np
from mpi4py import MPI

#Constants
H0=70 #km/s/Mpc

os.nice(20)
rank=MPI.COMM_WORLD.rank
size=MPI.COMM_WORLD.size
t=time.time()
#Read in the observations input file.  
input='sed_objects.txt'
flux,e_flux,limit,id=read_objects(input)

#Read in the luminosity distances
[zt,dlt]=pickle.load(open('z.save','rb'))
z_dl=[]
dl=[]
z_dl.extend(zt)
dl.extend(dlt)

#Info about the filters
allfilters=['u','g','r','i','z','J','H','Ks','W1','W2','W3','W4']
filt_lambda=[.3551,.4686,.6165,.7481,.8931,1.25,1.65,2.17,3.4,4.6,12.0,22.0] #microns

#Create dictionary of previously calculated magnitudes
master_optmags={}
master_irmags={}
zlist=np.linspace(0,5,51)
zlist[0]=.0001

fmu_ir=[]
fmu_opt=[]
ldust_filt=[]

for objindex, obj in enumerate(id):
    if rank==0:
        print 'Fitting for '+str(obj)
    #Separate fluxes into limits and actual detections
    det_flux=[] #Jy
    lim_flux=[]
    lam_detflux=[] #Microns
    lam_limflux=[]
    det_eflux=[]
    for index,fl in enumerate(flux[objindex]):
        if limit[objindex][index]==1:
            det_flux.append(fl)
            det_eflux.append(e_flux[objindex][index])
            lam_detflux.append(filt_lambda[index])
        if limit[objindex][index]==2:
            lim_flux.append(fl)
            lam_limflux.append(filt_lambda[index])


    filters=lam_detflux+lam_limflux
    filters.sort()

    if rank==0:
        
        #Open a file for this object
        f=open('results/pdf_info_'+str(obj)+'.txt','wb')
        f.write('# Chi_squared dof optical_model_index IR_model_index a redshift ldust\n')

    for zindex,z in enumerate(zlist): #Don't set z=0.  Get a zero distance.  Screws up later calculations!
        if rank==0:
            print 'z = '+str(z)
            print 'time per z '+str(time.time()-t)
        #Find closest z in z_dl
        ind=min(range(len(z_dl)),key=lambda k:abs(z_dl[k]-z))
        if abs(z_dl[ind]-z)< .0001:
            dist=dl[ind]*3.086e24 #cm
        else:
            print 'Error, '+str(z)+' not in z.save!'
            continue
        #Check how many filters with detections or upper limits will be optical and how many will be infrared.
        ninf=0
        nopt=0
        #Convert those filters that'll be used from Jy to Lsolar/Hz and calculate w (1/error^2).
        det_flux_z=[]
        lam_detflux_z=[]
        lim_flux_z=[]
        lam_limflux_z=[]
        w=[]
        for filt in filters:
            if filt/(1+z)<10 and filt/(1+z)>.005: #Filter will be optical at this redshift
                nopt+=1
            if filt/(1+z)>2.5: #Filter will be ir at this redshift
                ninf+=1
            if filt/(1+z)>.005: #Filter will be used at this redshift.  If filt<.005 microns, models will not make a significant contribution to the flux.
                if filt in lam_detflux: #Filter is a definite detection
                    lam_detflux_z.append(filt) #Observer frame
                    det_flux_z.append((det_flux[lam_detflux.index(filt)]*10**-23*2.604845e-34*4*pi*dist**2)/(1+z)) #Convert flux from Jy to Lsolar/Hz.  Emitted frame
                    w.append(1./((det_eflux[lam_detflux.index(filt)]*10**-23*2.604845e-34*4*pi*dist**2)**2)/(1+z))
                else: #Filter is an upper limit
                    lim_flux_z.append((lim_flux[lam_limflux.index(filt)]*10**-23*2.604845e-34*4*pi*dist**2)/(1+z))
                    lam_limflux_z.append(filt)

        
        opt_filters=[]
        ir_filters=[]
        
        #Read in the file for this redshift
        fname='modelflux-'+str(z)+'.tbl'
        if fname not in os.listdir('modelflux') and rank==0: #Check to make sure relavent file exists
            print 'ERROR: no model flux information for this redshift ('+str(z)+')!'
            continue
        reader=csv.reader(open('modelflux/'+fname,'rb'))
        data=[]
        data.extend(reader)
        data_length=len(data[3])
        if rank==0:
            print 'Fitting data for '+str(data_length)+' models'
            
        opt_mags=np.zeros(shape=(nopt,data_length))  
        ir_mags=np.zeros(shape=(ninf,data_length))   
        
        if round(float(data[1][0]),2)!=round(float(z),2):
            print 'ERROR: redshift in file '+fname+' doesnt match current program redshift of '+str(z)
            sys.exit()
        if len(fmu_opt)==0: #Next 3 arrays don't change with redshift so just read them in from the file once.
            fmu_opt=data[3]
        if len(fmu_ir)==0:
            fmu_ir=data[5]
        if len(ldust_filt)==0:
            ldust_filt=data[7]
        #Check to make sure system isn't overloaded.:
        m,n,o=os.getloadavg()
        if m>=20.0: #system is getting overloaded.  Kill program!
            print 'ERROR: System overload'
            sys.exit()
        
        curr_ir=0
        curr_opt=0
            
        for rindex,row in enumerate(data[9:]):
            #Check to make sure system isn't overloaded.:
            m,n,o=os.getloadavg()
            if m>=20.0: #system is getting overloaded.  Kill program!
                print 'ERROR: System overload'
                sys.exit()
            if 'filter' in row[0]:
                line=row[0].split('=')
                if float(line[1]) in lam_detflux_z or float(line[1]) in lam_limflux_z:
                    if row[1]=='opt':
                        opt_filters.append(float(line[1]))
                        opt_mags[curr_opt]=data[9+rindex+1]
                        curr_opt+=1
                    else:
                        ir_filters.append(float(line[1]))
                        ir_mags[curr_ir]=data[9+rindex+1]   
                        curr_ir+=1
                                  
        #Clear reader and data array to save space
        data=[]
        reader=[]

            
        #if rank==0:
            #print 'Optical and IR Filters: '
            #print opt_filters, ir_filters

        ###Chi squared portion
        const=[0]*6
        #DOF = number of detections - 1 b/c each model was fit to data using 1 parameter (a).

        dof=len(det_flux_z)-1
        for oindex in range(data_length): 
            if mod(oindex,10)==0:
                #Check to make sure system isn't overloaded.:
                m,n,o=os.getloadavg()
                if m>=20.0: #system is getting overloaded.  Kill program!
                    print 'ERROR: System overload'
                    sys.exit()
            #Check if this model had an ok formation age
            if opt_mags[0,oindex]!= -99:
                if len(ir_filters)!=0: #Make sure the object actually contains filters where IR models matter
                    counter=0
                    for iindex in range(data_length): 
                        if mod(iindex,10)==0:
                            #Check to make sure system isn't overloaded.:
                            m,n,o=os.getloadavg()
                            if m>=20.0: #system is getting overloaded.  Kill program!
                                print 'ERROR: System overload'
                                sys.exit()
                        #If stellar and ir models have ~ the same amount of dust luminosity from the ISM and stellar birth clouds,
                        #then the 2 models are compatible.
                        if counter == size-1:
                            counter=0
                            MPI.COMM_WORLD.Barrier()
                            if rank==0:
                                for j in range(size-1):
                                    const=MPI.COMM_WORLD.recv(source=j+1,tag=11)
                                    #Only record if it isn't a dummy array
                                    if const != []:
                                        f.write(str(const[0])+' '+str(const[1])+' '+str(const[2])+' '+str(const[3])+' '+str(const[4])+' '+str(const[5])+' '+str(const[6])+'\n')
                            
                        counter+=1
                        if iindex % (size-1) != rank-1:
                            continue
                        if math.fabs(float(fmu_opt[oindex])-float(fmu_ir[iindex])) <= 0.15:
                            #Scale ir model fluxes to have some total dust luminosity as predicted by
                            #the stellar model -> satisfies energy balance.
                            #Also add together stellar and dust fluxes for filters that have both ir and optical
                            #contributions (redshifted lambda is btw 2.5 and 10 microns)
                            #Separate out the added filters into those that are for filters with definite detections and those that are for filters that serve as upper limits.
                            added_filts_obs,added_filts_lim=add_filters(filters,opt_filters,ir_filters,lam_detflux_z,float(ldust_filt[oindex]),opt_mags[0:,oindex],ir_mags[0:,iindex],1)
                            #Compute the scaling factor
                            a=comp_a(det_flux_z,added_filts_obs,w)
                            #Check is the added fluxes in each band for this model are below the limits
                            good=0
                            for ind,added in enumerate(added_filts_lim):
                                    if lim_flux_z[ind]<a*added:
                                        good=1 #this model is no good
                                        #Send dummy array
                                        MPI.COMM_WORLD.send([],dest=0,tag=11)
                                        break
                            if good==0: #model is good.  compute chisq
                                cs=comp_chisq(det_flux_z,added_filts_obs,w,a)
                                if cs < 50: #If cs is higher, the model clearly doesn't fit and it wont contribute to the pdfs
                                    const=[cs,dof,oindex,iindex,a,z,ldust_filt[oindex]]
                                    MPI.COMM_WORLD.send(const,dest=0,tag=11)

                                else:
                                    #Send dummy array
                                    MPI.COMM_WORLD.send([],dest=0,tag=11)
                        else:
                            #send dummy array
                            MPI.COMM_WORLD.send([],dest=0,tag=11)

                            
                    #Wait for all processors to finish this iindex loop before moving on to the next one
                    MPI.COMM_WORLD.Barrier()
                    if rank==0:
                        N=data_length
                        if N%(size-1)==0:
                            r = size-1
                        else:
                            r = N%(size-1)
                        for j in range(r):
                            const=MPI.COMM_WORLD.recv(source=j+1,tag=11)
                            if const !=[]:
                                f.write(str(const[0])+' '+str(const[1])+' '+str(const[2])+' '+str(const[3])+' '+str(const[4])+' '+str(const[5])+' '+str(const[6])+'\n')
                                
                else: #Otherwise, just use the optical fluxes. Much smaller process, only use processor 0
                    if rank==0:
                        added_filts_obs,added_filts_lim=add_filters(filters,opt_filters,ir_filters,lam_detflux_z,0.0,opt_mags[0:,oindex],[],0)

                        #Calculate the scaling factor
                        a=comp_a(det_flux_z,added_filts_obs,w)
                        #Check is the added fluxes in each band for this model are below the limits
                        good=0
                        for ind,added in enumerate(added_filts_lim):
                            if lim_flux_z[ind]<a*added:
                                good=1 #this model is no good
                                break
                        if good==0: #model is good.  compute reduced chisq
                            cs=comp_chisq(det_flux_z,added_filts_obs,w,a)
                            if cs < 50: #If cs is higher, the model clearly doesn't fit and it wont contribute to the pdfs
                                const=[cs,dof,oindex,-99,a,z,ldust_filt[oindex]] #Put in a dummy number for iindex to indicate no ir model was used.
                                f.write(str(const[0])+' '+str(const[1])+' '+str(const[2])+' '+str(const[3])+' '+str(const[4])+' '+str(const[5])+' '+str(const[6])+'\n')
                                
    if rank==0:
        f.close()

MPI.COMM_WORLD.Barrier()

#####Convert each of the files to fits
#Check to make sure system isn't overloaded.:
m,n,o=os.getloadavg()
if m>=20.0: #system is getting overloaded.  Kill program!
    print 'ERROR: System overload'
    sys.exit()

for object in id:
    if rank==0:
        #Convert model to fits
        os.system('java -jar $HOME/bin/stilts/stilts.jar tpipe \
                 cmd="sort Chi_squared" ifmt=ascii ofmt=fits \
                 in=results/pdf_info_'+str(object)+'.txt \
                 out=results/pdf_info_'+str(object)+'.fits')
        


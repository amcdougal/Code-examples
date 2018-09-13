"""
These functions are used by sed_model.py, selectandplot_sed.py, and plot_sed.py.  This file must be in the same directory as those program files for the programs to run.  See the main program files for information about the programs. 

Created by: Alex Truebenbach
Last Editted: June 2014
"""
import csv
from scipy import interpolate, integrate
from pylab import *

#Constants
H0=70 #km/s/Mpc
    
def read_objects(f):
    """
    Reads in the file of objects that are to be fit by the program.  Fluxes are in Jy.
    Input: file name
    Outputs:
    ---------------------
    flux: flux[object][fluxes]
        -first dimension of array = object
        -2nd dimension: flux for each filter for that object.  list of 12 fluxes.
    e_flux: same structure as flux, except a list of the errors associated with each flux above.  If the flux is an upper limit or not included, then the error is 0.
    limits: same structure.  For each flux, lists either 0, 1, or 2.  0=don't include flux in fit, 1=use flux in the fit, 2=flux is an upper limit
    id: array with one element for each object. Lists the names of the objects.
    """
    reader=csv.reader(open(f,'rb'))
    data=[]
    data.extend(reader)
    data=data[7:]
    id=[]
    flux=[]
    e_flux=[]
    limit=[]
    mag=[]
    emag=[]
    
    for row in data:
        line=row[0]
        line=line.split()
        if line[0]=='1' or line[0]=='0' or line[0]=='2':
            l=[float(line[0]),float(line[1]),float(line[2]),float(line[3]),float(line[4]),float(line[5]),float(line[6]),float(line[7]),float(line[8]),float(line[9]),float(line[10]),float(line[11])]
            limit.append(l)
        else:
            id.append(line[0])
            f=[float(line[1]),float(line[3]),float(line[5]),float(line[7]),float(line[9]),float(line[11]),float(line[13]),float(line[15]),float(line[17]),float(line[19]),float(line[21]),float(line[23])] #Jy 
            ef=[float(line[2]),float(line[4]),float(line[6]),float(line[8]),float(line[10]),float(line[12]),float(line[14]),float(line[16]),float(line[18]),float(line[20]),float(line[22]),float(line[24])] #Jy.  1sigma errors
            flux.append(f)
            e_flux.append(ef)

    return flux,e_flux,limit,id


def comp_a(obs_flux,model_flux,w):
    """
    Compute the scaling factor for a model
    Inputs
    -------
    obs_flux=array of observed fluxes
    model_flux=array of model fluxes
    w=array of 1/error^2 for each observed flux.
    Output
    ------
    a=float of the scaling for the model.  Equivalent to the stellar mass of the model.
    Notes
    ------
    All three input arrays must be the same size.
    """
    if len(obs_flux)!=len(model_flux) or len(model_flux)!=len(w) or len(obs_flux)!=len(w):
        print 'ERROR! Arrays input in comp_chisq are different lengths!'
    num=0
    den=0
    for findex in range(len(obs_flux)):
        num+=model_flux[findex]*obs_flux[findex]*w[findex]
        den+=model_flux[findex]**2*w[findex]
    if den==0:
        a=0
    else:
        a=num/den

    return a

def comp_chisq(obs_flux,model_flux,w,a):
    """
    Calculate chi squared for a model, given a set of data and errors.
    Inputs
    ------
    model_flux=array of model fluxes
    obs_flux=array of observed fluxes
    w=array of 1/error^2 for each observed flux
    a=float, scaling factor by which each flux in model_fluc will be multiplied.
    Output
    -----
    cs=float, the chi squared value for this model.
    Notes
    ------
    model_flux,obs_flux, and w must all be the same length.
    """
    if len(obs_flux)!=len(model_flux) or len(model_flux)!=len(w) or len(obs_flux)!=len(w):
        print 'ERROR! Arrays input in comp_chisq are different lengths!'
    cs=0
    for findex in range(len(model_flux)):
        cs+=((obs_flux[findex]-a*model_flux[findex])**2*w[findex])

    return cs

def add_filters(filters,opt_filters,ir_filters,lam_detflux_z,ldust,opt_mags,ir_mags,ir_flag):
    #Scale ir model fluxes to have some total dust luminosity as predicted by
    #the stellar model -> satisfies energy balance.
    #Also add together stellar and dust fluxes for filters that have both ir and optical
    #contributions (redshifted lambda is btw 2.5 and 10 microns)
    #Separate out the added filters into those that are for filters with definite detections and those that are for filters that serve as upper limits.
    #IR model fluxes are normalized to 1 Lsolar so units are actually [per angstrom]

    added_filts_lim=[]
    added_filts_obs=[]
    if ir_flag==1: #IR mags are also given.
        for filt in filters:
            if filt in opt_filters and filt in ir_filters:
                if filt in lam_detflux_z:
                    added_filts_obs.append(float(opt_mags[opt_filters.index(filt)])+ldust*float(ir_mags[ir_filters.index(filt)]))
                else:
                    added_filts_lim.append(float(opt_mags[opt_filters.index(filt)])+ldust*float(ir_mags[ir_filters.index(filt)]))
            elif filt in opt_filters:
                if filt in lam_detflux_z:
                    added_filts_obs.append(float(opt_mags[opt_filters.index(filt)]))
                else:
                    added_filts_lim.append(float(opt_mags[opt_filters.index(filt)]))
            elif filt in ir_filters:
                if filt in lam_detflux_z:
                    added_filts_obs.append(ldust*float(ir_mags[ir_filters.index(filt)]))
                else:
                    added_filts_lim.append(ldust*float(ir_mags[ir_filters.index(filt)]))
    else: #Just use the optical mags
        for filt in filters:
            if filt in lam_detflux_z:
                added_filts_obs.append(float(opt_mags[opt_filters.index(filt)]))
            else: #Filter must be a limit
                added_filts_lim.append(float(opt_mags[opt_filters.index(filt)]))
    return added_filts_obs,added_filts_lim

def plot_model(z,dist,iindex,a,fprop,wl,fprop_ir,wl_ir,mstr,ldust):
    """
    Plot the given model.  Program converts both optical and IR models to Jy in the observed frame and combines the two templates.  It plots these models vs wavlength in microns in the frame of the observations.  (i.e., it redshifts the wavelengths based on the given redshift).
    Inputs
    --------
    z=float, redshift of the model.  The model will be shifted from this redshift to the frame of the observer.
    dist= float, luminosity distance of this redshfit in cm.
    iindex=integer, index of the infrared model to be plotted.  -99 if no ir model is used.
    a=float, scaling of the models.  Both the optical and ir models will be scaled by this number.
    fprop=array, flux of the optical model. Units=Lsolar/Angstrom
    wl=array, wavelengths of the optical model. Units=Angstrom.  In rest frame of the model.
    fprop_ir=array, flux of the IR model.  Units=Lsolar/Angstrom.  If no IR model was used, then this array is empty.
    wl_ir=array, wavelengths of the IR model.  Units=Angstrom.  In rest frame of the model.
    mstr=float, total stellar mass of the optical model.  Optical model will be normalized to 1 solar mass using this number.
    ldust=float, fraction of luminosity due to dust in the optical model, normalized to 1 solar mass.  The IR model is divided by this amount to preserve the energy balance between the optical and IR models.
    """
    #Convert Lsolar/Ang to Jy.  In rest frame
    fprop=[(i*(wl[j]**2)*a*3.846e33*10**23*(1+z))/(mstr*2.997925e18*4*pi*dist**2) for j,i in enumerate(fprop)]
    if iindex !=-99:
        fprop_ir=[(i*(wl[j]**2)*ldust*a*3.846e33*10**23*(1+z))/(2.997925e18*4*pi*dist**2) for j,i in enumerate(fprop_ir)]
        #Combine the two templates together
        wl_z=[(i*10**-4)*(1.+z) for i in wl]
        wlir_z=[(i*10**-4)*(1.+z) for i in wl_ir]
        
        #Plot the ir template by itself
        plot(wlir_z,fprop_ir)
        #Plot the optical template by itself
        plot(wl_z,fprop)
        
        #Combine the two wavelength ranges
        overlap=min(range(len(wl_z)),key=lambda i:abs(wl_z[i]-min(wlir_z)))
        overlap_u=min(range(len(wlir_z)), key=lambda i:abs(wlir_z[i]-max(wl_z)))
        if wl_z[overlap]>= min(wlir_z):
            overlap-=1
        
        wl_combo=np.concatenate((wl_z[0:overlap],wlir_z[0:overlap_u]))
       
        fir_new=np.concatenate(([0]*(overlap),fprop_ir[0:overlap_u]))
        fopt_new=np.concatenate((fprop,[0]*(len(wl_combo)-len(fprop))))
        #Interpolate the optical SED in this new wavelength grid
        x=interpolate.interp1d(wl_z,fprop,fill_value=0,bounds_error=False)

        #Add together the templates
        f_combo=[fir_new[j]+x(wl_combo[j]) for j in range(len(wl_combo))]
        
    else: 
        f_combo=fprop
        wl_combo=[i*10**-4*(1+z) for i in wl] #microns

    plot(wl_combo,f_combo,zorder=-100) 


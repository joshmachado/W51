### Recalculating core masses with a new flat temperature assumption derived from mean / median NH3 temps ###

from astropy.table import Table
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
from astropy import constants as const
import numpy as np


fp = '/Users/josh/GitHub/W51/'
###Retrieve source data###

# Dendrogram Catalog of Sources from Ginsburg et. al 2016
t = Table.read(fp+'data/dendro_merge_continuum_and_line.ipac', format='ascii.ipac')
# Temperature map derived from VLA D config NH3 observations and pyspeckit fitting
temp_cube = fits.open(fp+'data/par_maps.fits')[0]
wcs = WCS(temp_cube.header)
flux_cube = fits.open(fp+'data/W51_te_continuum_best_noise.fits')[0]
flux_wcs = WCS(flux_cube.header)



#Housekeeping
c = const.c #Speed of light m/s
c = c.to(u.cm/u.s)
k = const.k_B #Botlzmann constant J/K
beta = 1.75
kgas = 0.0114 * (u.cm)**2/u.g #Kappa
nu = 226*u.GHz
nu0 = 271.1*u.GHz
KTemp = np.zeros(len(t))
Sigma_g = [None] * len(t)
mass = np.zeros(len(t))
mass_uncertainty = np.zeros(len(t))
temp_uncertainty = np.zeros(len(t))
flux_uncertainty = np.zeros(len(t))
i = 0

#Determining mean & median ktemps
cat = Table.read(fp+'data/dendro_catalog.tex')
ktemps = np.array(cat['KTemps'])
ktemps = ktemps[1:len(ktemps)].astype('float64')
ktemps = ktemps[ktemps>2.9]
mean = np.mean(ktemps)
median = np.median(ktemps)

flux_x = np.zeros(len(t))
flux_y = np.zeros(len(t))

#Determining mass & uncertainties based off of MEAN NH3 derived temperature

while i < len(t)-1:
    #Grab coords / vals

    ra = float(t[i]['x_cen'])
    dec = float(t[i]['y_cen'])
        

    #Convert coordinates

    coord = SkyCoord(ra=ra*u.degree, dec=dec*u.degree, frame='fk5')
    pixcrd = wcs.celestial.wcs_world2pix(coord.ra.deg,coord.dec.deg,0)

    flux_x[i] = int(pixcrd[0])
    flux_y[i] = int(pixcrd[1])



    ### Some sources in dendrocat lie outside the temperature map. This eliminates sources outside VLA observations F.O.V. 
    if (flux_x[i] <300 and flux_y[i] < 300):


        ktemps = temp_cube.data[0, int(pixcrd[1]), int(pixcrd[0])]
    else:
        ktemps = 0.0

    #Units

    intensity = t[i]['peak_cont_flux'] * u.Jy
    intensity = intensity.to(u.erg/(u.cm**2))
    #intensity = intensity * u.erg / (u.cm)**2

    
    #Populate Table with corresponding source KTemp MEAN
    KTemp[i] = mean

    if (float(ktemps) != 0.0)  and (float(ktemps) != 2.7315):

        #Units
        tempNH3 = mean * u.K

        #Compute surface density & mass
            #Sigma_g = intensity * c^2/(2*k*T*kgas) * (nu^-2) * (nu/nu0)^-beta

        Sigma_g[i] = ((intensity)*(c**2)*((2*k.cgs*tempNH3*kgas)**(-1))*((nu.cgs)**(-2))*((nu/nu0)**(-beta))).value

            #mass  = Sigma_g * d^2 / mass_uncertainty
            
        mass[i] = (Sigma_g[i] * (5.41*u.kpc.to(u.cm))**2 / (u.Msun.to(u.g)))
        #5.41pc is current distant estimate to W51

        #Compute uncertainties

        temp_uncertainty[i] = temp_cube.data[6, int(pixcrd[1]), int(pixcrd[0])]

        flux_rms = flux_cube.data[int(flux_x[i]),int(flux_y[i])]*u.Jy
        flux_uncertainty[i] = (flux_rms.cgs).value


        #From variance formula. dmass/dintensity (dmdi), dmass/dtempNH3 (dmdT)

        #dSigma_gdi = c^2/(2*k*T*kgas) * (nu^-2) * (nu/nu0)^-beta
        dSigma_gdi = ((c**2)*((2*k.cgs*tempNH3*kgas)**(-1))*((nu.cgs)**(-2))*((nu/nu0)**(-beta))).value #Sigma_g temp component w/ uncertainty

        
        #dmdi = dSigma_gdi * d^2 / Msun
        #dmdi = dSigma_gdi * (float((5.41*u.kpc.to(u.cm)))**2 / (u.Msun.to(u.g)) #converting from mass density to mass
        dmdi = (dSigma_gdi) * float(5.41*u.kpc.to(u.cm))**2 / (u.Msun.to(u.g))

        
        #dSigma_gdT =   c^2/(2*k*kgas) * (nu^-2) * (nu/nu0)^-beta * T^-2 
        dSigma_gdT = ((intensity)* (c**2)*((2*k.cgs*kgas)**(-1)) *((nu.cgs)**(-2))*((nu/nu0)**(-beta)) * (1/(tempNH3**2))).value #Sigma_ g flux component w/ uncertainty

        #dmdT = dSigma_gdT * d^2 / Msun
        dmdT = (dSigma_gdT * (5.41*u.kpc.to(u.cm))**2 / (u.Msun.to(u.g))) #converting from mass density to mass


        #mass_uncertainty = sqrt(dmdi^2 * flux_uncertainty^2 + dmdT^2 * temp_uncertainty^2)
        mass_uncertainty[i] = (((dmdT**2)*temp_uncertainty[i+1]**2 + (dmdi**2)*flux_uncertainty[i]**2))**(0.5)
    else: #For no data
        Sigma_g[i] = 0
        mass[i] = 0
        mass_uncertainty[i] = 0
    i += 1

#Update table
t['mean nh3 mass'] = mass
i=0
#Determining mass & uncertainties based off of MEDIAN NH3 derived temperature
mass = [None] * len(t)
while i < len(t)-1:
    #Grab coords / vals

    ra = float(t[i]['x_cen'])
    dec = float(t[i]['y_cen'])
        

    #Convert coordinates

    coord = SkyCoord(ra=ra*u.degree, dec=dec*u.degree, frame='fk5')
    pixcrd = wcs.celestial.wcs_world2pix(coord.ra.deg,coord.dec.deg,0)

    flux_x[i] = int(pixcrd[0])
    flux_y[i] = int(pixcrd[1])



    ### Some sources in dendrocat lie outside the temperature map. This eliminates sources outside VLA observations F.O.V. 
    if (flux_x[i] <300 and flux_y[i] < 300):


        ktemps = temp_cube.data[0, int(pixcrd[1]), int(pixcrd[0])]
    else:
        ktemps = 0.0

    #Units

    intensity = t[i]['peak_cont_flux'] * u.Jy
    intensity = intensity.to(u.erg/(u.cm**2))
    #intensity = intensity * u.erg / (u.cm)**2

    
    #Populate Table with corresponding source KTemp MEDIAN
    KTemp[i] = median

    if (float(ktemps) != 0.0)  and (float(ktemps) != 2.7315):

        #Units
        tempNH3 = median * u.K

        #Compute surface density & mass
            #Sigma_g = intensity * c^2/(2*k*T*kgas) * (nu^-2) * (nu/nu0)^-beta

        Sigma_g[i] = ((intensity)*(c**2)*((2*k.cgs*tempNH3*kgas)**(-1))*((nu.cgs)**(-2))*((nu/nu0)**(-beta))).value

            #mass  = Sigma_g * d^2 / mass_uncertainty
            
        mass[i] = (Sigma_g[i] * (5.41*u.kpc.to(u.cm))**2 / (u.Msun.to(u.g)))
        #5.41pc is current distant estimate to W51

        #Compute uncertainties

        temp_uncertainty[i] = temp_cube.data[6, int(pixcrd[1]), int(pixcrd[0])]

        flux_rms = flux_cube.data[int(flux_x[i]),int(flux_y[i])]*u.Jy
        flux_uncertainty[i] = (flux_rms.cgs).value


        #From variance formula. dmass/dintensity (dmdi), dmass/dtempNH3 (dmdT)

        #dSigma_gdi = c^2/(2*k*T*kgas) * (nu^-2) * (nu/nu0)^-beta
        dSigma_gdi = ((c**2)*((2*k.cgs*tempNH3*kgas)**(-1))*((nu.cgs)**(-2))*((nu/nu0)**(-beta))).value #Sigma_g temp component w/ uncertainty

        
        #dmdi = dSigma_gdi * d^2 / Msun
        #dmdi = dSigma_gdi * (float((5.41*u.kpc.to(u.cm)))**2 / (u.Msun.to(u.g)) #converting from mass density to mass
        dmdi = (dSigma_gdi) * float(5.41*u.kpc.to(u.cm))**2 / (u.Msun.to(u.g))

        
        #dSigma_gdT =   c^2/(2*k*kgas) * (nu^-2) * (nu/nu0)^-beta * T^-2 
        dSigma_gdT = ((intensity)* (c**2)*((2*k.cgs*kgas)**(-1)) *((nu.cgs)**(-2))*((nu/nu0)**(-beta)) * (1/(tempNH3**2))).value #Sigma_ g flux component w/ uncertainty

        #dmdT = dSigma_gdT * d^2 / Msun
        dmdT = (dSigma_gdT * (5.41*u.kpc.to(u.cm))**2 / (u.Msun.to(u.g))) #converting from mass density to mass


        #mass_uncertainty = sqrt(dmdi^2 * flux_uncertainty^2 + dmdT^2 * temp_uncertainty^2)
        mass_uncertainty[i] = (((dmdT**2)*temp_uncertainty[i+1]**2 + (dmdi**2)*flux_uncertainty[i]**2))**(0.5)
    else: #For no data
        Sigma_g[i] = 0
        mass[i] = 0
        mass_uncertainty[i] = 0
    i += 1

#Update table
t['median nh3 mass'] = mass



t.write(fp+'data/dendro_catalog_mean_med.tex', format='latex', overwrite=True)



### Recalculating core masses with 20K flat temperature assumption to see if it yields the same results as the original catalog ###

from astropy.table import Table
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
from astropy import constants as const
import numpy as np

fp = '/Users/josh/GitHub/W51/'

#Retrieve source data
t = Table.read(fp+'data/byeye_catalog.tex')
temp_cube = fits.open(fp+'data/par_maps.fits')[0]
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
Sigma_g = [None] * len(t)
mass = [None] * len(t)
mass_uncertainty = [None] * len(t)
temp_uncertainty = np.zeros(len(t))
flux_uncertainty = np.zeros(len(t))
i = 0

tempNH3 = 20 * u.K

#Determining mass & uncertainties based off of MEAN NH3 derived temperature

while i < len(t)-1:
    if float(t['KTemp'][i+1]) > 2.9:

        #Grab coords / vals

        ra = float(t[i+1][1])
        dec = float(t[i+1][2])

        coord = SkyCoord(ra=ra*u.degree, dec=dec*u.degree, frame='fk5')
        pixcrd = flux_wcs.celestial.wcs_world2pix(coord.ra.deg,coord.dec.deg,0)

        flux_x = int(pixcrd[0])
        flux_y = int(pixcrd[1])

        temp_x = int(t[i+1]['xval'])
        temp_y = int(t[i+1]['yval'])

        #Units

        intensity = t[i+1]['peak'] * u.Jy
        intensity = intensity.to(u.erg/(u.cm**2))
        intensity = intensity * u.erg / (u.cm)**2
    

        #Compute surface density & mass
            #Sigma_g = intensity * c^2/(2*k*T*kgas) * (nu^-2) * (nu/nu0)^-beta

        Sigma_g[i+1] = float(t['beam_area'][i+1])**-1 * (intensity)*(c**2)*((2*k.cgs*tempNH3*kgas)**(-1))*((nu.cgs)**(-2))*((nu/nu0)**(-beta))

            #mass  = Sigma_g * beam_area * d^2 / Msun
            
        mass[i+1] = (Sigma_g[i+1] * float(t['beam_area'][i+1]) * (5.41*u.kpc.to(u.cm))**2 / (u.Msun.to(u.g))).value

        #Compute uncertainties

        temp_uncertainty[i+1] = temp_cube.data[6, temp_x, temp_y]

        flux_rms = flux_cube.data[flux_x,flux_y]*u.Jy
        flux_uncertainty[i+1] = (flux_rms.cgs).value


        #From variance formula. dmass/dintensity (dmdi), dmass/dtempNH3 (dmdT)

        #dSigma_gdi = 1/beam_area * c^2/(2*k*T*kgas) * (nu^-2) * (nu/nu0)^-beta
        dSigma_gdi = (float(t['beam_area'][i+1])**-1 * (c**2)*((2*k.cgs*tempNH3*kgas)**(-1))*((nu.cgs)**(-2))*((nu/nu0)**(-beta))).value #Sigma_g temp component w/ uncertainty

        
        #dmdi = dSigma_gdi * beam_area * d^2 / Msun
        dmdi = (dSigma_gdi * float(t['beam_area'][i+1]) * (5.41*u.kpc.to(u.cm))**2 / (u.Msun.to(u.g))) #converting from mass density to mass

        
        #dSigma_gdT = -1/beam_area * c^2/(2*k*kgas) * (nu^-2) * (nu/nu0)^-beta * T^-2 
        dSigma_gdT = ((1/float(t['beam_area'][i+1])) *(intensity)* (c**2)*((2*k.cgs*kgas)**(-1)) *((nu.cgs)**(-2))*((nu/nu0)**(-beta)) * (1/(tempNH3**2))).value #Sigma_ g flux component w/ uncertainty

        #dmdT = dSigma_gdT * beam_area * d^2 / Msun
        dmdT = (dSigma_gdT * float(t['beam_area'][i+1]) * (5.41*u.kpc.to(u.cm))**2 / (u.Msun.to(u.g))) #converting from mass density to mass


        #mass_uncertainty = sqrt(dmdi^2 * flux_uncertainty^2 + dmdT^2 * temp_uncertainty^2)
        mass_uncertainty[i+1] = (((dmdT**2)*temp_uncertainty[i+1]**2 + (dmdi**2)*flux_uncertainty[i+1]**2))**(0.5)
    else: #For no data
        Sigma_g[i+1] = 0
        mass[i+1] = 0
        mass_uncertainty[i+1] = 0
    i += 1

#Update table
t['20K nh3 mass'] = mass

t.write(fp+'data/byeye_catalog.tex', format='latex', overwrite=True)


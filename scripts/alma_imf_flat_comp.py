###Imports###
from astropy.table import Table
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import ascii, fits
from astropy.wcs import WCS
from astropy import wcs
from astropy import constants as const
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import powerlaw

fp = '/Users/josh/GitHub/W51/'
###Retrieve source data###


#RUN ALMA_IMF_CAT.PY FIRST
cat_w51e = ascii.read('/Users/josh/GitHub/W51/data/W51-E.sw.sources.fin.ok_dec2020.cat', data_start=0, format='commented_header', header_start=110,  comment="!")
cat_w51irs2 = ascii.read('/Users/josh/GitHub/W51/data/W51-IRS2.sw.sources.fin.ok_dec2020.cat', data_start=0, format='commented_header', header_start=110, comment="!")
parmap = fits.open('/Users/josh/GitHub/W51/data/par_maps.fits')
ww = wcs.WCS(parmap[0].header).celestial
xpix, ypix = ww.wcs_world2pix(cat_w51e['WCS_ACOOR'], cat_w51e['WCS_DCOOR'], 0)
temperatures_w51e = parmap[0].data[0][ypix.astype('int'), xpix.astype('int')]
xpix, ypix = ww.wcs_world2pix(cat_w51irs2['WCS_ACOOR'], cat_w51irs2['WCS_DCOOR'], 0)
temperatures_irs2 = parmap[0].data[0][ypix.astype('int'), xpix.astype('int')]


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
dist = 5.41*u.kpc
k_nu = kgas * (nu/nu0)**beta
Sigma_g_w51e = [None] * len(cat_w51e)
Sigma_g_irs2 = [None] * len(cat_w51irs2)
mass_w51e = np.zeros(len(cat_w51e))
mass_irs2 = np.zeros(len(cat_w51irs2))


###NEED FLUXES, TEMPS, COORDS
temps = np.append(temperatures_w51e, temperatures_irs2)
flux_x = np.zeros(len(temps))
flux_y = np.zeros(len(temps))


#CALC FOR w51e
for i in range(len(temperatures_w51e)):
    ra = cat_w51e['WCS_ACOOR'][i]
    dec = cat_w51e['WCS_DCOOR'][i]
    #Convert coordinates
    coord = SkyCoord(ra=ra*u.degree, dec=dec*u.degree, frame='fk5')
    pixcrd = wcs.celestial.wcs_world2pix(coord.ra.deg,coord.dec.deg,0)
    #corresponding pixels on flux map
    flux_x[i] = int(pixcrd[0])
    flux_y[i] = int(pixcrd[1])
    #corresponding pixel on temp map
    xpix, ypix = ww.wcs_world2pix(cat_w51e['WCS_ACOOR'][i], cat_w51e['WCS_DCOOR'][i], 0)
    ### Some sources in catalog lie outside the temperature map. This eliminates sources outside VLA observations F.O.V.
    if (flux_x[i] <300 and flux_y[i] < 300):


        ktemps = temp_cube.data[0, int(pixcrd[1]), int(pixcrd[0])]

    else:
        ktemps = 0.0

    #Units
    intensity = cat_w51e['FXP_BST02'][i] *u.Jy
    intensity = intensity.to(u.erg/(u.cm)**2)
    tempNH3 = temp_cube.data[0, int(xpix), int(ypix)] * u.K

    if np.isfinite(tempNH3):
        mass_w51e[i] = ((intensity * (dist.to(u.cm))**2 * c**2) / (2* k_nu * nu **2 * k.cgs * tempNH3)).to(u.Msun).value
    else:
        mass_w51e[i] = 0
#CALC FOR IRS2
for i in range(len(temperatures_irs2)):
    ra = cat_w51irs2['WCS_ACOOR'][i]
    dec = cat_w51irs2['WCS_DCOOR'][i]
    #Convert coordinates
    coord = SkyCoord(ra=ra*u.degree, dec=dec*u.degree, frame='fk5')
    pixcrd = wcs.celestial.wcs_world2pix(coord.ra.deg,coord.dec.deg,0)
    #corresponding pixels on flux map
    flux_x[i] = int(pixcrd[0])
    flux_y[i] = int(pixcrd[1])
    #corresponding pixel on temp map
    xpix, ypix = ww.wcs_world2pix(cat_w51irs2['WCS_ACOOR'][i], cat_w51irs2['WCS_DCOOR'][i], 0)
    ### Some sources in catalog lie outside the temperature map. This eliminates sources outside VLA observations F.O.V.
    if (flux_x[i] <300 and flux_y[i] < 300):


        ktemps = temp_cube.data[0, int(pixcrd[1]), int(pixcrd[0])]

    else:
        ktemps = 0.0

    #Units
    intensity = cat_w51irs2['FXP_BST02'][i] *u.Jy
    intensity = intensity.to(u.erg/(u.cm)**2)
    tempNH3 = temp_cube.data[0, int(xpix), int(ypix)] * u.K

    if np.isfinite(tempNH3):
        mass_irs2[i] = ((intensity * (dist.to(u.cm))**2 * c**2) / (2* k_nu * nu **2 * k.cgs * tempNH3)).to(u.Msun).value
    else:
        mass_irs2[i] = 0

masses = np.append(mass_w51e, mass_irs2)
masses = np.sort(masses)
masses = np.trim_zeros(masses)


#Flat temp assumption
temp = 20 * u.K
for i in range(len(temperatures_w51e)):
    ra = cat_w51e['WCS_ACOOR'][i]
    dec = cat_w51e['WCS_DCOOR'][i]
    #Convert coordinates
    coord = SkyCoord(ra=ra*u.degree, dec=dec*u.degree, frame='fk5')
    pixcrd = wcs.celestial.wcs_world2pix(coord.ra.deg,coord.dec.deg,0)
    #corresponding pixels on flux map
    flux_x[i] = int(pixcrd[0])
    flux_y[i] = int(pixcrd[1])
    #corresponding pixel on temp map
    xpix, ypix = ww.wcs_world2pix(cat_w51e['WCS_ACOOR'][i], cat_w51e['WCS_DCOOR'][i], 0)
    ### Some sources in catalog lie outside the temperature map. This eliminates sources outside VLA observations F.O.V.
    if (flux_x[i] <300 and flux_y[i] < 300):


        ktemps = temp_cube.data[0, int(pixcrd[1]), int(pixcrd[0])]

    else:
        ktemps = 0.0

    #Units
    intensity = cat_w51e['FXP_BST02'][i] *u.Jy
    intensity = intensity.to(u.erg/(u.cm)**2)
    tempNH3 = temp_cube.data[0, int(xpix), int(ypix)] * u.K

    if np.isfinite(tempNH3) and tempNH3 != 0:
        mass_w51e[i] = ((intensity * (dist.to(u.cm))**2 * c**2) / (2* k_nu * nu **2 * k.cgs * temp)).to(u.Msun).value
    else:
        mass_w51e[i] = 0
for i in range(len(temperatures_irs2)):
    ra = cat_w51irs2['WCS_ACOOR'][i]
    dec = cat_w51irs2['WCS_DCOOR'][i]
    #Convert coordinates
    coord = SkyCoord(ra=ra*u.degree, dec=dec*u.degree, frame='fk5')
    pixcrd = wcs.celestial.wcs_world2pix(coord.ra.deg,coord.dec.deg,0)
    #corresponding pixels on flux map
    flux_x[i] = int(pixcrd[0])
    flux_y[i] = int(pixcrd[1])
    #corresponding pixel on temp map
    xpix, ypix = ww.wcs_world2pix(cat_w51irs2['WCS_ACOOR'][i], cat_w51irs2['WCS_DCOOR'][i], 0)
    ### Some sources in catalog lie outside the temperature map. This eliminates sources outside VLA observations F.O.V.
    if (flux_x[i] <300 and flux_y[i] < 300):


        ktemps = temp_cube.data[0, int(pixcrd[1]), int(pixcrd[0])]

    else:
        ktemps = 0.0

    #Units
    intensity = cat_w51irs2['FXP_BST02'][i] *u.Jy
    intensity = intensity.to(u.erg/(u.cm)**2)
    tempNH3 = temp_cube.data[0, int(xpix), int(ypix)] * u.K

    if np.isfinite(tempNH3) and tempNH3 != 0:
        mass_irs2[i] = ((intensity * (dist.to(u.cm))**2 * c**2) / (2* k_nu * nu **2 * k.cgs * temp)).to(u.Msun).value
    else:
        mass_irs2[i] = 0

flat_masses = np.append(mass_w51e, mass_irs2)
flat_masses = np.sort(flat_masses)
flat_masses = np.trim_zeros(flat_masses)

plt.hist(masses, bins=np.logspace(np.log(0.1),np.log(1000), 25), alpha=0.4,
label='Ammonia Derived Temp. (N = '+str(len(masses))+')')
plt.hist(flat_masses, bins=np.logspace(np.log(0.1),np.log(1000), 25), alpha=0.4,
label='Flat 20K (N = '+str(len(masses))+')')

plt.xlim(0.1,1000)
plt.xscale('log')

plt.xlabel('Core Mass '+r'(M$_{\odot}$)')
plt.ylabel('Counts')
plt.legend()
plt.title('ALMA-IMF Core Mass Distribution')
plt.show()

def linearpoly(x,a,b):
    return b + a*x

def fitalinetocdf(masses):
## fits a line to the loglog scale cdf of a function returns the alpha calculated
    popt,pcov = curve_fit(linearpoly, np.log10(np.sort(masses)),np.log10(np.linspace(1,0,len(masses),endpoint=False)))
    gamma = popt[0]
    alpha_est = -gamma +  1
    return  alpha_est


###Fit Line to CDF alpha & gamma values
flat_alpha = fitalinetocdf(flat_masses)
nh3_alpha = fitalinetocdf(masses)

flat_gamma = -flat_alpha+1
nh3_gamma = -nh3_alpha+1

###Powerlaw alpha values
fit1 = powerlaw.Fit(flat_masses)
fit2 = powerlaw.Fit(masses)

pl_flat_alpha = fit1.power_law.alpha
pl_nh3_alpha = fit2.power_law.alpha

###Plot CDF
plt.plot(np.log10(np.sort(flat_masses)),np.log10(np.linspace(1,0,len(flat_masses), endpoint=False)), label=r'Flat Temp., $\alpha$ = -{:0.4f}'.format(pl_flat_alpha))
plt.plot(np.log10(np.sort(masses)),np.log10(np.linspace(1,0,len(masses), endpoint=False)), label=r'NH$_3$ Temp., $\alpha$ = -{:0.4f}'.format(pl_nh3_alpha))


plt.legend()
plt.title('CDF of ALMA-IMF Catalog for Flat v. Dynamic Temp')
plt.show()

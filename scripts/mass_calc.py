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
mass_uncertainty_w51e = np.zeros(len(cat_w51e))
mass_uncertainty_irs2 = np.zeros(len(cat_w51irs2))
temp_uncertainty_w51e = np.zeros(len(cat_w51e))
temp_uncertainty_irs2 = np.zeros(len(cat_w51irs2))
flux_uncertainty_w51e = np.zeros(len(cat_w51e))
flux_uncertainty_irs2 = np.zeros(len(cat_w51irs2))


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


fp = '/Users/josh/GitHub/W51/'
catalog = Table.read(fp+'data/coldnh3_catalog.tex')
dendrocat = Table.read(fp+'data/dendro_catalog.tex')

###Constant temperature derived masses
const_mass = np.array(catalog['peak_mass'])
const_mass = const_mass[1:len(const_mass)].astype('float64')

###Ammonia temp derived masses
mass = np.array(catalog['mass'])
mass = mass[1:len(mass)].astype('float64')

###Dendrocat constant temperature derived masses
const_dendro = np.array(dendrocat['peak_cont_mass'])
const_dendro = const_dendro[1:len(const_dendro)].astype('float64')

###Dendrocat ammonia temp derived masses
ammon_dendro = np.array(dendrocat['mass'])
ammon_dendro =  ammon_dendro[1:len(ammon_dendro)].astype('float64')

plt.hist(mass, bins=np.logspace(np.log(0.1),np.log(1000), 25), fc=(0, 0, 0, 0.35), ls='dashed', lw=3, histtype='step', hatch='/', fill=True,
label='By-Eye Catalog (N = '+str(len(mass))+')')
plt.hist(ammon_dendro, bins=np.logspace(np.log(0.1),np.log(1000), 25), ls='dotted', lw=3, fc=(0, 0, 0, 0.5), histtype='step', fill=True,
label='Dendrogram Catalog (N = '+str(len(ammon_dendro))+')')
plt.hist(masses, bins=np.logspace(np.log(0.1),np.log(1000), 25),ls='solid', lw=3, fc=(0, 0, 0, 0.8), histtype='step', fill=True,
label='ALMA-IMF Catalog (N = '+str(len(masses))+')')

plt.xlim(0.1,1000)
plt.xscale('log')

plt.xlabel('Core Mass '+r'(M$_{\odot}$)')
plt.ylabel('Counts')
plt.legend()
plt.title('Comparing Catalog Core Mass Distribution')
plt.show()
#pl.savefig('/Users/josh/GitHub/W51/fig_products/dendro_mass_hist.pdf')
#pl.close()




plt.show()

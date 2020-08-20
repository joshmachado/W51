#### NH3 Temp vs. 1.3mm Flux####

from astropy.table import Table
import numpy as np
import pylab as pl
from astropy.io import fits
from astropy.wcs import WCS

fp = '/Users/josh/GitHub/W51/'

### Catalogs for Source Positions & Temperatures ###
catalog = Table.read(fp+'data/byeye_catalog.tex')
dendrocat = Table.read(fp+'data/dendro_catalog.tex')



## By - Eye Catalog ##
temps = np.array(catalog['KTemp'])
temps = temps[1:len(temps)].astype(float)
fluxes = np.array(catalog['peak'])
fluxes = fluxes[1:len(fluxes)].astype(float)

##SOURCE UNCERTAINTIES##
temp_unc = np.array(catalog['KTemp uncertainty'])
temp_unc = temp_unc[1:len(temp_unc)].astype(float)
flux_unc = np.array(catalog['flux uncertainty'])
flux_unc = flux_unc[1:len(flux_unc)].astype(float)
#flux_unc = np.array()

## Need to remove sources w/o NH3 data ##
temp_unc = temp_unc[temps != 0]
flux_unc = flux_unc[temps != 0]
fluxes = fluxes[temps != 0]
temps = temps[temps != 0]



## Dendrogram Catalog ##
den_temps = np.array(dendrocat['KTemps'])
den_temps = den_temps[1:len(den_temps)].astype(float)
den_fluxes = np.array(dendrocat['peak_cont_flux'])
den_fluxes = den_fluxes[1:len(den_fluxes)].astype(float)

##SOURCE UNCERTAINTIES##
dendro_temp_unc = np.array(dendrocat['KTemp uncertainty'])
dendro_temp_unc = dendro_temp_unc[1:len(dendro_temp_unc)].astype(float)
dendro_flux_unc = np.array(dendrocat['noise'])
dendro_flux_unc = dendro_flux_unc[1:len(dendro_flux_unc)].astype(float)

## Need to remove sources w/o NH3 data ##
dendro_temp_unc = dendro_temp_unc[den_temps != 0]
dendro_flux_unc = dendro_flux_unc[den_temps != 0]
den_fluxes = den_fluxes[den_temps != 0]
den_temps = den_temps[den_temps != 0]

pl.plot(fluxes, temps, '.', alpha = 0.5, label='By-Eye Catalog')
pl.errorbar(fluxes, temps, yerr=0.015, xerr=flux_unc, ls='none', alpha = 0.5)
pl.plot(den_fluxes, den_temps, '.', alpha = 0.5, label='Dendrogram Catalog')
pl.errorbar(den_fluxes, den_temps, yerr=dendro_temp_unc, xerr=dendro_flux_unc, capsize = 5, ls='none', alpha = 0.5)
pl.xlabel('ALMA 1.3mm Flux (Jy)')
pl.ylabel('NH3 Temperature (K)')
pl.xscale('log')
pl.legend()
pl.show()
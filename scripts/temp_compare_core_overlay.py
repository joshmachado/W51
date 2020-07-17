from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
import pylab as pl
from astropy.table import Table

t = Table.read('/Users/Josh/W51/data/coldnh3_catalog.tex')
t2 = Table.read('/Users/Josh/W51/scripts/coldnh3_catalog_appended.tex')

#FOUR FOUR KTEMPS
ff_file = fits.open('/Users/Josh/W51/data/par_maps.fits')
ff_data = ff_file[0].data[0,:,:]
#FOUR FOUR KTEMPS OF CORES
ffcores = np.array(t['KTemp'])
ffcores = ffcores[1:len(ffcores)].astype('float64')

#COLD AMMONIA KTEMPS
cold_file = fits.open('/Users/Josh/W51/data/cold_par_maps.fits')
cold_data = cold_file[0].data[0,:,:]
#COLD AMMONIA KTEMPS OF CORES
coldcores = np.array(t2['cold_temp'])
coldcores[coldcores=='None'] = 0
coldcores = coldcores[1:len(coldcores)].astype('float64')


pl.plot(ff_data, cold_data, '.', alpha=0.25)
pl.plot(ffcores, coldcores, '*', color='red', alpha=0.5, label='Cores')
pl.xlabel('Four Four KTemp')
pl.ylabel('Cold Ammonia KTemp')
pl.legend()
pl.show()

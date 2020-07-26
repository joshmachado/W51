from astropy import units as u
from astropy.wcs import WCS
import pylab as pl
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
import numpy as np
from matplotlib_scalebar.scalebar import ScaleBar

#Establishing file names & paths
fp = '/Users/josh/GitHub/W51/'

nh3_11 = fp+"data/nh3_11_m0.fits"
par_maps = fp+"data/par_maps.fits"
################################################


#base = get_pkg_data_filename(nh3_11)
base_hdu = fits.open(nh3_11)[0]
base_wcs = WCS(base_hdu.header)
#cube = get_pkg_data_filename(par_maps)
cube_hdu = fits.open(par_maps)[0]
cube_wcs = WCS(cube_hdu.header)

ax = pl.subplot(projection=base_wcs.celestial)
ax.set(xlim=(30,265), ylim=(50,250))
ax.coords[0].set_ticks_visible(False)
ax.coords[1].set_ticks_visible(False)
ax.coords[0].set_ticklabel_visible(False)
ax.coords[1].set_ticklabel_visible(False)
pl.imshow(base_hdu.data, origin='lower')
ax.contour(cube_hdu.data[2,:,:], levels=[15.5,16,16.5,17,17.5,18], colors='white', alpha=0.8)

pl.show()

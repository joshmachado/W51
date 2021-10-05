### Displays regions where 4-4 emission was used

from astropy import units as u
from astropy.wcs import WCS
import pylab as plt
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
import numpy as np
from matplotlib_scalebar.scalebar import ScaleBar

#Establishing file names & paths
fp = '/Users/josh/GitHub/W51/'

m0_11 = fp+'data/m0_11.fits'
nh3_44 = fp+"data/nh3_44.fits"
par_maps = fp+"data/par_maps.fits"
m0_44 = fp+"data/m0_44.fits"

m0_11 = fits.open(m0_11)
m0_44 = fits.open(m0_44)

plt.imshow(m0_11[0].data, origin='lower')
plt.contour(m0_44[0].data, levels=[16.5], origin='lower')
plt.show()

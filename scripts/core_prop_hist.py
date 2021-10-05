from astropy.table import Table
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
fp = '/Users/josh/GitHub/W51/'

#Load Core Catalogs
catalog = Table.read(fp+'data/byeye_catalog.tex')
dendrocat = Table.read(fp+'data/dendro_catalog.tex')

coords = np.array([catalog['PeakRA'][1:], catalog['PeakDec'][1:]]).astype(float)
dendro_coords = np.array([dendrocat['x_cen'], dendrocat['y_cen']])

#Load Parameter Maps
par_map = fits.open(fp+'/data/par_maps.fits')
wcs = WCS(par_map[0].header)
#Tkin
tkin = par_map[0].data[0]
#Column Density
col_den = par_map[0].data[2]
#Line Width
lw = par_map[0].data[3]
#Centroid Velocity
centroid = par_map[0].data[4]

#COORDS
coords = SkyCoord(ra=coords[0]*u.degree, dec=coords[1]*u.degree, frame='fk5')
x,y = wcs.celestial.world_to_pixel(coords)
x,y = x.astype(int), y.astype(int)
x1,y1 = np.array(catalog['xval'][1:]).astype(int), np.array(catalog['yval'][1:]).astype(int)
#DENDROCOORDS
dendro_coords = SkyCoord(ra=dendro_coords[0]*u.degree, dec=dendro_coords[1]*u.degree, frame='fk5')
deny,denx = wcs.celestial.world_to_pixel(dendro_coords)
deny,denx = deny.astype(int), denx.astype(int)
deny[0] = 0
denx[0] = 0

#Trimming out zeros & nans
eye_tkin = tkin[x1,y1][np.isfinite(tkin[x1,y1])]
dendro_tkin = tkin[denx,deny][np.isfinite(tkin[denx,deny])]
dendro_tkin = dendro_tkin[dendro_tkin != 0]

eye_col = col_den[x1,y1][np.isfinite(col_den[x1,y1])]
dendro_col = col_den[denx,deny][np.isfinite(col_den[denx,deny])]
dendro_col = dendro_col[dendro_col != 0]

eye_lw = lw[x1,y1][np.isfinite(lw[x1,y1])]
dendro_lw = lw[denx,deny][np.isfinite(lw[denx,deny])]
dendro_lw = dendro_lw[dendro_lw != 0]

eye_cen = centroid[x1,y1][np.isfinite(centroid[x1,y1])]
dendro_cen = centroid[denx,deny][np.isfinite(centroid[denx,deny])]
dendro_cen = dendro_cen[dendro_cen != 0]

#Starting plotting
figsize = (15,8)
fig, ax = plt.subplots(2, 2, figsize=figsize)
#fig.suptitle('Core Property Distributions')

#TKIN
tkin_bins = np.linspace(0,160,9)
ax[0,0].hist(eye_tkin, bins=tkin_bins, alpha=0.4, label='By-Eye Catalog')
ax[0,0].hist(dendro_tkin, bins=tkin_bins, alpha=0.4, label='Dendrogram Catalog')
ax[0,0].set_xlabel(r'T$_{kin}$ (K)')
ax[0,0].set_ylabel('Counts')
ax[0,0].tick_params(direction='in')
ax[0,0].legend()

#COLUMN DENSITY
col_bins = np.linspace(15,18,7)
ax[1,0].hist(eye_col, bins=col_bins, alpha=0.4, label='By-Eye Catalog')
ax[1,0].hist(dendro_col, bins=col_bins, alpha=0.4, label='Dendrogram Catalog')
ax[1,0].set_xlabel(r'log$_{10}$(N$_{NH_3}$)')
ax[1,0].set_ylabel('Counts')
ax[1,0].tick_params(direction='in')
ax[1,0].legend()

#Line Width
line_bins = np.linspace(0,5,11)
ax[0,1].hist(eye_lw, bins=line_bins, alpha=0.4, label='By-Eye Catalog')
ax[0,1].hist(dendro_lw, bins=line_bins, alpha=0.4, label='Dendrogram Catalog')
ax[0,1].set_xlabel(r'Line Width (km s$^{-1}$)')
ax[0,1].set_ylabel('Counts')
ax[0,1].tick_params(direction='in')
ax[0,1].legend()

#Centroid
cen_bins = np.linspace(47.5,65,8)
ax[1,1].hist(eye_cen, bins=cen_bins, alpha=0.4, label='By-Eye Catalog')
ax[1,1].hist(dendro_cen, bins=cen_bins, alpha=0.4, label='Dendrogram Catalog')
ax[1,1].set_xlabel(r'Centroid Velocity (km s$^{-1}$)')
ax[1,1].set_ylabel('Counts')
ax[1,1].tick_params(direction='in')
ax[1,1].legend()


plt.show()

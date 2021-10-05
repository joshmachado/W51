from matplotlib import cbook
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from mpl_toolkits.axes_grid1 import make_axes_locatable

fp = '/Users/josh/GitHub/W51/'
fig, ax = plt.subplots(figsize=[8, 4], frameon=False)

par_map = fits.open(fp+'/data/par_maps.fits')
lw = par_map[0].data[3]

ax.imshow(lw, origin="lower")
ax.set_xlim(50,275)
ax.set_ylim(50,250)
ax.axis('off')
# inset axes....
axins = ax.inset_axes([0.5, 0.3, 1.55, 0.47])
axins.imshow(lw, origin="lower", vmin=0.1, vmax=0.6, cmap='plasma')
# sub region of the original image
x1, x2, y1, y2 = 182, 205, 125, 160
axins.set_xlim(x1, x2)
axins.set_ylim(y1, y2)
axins.set_xticklabels('')
axins.set_yticklabels('')

ax.indicate_inset_zoom(axins)
divider = make_axes_locatable(axins)
cax = divider.append_axes("right", size="5%", pad=0.05)
fig.colorbar(axins.imshow(lw, origin="lower", vmin=0.1, vmax=0.6, cmap='plasma'), cax=cax)
plt.show()

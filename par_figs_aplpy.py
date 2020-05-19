##############################################
## CREATING PARAMETER MAP FIGURES FOR PAPER ##
## TRYING THIS WITH APLPY TO FIX SCALE BAR  ##
##############################################


import aplpy
from astropy.io import fits
from astropy import units as u
from astropy.wcs import WCS
import pylab as pl
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
import numpy as np
from matplotlib_scalebar.scalebar import ScaleBar



pars = 'par_maps.fits'
cube = fits.open(pars)[0]



#tkin
tkin = aplpy.FITSFigure(cube.data[0,:,:], dimensions=[0,:,:])
tkin.show_colorscale(cmap='plasma', vmin=20,vmax=140)
tkin.add_colorbar()
tkin.colorbar.set_axis_label_font(size=16)
tkin.colorbar.set_axis_label_text('K')
tkin.colorbar.set_font(size=14)
tkin.add_scalebar(0.015)
tkin.scalebar.set_label('80pc')
tkin.add_beam()
tkin.set_title('Kinematic Temperature')
tkin.axis_labels.set_font(size=19)
tkin.tick_labels.set_font(size=14)



#tkin
ax = pl.subplot(projection=cube_wcs.celestial)
ax.set(xlim=(50,270), ylim=(60,240))
ax.coords[0].set_ticks_visible(False)
ax.coords[1].set_ticks_visible(False)
ax.coords[0].set_ticklabel_visible(False)
ax.coords[1].set_ticklabel_visible(False)
ax.set_axis_off
pl.imshow(cube_hdu.data[0,:,:], origin='lower', cmap='plasma', vmin=20, vmax=140)
cbar = pl.colorbar()
#cbar.ax.set_ylabel('K', rotation=0)
pl.box(on=None)
scalebar = ScaleBar(1.53)
pl.gca().add_artist(scalebar)
pl.title('Kinematic Temperature (K)', fontsize=18)
pl.savefig('tkin.pdf')
pl.close()

#tex
ax = pl.subplot(projection=cube_wcs.celestial)
ax.set(xlim=(50,270), ylim=(60,240))
ax.coords[0].set_ticks_visible(False)
ax.coords[1].set_ticks_visible(False)
ax.coords[0].set_ticklabel_visible(False)
ax.coords[1].set_ticklabel_visible(False)
ax.set_axis_off
pl.imshow(cube_hdu.data[1,:,:], origin='lower', cmap='bone')
cbar = pl.colorbar()
#cbar.ax.set_ylabel('K', rotation=0)
pl.box(on=None)
scalebar = ScaleBar(1.53)
pl.gca().add_artist(scalebar)
pl.title('Excitation Temperature (K)', fontsize=18)
pl.savefig('tex.pdf')
pl.close()

#column density
ax = pl.subplot(projection=cube_wcs.celestial)
ax.set(xlim=(50,270), ylim=(60,240))
ax.coords[0].set_ticks_visible(False)
ax.coords[1].set_ticks_visible(False)
ax.coords[0].set_ticklabel_visible(False)
ax.coords[1].set_ticklabel_visible(False)
ax.set_axis_off
pl.imshow(cube_hdu.data[2,:,:], origin='lower', cmap='bone', vmax=17.5)
cbar = pl.colorbar()
#cbar.ax.set_ylabel('log$_{10}$gm$^{-3}$', rotation=0)
pl.box(on=None)
scalebar = ScaleBar(1.53)
pl.gca().add_artist(scalebar)
pl.title('Column Density (log$_{10}$gm$^{-3}$)', fontsize=18)
pl.savefig('col_density.pdf')
pl.close()

#sigma
ax = pl.subplot(projection=cube_wcs.celestial)
ax.set(xlim=(50,270), ylim=(60,240))
ax.coords[0].set_ticks_visible(False)
ax.coords[1].set_ticks_visible(False)
ax.coords[0].set_ticklabel_visible(False)
ax.coords[1].set_ticklabel_visible(False)
ax.set_axis_off
pl.imshow(cube_hdu.data[3,:,:], origin='lower', cmap='summer', vmax=4)
cbar = pl.colorbar()
#cbar.ax.set_ylabel('km s$^{-1}$', rotation=0)
pl.box(on=None)
scalebar = ScaleBar(1.53)
pl.gca().add_artist(scalebar)
pl.title('Line Width (km s$^{-1}$)', fontsize=18)
pl.savefig('sigma.pdf')
pl.close()


#centroid
ax = pl.subplot(projection=cube_wcs.celestial)
ax.set(xlim=(50,270), ylim=(60,240))
ax.coords[0].set_ticks_visible(False)
ax.coords[1].set_ticks_visible(False)
ax.coords[0].set_ticklabel_visible(False)
ax.coords[1].set_ticklabel_visible(False)
ax.set_axis_off
pl.imshow(cube_hdu.data[4,:,:], origin='lower', cmap='winter', vmin=40, vmax=72)
cbar = pl.colorbar()
#cbar.ax.set_ylabel('km s$^{-1}$', rotation=0)
pl.box(on=None)
scalebar = ScaleBar(1.53)
pl.gca().add_artist(scalebar)
pl.title('Centroid Velocity (km s$^{-1}$)', fontsize=18)
pl.savefig('centroid.pdf')
pl.close()

#######################################################
## CREATING PARAMETER DISTRIBUTION FIGURES FOR PAPER ##
#######################################################


import scipy as sp
import pylab as pl
import numpy as np
from astropy.io import fits


#Establishing file names & paths
par_maps = "/Users/Josh/W51/data/par_maps.fits"
################################################

fitcube = fits.open(par_maps)[0]

#### Kinetic Temperature Histogram and Fitting ####
tex = fitcube.data[0,:,:]
tex = tex.ravel()
tex = np.sort(tex)
tex = np.trim_zeros(tex)
#tex = tex[tex<100]
#Here I effectively mask out the poor regions by removing temps above 100K.
#The poorly fit multi-component regions provided values up to 150K.
#There are only a few outlying pixels anyways. It doesn't really change
#the shape of the distribution.


pl.hist(tex, bins=50)
pl.xlabel('Kinetic Temperature (K)', fontsize=14)
pl.ylabel('Counts', fontsize=14)
pl.title('Kinetic Temperature Distribution', fontsize=18)
pl.savefig('ktemp_hist.pdf')

#### Rotational Temperature Histogram and Fitting ####
trot = fitcube.data[1,:,:]
trot = trot.ravel()
trot = np.sort(trot)
trot = np.trim_zeros(trot)
trot = trot[trot<29]
#Same crude masking technique descrbied above. I remove all rotational temp
#values above 29K. All poorly fit regions returned a value of 30K for this
#parameter.

pl.hist(trot, bins=50)
pl.xlabel('Rotational Temperature (K)', fontsize=14)
pl.ylabel('Counts', fontsize=14)
pl.title('Rotational Temperature Distribution', fontsize=18)
pl.savefig('rtemp_hist.pdf')

#### NH3 Column Density Histogram and Fitting ####
colden = fitcube.data[2,:,:]
colden = colden.ravel()
colden = np.sort(colden)
colden = np.trim_zeros(colden)
colden = colden[colden<18]
#I didn't provide any restrictions on which values were included in the 
#distribution. All values (as far as I understand) are within reason.

pl.hist(colden, bins=50)
pl.xlabel('NH3 Column Density (g/cm^2)', fontsize=14)
pl.ylabel('Counts', fontsize=14)
pl.title('NH$_3$ Column Density Distribution', fontsize=18)
pl.savefig('colden_hist.pdf')

#### Line Width Histogram and Fitting ####
linewidth = fitcube.data[3,:,:]
linewidth = linewidth.ravel()
linewidth = np.sort(linewidth)
linewidth = np.trim_zeros(linewidth)
linewidth = linewidth[linewidth<3.5]
#Another crude masking attempt to remove poorly fit multi-emission regions.
#All poorly fit multi-component emission regions had large line widths.

pl.hist(linewidth, bins=50)
pl.xlabel('Line Width (km/s)', fontsize=14)
pl.ylabel('Counts', fontsize=14)
pl.title('Line Width Distribution', fontsize=18)
pl.savefig('sigma_hist.pdf')

#### Centroid Velocity Histogram and Fitting ####
centroid = fitcube.data[4,:,:]
centroid = centroid.ravel()
centroid = np.sort(centroid)
centroid = np.trim_zeros(centroid)
#Centroid values were generally not effected by multi-component regions.
#The values were probably the easiest for the fitter to guess.

pl.hist(centroid, bins=50)
pl.xlabel('Centroid Velocity (km/s)', fontsize=14)
pl.ylabel('Counts', fontsize=14)
pl.title('Centroid Velocity Distribution', fontsize=18)
pl.savefig('vel_hist.pdf')


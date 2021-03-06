import pyspeckit
import astropy
try:
    import astropy.io.fits as pyfits

except ImportError:
    import pyfits
import numpy as np
import os
from astropy.convolution import convolve_fft,Gaussian2DKernel

# set up CASA-like shortcuts
F=False; T=True

# Some optional parameters for the script
# (if False, it will try to load an already-stored version 
# of the file)
fitcube = True

# Mask out low S/N pixels (to speed things up)
mask = pyfits.getdata('badpix_mask.fits')
mask = np.isfinite(mask) * (mask > 0)

# Load the data using a mask
# Then calibrate the data (the data we're loading in this case are in Janskys,
# but we want surface brightness in Kelvin for the fitting process)
cube11 = pyspeckit.Cube('unpb_mosaic_11.fits', maskmap=mask)
cube11.cube *= (13.6 * (300.0 /
    (pyspeckit.spectrum.models.ammonia.freq_dict['oneone']/1e9))**2 *
    1./cube11.header.get('BMAJ')/3600. * 1./cube11.header.get('BMIN')/3600. )
cube11.unit = "K"
cube22 = pyspeckit.Cube('unpb_mosaic_22.fits', maskmap=mask)
cube22.cube *= (13.6 * (300.0 /
        (pyspeckit.spectrum.models.ammonia.freq_dict['twotwo']/1e9))**2 *
        1./cube22.header.get('BMAJ')/3600. * 1./cube22.header.get('BMIN')/3600. )
cube22.unit = "K"
cube44 = pyspeckit.Cube('unpb_mosaic_44.fits', maskmap=mask)
cube44.cube *= (13.6 * (300.0 /
        (pyspeckit.spectrum.models.ammonia.freq_dict['fourfour']/1e9))**2 *
        1./cube44.header.get('BMAJ')/3600. * 1./cube44.header.get('BMIN')/3600. )
cube44.unit = "K"


# Compute an error map.  We use the 1-1 errors for all 3 because they're
# essentially the same, but you could use a different error map for each
# frequency
oneonemomentfn = '11_err_mosaic.fits'
errmap11 = (pyfits.getdata(oneonemomentfn).squeeze() * 13.6 *
            (300.0 /
             (pyspeckit.spectrum.models.ammonia.freq_dict['oneone']/1e9))**2
            * 1./cube11.header.get('BMAJ')/3600. *
            1./cube11.header.get('BMIN')/3600.)
# Interpolate errors across NaN pixels
errmap11[errmap11 != errmap11] = convolve_fft(errmap11,
                                              Gaussian2DKernel(3),
                                              nan_treatment='interpolate')[errmap11 != errmap11]

# Stack the cubes into one big cube.  The X-axis is no longer linear: there
# will be jumps from 1-1 to 2-2 to 4-4.  
cubes = pyspeckit.CubeStack([cube11,cube22,cube44], maskmap=mask)
cubes.unit = "K"

# Make a "moment map" to contain the initial guesses
# If you've already fit the cube, just re-load the saved version
if os.path.exists('mosaic_momentcube.fits'):
    momentcubefile = pyfits.open('mosaic_momentcube.fits')
    momentcube = momentcubefile[0].data
else:
    cube11.mapplot()
    # compute the moment at each pixel
    cube11.momenteach()
    momentcube = cube11.momentcube
    momentcubefile = pyfits.PrimaryHDU(data=momentcube, header=cube11.header)
if astropy.version.major >= 2 or (astropy.version.major==1 and astropy.version.minor>=3):
    momentcubefile.writeto('mosaic_momentcube.fits',overwrite=True)
else:
    momentcubefile.writeto('mosaic_momentcube.fits',clobber=True)

# Create a "guess cube".  Because we're fitting physical parameters in this
# case, we want to make the initial guesses somewhat reasonable
# As above, we'll just reload the saved version if it exists
guessfn = 'mosaic_stacked_guesscube_02.fits'
if os.path.exists(guessfn):
    guesscube = pyfits.open(guessfn)
    guesses = guesscube[0].data
else:
    guesses = np.zeros((6,)+cubes.cube.shape[1:])
    guesses[0,:,:] = 30                     # Kinetic temperature 
    guesses[1,:,:] = 11                     # Excitation  Temp
    guesses[2,:,:] = 15.5                   # log(column)
    guesses[3,:,:] = 3                      # Line width / 5 (the NH3 moment overestimates linewidth)                  
    guesses[4,:,:] = momentcube[1,:,:]      # Line centroid              
    guesses[5,:,:] = 0.5                    # F(ortho) - ortho NH3 fraction (fixed)

    guesscube = pyfits.PrimaryHDU(data=guesses, header=cube11.header)
    if astropy.version.major >= 2 or (astropy.version.major==1 and astropy.version.minor>=3):
        guesscube.writeto(guessfn, overwrite=True)
    else:
        guesscube.writeto(guessfn, clobber=True)

# This bit doesn't need to be in an if statment
if os.path.exists('mosaic_stacked_cube_03.fits'):
    cubes.load_model_fit('mosaic_stacked_cube_03.fits', 6, fittype='ammonia', _temp_fit_loc=(94,250))
    cubes.specfit.parinfo[5]['fixed'] = True
else:
    # excise guesses that fall out of the "good" range
    guesses[0,:,:][guesses[0,:,:] > 150] = 150
    guesses[0,:,:][guesses[0,:,:] < 2.7315] = 2.7315
    guesses[1,:,:][guesses[1,:,:] > 30] = 30
    guesses[1,:,:][guesses[1,:,:] < 2.7315] = 2.7315
    guesses[2,:,:][guesses[2,:,:] > 35] = 35
    guesses[2,:,:][guesses[2,:,:] < 5] = 5
    guesses[4,:,:][guesses[4,:,:] > 75] = 75
    guesses[4,:,:][guesses[4,:,:] < 40] = 40                         

    # do the fits
    # signal_cut means ignore any pixel with peak S/N less than this number
    # In this fit, many of the parameters are limited
    # start_from_point selects the pixel coordinates to start from
    # use_nearest_as_guess says that, at each pixel, the input guesses will be
    # set by the fitted parameters from the nearest pixel with a good fit
    # HOWEVER, because this fitting is done in parallel (multicore=12 means
    # 12 parallel fitting processes will run), this actually means that EACH
    # core will have its own sub-set of the cube that it will search for good 
    # fits. So if you REALLY want consistency, you need to do the fit in serial.
    cubes.fiteach(fittype='ammonia', skip_failed_fits=False, guesses=guesses, maskmap=mask,
            integral=False, verbose_level=3, fixed=[F,F,F,F,F,T], signal_cut=3,
            limitedmax=[T,T,T,T,T,T],
            maxpars=[150,30,25,15,75,1],
            limitedmin=[T,T,T,T,T,T],
            minpars=[2.7315,2.7315,5,0,40,0],
            use_nearest_as_guess=False, start_from_point=(128,143),
            multicore=12,
            errmap=errmap11)

    # Save the fitted parameters in a data cube
    fitcubefile = pyfits.PrimaryHDU(data=np.concatenate([cubes.parcube,cubes.errcube]), header=cubes.header)
    fitcubefile.header['PLANE1'] = 'TKIN'
    fitcubefile.header['PLANE2'] = 'TEX'
    fitcubefile.header['PLANE3'] = 'COLUMN'
    fitcubefile.header['PLANE4'] = 'SIGMA'
    fitcubefile.header['PLANE5'] = 'VELOCITY'
    fitcubefile.header['PLANE6'] = 'FORTHO'
    fitcubefile.header['PLANE7'] = 'eTKIN'
    fitcubefile.header['PLANE8'] = 'eTEX'
    fitcubefile.header['PLANE9'] = 'eCOLUMN'
    fitcubefile.header['PLANE10'] = 'eSIGMA'
    fitcubefile.header['PLANE11'] = 'eVELOCITY'
    fitcubefile.header['PLANE12'] = 'eFORTHO'
    fitcubefile.header['CDELT3'] = 1
    fitcubefile.header['CTYPE3'] = 'FITPAR'
    fitcubefile.header['CRVAL3'] = 0
    fitcubefile.header['CRPIX3'] = 1
    fitcubefile.writeto("mosaic_stacked_cube_03.fits")


# Now do some plotting things
import pylab as pl

# Set the map-to-plot to be the line centroid
cubes.mapplot.plane = cubes.parcube[4,:,:]
cubes.mapplot(estimator=None,vmin=50,vmax=76)

# Set the reference frequency to be the 1-1 line frequency
cubes.xarr.refX = pyspeckit.spectrum.models.ammonia.freq_dict['oneone']
cubes.xarr.refX_unit='Hz'

# If you wanted to view the spectra in velocity units, use this:
#cubes.xarr.convert_to_unit('km/s')
#cubes.plotter.xmin=55
#cubes.plotter.xmax=135

# Now replace the cube's plotter with a "special" plotter
# The "special" plotter puts the 1-1, 2-2, and 4-4 lines in their own separate
# windows

cubes.plot_special = pyspeckit.wrappers.fitnh3.plotter_override
cubes.plot_special_kwargs = {'fignum':3, 'vrange':[25,95]}
#cubes.plot_spectrum(160,99)

# make interactive
pl.ion()
pl.show()

import astropy.units as u
from astropy.utils import data
from spectral_cube import SpectralCube
import aplpy

data = 'cubename.fits' #data cube used to create moment map
cube =  SpectralCube.read(data)

lower = 0.0 #lower velocity bound for spectral slice
upper = 1.0 #upper velocity bound for spectral slice
cube_slice = cube.with_spectral_unit(u.km/u.s, velocity_convention='radio').spectral_slab(lower*u.km/u.s, upper*u.km/u.s)

### MOMENT MAPS ###

moment_0_slice = cube_slice.moment(order=0) #moment map of specified spectral slice
moment_0 = cube.moment(order=0) #moment map of full cube

moment_0.hdu
moment_0_slice.hdu

#generating full cube moment map image
full_name = 'Full_Moment_Map.png' #name of output image
f = aplpy.FITSFigure(moment_0.hdu)
f.show_colorscale()
f.show_colorbar()
f.save(full_name)

#generating spectral slice cube moment map image
slice_name = 'Slice_Moment_Map.png' #name of output image
s = aplpy.FITSFigure(moment_0_slice.hdu)
s.show_colorscale()
s.show_colorbar()
s.save(slice_name)

### PEAK INTENSITY MAPS ###

full_peak = cube.max(axis=0)
slice_peak = cube_slice.max(axis=0)

full_peak.hdu
slice_peak.hdu

#generating full cube peak intensity image
fpeak_name = 'Full_Peak.png' #name of output image
fp = aplpy.FITSFigure(full_peak.hdu)
fp.show_colorscale()
fp.show_colorbar()
fp.save(fpeak_name)

#generating slice cube peak intensity image
speak_name = 'Slice_Peak.png' #name of output image
sp = aplpy.FITSFigure(slice_peak.hdu)
sp.show_colorscale()
sp.show_colorbar()
sp.save(speak_name)

exit
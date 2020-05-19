import aplpy
from astropy.io import fits

#Defining file names

nh3_11 = 'm0_11.fits'
data_11 = fits.open(nh3_11)[0].data

nh3_22 = 'm0_22.fits'
data_22 = fits.open(nh3_22)[0].data

nh3_33 = 'm0_33.fits'
data_33 = fits.open(nh3_33)[0].data

nh3_44 = 'm0_44.fits'
data_44 = fits.open(nh3_44)[0].data



#Making and saving m0 maps

m11 = aplpy.FITSFigure(nh3_11)
m11.show_colorscale(cmap='magma')
m11.add_colorbar()
m11.colorbar.set_axis_label_font(size=16)
m11.colorbar.set_axis_label_text('K km s$^{-1}$')
m11.colorbar.set_font(size=14)
m11.add_scalebar(0.015)
m11.scalebar.set_label('80pc')
#ZOOM IN
center = m11.pixel2world(147,149)
m11.recenter(center[0], center[1], width=0.07, height=0.053)
m11.add_beam()
m11.show_contour(data_11, levels=(-175, -150, -100, 400,800,1000), colors=('white', 'white', 'white','black', 'black', 'black'))
m11.set_title('NH$_3$ (1,1)', fontsize=24)
m11.axis_labels.set_font(fontsize=19)
m11.tick_labels.set_font(size=14)
m11.save('m0_11.pdf')


m22 = aplpy.FITSFigure(nh3_22)
m22.show_colorscale(cmap='magma')
m22.add_colorbar()
m22.colorbar.set_font(size=14)
m22.colorbar.set_axis_label_font(size=16)
m22.colorbar.set_axis_label_text('K km s$^{-1}$')
m22.add_scalebar(0.015)
m22.scalebar.set_label('80pc')
#ZOOM IN
center = m22.pixel2world(147,149)
m22.recenter(center[0], center[1], width=0.07, height=0.053)
m22.add_beam()
m22.show_contour(data_22, levels=(400,800,1000), colors=('black', 'black', 'black'))
m22.set_title('NH$_3$ (2,2)', fontsize=24)
m22.axis_labels.set_font(fontsize=19)
m22.tick_labels.set_font(size=14)
m22.save('m0_22.pdf')


m33 = aplpy.FITSFigure(nh3_33)
m33.show_colorscale(cmap='magma')
m33.add_colorbar()
m33.colorbar.set_font(size=14)
m33.colorbar.set_axis_label_font(size=16)
m33.colorbar.set_axis_label_text('K km s$^{-1}$')
m33.add_scalebar(0.015)
m33.scalebar.set_label('80pc')
#ZOOM IN
center = m33.pixel2world(147,149)
m33.recenter(center[0], center[1], width=0.07, height=0.053)
m33.add_beam()
m33.show_contour(data_33, levels=(400,800,1000), colors=('black', 'black', 'black'))
m33.set_title('NH$_3$ (3,3)', fontsize=24)
m33.axis_labels.set_font(fontsize=19)
m33.tick_labels.set_font(size=14)
m33.save('m0_33.pdf')


m44 = aplpy.FITSFigure(nh3_44)
m44.show_colorscale(cmap='magma')
m44.add_colorbar()
m44.colorbar.set_font(size=14)
m44.colorbar.set_axis_label_font(size=16)
m44.colorbar.set_axis_label_text('K km s$^{-1}$')
m44.add_scalebar(0.015)
m44.scalebar.set_label('80pc')
#ZOOM IN
center = m44.pixel2world(147,149)
m44.recenter(center[0], center[1], width=0.07, height=0.053)
m44.add_beam()
m44.show_contour(data_44, levels=(400,800,1000), colors=('black', 'black', 'black'))
m44.set_title('NH$_3$ (4,4)', fontsize=24)
m44.axis_labels.set_font(fontsize=19)
m44.tick_labels.set_font(size=14)
m44.save('m0_44.pdf')

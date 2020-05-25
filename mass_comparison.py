from astropy.table import Table
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
from matplotlib.ticker import NullFormatter
import pylab as pl


#data
t = Table.read('/Users/Josh/W51/data/coldnh3_catalog.tex')
new_mass = t['mass']
old_mass = t['peak_mass']

nullfmt = NullFormatter()

#axes definition
left, width = 0.1, 0.65
bottom, height = 0.1, 0.65
bottom_h = left_h = left + width + 0.02

rect_scatter = [left, bottom, width, height]
rect_histx = [left, bottom_h, width, 0.2] #new masses
rect_histy = [left_h, bottom, 0.2, height] #old masses

#start with scatter
pl.figure(1, figsize=(8,8))
axScatter = pl.axes(rect_scatter)
axHistx = pl.axes(rect_histx)
axHisty = pl.axes(rect_histy)

axHistx.xaxis.set_major_formatter(nullfmt)
axHisty.yaxis.set_major_formatter(nullfmt)

#scatter plot
axScatter.scatter(new_mass, old_mass)

#limits

binwidth = 10
xymax = 10
lim = (int(xymax/binwidth) + 1) * binwidth

axScatter.set_xlim((-lim,lim))
axScatter.set_ylim((-lim,lim))


bins = np.arange(-lim,lim+binwidth, binwidth)
axHistx.hist(new_mass, bins=bins)
axHisty.hist(old_mass, bins=bins, orientation='horizontal')

axHistx.set_xlim(axScatter.get_xlim())
axHisty.set_ylim(axScatter.get_ylim())


pl.plot(t['mass'],t['peak_mass'],'.')
pl.xscale('log')
pl.yscale('log')
pl.xlabel('Ammonia Constrained Masses', fontsize='15')
pl.ylabel('Ginsburg et al. 2017 Masses', fontsize='15')
pl.title('Mass Comparison', fontsize='18')
x = np.logspace(0,2.1,base=10)
pl.plot(x, x**1)

pl.show()

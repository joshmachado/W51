from astropy.table import Table
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
from matplotlib.ticker import NullFormatter
import pylab as pl


#data
t = Table.read('/Users/adam/work/students/JoshMachado2019/W51/data/coldnh3_catalog.tex', data_start=4)
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
fig1 = pl.figure(1, figsize=(8,8))
fig1.clf()
axScatter = pl.axes(rect_scatter)
axHistx = pl.axes(rect_histx)
axHisty = pl.axes(rect_histy)

#axHistx.xaxis.set_major_formatter(nullfmt)
#axHisty.yaxis.set_major_formatter(nullfmt)

#scatter plot
axScatter.scatter(new_mass, old_mass)

axScatter.set_xscale('log')
axScatter.set_yscale('log')

#limits

binwidth = 10
xymax = 10
#lim = (int(xymax/binwidth) + 1) * binwidth


#axScatter.set_xlim((-lim,lim))
#axScatter.set_ylim((-lim,lim))

xmin, xmax = axScatter.get_xlim()
xmax = 500
xmin=0.5
axScatter.set_xlim(xmin, xmax)
axScatter.set_ylim(xmin, xmax)
axScatter.plot([xmin,xmax], [xmin,xmax], 'k--')
axScatter.plot(np.array([xmin,xmax])*0.25, np.array([xmin,xmax]), 'k:')

bins = np.arange(xmin, xmax, binwidth)
bins = np.logspace(np.log10(xmin), np.log10(xmax), 20)
axHistx.hist(new_mass, bins=bins)
axHisty.hist(old_mass, bins=bins, orientation='horizontal')
axHistx.set_xscale('log')
axHisty.set_yscale('log')

axHistx.set_xticklabels([])
axHistx.tick_params(direction='in')
axHisty.set_yticklabels([])
axHisty.tick_params(direction='in')

axHistx.set_xlim(axScatter.get_xlim())
axHisty.set_ylim(axScatter.get_ylim())


# pl.plot(t['mass'],t['peak_mass'],'.')
# pl.xscale('log')
# pl.yscale('log')
axScatter.set_xlabel('Ammonia Constrained Masses', fontsize='15')
axScatter.set_ylabel('Ginsburg et al. 2017 Masses', fontsize='15')
# pl.title('Mass Comparison', fontsize='18')
# x = np.logspace(0,2.1,base=10)
# pl.plot(x, x**1)

pl.show()
pl.savefig('../fig_products/mass_vs_ammoniamass.png', bbox_inches='tight')

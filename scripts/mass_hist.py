#### MASS HISTOGRAMS ####

from astropy.table import Table
import numpy as np
import pylab as pl


catalog = Table.read('/Users/Josh/W51/data/coldnh3_catalog.tex')

###Constant temperature derived masses
const_mass = np.array(catalog['peak_mass'])
const_mass = const_mass[1:len(const_mass)].astype('float64')

###Ammonia temp derived masses
mass = np.array(catalog['mass'])
mass = mass[1:len(mass)].astype('float64')

pl.hist(const_mass, bins=np.logspace(np.log(0.1),np.log(1000), 25), alpha=0.5, label='Flat Temp. Assumption')
pl.hist(mass, bins=np.logspace(np.log(0.1),np.log(1000), 25), alpha=0.5, label='Ammonia Derived Temp.')
pl.xscale('log')
pl.xlim(0.1,1000)
pl.xlabel('Core Mass')
pl.ylabel('Counts')
pl.legend()
pl.savefig('/Users/Josh/W51/fig_products/mass_hist.pdf')
#### MASS HISTOGRAMS ####

from astropy.table import Table
import numpy as np
import pylab as pl

fp = '/Users/josh/GitHub/W51/'
catalog = Table.read(fp+'data/coldnh3_catalog.tex')
dendrocat = Table.read(fp+'data/dendro_catalog.tex')

###Constant temperature derived masses
const_mass = np.array(catalog['peak_mass'])
const_mass = const_mass[1:len(const_mass)].astype('float64')

###Ammonia temp derived masses
mass = np.array(catalog['mass'])
mass = mass[1:len(mass)].astype('float64')

###Dendrocat constant temperature derived masses
const_dendro = np.array(dendrocat['peak_cont_mass'])
const_dendro = const_dendro[1:len(const_dendro)].astype('float64')

###Dendrocat ammonia temp derived masses
ammon_dendro = np.array(dendrocat['mass'])
ammon_dendro =  ammon_dendro[1:len(ammon_dendro)].astype('float64')

pl.hist(const_mass, bins=np.logspace(np.log(0.1),np.log(1000), 25), alpha=0.5, label='Flat Temp. Assumption')
pl.hist(mass, bins=np.logspace(np.log(0.1),np.log(1000), 25), alpha=0.5, label='Ammonia Derived Temp.')

#pl.hist(const_dendro, bins=np.logspace(np.log(0.1),np.log(1000), 25), alpha=0.5, label='Dendrocat Flat Temp. Assumption')
#pl.hist(ammon_dendro, bins=np.logspace(np.log(0.1),np.log(1000), 25), alpha=0.5, label='Dendrocat Ammonia Derived Temp.')

pl.xscale('log')
pl.xlim(0.1,1000)
pl.xlabel('Core Mass')
pl.ylabel('Counts')
pl.legend()
pl.show()
#pl.savefig('/Users/Josh/W51/fig_products/mass_hist.pdf')
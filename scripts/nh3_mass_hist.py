#### AMMONIA DERIVED MASS HISTOGRAMS ####

from astropy.table import Table
import numpy as np
import pylab as pl

fp = '/Users/josh/GitHub/W51/'
catalog = Table.read(fp+'data/coldnh3_catalog.tex')
dendrocat = Table.read(fp+'data/dendro_catalog.tex')


###Ammonia temp derived masses for 'by-eye' catalog
mass = np.array(catalog['mass'])
mass = mass[1:len(mass)].astype('float64')

###Dendrocat ammonia temp derived masses
ammon_dendro = np.array(dendrocat['mass'])
ammon_dendro =  ammon_dendro[1:len(ammon_dendro)].astype('float64')

pl.hist(mass, bins=np.logspace(np.log(0.1),np.log(1000), 25), alpha=0.5, label='By Eye Ammonia Derived Masses')
pl.hist(ammon_dendro, bins=np.logspace(np.log(0.1),np.log(1000), 25), alpha=0.5, label='Dendrogram Ammonia Derived Masses')


pl.xscale('log')
pl.xlim(0.1,1000)
pl.xlabel('Core Mass')
pl.ylabel('Counts')
pl.legend()
pl.show()
#pl.savefig('/Users/Josh/W51/fig_products/mass_hist.pdf')
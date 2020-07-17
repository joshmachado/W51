### MASS V MASS ###

from astropy.table import Table
import numpy as np
import pylab as pl


###Comparing mass derived from constant temperature assumption

catalog = Table.read('/Users/Josh/W51/data/coldnh3_catalog.tex')
dendrocat = Table.read('/Users/Josh/W51/data/dendro_catalog.tex')


###Constant temperature derived masses
const_mass = np.array(catalog['peak_mass'])
const_mass = const_mass[1:len(const_mass)].astype('float64')

###Ammonia temp derived masses
mass = np.array(catalog['mass'])
mass = mass[1:len(mass)].astype('float64')

###Dendrogram catalog original masses
dendromass = np.array(dendrocat['peak_cont_mass'])
dendromass = dendromass[1:len(dendromass)].astype('float64')

###Dendrogram catalog original masses (T corrected)
tdendro = np.array(dendrocat['T_corrected_mass'])
tdendro = tdendro[1:len(tdendro)].astype('float64')

###Ammonia temp derived DENDRO masses
ammon_dendro = np.array(dendrocat['mass'])
ammon_dendro = ammon_dendro[1:len(ammon_dendro)].astype('float64')

#One to One line
x = np.linspace(0,10000,10000)
y = x

pl.plot(x,y, '-r')
pl.plot(ammon_dendro, dendromass, '.')
pl.xlabel('Ammonia Temp Derived Core Mass')
pl.ylabel('Flat Temp Derived Core Mass')
pl.xscale('log')
pl.yscale('log')
pl.show()
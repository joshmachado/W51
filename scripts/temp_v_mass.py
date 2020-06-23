###COMPARING MASS VS. AMMONIA DERIVED KINETIC TEMPERATURES###

from astropy.table import Table
import numpy as np
import pylab as pl


###Comparing mass derived from constant temperature assumption

catalog = Table.read('/Users/Josh/W51/data/coldnh3_catalog.tex')
dendrocat = Table.read('/Users/Josh/W51/data/dendro_merge_continuum_and_line.ipac', format='ascii.ipac')
ktemps = np.array(catalog['KTemp'])
ktemps = ktemps[1:len(ktemps)].astype('float64')

###Constant temperature derived masses
const_mass = np.array(catalog['peak_mass'])
const_mass = const_mass[1:len(const_mass)].astype('float64')

###Ammonia temp derived masses
mass = np.array(catalog['mass'])
mass = mass[1:len(mass)].astype('float64')

###Dendrogram catalog original masses
dendromass = np.array(dendrocat['peak_cont_mass']).astype('float64')

###Dendrogram catalog original masses (T corrected)
tdendro = np.array(dendrocat['T_corrected_mass']).astype('float64')

#pl.plot(ktemps, const_mass, '.', color='blue', label='Const. Temp Mass')
#pl.plot(ktemps, mass, '.', color='green', label='Ammonia Temp Mass')
#pl.xlabel('Kinetic Temperature')
#pl.ylabel('Core Mass')
#pl.show()

pl.plot(ktemps, dendromass, '.', color='red', label='Dendrogram Peak Cont Mass')
pl.plot(ktemps. tdendro, '.', color='purple', label='Dendrogram T Corrected Mass')
pl.xlabel('Kinetic Temperature')
pl.ylabel('Core Mass')
pl.show()

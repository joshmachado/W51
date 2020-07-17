###COMPARING MASS VS. AMMONIA DERIVED KINETIC TEMPERATURES###

from astropy.table import Table
import numpy as np
import pylab as pl


###Comparing mass derived from constant temperature assumption

dendrocat = Table.read('/Users/Josh/W51/data/dendro_catalog.tex')
ktemps = np.array(dendrocat['KTemps']).astype('float64')
ktemps = ktemps[1:len(ktemps)]



###Ammonia temp derived masses
mass = np.array(dendrocat['mass']).astype('float64')
mass = mass[1:len(mass)]

###Dendrogram catalog original masses
dendromass = np.array(dendrocat['peak_cont_mass']).astype('float64')
dendromass = dendromass[1:len(dendromass)]

###Dendrogram catalog original masses (T corrected)
tdendro = np.array(dendrocat['T_corrected_mass'])
tdendro = tdendro[1:len(dendrocat)].astype('float64')

pl.plot(ktemps, const_mass, '.', color='blue', label='Const. Temp Mass')
pl.plot(ktemps, mass, '.', color='green', label='Ammonia Temp Mass')
pl.xlabel('Kinetic Temperature')
pl.ylabel('Core Mass')
pl.show()

#pl.plot(ktemps, mass, '.', color='blue', label='Ammonia Temp Derived Masses')
pl.plot(ktemps, dendromass, '.', color='red', label='Dendrogram Peak Cont Mass')
pl.plot(ktemps, tdendro, '.', color='purple', label='Dendrogram T Corrected Mass')
pl.xlabel('Kinetic Temperature')
pl.ylabel('Core Mass')
pl.legend()
pl.show()

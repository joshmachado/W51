###COMPARING MASS VS. AMMONIA DERIVED KINETIC TEMPERATURES###

from astropy.table import Table
import numpy as np
import pylab as pl


###Comparing mass derived from constant temperature assumption
fp = '/Users/josh/GitHub/W51/'

catalog = Table.read(fp+'data/coldnh3_catalog.tex')
dendrocat = Table.read(fp+'data/dendro_catalog.tex')

###Non dendro ktemps
ktemps = np.array(catalog['KTemp'])
ktemps = ktemps[1:len(ktemps)].astype('float64')

###Dendro KTemps
dendrokt = np.array(dendrocat['KTemps'])
dendrokt = dendrokt[1:len(dendrokt)].astype('float64')

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

pl.plot(dendrokt, dendromass, '.', color='blue', label='Const. Temp Mass')
pl.plot(dendrokt, ammon_dendro, '.', color='green', label='Ammonia Temp Mass')
pl.xlabel('Kinetic Temperature')
pl.ylabel('Core Mass')
#pl.xscale('log')
pl.yscale('log')
pl.legend()
pl.show()

#pl.plot(ktemps, dendromass, '.', color='red', label='Dendrogram Peak Cont Mass')
#pl.plot(ktemps. tdendro, '.', color='purple', label='Dendrogram T Corrected Mass')
#pl.xlabel('Kinetic Temperature')
#pl.ylabel('Core Mass')
#pl.show()

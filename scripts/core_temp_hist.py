#### CORE TEMP HISTOGRAMS ####

from astropy.table import Table
import numpy as np
import pylab as pl

fp = '/Users/josh/GitHub/W51/'
catalog = Table.read(fp+'data/byeye_catalog.tex')
dendrocat = Table.read(fp+'data/dendro_catalog.tex')

###By-eye Core Temps
byeye = np.array(catalog['KTemp'])
byeye = byeye[1:len(byeye)].astype('float64')


###Dendrocat Core Temps
dendro = np.array(dendrocat['KTemps'])
dendro = dendro[1:len(dendro)].astype('float64')

binlist = [40, 80, 120]
pl.hist(dendro,  bins = binlist, alpha=0.5, label='Dendrogram Sources')
pl.hist(byeye,  bins = binlist, alpha=0.5, label='By-Eye Sources')


pl.xlabel('Core Temperature (K)')
pl.ylabel('Counts')
pl.legend()
pl.title('Core Temperatures')
pl.show()
#pl.savefig('/Users/josh/GitHub/W51/fig_products/dendro_mass_hist.pdf')
#pl.close()
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
byeye = np.trim_zeros(np.sort(byeye))
byeye = byeye[byeye > 3]
###Dendrocat Core Temps
dendro = np.array(dendrocat['KTemps'])
dendro = dendro[1:len(dendro)].astype('float64')
dendro = np.trim_zeros(np.sort(dendro))
dendro = dendro[np.where(np.isnan(dendro) == False)]
dendro = dendro[dendro > 3]
binlist = np.array([np.arange(start=0, stop=140, step=10)])
pl.hist(dendro,  bins = binlist[0], alpha=0.5, label='Dendrogram Sources')
pl.hist(byeye,  bins = binlist[0], alpha=0.5, label='By-Eye Sources')


pl.xlabel('Core Temperature (K)')
pl.ylabel('Counts')
pl.legend()
pl.title('Core Temperatures')
pl.show()
#pl.savefig('/Users/josh/GitHub/W51/fig_products/dendro_mass_hist.pdf')
#pl.close()

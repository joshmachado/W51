#### NH3 Temp vs. 1.3mm Flux####

from astropy.table import Table
import numpy as np
import pylab as pl

fp = '/Users/josh/GitHub/W51/'
catalog = Table.read(fp+'data/byeye_catalog.tex')
dendrocat = Table.read(fp+'data/dendro_catalog.tex')
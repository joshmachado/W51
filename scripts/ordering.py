#### Sorting sources by flux ####
#### Does order / index of sources change when sorting by 1.3mm flux vs NH3 derived mass? ####


#Sorting sources by 1.3mm flux#

from astropy.table import Table 
import numpy as np

###By-eye catalog
cat = Table.read('/Users/josh/GitHub/W51/data/coldnh3_catalog.tex')
#Fluxes come from cat['peak']

#Catalog is already sorted by fluxes, low to high
fluxes = np.array(cat['peak'])
fluxes = fluxes[1:len(fluxes)]
fluxes = fluxes.astype('float64')

temps = np.array(cat['KTemp'])
temps = temps[1:len(temps)]
temps = temps.astype('float64')

masses = np.array(cat['mass'])
masses = masses[1:len(masses)]
masses = masses.astype('float64')

index = np.array(np.arange(len(fluxes)).astype('int'))

x = np.ndarray(shape=(4, int(len(index))))
x[0] = index
x[1] = fluxes
x[2] = temps
x[3] = masses


##Sorting by mass instead of flux
y = np.sort(x)

##Check to make sure sorting worked (mass arrays should be different - false)
comp = masses == y[3]
eq_arr = comp.all()
print(eq_arr) #Want false! That means arrays are not equal, and have been sorted

##See if indecies changes
comp = index == y[0]
eq_arr = comp.all()
print(eq_arr) #True means that the index does not change when sorted by flux or mass
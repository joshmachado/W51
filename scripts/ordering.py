#### Sorting sources by flux ####
#### Does order / index of sources change when sorting by 1.3mm flux vs NH3 derived mass? ####


#Sorting sources by 1.3mm flux#

from astropy.table import Table 
import numpy as np
import pylab as pl

###By-eye catalog
cat = Table.read('/Users/josh/GitHub/W51/data/byeye_catalog.tex')
#Fluxes come from cat['peak']

#Catalog is already sorted by fluxes, low to high

masses = np.array(cat['mass'])
masses = masses[1:len(masses)]
masses = masses.astype('float64')

fluxes = np.array(cat['peak'])
fluxes = fluxes[1:len(fluxes)]
fluxes = fluxes.astype('float64')
#TRIM ZERO MASSES
fluxes = fluxes[masses!=0]

temps = np.array(cat['KTemp'])
temps = temps[1:len(temps)]
temps = temps.astype('float64')
#TRIM ZERO MASSES
temps = temps[masses!=0]

#TRIM ZERO MASSES
masses = masses[masses!=0]

index = np.array(np.arange(len(fluxes)).astype('int'))

x = np.ndarray(shape=(4, int(len(index))))
x[0] = index
x[1] = fluxes
x[2] = temps
x[3] = masses

z = np.argsort(x)
# z[3] # this is the index order of masses low to high

##See if indecies changes
comp = index == z[3]
eq_arr = comp.all()
print(eq_arr) #True means that the index does not change when sorted by flux or mass
# FALSE. z[3] is the index order when sorted by mass

##Quantifying Mass Change##

reorder = masses[z[3]]
diff = masses - reorder
change = np.zeros(len(diff))
i = 0
while i < len(diff):
	if diff[i] < -0.1:
		change[i] = -1
	if diff[i] > 0.1:
		change[i] = 1
	i+=1

pl.hist(change, bins= [-1,0,1,2])
#pl.xaxis.set_visible(False)
pl.xlabel('Decrease | No Change | Increase')
pl.ylabel('Counts')
pl.title('Mass Change in By-Eye Catalog')
pl.show()


#### Sorting sources by flux ####
#### Does order / index of sources change when sorting by 1.3mm flux vs NH3 derived mass? ####


#Sorting sources by 1.3mm flux#

from astropy.table import Table 
import numpy as np
import pylab as pl

###By-eye catalog
cat = Table.read('/Users/josh/GitHub/W51/data/dendro_catalog.tex')
#Fluxes come from cat['peak']

#Catalog is already sorted by fluxes, low to high

masses = np.array(cat['mass'])
masses = masses[1:len(masses)]
masses = masses.astype('float64')
nans = np.isnan(masses)
fluxes = np.array(cat['peak_cont_flux'])
fluxes = fluxes[1:len(fluxes)]
fluxes = fluxes.astype('float64')
#TRIM ZERO MASSES
i = 0
while i < len(fluxes):
	if masses[i] == 0:
		fluxes[i] = 0
	if nans[i] == True:
		fluxes[i] = 0
	i+=1

fluxes = fluxes[fluxes!=0]

temps = np.array(cat['KTemps'])
temps = temps[1:len(temps)]
temps = temps.astype('float64')

i = 0
while i < len(temps):
	if temps[i] == 0:
		temps[i] = 0
	if nans[i] == True:
		temps[i] = 0
	i+=1
#TRIM ZERO MASSES
temps = temps[temps!=0]
#temps = temps[~np.isnan(temps)]

#TRIM ZERO MASSES
masses = masses[~np.isnan(masses)]
masses = masses[masses!=0]

index = np.array(np.arange(len(fluxes)).astype('int'))

x = np.ndarray(shape=(3, int(len(index))))
x[0] = index
x[1] = fluxes
x[2] = masses

z = np.argsort(x)
# z[3] # this is the index order of masses low to high

##See if indecies changes
comp = index == z[2]
eq_arr = comp.all()
print(eq_arr) #True means that the index does not change when sorted by flux or mass
# FALSE. z[3] is the index order when sorted by mass

##Quantifying Mass Change##

reorder = masses[z[2]]
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
pl.title('Mass Change in Dendrogram Catalog')
pl.show()

absdiff = [abs(value) for value in diff]
print('Average Difference: ' + str(np.mean(absdiff)))
print('Greatest Variation: ' + str(np.max(absdiff)))
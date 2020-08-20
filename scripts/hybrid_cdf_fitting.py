from astropy.table import Table
from astropy import units as u
import powerlaw
import pylab as pl
import numpy as np
from astropy import stats
import astropy
from scipy.optimize import curve_fit

###
#Identifying Catalogs
###

###By-eye catalog
cat = Table.read('/Users/josh/GitHub/W51/data/byeye_catalog.tex')
###Dendrogram catalog
dendrocat = Table.read('/Users/josh/GitHub/W51/data/dendro_catalog.tex')

### Using Three Mass Sets Here:


###### 1. Masses calculated with NH3 temp
#By-eye Ammmonia Temperature Derived Masses
nh3mass = np.asarray(cat['mass'])
nh3mass = nh3mass[1:len(cat['mass'])].astype(float)
nh3mass = np.sort(nh3mass[nh3mass != 0])
#Dendrogram Ammonia Temperature Derived Masses
dendronh3mass = np.asarray(dendrocat['mass'])
dendronh3mass = dendronh3mass[1:len(dendrocat['mass'])].astype(float)
dendronh3mass = dendronh3mass[~np.isnan(dendronh3mass)]
dendronh3mass = np.sort(dendronh3mass)

###### 2. Masses calculated with flat temp (mean)
#By-eye Flat Temperature Assumption (Mean Temp)
flat_mass = np.asarray(cat['mean nh3 mass'])
flat_mass = flat_mass[1:len(cat['mean nh3 mass'])].astype(float)
flat_mass = flat_mass[np.array(cat['mass'])[1:len(cat['mass'])].astype(float) != 0] #Removes all by-eye sources that don't have NH3 data
flat_mass = np.sort(np.trim_zeros(flat_mass))
#Dendrogram Flat Temperature Assumption (Mean Temp)
flat_dendro = np.asarray(dendrocat['mean nh3 mass'])
flat_dendro = flat_dendro[1:len(dendrocat['mean nh3 mass'])]
selection = np.logical_and(np.array(dendrocat['mass'])[1:len(dendrocat['mass'])].astype(float) != 0, np.array(dendrocat['mass'])[1:len(dendrocat['mass'])].astype(float) != np.nan)
flat_dendro = flat_dendro[selection] #Removes all dendrocat sources that don't have NH3 data
flat_mass = np.sort(np.trim_zeros(flat_mass))
flat_dendro = flat_dendro[~np.isnan(flat_dendro)].astype(float)
flat_dendro = np.sort(np.trim_zeros(flat_dendro))


###### 3. Hybrid, with M(T<80K) from flat temp (mean) and M(T>80K) from NH3 temp
#By-eye Hybrid Mass Set
nh3temp = np.array(cat['KTemp'])[1:len(cat['KTemp'])].astype(float)
nh3cores = np.array(cat['mass'][1:len(cat['mass'])].astype(float))
hottest_cores = nh3cores[nh3temp > 80]
flat_cores = np.array(cat['mean nh3 mass'][1:len(cat['mean nh3 mass'])].astype(float))
flat_cores = flat_cores[np.array(cat['mass'])[1:len(cat['mass'])].astype(float) != 0] #Removes all by-eye sources that don't have NH3 data (MASSES)
nh3temp = nh3temp[np.array(cat['mass'])[1:len(cat['mass'])].astype(float) != 0] #Removes all by-eye sources that don't have NH3 data (TEMPERATURES)
flat_cores = flat_cores[nh3temp < 80]
byeye_hybrid = np.append(flat_cores, hottest_cores)
#Dendrogram Hybrid Mass Set
dendro_nh3temp = np.array(dendrocat['KTemps'])[1:len(dendrocat['KTemps'])].astype(float)
dendro_nh3cores = np.array(dendrocat['mass'][1:len(dendrocat['mass'])].astype(float))
dendro_hottest_cores = dendro_nh3cores[dendro_nh3temp > 80]
dendro_flat_cores = np.array(dendrocat['mean nh3 mass'][1:len(dendrocat['mean nh3 mass'])].astype(float))
dendro_flat_cores = dendro_flat_cores[np.array(dendrocat['mass'])[1:len(dendrocat['mass'])].astype(float) != 0] #Removes all dendrocat sources that don't have NH3 data (MASSES)
dendro_nh3temp = dendro_nh3temp[np.array(dendrocat['mass'])[1:len(dendrocat['mass'])].astype(float) != 0] #Removes all dendrocat sources that don't have NH3 data (TEMPERATURES)
dendro_flat_cores = dendro_flat_cores[dendro_nh3temp < 80]
dendro_hybrid = np.append(dendro_flat_cores, dendro_hottest_cores)


###Powerlaw alpha values
fit1 = powerlaw.Fit(nh3mass)
fit2 = powerlaw.Fit(dendronh3mass)
fit3 = powerlaw.Fit(flat_mass)
fit4 = powerlaw.Fit(flat_dendro)
fit5 = powerlaw.Fit(byeye_hybrid)
fit6 = powerlaw.Fit(dendro_hybrid)

nh3_alpha_eye = fit1.power_law.alpha
nh3_alpha_dendro = fit2.power_law.alpha
flat_alpha_eye = fit3.power_law.alpha
flat_alpha_dendro = fit4.power_law.alpha
hybrid_alpha_eye = fit5.power_law.alpha
hybrid_alpha_dendro = fit6.power_law.alpha


###Plot CDF
#pl.plot(np.log10(np.sort(nh3mass)),np.log10(np.linspace(1,0,len(nh3mass), endpoint=False)), label=r'By-eye NH3 Temp., $\alpha$ = -{:0.4f}'.format(nh3_alpha_eye))
pl.plot(np.log10(np.sort(dendronh3mass)),np.log10(np.linspace(1,0,len(dendronh3mass), endpoint=False)), label=r'Dendrogram NH3 Temp., $\alpha$ = -{:0.4f}'.format(nh3_alpha_dendro))
#pl.plot(np.log10(np.sort(flat_mass)),np.log10(np.linspace(1,0,len(flat_mass), endpoint=False)), label=r'By-eye Flat Mean Temp., $\alpha$ = -{:0.4f}'.format(flat_alpha_eye))
pl.plot(np.log10(np.sort(flat_dendro)),np.log10(np.linspace(1,0,len(flat_dendro), endpoint=False)), label=r'Dendrogram Flat Mean Temp., $\alpha$ = -{:0.4f}'.format(flat_alpha_dendro))
#pl.plot(np.log10(np.sort(byeye_hybrid)),np.log10(np.linspace(1,0,len(byeye_hybrid), endpoint=False)), label=r'By-eye "Hybrid" Temp., $\alpha$ = -{:0.4f}'.format(hybrid_alpha_eye))
pl.plot(np.log10(np.sort(dendro_hybrid)),np.log10(np.linspace(1,0,len(dendro_hybrid), endpoint=False)), label=r'Dendrogram "Hybrid" Temp., $\alpha$ = -{:0.4f}'.format(hybrid_alpha_dendro))


pl.legend()
pl.title('CDF of Dendrogram for Mean Temp, NH3 Temp & "Hybrid"')
pl.show()



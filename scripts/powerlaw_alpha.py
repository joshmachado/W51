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
cat = Table.read('/Users/josh/GitHub/W51/data/coldnh3_catalog.tex')
###Dendrogram catalog
dendrocat = Table.read('/Users/josh/GitHub/W51/data/dendro_catalog.tex')

###
#Identifying datasets
###

###By-eye Flat Temperature Assumption (20K)
const_mass = cat['peak_mass']
const_mass = const_mass[1:len(cat['peak_mass'])]
const_mass = np.asarray(const_mass).astype(float)
const_mass = np.sort(const_mass)
const_mass = np.trim_zeros(const_mass)

###By-eye Ammmonia Temperature Derived Masses
nh3mass = cat['mass']
nh3mass = nh3mass[1:len(cat['mass'])]
nh3mass = np.asarray(nh3mass).astype(float)
nh3mass = np.sort(nh3mass)
nh3mass = np.trim_zeros(nh3mass)

###Dendrogram Flat Temperature Assumption (20K)
const_dendro = dendrocat['peak_cont_mass']
const_dendro = const_dendro[1:len(dendrocat['peak_cont_mass'])]
const_dendro = np.asarray(const_dendro).astype(float)
const_dendro = np.sort(const_dendro)
const_dendro = np.trim_zeros(const_dendro)
const_dendro = const_dendro[~np.isnan(const_dendro)]

###Dendrogram Ammonia Temperature Derived Masses
nh3dendro = dendrocat['mass']
nh3dendro = nh3dendro[1:len(dendrocat['mass'])]
nh3dendro = np.asarray(nh3dendro).astype(float)
nh3dendro = np.sort(nh3dendro)
nh3dendro = np.trim_zeros(nh3dendro)
nh3dendro = nh3dendro[~np.isnan(nh3dendro)]


fit1 = powerlaw.Fit(const_mass)
fit2 = powerlaw.Fit(nh3mass)
fit3 = powerlaw.Fit(const_dendro)
fit4 = powerlaw.Fit(nh3dendro)

alpha1 = fit1.power_law.alpha
alpha2 = fit2.power_law.alpha
alpha3 = fit3.power_law.alpha
alpha4 = fit4.power_law.alpha


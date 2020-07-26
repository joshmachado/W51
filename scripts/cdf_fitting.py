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



def linearpoly(x,a,b):
    return b + a*x

def fitalinetocdf(masses):
## fits a line to the loglog scale cdf of a function returns the alpha calculated
    popt,pcov = curve_fit(linearpoly, np.log10(np.sort(masses)),np.log10(np.linspace(1,0,len(masses),endpoint=False)))
    gamma = popt[0]
    alpha_est = -gamma +  1
    return  alpha_est    

flat_alpha = fitalinetocdf(const_mass)
nh3_eye_alpha = fitalinetocdf(nh3mass)
flat_dendro_alpha = fitalinetocdf(const_dendro)
nh3_dendro_alpha = fitalinetocdf(nh3dendro)

flat_gamma = -flat_alpha+1
nh3_eye_gamma = -nh3_eye_alpha+1
flat_dendro_gamma = -flat_dendro_alpha+1
nh3_dendro_gamma = -nh3_dendro_alpha+1


###Plot CDF
pl.plot(np.log10(np.sort(const_mass)),np.log10(np.linspace(1,0,len(const_mass), endpoint=False)), label=r'By-eye Flat Temp., $\alpha$ = -1.6432')
pl.plot(np.log10(np.sort(nh3mass)),np.log10(np.linspace(1,0,len(nh3mass), endpoint=False)), label=r'By-eye NH3 Temp., $\alpha$ = -1.6752')
pl.plot(np.log10(np.sort(const_dendro)),np.log10(np.linspace(1,0,len(const_dendro), endpoint=False)), label=r'Dendrogram Flat Temp., $\alpha$ = -1.9365')
pl.plot(np.log10(np.sort(nh3dendro)),np.log10(np.linspace(1,0,len(nh3dendro), endpoint=False)), label=r'Dendrogram NH3 Temp., $\alpha$ = -2.0763')


pl.legend()
pl.title('CDF of By-eye vs. Dendrogram for Flat & Dynamic Temp')
pl.show()
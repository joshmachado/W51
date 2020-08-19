### Reviewing how DENDRO CDF changes between default catalog and recalculated 20K catalog ###


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
cat = Table.read('/Users/josh/GitHub/W51/data/dendro_catalog.tex')

###
#Identifying datasets
###

###Dendro Ammmonia Temperature Derived Masses
nh3mass = cat['mass']
nh3mass = nh3mass[1:len(cat['mass'])]
nh3mass = np.asarray(nh3mass).astype(float)


###Dendro Flat Temperature Assumption (20K) DEFAULT
const_mass = cat['peak_cont_mass']
const_mass = const_mass[1:len(cat['peak_cont_mass'])]
const_mass = np.array(const_mass[nh3mass != 0]) ##Removes all by-eye sources that don't have NH3 data

const_mass = np.asarray(const_mass).astype(float)
const_mass = np.sort(const_mass)
const_mass = np.trim_zeros(const_mass)


###By-eye Flat Temperature Assumption (20K) RECALCULATED
re_mass = cat['20K nh3 mass']
re_mass = re_mass[1:len(cat['20K nh3 mass'])]
re_mass = re_mass[nh3mass != 0] ##Removes all by-eye sources that don't have NH3 data
re_mass = re_mass[re_mass != 'None']
re_mass = np.asarray(re_mass).astype(float)
re_mass = np.sort(re_mass)
re_mass = np.trim_zeros(re_mass)




def linearpoly(x,a,b):
    return b + a*x

def fitalinetocdf(masses):
## fits a line to the loglog scale cdf of a function returns the alpha calculated
    popt,pcov = curve_fit(linearpoly, np.log10(np.sort(masses)),np.log10(np.linspace(1,0,len(masses),endpoint=False)))
    gamma = popt[0]
    alpha_est = -gamma +  1
    return  alpha_est    

###Fit Line to CDF alpha & gamma values
flat_alpha = fitalinetocdf(const_mass)
recalc_alpha = fitalinetocdf(re_mass)

flat_gamma = -flat_alpha+1
recalc_gamma = -recalc_alpha+1


###Powerlaw alpha values
fit1 = powerlaw.Fit(const_mass)
fit2 = powerlaw.Fit(re_mass)

pl_flat_alpha = fit1.power_law.alpha
pl_recalc_alpha = fit2.power_law.alpha

###Plot CDF
pl.plot(np.log10(np.sort(const_mass)),np.log10(np.linspace(1,0,len(const_mass), endpoint=False)), label=r'By-eye Flat Temp., $\alpha$ = -{:0.4f}'.format(pl_flat_alpha))
pl.plot(np.log10(np.sort(re_mass)),np.log10(np.linspace(1,0,len(re_mass), endpoint=False)), label=r'By-eye Flat Temp. RECALCULATED., $\alpha$ = -{:0.4f}'.format(pl_recalc_alpha))

#Fit lines to cdf
#intercept = 0.34281145
#x = np.logspace(-1.5,.5, base=10)
#gamma2 = -1.7934017965102476+1



#pl.plot(x,flat_gamma*x+intercept, label=r'$\alpha$ = -1.7131')

#pl.plot(x,gamma2*x+intercept+0.125, label=r'Bootstrapped $\alpha$ = -1.7934017965102476 ')
#pl.plot(x, (intercept)*x**(-gamma))

pl.legend()
pl.title('CDF of Default Flat vs Recalc 20K Temp - Dendrogram')
pl.show()
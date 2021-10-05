##############
## imports ##
##############
from astropy.table import Table
from astropy import units as u
from astropy.io import fits
from astropy import constants as const
import numpy as np
import powerlaw
import pylab as plt
from astropy import stats
import astropy
from scipy.optimize import curve_fit
##############################
## file paths & definitions ##
##############################
fp = '/Users/josh/GitHub/W51/'
catalog = Table.read(fp+'data/byeye_catalog.tex')
nh3_cat = np.array(catalog['KTemp'][1:])
dendrocat = Table.read(fp+'data/dendro_catalog.tex')
nh3_dendro = np.array(dendrocat['KTemps'][1:])
c = const.c #Speed of light m/s
c = c.to(u.cm/u.s)
k = const.k_B #Botlzmann constant J/K
beta = 1.75
kgas = 0.0114 * (u.cm)**2/u.g #Kappa
nu = 226*u.GHz
nu0 = 271.1*u.GHz
k_nu = kgas * (nu/nu0)**beta
dist = 5.41*u.kpc
cat_mass = [None] * len(catalog)
dendro_mass = [None] * len(dendrocat)
cutoff = 40
##########################
## grab cores above 30K ##
##########################
warm_cat = np.array(catalog['KTemp'][1:])
warm_cat = np.array([float(i) for i in warm_cat])
warm_cat = warm_cat[warm_cat > cutoff]
warm_dendro = dendrocat['KTemps'][1:]
warm_dendro = np.array(warm_dendro[warm_dendro > cutoff])

#######################################
## determine mean temp of 30K+ cores ##
#######################################
cat_temp = np.mean(warm_cat)
dendro_temp = np.mean(warm_dendro)

#######################################################
## use new mean to compute new masses for 30K+ cores ##
#######################################################
for i in range(len(nh3_cat)):
    if float(nh3_cat[i]) > cutoff:
        temp = cat_temp * u.K
    else:
        temp = 20 * u.K
    intensity = catalog['peak'][i+1] * u.Jy
    intensity = intensity.to(u.erg/((u.cm)**2))
    intensity = intensity * u.erg / (u.cm)**2
    cat_mass[i] = ((intensity * (dist.to(u.cm))**2 * c**2) / (2* k_nu * nu **2 * k.cgs * temp)).to(u.Msun).value

for i in range(len(nh3_dendro)):
    if float(nh3_dendro[i]) > cutoff:
        temp = dendro_temp * u.K
    else:
        temp = 20 * u.K
    intensity = dendrocat['peak_cont_flux'][i+1] * u.Jy
    intensity = intensity.to(u.erg/((u.cm)**2))
    dendro_mass[i] = ((intensity * (dist.to(u.cm))**2 * c**2) / (2* k_nu * nu **2 * k.cgs * temp)).to(u.Msun).value

##################################################################
## combine sub 30K and 30K+ flat assumptions to recalculate CDF ##
##################################################################
def linearpoly(x,a,b):
    return b + a*x

def fitalinetocdf(masses):
## fits a line to the loglog scale cdf of a function returns the alpha calculated
    popt,pcov = curve_fit(linearpoly, np.log10(np.sort(masses)),np.log10(np.linspace(1,0,len(masses),endpoint=False)))
    gamma = popt[0]
    alpha_est = -gamma +  1
    return  alpha_est


cat_mass = cat_mass[:-1]
dendro_mass = dendro_mass[:-1]
###By-eye Masses
nh3mass = catalog['mass']
nh3mass = nh3mass[1:len(catalog['mass'])]
nh3mass = np.asarray(nh3mass).astype(float)
const_mass = catalog['20K nh3 mass']
const_mass = const_mass[1:len(catalog['20K nh3 mass'])]
const_mass = np.asarray(const_mass).astype(float)
const_mass = const_mass[nh3mass != 0] ##Removes all by-eye sources that don't have NH3 data
const_mass = np.sort(const_mass)
const_mass = np.trim_zeros(const_mass)
nh3mass = np.sort(nh3mass)
nh3mass = np.trim_zeros(nh3mass)


###Dendrogram Masses
nh3dendro = dendrocat['mass']
nh3dendro = nh3dendro[1:len(dendrocat['mass'])]
nh3dendro = np.asarray(nh3dendro).astype(float)
const_dendro = dendrocat['20K nh3 mass']
const_dendro = const_dendro[1:len(dendrocat['20K nh3 mass'])]
const_dendro = const_dendro[nh3dendro != 0] ##Removes all dendrocat sources that don't have NH3 data
const_dendro = const_dendro[const_dendro != 'None']
const_dendro = np.asarray(const_dendro).astype(float)
const_dendro = np.sort(const_dendro)
const_dendro = np.trim_zeros(const_dendro)
const_dendro = const_dendro[~np.isnan(const_dendro)]
nh3dendro = np.sort(nh3dendro)
nh3dendro = np.trim_zeros(nh3dendro)
nh3dendro = nh3dendro[~np.isnan(nh3dendro)]

###Fit Line to CDF alpha & gamma values
flat_alpha = fitalinetocdf(const_mass)
nh3_eye_alpha = fitalinetocdf(nh3mass)
flat_dendro_alpha = fitalinetocdf(const_dendro)
nh3_dendro_alpha = fitalinetocdf(nh3dendro)
two_tier_byeye_alpha = fitalinetocdf(cat_mass)
two_tier_dendro_alpha = fitalinetocdf(dendro_mass)

flat_gamma = -flat_alpha+1
nh3_eye_gamma = -nh3_eye_alpha+1
flat_dendro_gamma = -flat_dendro_alpha+1
nh3_dendro_gamma = -nh3_dendro_alpha+1
two_tier_byeye_gamma = -two_tier_byeye_alpha+1
two_tier_dendro_gamma = -two_tier_dendro_alpha+1

###Powerlaw alpha values
fit1 = powerlaw.Fit(const_mass)
fit2 = powerlaw.Fit(nh3mass)
fit3 = powerlaw.Fit(const_dendro)
fit4 = powerlaw.Fit(nh3dendro)
fit5 = powerlaw.Fit(cat_mass)
fit6 = powerlaw.Fit(dendro_mass)

pl_flat_alpha = fit1.power_law.alpha
pl_nh3_alpha = fit2.power_law.alpha
pl_flat_den_alpha = fit3.power_law.alpha
pl_nh3_den_alpha = fit4.power_law.alpha
pl_tt_eye_alpha = fit5.power_law.alpha
pl_tt_den_alpha = fit6.power_law.alpha

###Plot CDF
plt.plot(np.log10(np.sort(const_mass)),np.log10(np.linspace(1,0,len(const_mass), endpoint=False)), label=r'By-eye Flat Temp., $\alpha$ = -{:0.4f}'.format(pl_flat_alpha))
plt.plot(np.log10(np.sort(nh3mass)),np.log10(np.linspace(1,0,len(nh3mass), endpoint=False)), label=r'By-eye NH3 Temp., $\alpha$ = -{:0.4f}'.format(pl_nh3_alpha))
plt.plot(np.log10(np.sort(const_dendro)),np.log10(np.linspace(1,0,len(const_dendro), endpoint=False)), label=r'Dendrogram Flat Temp., $\alpha$ = -{:0.4f}'.format(pl_flat_den_alpha))
plt.plot(np.log10(np.sort(nh3dendro)),np.log10(np.linspace(1,0,len(nh3dendro), endpoint=False)), label=r'Dendrogram NH3 Temp., $\alpha$ = -{:0.4f}'.format(pl_nh3_den_alpha))
plt.plot(np.log10(np.sort(cat_mass)),np.log10(np.linspace(1,0,len(cat_mass), endpoint=False)), label=r'By-Eye Two-Tier Flat Temp., $\alpha$ = -{:0.4f}'.format(pl_tt_eye_alpha))
plt.plot(np.log10(np.sort(dendro_mass)),np.log10(np.linspace(1,0,len(dendro_mass), endpoint=False)), label=r'Dendrogram Two-Tier Flat Temp., $\alpha$ = -{:0.4f}'.format(pl_tt_den_alpha))


plt.legend()
plt.title('CDF of By-eye vs. Dendrogram for Flat 20K, Dynamic Temp & Two-Tier Flat')
plt.show()

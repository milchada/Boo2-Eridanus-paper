#calculate velocity dispersions using method of Walker+ 2006
import numpy as np
from astropy.io import fits
import scipy.optimize as opt 
import matplotlib.pylab as plt

keckfile_boo2 = 'alldata_Boo2.fits'#str(input('Keck file name: '))
#read keck file
boo2_fits = fits.open(keckfile_boo2)[1]
mbr_boo2 = boo2_fits.data['MEMBER'] 
v_boo2 = boo2_fits.data['VCORR'][mbr_boo2==1][1:] #removes spec binary
sigma_boo2 = boo2_fits.data['MCERR_VCORR'][mbr_boo2==1][1:]

keckfile_eri = 'alldata_Eri.fits'

eri_fits = fits.open(keckfile_eri)[1]
mbr_eri = eri_fits.data['MEMBER']
v_eri = eri_fits.data['VCORR'][mbr_eri==1]
sigma_eri = eri_fits.data['MCERR_VCORR'][mbr_eri==1]

methods = ['Nelder-Mead', 'Powell', 'CG', 'BFGS', 'L-BFGS-B', 'TNC', 'COBYLA', 'SLSQP'] 

def log_p(v_p_s_p,v_obs,s_obs):
	v_p = v_p_s_p[0]
	sigma_p = np.exp(v_p_s_p[1])
	N=len(v_obs)
	return -1*(-1./2*sum(np.log(s_obs**2+sigma_p**2))-1./2*sum((v_obs-v_p)**2/(s_obs**2+sigma_p**2))-N/2.*np.log(2*np.pi))
# log probability is negative, we want to minimize magnitude, #hence multiplied by -1.
	
def dispersion_mean(v_obs, s_obs,opt_method='Nelder-Mead'):	
	#input guess for mean velocity and dispersion 
	first_guess_vp = np.mean(v_obs)
	first_guess_sp = np.log(np.std(v_obs))

	result = opt.minimize(log_p, x0=np.array([first_guess_vp, first_guess_sp]), args = (v_obs, s_obs), method=opt_method)
	return result['x'] # #prints values of v_peak and sigma that minimize log(p)

#now see how dispersion changes on removing any individual stars 
def excluding_stars(method, v_obs = v_boo2, s_obs = sigma_boo2):
	vp_disp_list = []
	for star in xrange(len(v_obs)):
		vp_disp_list.append(dispersion_mean(np.delete(v_obs, star),np.delete(s_obs, star),method))

	return np.array(vp_disp_list)

def exclusion_results(v_obs = v_boo2, s_obs = sigma_boo2):
	results = {}
	for method in methods:
		results[method] = np.append(excluding_stars(method), dispersion_mean(v_obs, s_obs, method))
	return results

#calculate errors using the elements of covariance matrix, eq. 10
#check if star that causes largest change in sigma is the binary mentioned in the paper - SDSS name is ra, dec
###checked - it is! yay!
def covariance_matrix(v_obs = v_boo2, s_obs = sigma_boo2):
	vmax, log_smax = dispersion_mean(v_obs, s_obs)
	smax = np.exp(log_smax)
	o2 = s_obs**2+smax**2
	term1 = 1*sum(1./o2)
	term2 = sum(-1./o2 + (2*(smax**2)-(v_obs-vmax)**2)/(o2**2)+4*(v_obs-vmax)**2*(smax**2)/(o2**3))
	term3 = 2*sum((vmax - v_obs)*smax/(o2**2))
	
	detA = term1*term2-(term3**2)
	return np.sqrt(abs(term1/detA)), np.sqrt(term2/detA)

#how to use MCMC to calculate errorbars, relaxing condition of Gaussianity?
import emcee
import corner

"""
Note: First I sampled the mean and std, but this is wrong.
We should sample individual velocities. That takes into account all measurement errors.
So we have as many dimensions as there are stars
Then likelihood (v_array, v_p_s_p) = log_p(v_array, s_obs, )
"""

def mcmc_allstars(v_obs, s_obs, filename=None, nwalkers = 100, ntrials = 1000, method='Nelder-Mead'):
	vmax, log_smax = dispersion_mean(v_obs, s_obs)
	smax = np.exp(log_smax)
	def lnprior_all(vp):
		if sum(vp[i]>(vmax-5*smax) for i in xrange(len(vp))) and sum(vp[i]<(vmax+5*smax) for i in xrange(len(vp))):
			return 0.0
		else:
			return -np.inf

	def lnprob_all(v_sample):
		lp = lnprior_all(v_sample)
		if not np.isfinite(lp):
			return -np.inf
		return lp - log_p([vmax, log_smax], v_sample, s_obs) #subtracted because recall log_p above is NEGATIVE of log of likelihood
        #s_obs is measured errors - not sampling anything here

	ndim = len(v_obs)
	pos = [v_obs + 0.1*np.random.random(ndim) for i in range(int(nwalkers))]
	sampler = emcee.EnsembleSampler(int(nwalkers), ndim, lnprob_all)
	sampler.run_mcmc(pos, int(ntrials))
	samples = sampler.chain[:, 50:, :].reshape((-1, ndim)) #removes first 50 TRIALS
	mean_logstd = [dispersion_mean(samples[row], s_obs, opt_method = method) for row in xrange(len(samples))]
	mean_std = [np.array([mean_logstd[row][0],np.exp(mean_logstd[row][1])]) for row in xrange(len(mean_logstd))]
	if filename != None:
		fig = corner.corner(mean_std, labels=["$v_p$", "$s_p$"], truths=[vmax, smax], truth_color='r', range=((vmax-4,vmax+4),(0,max(3*smax,3))))
		fig.savefig(filename)
	return samples, np.array(mean_std)

def non_gaussian_errors(sample, binsize, confidence = 0.68):
	hist, bins = np.histogram(sample, range=(min(sample),max(sample)), bins = (max(sample)-min(sample))/binsize)
	c=np.cumsum(hist)
	s=sum(hist)
	minbin = s*(1-confidence)/2.
	maxbin = s*(1+confidence)/2.
	return min(bins[c>minbin]), max(bins[c<maxbin])

def print_kinematics(v_obs, s_obs, filename=None,binsize=0.01, confidence = 0.68, nwalkers=1e2, ntrials=1e2, method = 'Nelder-Mead'):
    dict = {}
    s, ms = mcmc_allstars(v_obs, s_obs, filename, nwalkers, ntrials, method)
    maxlike = dispersion_mean(v_obs, s_obs)
    cov_err = covariance_matrix(v_obs, s_obs)
    minmax_v = non_gaussian_errors(ms[:,0], binsize, confidence)
    minmax_s = non_gaussian_errors(ms[:,1], binsize, confidence)
    v_hist, v_bins = np.histogram(ms[:,0], bins = (max(ms[:,0])-min(ms[:,0]))/binsize)
    s_hist, s_bins = np.histogram(ms[:,1], bins = (max(ms[:,1])-min(ms[:,1]))/binsize)
    v_max = v_bins[v_hist==max(v_hist)]
    s_max = s_bins[s_hist==max(s_hist)]
    mcmc_sigma = (s_max, minmax_s[1]-s_max, s_max-minmax_s[0])
    mcmc_v = (v_max, minmax_v[1]-v_max, v_max-minmax_v[0])
    ml_sigma = (np.exp(maxlike[1]), cov_err[1], cov_err[1])
    ml_v = (maxlike[0], cov_err[0], cov_err[0])
    print 'ML: mean velocity = %0.02f +%0.02f -%0.02f ' % ml_v
    print 'ML: dispersion = %0.02f +%0.02f -%0.02f ' % ml_sigma
    print 'MCMC: mean velocity = %0.02f +%0.02f -%0.02f ' % mcmc_v
    print 'MCMC: dispersion = %0.02f +%0.02f -%0.02f' % mcmc_sigma
    dict['ml_v']=ml_v
    dict['ml_s']=ml_sigma
    dict['mcmc_v']=mcmc_v
    dict['mcmc_s']=mcmc_sigma
    return dict
	
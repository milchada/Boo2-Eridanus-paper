import pandas as pd
import emcee
import numpy as np


df=pd.read_csv("Peg3_velocities.csv")

def modelemcee(x, v, v_err, N):
	logL= -((N*np.log(2*np.pi)/2)+0.5*(np.log(v_err**2.+(x[1])**2.).sum())+0.5*((v-x[0])**2./(v_err**2.+(x[1])**2.)).sum())
	if x[0] < -250 or x[0] > -200:
		return -np.inf
	if x[1] > 30 or x[1] < 0:
		return -np.inf
	return logL


import scipy.optimize as op
nll = lambda *args: -modelemcee(*args)
result = op.minimize(nll, [np.median(df.v),np.std(df.v)], args=(df.v, df.v_err, len(df)))
print result #Maximum Likelihood
#ok so his import velocity dispersion is the mean of the observed uncertainties, whereas I use the std of the measured velocities. editing..

ndim, nwalkers = int(1e2), int(1e2) #he had 200, i.e. not many realizations.. tweaking
pos = np.array([result["x"] + 1e-4*np.random.randn(ndim) for i in range(nwalkers)])
#ah ha! he is also doing the thing where the only dimensions being sampled are v_mean and std.
sampler = emcee.EnsembleSampler(nwalkers, ndim, modelemcee, args=(df.v, df.v_err, len(df)), threads = 4)
sampler.run_mcmc(pos,10000)


import corner
samples = sampler.chain[50:, 50:].reshape((-1, ndim))
print np.percentile(samples,[16,50,84],axis=0) # 16, 50, 84 percentiles
fig = corner.corner(samples,bins=200,labels=[ '$<V>\/(km/s)$', '$\sigma_v\/(km/s)$'],quantiles=[0.16, 0.5, 0.84],max_n_ticks=6)
fig.savefig('test.png')


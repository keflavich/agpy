import numpy as np
old_errsettings = np.geterr()
import pymc # pymc breaks np error settings
np.seterr(**old_errsettings)
import itertools
from . import blackbody


def fit_blackbody_montecarlo(frequency, seds, errors=None,
        temperature_guess=10, beta_guess=None, scale_guess=None,
        blackbody_function=blackbody, quiet=True, return_MC=True,
        nsamples=5000, burn=1000, min_temperature=0, max_temperature=100,
        scale_keyword='scale', max_scale=1e60,
        multivariate=False,
        **kwargs):
    """
    Parameters
    ----------
    frequency : array
        Array of frequency values
    flux : array
        array of flux values
    err : array (optional)
        Array of error values (1-sigma, normal)
    temperature_guess : float
        Input / starting point for temperature
    min_temperature : float
    max_temperature : float
        Lower/Upper limits on fitted temperature
    beta_guess : float (optional)
        Opacity beta value
    scale_guess : float
        Arbitrary scale value to apply to model to get correct answer
    blackbody_function: function
        Must take x-axis (e.g. frequency), temperature, then scale and beta
        keywords (dependence on beta can be none)
    return_MC : bool
        Return the pymc.MCMC object?
    nsamples : int
        Number of samples to use in determining the posterior distribution
        (the answer)
    burn : int
        number of initial samples to ignore
    scale_keyword : ['scale','logscale','logN']
        What scale keyword to pass to the blackbody function to determine
        the amplitude
    kwargs : kwargs
        passed to blackbody function
    """

    d = {}

    d['temperature'] = pymc.distributions.Uniform('temperature',
            min_temperature, max_temperature, value=temperature_guess)
    d['scale'] = pymc.distributions.Uniform('scale',0,max_scale,
            value=scale_guess)
    if beta_guess is not None:
        d['beta'] = pymc.distributions.Uniform('beta',0,10,
                value=beta_guess)
    else:
        d['beta'] = pymc.distributions.Uniform('beta',0,10,
                value=1)

    covar_list = dict([((i,j),pymc.Uninformative('%s-%s' % (i,j),value=(i==j))) for i,j in
        itertools.combinations_with_replacement(('t','b','s'),2)])
    for i,j in itertools.permutations(('t','b','s'),2):
        if (i,j) in covar_list:
            covar_list[(j,i)] = covar_list[(i,j)]
    covar_grid = [[covar_list[(i,j)] for i in ('t','b','s')] for j in ('t','b','s')]
    d['tbcov'] = pymc.MvNormalCov('tbcov', mu=[d['temperature'],d['beta'],d['scale']],
            C=covar_grid, value=[d['temperature'],d['beta'],d['scale']])

    precision_list = dict([((i,j),pymc.Uninformative('%s-%s' % (i,j),value=(i==j))) for i,j in
        itertools.combinations_with_replacement(('t','b','s'),2)])
    for i,j in itertools.permutations(('t','b','s'),2):
        if (i,j) in precision_list:
            precision_list[(j,i)] = precision_list[(i,j)]
    precision_grid = [[precision_list[(i,j)] for i in ('t','b','s')] for j in ('t','b','s')]
    # need to force tau > 0...
    d['tbprec'] = pymc.MvNormalCov('tbprec', mu=[d['temperature'],d['beta'],d['scale']],
            C=precision_grid, value=[1,1,1])

    for ii,(sed,err) in enumerate(zip(seds,errors)):
        d['t_%i' % ii] = pymc.Normal('t_%i' % ii, mu=d['tbcov'][0],
                tau=d['tbprec'][0])
        d['b_%i' % ii] = pymc.Normal('b_%i' % ii, mu=d['tbcov'][1],
                tau=d['tbprec'][1])
        d['s_%i' % ii] = pymc.Normal('s_%i' % ii, mu=d['tbcov'][2],
                tau=d['tbprec'][2])


        def bb_model(temperature=d['t_%i' % ii],
                scale=d['s_%i' % ii],
                beta=d['b_%i' % ii]):
            kwargs[scale_keyword] = scale
            y = blackbody_function(frequency, temperature, 
                    beta=beta, normalize=False, **kwargs)
            #print kwargs,beta,temperature,(-((y-flux)**2)).sum()
            return y

        d['bb_model_%i' % ii] = pymc.Deterministic(
                eval=bb_model,
                name='bb_model_%i' % ii,
                parents={'temperature': d['t_%i' % ii], 
                    'scale':d['s_%i' % ii],
                    'beta':d['b_%i' % ii]},
                doc='Blackbody SED model.',
                trace=True,
                verbose=0,
                dtype=float,
                plot=False,
                cache_depth=2)
                

        if err is None:
            d['err_%i' % ii] = pymc.distributions.Uninformative('error_%i' %
                    ii, value=1.)
        else:
            d['err_%i' % ii] = pymc.distributions.Uninformative('error_%i' %
                    ii, value=err, observed=True)

        d['flux_%i' % ii] = pymc.distributions.Normal('flux_%i' % ii,
                mu=d['bb_model_%i' % ii],  tau=1./d['err_%i' % ii]**2,
                value=sed, observed=True)

    #print d.keys()
    MC = pymc.MCMC(d)
    
    if nsamples > 0:
        MC.sample(nsamples, burn=burn)
        if return_MC:
            return MC

        MCfit = pymc.MAP(MC)
        MCfit.fit()
        T = MCfit.temperature.value
        scale = MCfit.scale.value

        if beta_guess is not None:
            beta = MCfit.beta.value
            return T,scale,beta
        else:
            return T,scale

    return MC

import numpy as np
from math import gamma

def Bisect(f, a, b, args):

    '''
    Replace scipy.optimize.bisect.
    '''

    for i in range(54):

        m = (a + b) / 2.0

        if f(a, args) * f(m, args) > 0:

            a = m

        else:

            b = m

    return m

def NegativeBinomial(k, p, n):

    '''
    Replace scipy.stats.nbinom.pmf.
    '''

    def binom(n, k):

        return 1 if k == 0 or k == n else np.prod([(n + 1 - i) / i for i in range(1, k + 1)])

    return binom(k + n - 1, n - 1) * ((1 - p) ** k) * (p ** n)

def epgPerPerson(x, params):

    '''
    This function calculates the total eggs per gram as
    a function of the mean worm burden.

    Parameters
    ----------
    x: float
        array of mean worm burdens;

    params: dict
        dictionary containing the parameter names and values;
    '''

    return params['lambda'] * x * params['z'] * (1 + x * (1 - params['z']) / params['k']) ** (- params['k'] - 1)

def fertilityFunc(x, params):

    '''
    This function calculates the multiplicative fertility correction factor
    to be applied to the mean eggs per person function.

    Parameters
    ----------
    x: float
        array of mean worm burdens;

    params: dict
        dictionary containing the parameter names and values;
    '''

    a = 1 + x * (1 - params['z']) / params['k']
    b = 1 + 2 * x / params['k'] - params['z'] * x / params['k']

    return 1 - (a / b) ** (params['k'] + 1)

def monogFertilityConfig(params, N=30):

    '''
    This function calculates the monogamous fertility
    function parameters.

    Parameters
    ----------
    params: dict
        dictionary containing the parameter names and values;

    N: int
        resolution for the numerical integration
    '''

    return dict(c_k=gamma(params['k'] + 0.5) * (2 * params['k'] / np.pi) ** 0.5 / gamma(params['k'] + 1),
    cosTheta=np.cos(np.linspace(start=0, stop=2 * np.pi, num=N + 1)[:N]))

def monogFertilityFuncApprox(x, params):

    '''
    This function calculates the fertility factor for monogamously mating worms.

    Parameters
    ----------

    x: float
        mean worm burden;

    params: dict
        dictionary containing the parameter names and values;
    '''

    if x > 25 * params['k']:

        return 1 - params['monogParams']['c_k'] / np.sqrt(x)

    else:

        g = x / (x + params['k'])
        integrand = (1 - params['monogParams']['cosTheta']) * (1 + g * params['monogParams']['cosTheta']) ** (- 1 - params['k'])
        integral = np.mean(integrand)

        return 1 - (1 - g) ** (1 + params['k']) * integral

def epgMonog(x, params):

    '''
    This function calculates the generation of eggs with monogamous
    reproduction taken into account.

    Parameters
    ----------
    x: float
        array of mean worm burdens;

    params: dict
        dictionary containing the parameter names and values;
    '''

    return epgPerPerson(x, params) * np.vectorize(monogFertilityFuncApprox)(x, params)

def epgFertility(x, params):

    '''
    This function calculates the generation of eggs with
    sexual reproduction taken into account.

    Parameters
    ----------
    x: float
        array of mean worm burdens;

    params: dict
        dictionary containing the parameter names and values;
    '''

    return epgPerPerson(x, params) * fertilityFunc(x, params)

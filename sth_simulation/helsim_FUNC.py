import numpy as np
import pandas as pd
import copy
import pkg_resources

import sth_simulation.ParallelFuncs as ParallelFuncs

def setIntensityCount(paramFileName):

    '''
    This function defines the intensity thresholds and converts them to egg counts.
    Note that 24 is a conversion factor from epg (eggs per gram faeces) to egg counts
    on a microscopy slide of a faecal smear test.

    Ascaris: moderate intensity: 5,000, high intensity: 50,000.
    Trichuris: moderate intensity: 1,000, high intensity: 10,000.
    Hookworm: moderate intensity: 2,000, high intensity: 4,000.

    Parameters
    ----------
    paramFileName: str
        name of the input text file;

    Returns
    -------
    moderateIntensityCount: float
        moderate intensity threshold;

    highIntensityCount: float
        high intensity threshold;
    '''

    if paramFileName == 'AscarisParameters_moderate.txt' or paramFileName == 'AscarisParameters_high.txt':

        moderateIntensityCount = 5000 / 24
        highIntensityCount = 50000 / 24

    elif paramFileName == 'TrichurisParameters_moderate.txt' or paramFileName == 'TrichurisParameters_high.txt':

        moderateIntensityCount = 1000 / 24
        highIntensityCount = 10000 / 24

    elif paramFileName == 'HookwormParameters_moderate.txt' or paramFileName == 'HookwormParameters_high.txt':

        moderateIntensityCount = 2000 / 24
        highIntensityCount = 4000 / 24

    elif paramFileName == 'SCH_MansoniParameters.txt' or paramFileName == 'SCH_HaematobiumParameters.txt':

        moderateIntensityCount = 100 / 24
        highIntensityCount = 400 / 24

    return moderateIntensityCount, highIntensityCount

def readParam(fileName):

    '''
    This function extracts the parameter values stored
    in the input text files into a dictionary.

    Parameters
    ----------
    fileName: str
        name of the input text file;

    Returns
    -------
    params: dict
        dictionary containing the parameter names and values;
    '''

    DATA_PATH = pkg_resources.resource_filename('sth_simulation', 'data/')

    with open(DATA_PATH + fileName) as f:

        contents = f.readlines()

    params = []

    for content in contents:

        line = content.split('\t')

        if len(line) >= 2:

            try:

                line[1] = np.array([np.float(x) for x in line[1].split(' ')])

                if len(line[1]) == 1:

                    line[1] = line[1][0]

            except:

                pass

            params.append(line[:2])

    params = dict(params)

    return params

def readParams(paramFileName, demogFileName='Demographies.txt', demogName='Default'):

    '''
    This function organizes the model parameters and
    the demography parameters into a unique dictionary.

    Parameters
    ----------
    paramFileName: str
        name of the input text file with the model parameters;

    demogFileName: str
        name of the input text file with the demography parameters;

    demogName: str
        subset of demography parameters to be extracted;

    Returns
    -------
    params: dict
        dictionary containing the parameter names and values;
    '''

    demographies = readParam(demogFileName)
    parameters = readParam(paramFileName)

    params = {'numReps': np.int(parameters['repNum']),
              'maxTime': parameters['nYears'],
              'N': np.int(parameters['nHosts']),
              'R0': parameters['R0'],
              'lambda': parameters['lambda'],
              'gamma': parameters['gamma'],
              'k': parameters['k'],
              'sigma': parameters['sigma'],
              'LDecayRate': parameters['ReservoirDecayRate'],
              'DrugEfficacy': parameters['drugEff'],
              'contactAgeBreaks': parameters['contactAgeBreaks'],
              'treatmentAgeBreaks': parameters['treatmentBreaks'],
              'contactRates': parameters['betaValues'],
              'rho': parameters['rhoValues'],
              'outputFrequency': parameters['outputFrequency'],
              'outputOffset': parameters['outputOffset'],
              'highBurdenBreaks': parameters['highBurdenBreaks'],
              'highBurdenValues': parameters['highBurdenValues'],
              'k_epg': parameters['k_epg'],
              'nNodes': parameters['nNodes'],
              'maxStep': parameters['maxStep'],
              'seed': parameters['seed'],
              'demogType': demogName,
              'hostMuData': demographies[demogName + '_hostMuData'],
              'muBreaks': np.append(0, demographies[demogName + '_upperBoundData']),
              'SR': [True if parameters['StochSR'] == 'TRUE' else False][0],
              'reproFuncName': parameters['reproFuncName'],
              'unfertilized': parameters['unfertilized'] == "TRUE",
              'z': np.exp(-parameters['gamma']),
              'psi': 1.0}

    return params

def configure(params):

    '''
    This function defines a number of additional parameters.

    Parameters
    ----------
    params: dict
        dictionary containing the initial parameter names and values;

    Returns
    -------
    params: dict
        dictionary containing the updated parameter names and values;
    '''

    # level of discretization for the drawing of lifespans
    dT = 0.1

    # definition of the reproduction function
    params['reproFunc'] = getattr(ParallelFuncs, params['reproFuncName'])

    # max age cutoff point
    params['maxHostAge'] = np.min([np.max(params['muBreaks']), np.max(params['contactAgeBreaks'])])

    # full range of ages
    params['muAges'] = np.arange(start=0, stop=np.max(params['muBreaks']), step=dT) + 0.5 * dT

    params['hostMu'] = params['hostMuData'][pd.cut(x=params['muAges'], bins=params['muBreaks'],
    labels=np.arange(start=0, stop=len(params['hostMuData']))).to_numpy()]

    # probability of surviving
    params['hostSurvivalCurve'] = np.exp(-np.cumsum(params['hostMu']) * dT)

    # the index for the last age group before the cutoff in this discretization.
    maxAgeIndex = np.argmax([params['muAges'] > params['maxHostAge']]) - 1

    # cumulative probability of dying
    params['hostAgeCumulDistr'] = np.append(np.cumsum(dT * params['hostMu'] * np.append(1, params['hostSurvivalCurve'][: -1]))[:maxAgeIndex], 1)
    params['contactAgeGroupBreaks'] = np.append(params['contactAgeBreaks'][:-1], params['maxHostAge'])
    params['treatmentAgeGroupBreaks'] = np.append(params['treatmentAgeBreaks'][:-1], params['maxHostAge'] + dT)

    # monogamous fertility function parameters
    if params['reproFuncName'] == 'epgMonog':
        params['monogParams'] = ParallelFuncs.monogFertilityConfig(params)

    # time points included in simulation output
    params['outTimings'] = np.arange(start=params['outputOffset'], stop=params['maxTime'],
    step=1 / params['outputFrequency'])

    # make sure that the last time point is always included
    if params['outTimings'][-1] != params['maxTime']:
        params['outTimings'] = np.append(params['outTimings'], params['maxTime'])

    # if there are any time points smaller than zero, then reset those to zero
    params['outTimings'][params['outTimings'] < 0] = 0

    # make sure that there are no duplicate time points
    params['outTimings'] = np.unique(np.round(params['outTimings'], 4))

    return params

def setupSD(params):

    '''
    This function sets up the simulation to initial conditions
    based on analytical equilibria.

    Parameters
    ----------
    params: dict
        dictionary containing the parameter names and values;

    Returns
    -------
    SD: dict
        dictionary containing the equilibrium parameter settings;
    '''

    si = np.random.gamma(size=params['N'], scale=1 / params['k'], shape=params['k'])

    lifeSpans = getLifeSpans(params['N'], params)
    trialBirthDates = - lifeSpans * np.random.uniform(low=0, high=1, size=params['N'])
    trialDeathDates = trialBirthDates + lifeSpans

    communityBurnIn = 1000

    while np.min(trialDeathDates) < communityBurnIn:

        earlyDeath = np.where(trialDeathDates < communityBurnIn)[0]
        trialBirthDates[earlyDeath] = trialDeathDates[earlyDeath]
        trialDeathDates[earlyDeath] += getLifeSpans(len(earlyDeath), params)

    demography = {'birthDate': trialBirthDates - communityBurnIn, 'deathDate': trialDeathDates - communityBurnIn}

    contactAgeGroupIndices = pd.cut(x=-demography['birthDate'], bins=params['contactAgeGroupBreaks'],
    labels=np.arange(start=0, stop=len(params['contactAgeGroupBreaks']) - 1)).to_numpy()

    treatmentAgeGroupIndices = pd.cut(x=-demography['birthDate'], bins=params['treatmentAgeGroupBreaks'],
    labels=np.arange(start=0, stop=len(params['treatmentAgeGroupBreaks']) - 1)).to_numpy()

    meanBurdenIndex = pd.cut(x=-demography['birthDate'], bins=np.append(0, params['equiData']['ageValues']),
    labels=np.arange(start=0, stop=len(params['equiData']['ageValues']))).to_numpy()

    wTotal = np.random.poisson(lam=si * params['equiData']['hatProfile'][meanBurdenIndex] * 2, size=params['N'])
    worms = dict(total=wTotal, female=np.random.binomial(n=wTotal, p=0.5, size=params['N']))
    stableFreeLiving = params['equiData']['L_stable'] * 2

    SD = {'si': si,
          'worms': worms,
          'freeLiving': stableFreeLiving,
          'demography': demography,
          'contactAgeGroupIndices': contactAgeGroupIndices,
          'treatmentAgeGroupIndices': treatmentAgeGroupIndices}

    return SD

def calcRates(params, SD):

    '''
    This function calculates the event rates; the events are
    new worms and worms death.

    Parameters
    ----------
    params: dict
        dictionary containing the parameter names and values;

    SD: dict
        dictionary containing the equilibrium parameter values;

    Returns
    -------
    array of event rates;
    '''

    hostInfRates = SD['freeLiving'] * SD['si'] * params['contactRates'][SD['contactAgeGroupIndices']]
    deathRate = params['sigma'] * np.sum(SD['worms']['total'])

    return np.append(hostInfRates, deathRate)

def doFreeLive(params, SD, dt):

    '''
    This function updates the freeliving population deterministically.

    Parameters
    ----------
    params: dict
        dictionary containing the parameter names and values;

    SD: dict
        dictionary containing the initial equilibrium parameter values;

    dt: float
        time interval;

    Returns
    -------
    SD: dict
        dictionary containing the updated equilibrium parameter values;
    '''

    # polygamous reproduction; female worms produce fertilised eggs only if there's at least one male worm around
    if params['reproFuncName'] == 'epgFertility' and params['SR']:
        productivefemaleworms = np.where(SD['worms']['total'] == SD['worms']['female'], 0, SD['worms']['female'])

    elif params['reproFuncName'] == 'epgFertility' and not params['SR']:
        productivefemaleworms = SD['worms']['female']

    # monogamous reproduction; only pairs of worms produce eggs
    elif params['reproFuncName'] == 'epgMonog':
        productivefemaleworms = np.minimum(SD['worms']['total'] - SD['worms']['female'], SD['worms']['female'])

    eggOutputPerHost = params['lambda'] * productivefemaleworms * np.exp(-productivefemaleworms * params['gamma'])
    eggsProdRate = 2 * params['psi'] * np.sum(eggOutputPerHost * params['rho'][SD['contactAgeGroupIndices']]) / params['N']
    expFactor = np.exp(-params['LDecayRate'] * dt)

    SD['freeLiving'] = SD['freeLiving'] * expFactor + eggsProdRate * (1 - expFactor) / params['LDecayRate']

    return SD

def doDeath(params, SD, t):

    '''
    Death and aging function.

    Parameters
    ----------
    params: dict
        dictionary containing the parameter names and values;

    SD: dict
        dictionary containing the initial equilibrium parameter values;

    t: int
        time step;

    Returns
    -------
    SD: dict
        dictionary containing the updated equilibrium parameter values;
    '''

    # identify the indices of the dead
    theDead = np.where(SD['demography']['deathDate'] < t)[0]

    if len(theDead) != 0:

        # update the birth dates and death dates
        SD['demography']['birthDate'][theDead] = t - 0.001
        SD['demography']['deathDate'][theDead] = t + getLifeSpans(len(theDead), params)

        # they also need new force of infections (FOIs)
        SD['si'][theDead] = np.random.gamma(size=len(theDead), scale=1 / params['k'], shape=params['k'])

        # kill all their worms
        SD['worms']['total'][theDead] = 0
        SD['worms']['female'][theDead] = 0

    # update the contact age categories
    SD['contactAgeGroupIndices'] = pd.cut(x=t - SD['demography']['birthDate'], bins=params['contactAgeGroupBreaks'],
    labels=np.arange(0, len(params['contactAgeGroupBreaks']) - 1)).to_numpy()

    # update the treatment age categories
    SD['treatmentAgeGroupIndices'] = pd.cut(x=t - SD['demography']['birthDate'], bins=params['treatmentAgeGroupBreaks'],
    labels=np.arange(0, len(params['treatmentAgeGroupBreaks']) - 1)).to_numpy()

    return SD

def getAttendance(adherenceFactors, coverage):

    '''
    This function calculates the probability of attendance for each individual.

    Parameters
    ----------
    adherenceFactors: float
        array of adherence factors;

    coverage: float
        array of coverages for each age group;

    Returns
    -------
    boolean array with True / False attendance indicators;
    '''

    # attendance is undefined if coverage is NaN or zero
    attendance = []
    for i in range(len(coverage)):
        if pd.notna(coverage[i]) and coverage[i] != 0:
            p = adherenceFactors[i] ** ((1 - coverage[i]) / coverage[i])
            attendance.append((np.random.uniform(low=0, high=1, size=1) < p)[0])
        else:
            attendance.append(False)

    return attendance

def doChemo(params, SD, mdaRound):

    '''
    Chemotherapy function.

    Parameters
    ----------
    params: dict
        dictionary containing the parameter names and values;

    SD: dict
        dictionary containing the initial equilibrium parameter values;

    mdaRound: int;
        index of current chemotherapy time, used to extract the
        corresponding age-dependent coverage fractions.

    Returns
    -------
    SD: dict
        dictionary containing the updated equilibrium parameter values;
    '''

    # coverage is undefined if -birthDate is negative and treatmentAgeGroupIndices is NaN
    coverage = []
    for i in range(len(SD['treatmentAgeGroupIndices'])):
        if pd.notna(SD['treatmentAgeGroupIndices'][i]) and SD['treatmentAgeGroupIndices'][i] > 0:
            coverage.append(params['coverage'][mdaRound, :][np.int(SD['treatmentAgeGroupIndices'][i])])
        else:
            coverage.append(np.nan)

    # decide which individuals are treated, treatment is random
    attendance = getAttendance(np.random.uniform(low=0, high=1, size=params['N']), coverage)

    # calculate the number of dead worms
    femaleToDie = np.random.binomial(size=np.sum(attendance), n=SD['worms']['female'][attendance],
    p=params['DrugEfficacy'])

    maleToDie = np.random.binomial(size=np.sum(attendance), n=SD['worms']['total'][attendance] -
    SD['worms']['female'][attendance], p=params['DrugEfficacy'])

    SD['worms']['female'][attendance] -= femaleToDie
    SD['worms']['total'][attendance] -= (maleToDie + femaleToDie)

    return SD

def getPsi(params):

    '''
    This function calculates the psi parameter.

    Parameters
    ----------
    params: dict
        dictionary containing the parameter names and values;

    Returns
    -------
    value of the psi parameter;
    '''

    # higher resolution
    deltaT = 0.1

    # inteval-centered ages for the age intervals, midpoints from 0 to maxHostAge
    modelAges = np.arange(start=0, stop=params['maxHostAge'], step=deltaT) + 0.5 * deltaT

    # hostMu for the new age intervals
    hostMu = params['hostMuData'][pd.cut(x=modelAges, bins=params['muBreaks'], labels=np.arange(start=0,
    stop=len(params['hostMuData']))).to_numpy()]

    hostSurvivalCurve = np.exp(-np.cumsum(hostMu * deltaT))
    MeanLifespan = np.sum(hostSurvivalCurve[:len(modelAges)]) * deltaT

    modelAgeGroupCatIndex = pd.cut(x=modelAges, bins=params['contactAgeGroupBreaks'], labels=np.arange(start=0,
    stop=len(params['contactAgeGroupBreaks']) - 1)).to_numpy()

    betaAge = params['contactRates'][modelAgeGroupCatIndex]
    rhoAge = params['rho'][modelAgeGroupCatIndex]

    wSurvival = np.exp(-params['sigma'] * modelAges)

    B = np.array([np.sum(betaAge[: i] * np.flip(wSurvival[: i])) * deltaT for i in range(1, 1 + len(hostMu))])

    return params['R0'] * MeanLifespan * params['LDecayRate'] / (params['lambda'] * params['z'] *
           np.sum(rhoAge * hostSurvivalCurve * B) * deltaT)

def getLifeSpans(nSpans, params):

    '''
    This function draws the lifespans from the population survival curve.

    Parameters
    ----------
    nSpans: int
        number of drawings;

    params: dict
        dictionary containing the parameter names and values;

    Returns
    -------
    array containing the lifespan drawings;
    '''

    u = np.random.uniform(low=0, high=1, size=nSpans) * np.max(params['hostAgeCumulDistr'])
    spans = np.array([np.argmax(u[i] < params['hostAgeCumulDistr']) for i in range(nSpans)])

    return params['muAges'][spans]

def getEquilibrium(params):

    '''
    This function returns a dictionary containing the equilibrium worm burden
    with age and the reservoir value as well as the breakpoint reservoir value
    and other parameter settings.

    Parameters
    ----------
    params: dict
        dictionary containing the parameter names and values;

    Returns
    -------
    dictionary containing the equilibrium parameter settings;
    '''

    # higher resolution
    deltaT = 0.1

    # inteval-centered ages for the age intervals, midpoints from 0 to maxHostAge
    modelAges = np.arange(start=0, stop=params['maxHostAge'], step=deltaT) + 0.5 * deltaT

    # hostMu for the new age intervals
    hostMu = params['hostMuData'][pd.cut(x=modelAges, bins=params['muBreaks'], labels=np.arange(start=0,
    stop=len(params['hostMuData'])))]

    hostSurvivalCurve = np.exp(-np.cumsum(hostMu * deltaT))
    MeanLifespan = np.sum(hostSurvivalCurve[:len(modelAges)]) * deltaT

    modelAgeGroupCatIndex = pd.cut(x=modelAges, bins=params['contactAgeBreaks'], labels=np.arange(start=0,
    stop=len(params['contactAgeBreaks']) - 1)).to_numpy()

    betaAge = params['contactRates'][modelAgeGroupCatIndex]
    rhoAge = params['rho'][modelAgeGroupCatIndex]

    wSurvival = np.exp(-params['sigma'] * modelAges)

    # this variable times L is the equilibrium worm burden
    Q = np.array([np.sum(betaAge[: i] * np.flip(wSurvival[: i])) * deltaT for i in range(1, 1 + len(hostMu))])

    # converts L values into mean force of infection
    FOIMultiplier = np.sum(betaAge * hostSurvivalCurve) * deltaT / MeanLifespan

    # upper bound on L
    SRhoT = np.sum(hostSurvivalCurve * rhoAge) * deltaT
    R_power = 1 / (params['k'] + 1)
    L_hat = params['z'] * params['lambda'] * params['psi'] * SRhoT * params['k'] * (params['R0'] ** R_power - 1) / \
    (params['R0'] * MeanLifespan * params['LDecayRate'] * (1 - params['z']))

    # now evaluate the function K across a series of L values and find point near breakpoint;
    # L_minus is the value that gives an age-averaged worm burden of 1; negative growth should
    # exist somewhere below this
    L_minus = MeanLifespan / np.sum(Q * hostSurvivalCurve * deltaT)
    test_L = np.append(np.linspace(start=0, stop=L_minus, num=10), np.linspace(start=L_minus, stop=L_hat, num=20))

    def K_valueFunc(currentL, params):

        return params['psi'] * np.sum(params['reproFunc'](currentL * Q, params) * rhoAge * hostSurvivalCurve * deltaT) / \
               (MeanLifespan * params['LDecayRate']) - currentL

    K_values = np.vectorize(K_valueFunc)(currentL=test_L, params=params)

    # now find the maximum of K_values and use bisection to find critical Ls
    iMax = np.argmax(K_values)
    mid_L = test_L[iMax]

    if K_values[iMax] < 0:
        return dict(stableProfile=0 * Q,
                    ageValues=modelAges,
                    L_stable=0,
                    L_breakpoint=np.nan,
                    K_values=K_values,
                    L_values=test_L,
                    FOIMultiplier=FOIMultiplier)

    # find the top L
    L_stable = ParallelFuncs.Bisect(f=K_valueFunc, a=mid_L, b=4 * L_hat, args=(params))

    # find the unstable L
    L_break = test_L[1] / 50

    if K_valueFunc(L_break, params) < 0:  # if it is less than zero at this point, find the zero
        L_break = ParallelFuncs.Bisect(f=K_valueFunc, a=L_break, b=mid_L, args=(params))

    stableProfile = L_stable * Q
    hatProfile = L_hat * Q

    return dict(stableProfile=stableProfile,
                hatProfile=hatProfile,
                ageValues=modelAges,
                hostSurvival=hostSurvivalCurve,
                L_breakpoint=L_break,
                K_values=K_values,
                L_values=test_L,
                FOIMultiplier=FOIMultiplier,
                L_stable=L_stable,
                L_hat=L_hat)

def getSetOfEggCounts(total, female, params):

    '''
    This function returns a set of readings of egg counts from a vector of individuals,
    according to their reproductive biology.

    Parameters
    ----------
    total: int
        array of total worms;

    female: int
        array of female worms;

    params: dict
        dictionary containing the parameter names and values;

    Returns
    -------
    random set of egg count readings from a single sample;
    '''

    if params["unfertilized"]:

        meanCount = female * params['lambda'] * params['z'] ** female

    else:

        eggProducers = np.where(total == female, 0, female)
        meanCount = eggProducers * params['lambda'] * params['z'] ** eggProducers

    return np.random.negative_binomial(size=len(meanCount), p=params['k_epg'] / (meanCount + params['k_epg']), n=params['k_epg'])

def getVillageMeanCountsByHost(SD, params, nSamples=1):

    '''
    This function returns the mean egg count across readings by host.

    Parameters
    ----------
    SD: dict
        dictionary containing the equilibrium parameter values;

    params: dict
        dictionary containing the parameter names and values;

    nSamples: int
        number of samples;

    Returns
    -------
    array of mean egg counts;
    '''

    return getSetOfEggCounts(SD['worms']['total'], SD['worms']['female'], params)

def getWormCountsByVillage(SD, t, ageBand, params, nSamples=1, hostSampleSizeFrac=1.0):

    '''
    This function provides sampled, age-cat worm counts.

    Parameters
    ----------
    SD: dict
        dictionary containing the equilibrium parameter values;

    t: int
        time step;

    ageBand: int
        array with age group boundaries;

    params: dict
        dictionary containing the parameter names and values;

    nSamples: int
        number of samples;

    hostSampleSizeFrac: float;
        host sample size fraction;

    Returns
    -------
    dictionary with sampled mean egg count and village sample size;
    '''

    # get readings from the hosts
    meanEggCounts = getVillageMeanCountsByHost(SD, params, nSamples)

    # get ages, filter age group
    ageGroups = pd.cut(x=t - SD['demography']['birthDate'], bins=np.append(-10, np.append(ageBand, 150)),
    labels=np.array([1, 2, 3])).to_numpy()

    currentAgeGroupMeanEggCounts = meanEggCounts[ageGroups == 2]
    villageSampleSize = np.int(np.floor(len(currentAgeGroupMeanEggCounts) * hostSampleSizeFrac))
    meanEggCountSample = np.random.choice(a=currentAgeGroupMeanEggCounts, size=villageSampleSize, replace=False)

    return dict(meanEggCountSample=meanEggCountSample, villageSampleSize=villageSampleSize)

def getAgeCatSampledPrevByVillage(SD, t, ageBand, params, nSamples=1, hostSampleSizeFrac=1.0):

    '''
    This function provides sampled, age-cat worm prevalence.

    Parameters
    ----------
    SD: dict
        dictionary containing the equilibrium parameter values;

    t: int
        time step;

    ageBand: int
        array with age group boundaries;

    params: dict
        dictionary containing the parameter names and values;

    nSamples: int
        number of samples;

    hostSampleSizeFrac: float;
        host sample size fraction;

    Returns
    -------
    sampled worm prevalence;
    '''

    countData = getWormCountsByVillage(SD, t, ageBand, params, nSamples, hostSampleSizeFrac)

    return np.sum(countData['meanEggCountSample'] > 0.9) / countData['villageSampleSize']

def getMediumHeavyPrevalenceByVillage(SD, t, ageBand, eggCountThreshold, params, nSamples=1,
                                      hostSampleSizeFrac=1.0):

    '''
    This function calculates the prevalence by village.

    Parameters
    ----------
    SD: dict
        dictionary containing the equilibrium parameter values;

    t: int
        time step;

    ageBand: int
        array with age group boundaries;

    eggCountThreshold: float
        intensity threshold;

    params: dict
        dictionary containing the parameter names and values;

    nSamples: int
        number of samples;

    hostSampleSizeFrac: float;
        host sample size fraction;

    Returns
    -------
    prevalence by village;
    '''

    countData = getWormCountsByVillage(SD, t, ageBand, params, nSamples, hostSampleSizeFrac)

    return np.sum(countData['meanEggCountSample'] >= eggCountThreshold) / countData['villageSampleSize']

def doEvents(rates, SD, dt):

    '''
    This function enacts multiple events; the events are
    new worms and worms death.

    Parameters
    ----------
    rates: float
        array of event rates;

    SD: dict
        dictionary containing the initial equilibrium parameter values;

    dt: float
        time interval;

    Returns
    -------
    SD: dict
        dictionary containing the updated equilibrium parameter values;
    '''

    # make new worms
    L = len(rates) - 1

    females = np.random.poisson(lam=dt * rates[:-1] / 2, size=L)
    males = np.random.poisson(lam=dt * rates[:-1] / 2, size=L)

    SD['worms']['total'] += males + females
    SD['worms']['female'] += females

    total = np.sum(SD['worms']['total'])

    # kill some worms
    if total > 0:

        females = np.random.poisson(lam=dt * rates[-1] * SD['worms']['female'] / total, size=L)
        males = np.random.poisson(lam=dt * rates[-1] * (SD['worms']['total'] - SD['worms']['female']) / total, size=L)

        # don't kill worms that don't exist
        females = np.minimum(females, SD['worms']['female'])
        males = np.minimum(males, SD['worms']['total'] - SD['worms']['female'])

        SD['worms']['total'] -= males + females
        SD['worms']['female'] -= females

    return SD

def doRealization(params, seed):

    '''
    This function generates a single simulation path
    starting from time zero.

    Parameters
    ----------
    params: dict
        dictionary containing the parameter names and values;

    seed: int
        random seed;

    Returns
    -------
    simData: dict
        dictionary with simulation results;
    '''

    # set the random seed
    np.random.seed(seed)

    # initialize the simulation data
    simData = setupSD(params)

    # start time
    t = 0

    # end time
    maxTime = copy.deepcopy(params['maxTime'])

    # max step
    maxStep = copy.deepcopy(params['maxStep'])

    # time at which to update the freelive population
    freeliveTime = t

    # times at which data should be recorded
    outTimes = copy.deepcopy(params['outTimings'])

    # times at which chemotherapy is given
    chemoTimes = copy.deepcopy(params['chemoTimings'])

    # time when data should be recorded next
    nextOutIndex = np.argmin(outTimes)
    nextOutTime = outTimes[nextOutIndex]

    # time at which individuals' age is advanced next
    ageingInt = 1 / 52
    nextAgeTime = 1 / 52
    tau = 1 / 52

    # time at which individuals receive next chemotherapy
    nextChemoIndex = np.argmin(chemoTimes)
    nextChemoTime = chemoTimes[nextChemoIndex]

    # next event
    nextStep = np.min([nextOutTime, t + maxStep, nextAgeTime, nextChemoTime])

    # initialise empty list to store results
    prevKKSAC = []
    prevMHISAC = []

    # counter for chemotherapy treatments
    mdaRound = 0

    # run stochastic algorithm
    while t + 0.001 < maxTime:

        rates = calcRates(params, simData)
        sumRates = np.sum(rates)

        if sumRates > 0.001:

            dt = np.min([tau, nextStep - t, maxTime - t])
            simData = doEvents(rates, simData, dt)

        else:

            dt = nextStep - t

        if t + dt < nextStep:

            t += dt

        else:

            t = nextStep
            simData = doFreeLive(params, simData, nextStep - freeliveTime)
            freeliveTime = nextStep
            timeBarrier = nextStep + 0.001

            # ageing and death
            if timeBarrier > nextAgeTime:

                simData = doDeath(params, simData, t)
                nextAgeTime += ageingInt

            # chemotherapy
            if timeBarrier > nextChemoTime:

                simData = doDeath(params, simData, t)
                simData = doChemo(params, simData, mdaRound)
                chemoTimes[nextChemoIndex] = maxTime + 10
                nextChemoIndex = np.argmin(chemoTimes)
                nextChemoTime = chemoTimes[nextChemoIndex]
                mdaRound += 1

            # output
            if timeBarrier > nextOutTime:

                prevKKSAC.append(getAgeCatSampledPrevByVillage(simData, t, np.array([5, 14]), params))
                prevMHISAC.append(getMediumHeavyPrevalenceByVillage(simData, t, np.array([5, 14]), params['moderateIntensityCount'], params))

                outTimes[nextOutIndex] = maxTime + 10
                nextOutIndex = np.argmin(outTimes)
                nextOutTime = outTimes[nextOutIndex]

            nextStep = np.min([nextOutTime, t + maxStep, nextAgeTime])

    # save the state of the simulations
    simData['state'] = np.random.get_state()

    # save the simulation times
    simData['times'] = {'start_time': 0,
                        'end_time': t,
                        'maxTime': maxTime,
                        'freeliveTime': freeliveTime,
                        'nextStep': nextStep,
                        'nextAgeTime': nextAgeTime}

    # save the output times
    simData['outTimings'] = params['outTimings']

    # save the simulated prevalence
    simData['prevKKSAC'] = prevKKSAC
    simData['prevMHISAC'] = prevMHISAC

    return simData

def addRealization(params, simData, times, state):

    '''
    This function generates a single simulation path
    starting from a past time point different from zero.

    Parameters
    ----------
    params: dict
        dictionary containing the parameter names and values;

    simData: dict
        dictionary with previous simulation results;

    times: dict
        dictionary with previous simulation times;

    state: array
        random state;

    Returns
    -------
    simData: dict
        dictionary with new simulation results;
    '''

    # set the random state
    np.random.set_state(state)

    # start time
    t = copy.deepcopy(times['end_time'])

    # end time
    maxTime = copy.deepcopy(params['maxTime'])

    # max step
    maxStep = copy.deepcopy(params['maxStep'])

    # time at which to update the freelive population
    freeliveTime = copy.deepcopy(times['freeliveTime'])

    # times at which data should be recorded
    outTimes = copy.deepcopy(params['outTimings'])

    # times at which chemotherapy is given
    chemoTimes = copy.deepcopy(params['chemoTimings'])

    # time when data should be recorded next
    nextOutIndex = np.argmin(outTimes)
    nextOutTime = outTimes[nextOutIndex]

    # time at which individuals' age is advanced next
    ageingInt = 1 / 52
    nextAgeTime = copy.deepcopy(times['nextAgeTime'])
    tau = 1 / 52

    # time at which individuals receive next chemotherapy
    nextChemoIndex = np.argmin(chemoTimes)
    nextChemoTime = chemoTimes[nextChemoIndex]

    # next event
    nextStep = copy.deepcopy(times['nextStep'])

    # initialise empty list to store results
    prevKKSAC = []
    prevMHISAC = []

    # counter for chemotherapy treatments
    mdaRound = 0

    # run stochastic algorithm
    while t + 0.001 < maxTime:

        rates = calcRates(params, simData)
        sumRates = np.sum(rates)

        if sumRates > 0.001:

            dt = np.min([tau, nextStep - t, maxTime - t])
            simData = doEvents(rates, simData, dt)

        else:

            dt = nextStep - t

        if t + dt < nextStep:

            t += dt

        else:

            t = nextStep
            simData = doFreeLive(params, simData, nextStep - freeliveTime)
            freeliveTime = nextStep
            timeBarrier = nextStep + 0.001

            # ageing and death
            if timeBarrier > nextAgeTime:

                simData = doDeath(params, simData, t)
                nextAgeTime += ageingInt

            # chemotherapy
            if timeBarrier > nextChemoTime:

                simData = doDeath(params, simData, t)
                simData = doChemo(params, simData, mdaRound)
                chemoTimes[nextChemoIndex] = maxTime + 10
                nextChemoIndex = np.argmin(chemoTimes)
                nextChemoTime = chemoTimes[nextChemoIndex]

                mdaRound += 1

            # output
            if timeBarrier > nextOutTime:

                prevKKSAC.append(getAgeCatSampledPrevByVillage(simData, t, np.array([5, 14]), params))
                prevMHISAC.append(getMediumHeavyPrevalenceByVillage(simData, t, np.array([5, 14]), params['moderateIntensityCount'], params))

                outTimes[nextOutIndex] = maxTime + 10
                nextOutIndex = np.argmin(outTimes)
                nextOutTime = outTimes[nextOutIndex]

            nextStep = np.min([nextOutTime, t + maxStep, nextAgeTime])

    # save the state of the simulations
    simData['state'] = np.random.get_state()

    # save the simulation times
    simData['times'] = {'start_time': times['end_time'],
                        'end_time': t,
                        'maxTime': maxTime,
                        'freeliveTime': freeliveTime,
                        'nextStep': nextStep,
                        'nextAgeTime': nextAgeTime}

    # save the output times
    simData['outTimings'] = params['outTimings']

    # save the simulated prevalence
    simData['prevKKSAC'] = prevKKSAC
    simData['prevMHISAC'] = prevMHISAC

    return simData

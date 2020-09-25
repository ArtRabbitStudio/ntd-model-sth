from joblib import Parallel, delayed
import multiprocessing
import copy

from sth_simulation.helsim_FUNC import *

num_cores = multiprocessing.cpu_count()

def loadParameters(paramFileName, demogName):

    '''
    This function loads all the parameters from the input text
    files and organizes them in a dictionary.

    Parameters
    ----------
    paramFileName: str
        name of the input text file with the model parameters;

    demogName: str
        subset of demography parameters to be extracted;

    Returns
    -------
    params: dict
        dictionary containing the parameter names and values;
    '''

    # load the parameters
    params = readParams(paramFileName=paramFileName, demogName=demogName)

    # configure the parameters
    params = configure(params)

    # update the parameters
    params['chemoTimings'] = updateChemoTimings(params)
    params['psi'] = getPsi(params)
    params['equiData'] = getEquilibrium(params)
    params['moderateIntensityCount'], params['highIntensityCount'] = setIntensityCount(paramFileName)

    return params

def doRealization(params, i):

    '''
    This function generates a single simulation path.

    Parameters
    ----------
    params: dict
        dictionary containing the parameter names and values;

    i: int
        iteration number;

    Returns
    -------
    results: list
        list with simulation results;
    '''

    # setup simulation data for village
    simData = setupSD(params)

    # start time
    t = 0

    # end time
    maxTime = copy.deepcopy(params['maxTime'])

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
    maxStep = copy.deepcopy(params['maxStep'])
    
    # time at which individuals receive next chemotherapy
    nextChemoIndex = np.argmin(chemoTimes)
    nextChemoTime = chemoTimes[nextChemoIndex]

    # next event
    nextStep = np.min([nextOutTime, t + maxStep, nextChemoTime, nextAgeTime])

    results = list()  # initialise empty list to store results

    # run stochastic algorithm
    while t < maxTime:

        rates = calcRates(params, simData)
        sumRates = np.sum(rates)

        if sumRates > 0.001:

            dt = np.random.exponential(scale=1 / sumRates, size=1)[0]

        else:

            dt = maxTime

        if t + dt < nextStep:

            t += dt

            simData = doEvent(rates, simData)

        else:

            simData = doFreeLive(params, simData, nextStep - freeliveTime)

            t = nextStep
            freeliveTime = nextStep
            timeBarrier = nextStep + 0.001

            # ageing and death
            if timeBarrier > nextAgeTime:

                simData = doDeath(params, simData, t)

                nextAgeTime += ageingInt

            # chemotherapy
            if timeBarrier > nextChemoTime:

                simData = doDeath(params, simData, t)
                simData = doChemo(params, simData)

                chemoTimes[nextChemoIndex] = maxTime + 10
                nextChemoIndex = np.argmin(chemoTimes)
                nextChemoTime = chemoTimes[nextChemoIndex]

            if timeBarrier > nextOutTime:

                results.append(dict(
                    iteration=i,
                    time=t,
                    # worms=copy.deepcopy(simData['worms']),
                    # hosts=copy.deepcopy(simData['demography']),
                    # freeLiving=copy.deepcopy(simData['freeLiving']),
                    prevKKSAC=getAgeCatSampledPrevByVillage(simData, t, np.array([5, 14]), params),
                    prevKKAll=getAgeCatSampledPrevByVillage(simData, t,  np.array([0, 100]), params),
                    prevTrue=villageTruePrev(simData, params),
                    prevTrue2=villageTruePrev2(simData),
                    meanIntensitySAC=getMeanInfectionIntensity(simData, t, np.array([5, 14]), params)['mean'],
                    meanIntensityAll=getMeanInfectionIntensity(simData, t, np.array([0, 100]), params)['mean'],
                    prevHISAC=getMediumHeavyPrevalenceByVillage(simData, t, np.array([5, 14]), params['highIntensityCount'], params),
                    prevMHISAC=getMediumHeavyPrevalenceByVillage(simData, t, np.array([5, 14]), params['moderateIntensityCount'], params),
                ))

                outTimes[nextOutIndex] = maxTime + 10
                nextOutIndex = np.argmin(outTimes)
                nextOutTime = outTimes[nextOutIndex]

            nextStep = np.min([nextOutTime, t + maxStep, nextAgeTime])

    return results

def STH_Simulation(paramFileName, demogName, numReps=None):

    '''
    This function generates multiple simulation paths.

    Parameters
    ----------
    paramFileName: str
        name of the input text file with the model parameters;

    demogName: str
        subset of demography parameters to be extracted;

    numReps: int
        number of simulations;

    Returns
    -------
    df: data frame
        data frame with simulation results;
    '''

    # initialize the parameters
    params = loadParameters(paramFileName, demogName)

    # extract the number of simulations
    if numReps is None:
        numReps = params['numReps']

    print( 'run sims' )
    # run the simulations
    results = Parallel(n_jobs=num_cores)(delayed(doRealization)(params, i) for i in range(numReps))

	print( 'flatten output' )
    # flatten the output
    output = [item for sublist in results for item in sublist]

	print ( 'transform to data frame' )
    # transform the output to data frame
    df = pd.DataFrame(output)

    return df

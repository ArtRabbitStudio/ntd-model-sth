from joblib import Parallel, delayed
import multiprocessing
import pickle
import time
import functools
import datetime

from sth_simulation.helsim_FUNC import *

def timer(func):
    """Print the runtime of the decorated function"""
    @functools.wraps(func)
    def wrapper_timer(*args, **kwargs):
        start_time = time.perf_counter()    # 1
        print(f"-> Running {func.__name__!r}, starting at {datetime.datetime.now()}")
        value = func(*args, **kwargs)
        end_time = time.perf_counter()      # 2
        run_time = end_time - start_time    # 3
        print(f"=> Finished {func.__name__!r} in {run_time:.4f} secs\n\n")
        return value
    return wrapper_timer

@timer
def STH_Simulation(paramFileName, demogName, MDAFilePath, PrevFilePath, RkFilePath=None, nYears=None,
                   outputFrequency=None, numReps=None, SaveOutput=False, OutSimFilePath=None, InSimFilePath=None):

    '''
    Longitudinal simulations.

    Parameters:
    -----------
    paramFileName: str
        This is the name of the input text file with the
        model parameters.

    demogName: str
        This is the subset of demography parameters to be
        extracted from the Demographies.txt file.

    MDAFilePath: str
        This is the path to the input CSV file with the
        MDA times and with the respective coverage
        fractions for each age group.

    PrevFilePath: str
        This is the path where the output CSV file with
        the simulated prevalence will be saved.

    RkFilePath: str
        This is the path to the input CSV file with the
        random seed, R0, and k to be used for each simulation.
        If not provided, the code will use the parameters
        in the input text file.

    nYears: int
        This is the number of years for which the
        simulations will be performed.

    outputFrequency: int
        This is the number of time points saved down
        in the CSV file for each year.

    numReps: int
        This is the number of simulations. If not provided,
        the number of simulations will be inferred from the
        input text file or CSV file.

    SaveOutput: bool
        If True, the last state of the simulations will
        be saved in a pickle file. If False, the last
        state of the simulations will not be saved.

    OutSimFilePath: str
        This is the path where the output pickle file with
        the last state of the simulations will be saved. It
        is only required when SaveOutput = True.

    InSimFilePath: str
        This is the path where the input pickle file with
        the last state of the simulations has been saved.
        If this is provided, the code will resume the
        previous simulations from this state. If this is
        not provided, the code will start new simulations
        from scratch.

    Returns:
    -----------
    None.
    '''

    # make sure that the user has provided all the necessary inputs
    if '.csv' not in PrevFilePath:
        message = 'Please provide the directory to the output CSV file with the simulated prevalence.'

    elif RkFilePath is not None and '.csv' not in RkFilePath:
        message = 'Please provide the directory to the input CSV file with the simulation parameters.'

    elif MDAFilePath is not None and '.csv' not in MDAFilePath:
        message = 'Please provide the directory to the input CSV file with the MDA parameters.'

    elif SaveOutput and (OutSimFilePath is None or '.p' not in OutSimFilePath):
        message = 'Please provide the directory to the output pickle file.'

    elif InSimFilePath is not None and '.p' not in InSimFilePath:
        message = 'Please provide the directory to the input pickle file.'

    else:

        # load the parameters from the TXT files
        params = readParams(paramFileName=paramFileName, demogName=demogName)

        # load the MDA parameters from the CSV file
        timeparams = pd.read_csv(MDAFilePath)
        timeparams.columns = [s.replace(' ', '') for s in timeparams.columns]
        params['chemoTimings'] = timeparams.iloc[:, 0].values
        params['coverage'] = timeparams.iloc[:, 1:].values / 100

        # overwrite the output frequency
        if outputFrequency is not None:
            params['outputFrequency'] = outputFrequency

        # overwrite the number of years
        if nYears is not None:
            params['maxTime'] = nYears

        if RkFilePath is not None:

            # load the simulation parameters from the CSV file
            simparams = pd.read_csv(RkFilePath)
            simparams.columns = [s.replace(' ', '') for s in simparams.columns]

            # overwrite the number of simulations
            if numReps is not None:
                params['numReps'] = numReps
            else:
                params['numReps'] = simparams.shape[0]

            # define the lists of random seeds, R0 and k
            seed = simparams.iloc[:numReps, 0].tolist()
            R0 = simparams.iloc[:numReps, 1].tolist()
            k = simparams.iloc[:numReps, 2].tolist()

        else: # use the simulation parameters in the text file

            # overwrite the number of simulations
            if numReps is not None:
                params['numReps'] = numReps

            # define the lists of random seeds, R0 and k
            seed = [np.int(params['seed'] + i) for i in range(params['numReps'])]
            R0 = [params['R0']] * params['numReps']
            k = [params['k']] * params['numReps']

        if InSimFilePath is None:  # start new simulations

            def multiple_simulations(params, i):

                # update the parameters
                params['R0'] = R0[i]
                params['k'] = k[i]

                # configure the parameters
                params = configure(params)
                params['psi'] = getPsi(params)
                params['equiData'] = getEquilibrium(params)
                params['moderateIntensityCount'], params['highIntensityCount'] = setIntensityCount(paramFileName)

                # generate a simulation path
                return doRealization(params=params, seed=seed[i])

        else:  # continue previous simulations

            def multiple_simulations(params, i):

                # load the previous simulation results
                data = pickle.load(open(InSimFilePath, 'rb'))[i]

                # extract the previous simulation output
                keys = ['si', 'worms', 'freeLiving', 'demography', 'contactAgeGroupIndices', 'treatmentAgeGroupIndices']
                simData = dict((key, data[key]) for key in keys)

                # extract the previous random state
                state = data['state']

                # extract the previous simulation times
                times = data['times']

                # increment the simulation times
                params['maxTime'] += times['maxTime']
                params['outputOffset'] = times['end_time'] + 1 / params['outputFrequency']

                # update the parameters
                params['R0'] = R0[i]
                params['k'] = k[i]

                # configure the parameters
                params = configure(params)
                params['psi'] = getPsi(params)
                params['equiData'] = getEquilibrium(params)
                params['moderateIntensityCount'], params['highIntensityCount'] = setIntensityCount(paramFileName)

                # add a simulation path
                return addRealization(params=params, simData=simData, times=times, state=state)

        # run the simulations
        num_cores = multiprocessing.cpu_count()

        start_time = time.time()

        out = Parallel(n_jobs=num_cores)(delayed(multiple_simulations)(params, i) for i in range(params['numReps']))

        end_time = time.time()

        # save the simulated prevalence in a CSV file
        years = out[0]['outTimings'][:len(out[0]['prevKKSAC'])]
        years = [str(np.int(np.round(year))) if outputFrequency == 1 else str(np.round(year, 2)) for year in years]
        columns = ['Random Generator', 'R0', 'k'] + ['prevKKSAC year ' + year for year in years]

        df = pd.DataFrame(columns=columns)
        df['Random Generator'] = seed
        df['R0'] = R0
        df['k'] = k

        for i in range(len(out)):
            df.iloc[i, 3:] = out[i]['prevKKSAC']

        df.to_csv(PrevFilePath, index=None)

        # save the simulated data in a pickle file
        if SaveOutput:
            pickle.dump(out, open(OutSimFilePath, 'wb'))

        message = 'Running time: ' + format((end_time - start_time), '.0f') + ' seconds.'

    print(message)

    return None

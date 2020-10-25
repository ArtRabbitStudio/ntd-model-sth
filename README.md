# STH Simulation Model

To run the STH simulation model, import the `STH_Simulation()` function from the `helsim_RUN` module in the `sth_simulation` package.

The `STH_Simulation()` function requires the following inputs:

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
        
Some examples are included in the Jupyter notebook `sth_tests.ipynb` in the `tests` folder.

### How to run

Install [pipenv](https://drive.google.com/drive/folders/1Or6lUkymYd_p031xKGZLcnTV4GYf-oYb) according to the instructions for your OS, then `cd` to the project directory and run:

```
    $ pipenv install . # sets up per-project python environment ('env')
    $ pipenv shell # starts a per-project shell using that env
    (ntd-model-sth) $ python sth_run.py # runs the model
```
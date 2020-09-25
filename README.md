# STH Simulation Model

To run the STH simulation model, import the `STH_Simulation()` function from the `helsim_RUN` module in the `sth_simulation` package.

The `STH_Simulation()` function requires the following inputs:

- `paramFileName`: name of the input file with the model parameters. The available files are listed in the following table.

| paramFileName |
| :--- | 
| AscarisParameters_high.txt | 
| AscarisParameters_moderate.txt  | 
| HookwormParameters_high.txt | 
| HookwormParameters_moderate.txt  | 
| TrichurisParameters_high.txt | 
| TrichurisParameters_moderate.txt  | 

- `demogName`: name of the demography. The available demographies are listed in the following table.

| demogName | 
| :--- | 
| Default | 
| WHOGeneric  | 
| UgandaRural | 
| KenyaKDHS  | 
| Flat | 

- `numReps`: number of simulations. If not provided, the number of simulations is extracted from the parameters file.

The `STH_Simulation()` function returns a data frame with the following columns: `iteration`, `time`, `prevKKSAC`, 
`prevKKAll`, `prevTrue`, `prevTrue2`, `meanIntensitySAC`, `meanIntensityAll`, `prevHISAC`, and `prevMHISAC`.

The output data frame can be exported in several different formats; see `sth_results.json` for an example of the results in JSON format.

See also `sth_run.py` for an example of how to use the `STH_Simulation()` function.
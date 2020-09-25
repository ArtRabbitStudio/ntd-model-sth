from sth_simulation.helsim_RUN import STH_Simulation

df = STH_Simulation(paramFileName='AscarisParameters_high.txt', demogName='KenyaKDHS', numReps=100)

df.to_json('sth_results.json')

df.groupby(by='time')['prevKKSAC'].mean().plot(title='prevKKSAC')


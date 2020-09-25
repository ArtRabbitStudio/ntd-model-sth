from sth_simulation.helsim_RUN import STH_Simulation

print('pre-run')
df = STH_Simulation(paramFileName='AscarisParameters_high.txt', demogName='KenyaKDHS', numReps=100)
print('post-run')

df.to_json('sth_results.json')

df.groupby(by='time')['prevKKSAC'].mean().plot(title='prevKKSAC')


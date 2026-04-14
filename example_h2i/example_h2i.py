# Scripts and configuration files adapted from Example 29_wind_ard in the H2I github repository
# https://github.com/NatLabRockies/H2Integrate/tree/main/examples/29_wind_ard

#%%
from h2integrate.core.h2integrate_model import H2IntegrateModel

#%%
h2i = H2IntegrateModel("example_h2i.yaml")

# Adjust configuration to avoid bug with "shear"

# Run the model
h2i.run()

h2i.post_process()

# %%

import pandas as pd
import matplotlib.pyplot as plt

power_wind = h2i.prob.get_val("wind.ard_electricity_out")

power_solar = h2i.prob.get_val("solar.electricity_out")


wf_aep = sum(power_wind)
pv_aep = sum(power_solar)


# print(, flush=True)
print('Solar PV Farm AEP: {:.1f} GWh'.format(pv_aep*1e-6), flush=True)

fig, ax = plt.subplots(1,1, figsize = (10,2.5))
ax.plot(power_wind*1e-6)
ax.set_ylabel('Wind farm power [MW]')
ax.set_xlabel('Time [h]')
ax.annotate('Wind Farm AEP: {:.1f} GWh'.format(wf_aep*1e-9), xy = [2875, 225])
ax.set_ylim([-10, 250])
ax.set_xlim([2700,2940])

fig, ax = plt.subplots(1,1, figsize = (10,2.5))
ax.plot(power_solar*1e-3)
ax.set_ylabel('Solar PV farm power [MW]')
ax.set_xlabel('Time [h]')
ax.annotate('Solar PV Farm AEP: {:.1f} GWh'.format(pv_aep*1e-6), xy = [2870, 350])
ax.set_ylim([-10, 400])
ax.set_xlim([2700,2940])


#%% save to file

filename_power = '../power_timeseries.csv'

dt = pd.DataFrame({ 'power_wind': power_wind, 'power_solar': power_solar*1e3})

dt.to_csv(filename_power)
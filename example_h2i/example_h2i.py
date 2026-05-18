# Scripts and configuration files adapted from Example 29_wind_ard in the H2I github repository
# https://github.com/NatLabRockies/H2Integrate/tree/main/examples/29_wind_ard

# Run with H2I version 0.8.0

#%%
from h2integrate.core.h2integrate_model import H2IntegrateModel
import pandas as pd
import matplotlib.pyplot as plt

# Initialize model
h2i = H2IntegrateModel("example_h2i.yaml")

# Run the model
h2i.run()

h2i.post_process()

# %%
power_wind = h2i.prob.get_val("wind.ard_electricity_out")
power_solar = h2i.prob.get_val("solar.electricity_out")*1e3

wf_aep = h2i.prob.get_val("wind.annual_electricity_produced")[0] #in kWh
pv_aep = h2i.prob.get_val("solar.annual_electricity_produced")[0] #in kWh

#%%
fig, ax = plt.subplots(1,1, figsize = (10,2.5))
ax.plot(power_wind*1e-6)
ax.set_ylabel('Wind farm power [MW]')
ax.set_xlabel('Time [h]')
ax.annotate('Wind Farm AEP: {:.1f} GWh'.format(wf_aep*1e-6), xy = [2875, 225])
ax.set_ylim([-10, 250])
ax.set_xlim([2700,2940])

fig, ax = plt.subplots(1,1, figsize = (10,2.5))
ax.plot(power_solar*1e-6)
ax.set_ylabel('Solar PV farm power [MW]')
ax.set_xlabel('Time [h]')
ax.annotate('Solar PV Farm AEP: {:.1f} GWh'.format(pv_aep*1e-6), xy = [2870, 350])
ax.set_ylim([-10, 400])
ax.set_xlim([2700,2940])


#%% save to file

# dt = pd.DataFrame({ 'power_wind': power_wind})
# dt.to_csv('../power_timeseries_wind.csv')

# dt = pd.DataFrame({ 'power_solar': power_solar})
# dt.to_csv('../power_timeseries_solar.csv')
# %%

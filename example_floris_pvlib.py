'''
Example script to simulate the IEA Wind Reference HPP, using the software FLORIS and pvlib.
'''
#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime
import yaml
from windIO.yaml import load_yaml
# ARD/FLORIS imports
from floris import FlorisModel, TimeSeries

# pvlib imports
from pvlib.location import Location
from pvlib.modelchain import ModelChain
from pvlib.pvsystem import PVSystem
from pvlib.temperature import TEMPERATURE_MODEL_PARAMETERS

#%%
plot_figures = True
write_turbine = True # True of False. Specifies if turbine as specified in yaml files will be written to a file

# ------
# Load the power plant data 
# ------

file_path = 'Reference_HPP.yaml'

hpp_dat = load_yaml(file_path)

p_max = hpp_dat['grid_connection_capacity']  # Maximum power allowed to be send to the grid

# ------
# Load the wind farm and wind resource data
# ------

wind_farm_dat = hpp_dat['wind_farm']['wind_farm']
wind_resource_dat = hpp_dat['wind_farm']['site']['energy_resource']
solar_resource_dat = hpp_dat['solar_pv_farm']['site']['energy_resource']

# Get x and y positions of the wind farm layout
x = wind_farm_dat['layouts']['coordinates']['x'] 
y = wind_farm_dat['layouts']['coordinates']['y'] 
n_turbines = len(x)

rp = wind_farm_dat['turbine']['performance']['rated_power']  # wind turbine rated power [W]



# Load wind speed and wind direction time series
ws_time = wind_resource_dat['wind_resource']['wind_speed']
wd_time = wind_resource_dat['wind_resource']['wind_direction']
ti = wind_resource_dat['wind_resource']['turbulence_intensity']['data']


time_stamp = np.arange(len(wd_time))/24 


#%%
# ------
# Run the wind farm simulation
# ------

wind_timeseries = TimeSeries(
            wind_directions=np.array(wd_time),
            wind_speeds=np.array(ws_time),
            turbulence_intensities=np.array(ti),
        )


fmodel = FlorisModel('defaults')

fmodel.set(layout_x=x, layout_y=y, turbine_type=['nrel_5mw']*n_turbines, wind_data=wind_timeseries)

fmodel.set_operation_model("peak-shaving")
# Run the time simulation 
fmodel.run()

power_wind = fmodel.get_farm_power()
wf_aep = sum(power_wind)
wf_cf = wf_aep/(n_turbines*rp*8760)

print('Wind Farm AEP: {:.1f} GWh'.format(wf_aep*1e-9))
print('Wind Farm capacity factor: {:.1%}'.format(wf_cf))
#%%
if plot_figures:
      fig, ax = plt.subplots(1,1, figsize = (10,2.5))
      ax.plot(power_wind*1e-6)
      ax.set_ylabel('Power [MW]')
      ax.set_xlabel('Time [h]')
      ax.set_xlim([2700,2940])


#%%
# ------
# Load the solar pv farm and solar resource data
# ------

# Define location
latitude = hpp_dat['latitude']
longitude = hpp_dat['longitude']
tz = 'Etc/GMT+1'  # Timezone for the location

# Import solar resource data
time = [ datetime.strptime(x, '%Y-%m-%dT%H:%M:%S%zZ') for x in solar_resource_dat['solar_resource']['time']]
ghi = solar_resource_dat['solar_resource']['ghi']
dhi = solar_resource_dat['solar_resource']['dhi']
dni = solar_resource_dat['solar_resource']['dni']

# Import pv system and pv module parameters
pvsystem = hpp_dat['solar_pv_farm']['pvsystems']
pvmodule = hpp_dat['solar_pv_farm']['pvsystems']['modules']

v_oc = pvmodule['voltage_open_circuit'] # Open-circuit voltage (V)
v_mp = pvmodule['voltage_maximum_power'] # 44.36 #Voltage at maximum-power point (V)
i_mp = pvmodule['current_maximum_power'] # 14.43 #Current at the maximum-power point (A)
i_sc = pvmodule['current_short_circuit'] # 15.13 #Short-circuit current (A)
t_voc = pvmodule['temperature_coefficient_voc'] # Temperature coefficient of Voc (%/deg celcius)
p_dc = pvmodule['nominal_power_rating'] #Nameplate capacity

module_parameters =  {'Voco': v_oc,
                     'Impo': i_mp,
                     'Vmpo': v_mp,
                     'Isco': i_sc,
                     'Bvoco': t_voc*v_oc,
                     'pdc0': p_dc,
                     'gamma_pdc': -0.004, # module temperature coefficient (-/deg c.) Assumed value
                     }


ac_capacity = pvsystem['inverter']['ac_capacity']
modules_per_string = pvsystem['n_modules_per_string'] 
strings_per_inverter = pvsystem['n_strings_per_inverter'] 

inverter_parameters = {'pdc0': ac_capacity,
                      'eta_inv_nom': pvsystem['inverter']['efficiency']}

tilt = pvsystem['tilt'] 
surface_azimuth = pvsystem['surface_azimuth']

n_systems = hpp_dat['solar_pv_farm']['n_systems'] # number of PV systems in the farm

solar_pv_capacity = n_systems*ac_capacity 

# Retrieve temperature model parameters. Assumed value
temperature_model_parameters = TEMPERATURE_MODEL_PARAMETERS['sapm']['open_rack_glass_glass']

#%%
# ------
# Run the solar PV farm simulation with pvlib
# ------

# new_time = [t.timestamp() for t in time]

weather = pd.DataFrame(index = pd.DatetimeIndex(time, tz = 'Europe/Berlin'), data={'ghi': ghi, 'dhi': dhi, 'dni': dni})

location = Location(latitude, longitude, tz)

system = PVSystem(surface_tilt=tilt, surface_azimuth=surface_azimuth,
                  module_parameters=module_parameters,
                  inverter_parameters=inverter_parameters,
                  temperature_model_parameters=temperature_model_parameters,
                  strings_per_inverter=strings_per_inverter, modules_per_string = modules_per_string)

# Create a model chain object using the pvwatts model
model_chain = ModelChain.with_pvwatts(system, location)

# Run the simulation
model_chain.run_model(weather)

# Extract the AC power output for one pv system
power_per_pvsystem = model_chain.results.ac

power_solar = np.array(power_per_pvsystem.to_list())*n_systems

pv_aep = sum(power_solar)
pv_cf = np.mean(power_solar)/np.max(power_solar)

print('Solar PV farm AEP: {:.1f} GWh'.format(pv_aep*1e-9))
print('Solar PV farm capacity factor: {:.1%}'.format(pv_cf))

if plot_figures:
      fig, ax = plt.subplots(1,1, figsize = (10,3))
      ax.plot(power_solar)
      ax.set_ylabel('Power [W]')
      ax.set_xlabel('Time [days]')
      ax.set_xlim([2700, 2940])


#%%
# ------
# Load storage system data data 
# ------

batt_dat = hpp_dat['storage_system']

rt_efficiency_dc = batt_dat['battery_systems']['round_trip_efficiency_nominal']
power_conversion_efficiency = batt_dat['power_conditioning_unit']['efficiency']

rt_efficiency_ac = rt_efficiency_dc*power_conversion_efficiency**2

# power and energy capacity (in W and Wh)
n_battsystems = batt_dat['n_systems']
p_cap_w = n_battsystems * batt_dat['battery_systems']['power_capacity'] 
e_cap_w = n_battsystems * batt_dat['battery_systems']['energy_capacity']

#%% save to file

# filename_power = '../power_timeseries.csv'

# dt = pd.DataFrame({ 'power_wind': power_wind, 'power_solar': power_solar})

# dt.to_csv(filename_power)


'''
Example script to simulate the IEA Wind Reference HPP, using the software PyWake and pvlib.
'''
#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime
from windIO.yaml import load_yaml
# PyWake imports
import xarray as xr
from py_wake.site import XRSite
from py_wake.wind_turbines import WindTurbine
from py_wake.wind_turbines.power_ct_functions import PowerCtTabular
from py_wake.rotor_avg_models import RotorCenter
from py_wake.literature.gaussian_models import Bastankhah_PorteAgel_2014
# pvlib imports
from pvlib.location import Location
from pvlib.modelchain import ModelChain
from pvlib.pvsystem import PVSystem
from pvlib.temperature import TEMPERATURE_MODEL_PARAMETERS

#%%
plot_figures = True

# ------
# Load the power plant data 
# ------

file_path = 'hpp-reference-designs/Reference_HPP.yaml'

hpp_dat = load_yaml(file_path)

p_max = hpp_dat['grid_connection_capacity']  # Maximum power allowed to be send to the grid

# ------
# Load the wind farm and wind resource data
# ------

wind_farm_dat = hpp_dat['wind_farm']
resource_dat = hpp_dat['site']['energy_resource']

# Extract wind resource data (wind rose)
A = resource_dat['wind_resource']['weibull_a']
k = resource_dat['wind_resource']['weibull_k']
freq = resource_dat['wind_resource']['sector_probability']

# Get x and y positions of the wind farm layout
x = wind_farm_dat['layouts']['coordinates']['x'] 
y = wind_farm_dat['layouts']['coordinates']['y'] 

# Load turbine data - NREL 5 MW, extracted from report
wt_dat = wind_farm_dat['turbines'] 
hh = wt_dat['hub_height']
rd =wt_dat['rotor_diameter']  # rotor diameter
rp = wt_dat['performance']['rated_power'] # rated power

# Recreate power curve from power coefficient curve
factor = 0.5 * 1.225* np.pi *(rd/2)**2
p_ws = wt_dat['performance']['Cp_curve']['Cp_wind_speeds']
p = [cp * factor*u**3 for cp, u in zip(wt_dat['performance']['Cp_curve']['Cp_values'], p_ws)]

ct = wt_dat['performance']['Ct_curve']['Ct_values']
ct_ws = wt_dat['performance']['Ct_curve']['Ct_wind_speeds']

cut_in = min(p_ws) # cut-in wind speed
cut_out = max(p_ws) # cut-out wind speed

# Interpolate the data on a finer discretization
int_speeds = np.linspace(0, cut_out+5, 10000)  
ps_int = np.interp(int_speeds, p_ws, p)
cts_int = np.interp(int_speeds, ct_ws, ct)
ps_int[ (int_speeds < cut_in) | (int_speeds > cut_out) ] = 0 # Set power to zero outside of operating range

# Load wind speed and wind direction time series
ws_time = resource_dat['wind_resource']['wind_speed']
wd_time = resource_dat['wind_resource']['wind_direction']
ti = resource_dat['wind_resource']['turbulence_intensity']['data']

time_stamp = np.arange(len(wd_time))/24 

#%%
# ------
# Run the wind farm simulation
# ------

# Define turbines, site and wind farm objects in pywake
windTurbines = WindTurbine(name='NREL 5 MW', diameter=rd, hub_height=hh, 
                      powerCtFunction=PowerCtTabular(int_speeds, ps_int, power_unit='W', ct=cts_int))

# Recreate the wind speed and wind direction vectors corresponding to the wind rose
wd = np.linspace(0, 360, len(A['data'])+1)[:-1]
ws = p_ws

# Create the object representing the site and the wind farm model 
site = XRSite(ds=xr.Dataset(data_vars= dict(P=1)))
wf_model =  Bastankhah_PorteAgel_2014(site, windTurbines, turbulenceModel=None, k=0.05, rotorAvgModel=RotorCenter()) # Bastankah model

# Run the time simulation 
sim_res_time = wf_model(x, y, # wind turbine positions
                        wd=wd_time, # Wind direction time series
                        ws=ws_time, # Wind speed time series
                        time=time_stamp, # time stamps                         
                        TI=ti, # turbulence intensity time series
                  )

power_wind = sim_res_time.Power.sum(['wt']).values.tolist()
wf_aep = sum(power_wind)
wf_cf = np.mean(power_wind)/np.max(power_wind)


print('Wind Farm AEP: {:.1f} GWh'.format(wf_aep*1e-9))
print('Wind Farm capacity factor: {:.1%}'.format(wf_cf))


if plot_figures:
      fig, ax = plt.subplots(1,1, figsize = (10,3))
      sim_res_time.Power.sum(['wt']).plot(color = 'tab:green', ax=ax)
      ax.set_ylabel('Power [W]')
      ax.set_xlabel('Time [days]')

#%%
# ------
# Load the solar pv farm and solar resource data
# ------

# Define location
latitude = hpp_dat['site']['latitude']
longitude = hpp_dat['site']['longitude']
tz = 'Etc/GMT+1'  # Timezone for the location

# Import solar resource data
time = [ datetime.strptime(x, '%Y-%m-%dT%H:%M:%S+00:00Z')  for x in resource_dat['solar_resource']['time']]
ghi = resource_dat['solar_resource']['ghi']
dhi = resource_dat['solar_resource']['dhi']
dni = resource_dat['solar_resource']['dni']

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

weather = pd.DataFrame(index = time, data={'ghi': ghi, 'dhi': dhi, 'dni': dni})

location = Location(latitude, longitude, tz)

system = PVSystem(surface_tilt=tilt, surface_azimuth=surface_azimuth,
                  module_parameters=module_parameters,
                  inverter_parameters=inverter_parameters,
                  temperature_model_parameters=temperature_model_parameters,
                  strings_per_inverter=strings_per_inverter, modules_per_string = modules_per_string)

# Create a model chain object. Models chosen arbitrarily to work with the available data.
model_chain = ModelChain(system, location, dc_model='pvwatts', aoi_model='no_loss', ac_model='pvwatts')

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
      ax.plot(time, power_solar)
      ax.set_ylabel('Power [W]')
      ax.set_xlabel('Time [days]')

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

# IEA Wind Reference Hybrid Power Plant Design (v0.1)

This repository contains the data describing the IEA Wind reference hybrid power plant design (v0.1). The power plant is designed for a site in Denmark and consists of a 325 MW wind farm, 400 MW solar PV farm, and a 150 MW / 300 MWh battery storage system, co-located behind a grid connection with a capacity of 300 MW.

The wind farm consists of 65 NREL 5 MW wind turbines arranged in a grid layout. The wind farm description includes the electrical cable layout. The solar PV farm consists of 59 solar PV systems, each with a capacity of 6.8 MWp. The storage system is made of 60 Lithium-Ion-Phosphate units of capacity 2.5MW / 5.015 MW each. 

The power plant is described in YAML files adapted from the [WindIO](https://github.com/IEAWindSystems/windIO) format, with additional fields and files to describe the solar resource, solar PV farm, and battery system. The repository provides example scripts to load the data and generate time series of wind and solar power using [PyWake](https://topfarm.pages.windenergy.dtu.dk/PyWake/), [FLORIS](https://natlabrockies.github.io/floris/index.html), [pvlib](https://pvlib-python.readthedocs.io/en/stable/index.html) and [H2Integrate](https://h2integrate.readthedocs.io/).

The following data are available:

- Wind Farm
    - Wind farm boundaries
    - Coordinates of turbine positions using the EPSG:25831 coordinate reference system.
    - Coordinates of the substation using EPSG:25831 coordinate reference system.
    - Layout of the electrical collection network, including cable characteristics. 
    - Turbine characteristics (power coefficient and thrust curve, hub height, rotor diameter)
    - Wind resource:
        - Wind rose at hub height and shear
        - Time series of wind speed and direction at hub height for the year 2022
- Solar PV Farm
    - Characteristics of the solar PV module
    - Number of modules per strings
    - Number of strings per inverter
    - Inverter capacity and efficiency
    - Tilt and surface azimuth
    - Number of strings per combiner boxes and per MPP inputs
    - Solar resource, with the time series of GHI, DHI, DNI, and temperature at 2 m for the year 2022
- Battery
    - Round-trip efficiency (nominal and at 0.25C) at the DC connection
    - Cycling lifetime in full-load cycles
    - Power and energy capacity
    - Efficiency of the power conditioning unit
- Site
    - Latitude and longitude
    - Grid connection capacity

## Acknowledgment and related work

This repository documents the work conducted by work package 2 of the [IEA Wind Task TCP 50](https://iea-wind.org/task50) (Phase I) on reference designs for hybrid power plants. 

The design was generated with the help of the following open-source tools:
- [HyDesign](https://topfarm.pages.windenergy.dtu.dk/hydesign/): an open-source framework for design, control and optimization of utility-scale hybrid power plants, developed by the Technical University of Denmark.
- [H2Integrate](https://h2integrate.readthedocs.io/en/latest/intro.html): an open-source python package for modeling and designing hybrid energy systems producing electricity, hydrogen, ammonia, steel, and other products, developed by the National Laboratories of the Rockies.
- [ard](https://github.com/NLRWindSystems/Ard), a wind farm optimization suite for wind energy, built for modular, gradient-enabled multi-disciplinary and multi-fidelity optimizations, and developed by the National Laboratories of the Rockies.


## Contact

For questions about the data or repository, please reach out to: Jenna Iori - j.iori@tudelft.nl


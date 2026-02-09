This repository contains the code for the paper Guo, Zhou, Qi, Fan (2026) [*Conditional effects of cross-product substitution on systemic risk in multilayer food trade networks5*].

The function ‘shock_response’ simulates shock response and cascade starting with a change (decline) in production of a particular product in the target country

The function ‘sim_diagnostics’ takes the output of shock_response and tests whether it respects the mass balance equation

The script ‘multiNet_simShock_Script’ simulates shock cascades under different parameter configurations and network data

The script ‘multiNet_simAnalysis_script’ validates mass balance equations and calculates impacts (consumption deficits) of single/multiple shocks on single country-product, country, product layer and overall network

The input folder contains a list of country information, a list of crop information, as well as data for 2 two-layer networks and 2 three-layer networks.
The network data includes production, reserves, consumption, per capita GDP, and trade data spanning 31 years from 1993 to 2023, along with substitutability correlation matrices (where a value of 1 indicates that the row product can be used to substitute for the column product in that country).

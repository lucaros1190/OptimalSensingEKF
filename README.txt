# Project Title

**Author:** Benhamouche Ouassim and Luca Rossini

## Project Structure

- **matrix_generator**  
  Creates synthetic data from a past simulation with the addition of multiplicative noise.

- **create_data**  
  Adjusts the wind data to match the size and structure of the parcels.

- **ModelFunctions_opt3**  
  Contains the core model functions and optimization algorithms. All main computational logic is implemented here.

- **Main_ch_opt**  
  Primary script for running the main simulation.

- **Parallel_main**  
  Parallelized version of the main script used to speed up simulations.

- **ParametersInputs_opt**  
  Defines most of the initial parameters and configuration settings used in the simulations.

- **Vary_observation**  
  Alternative main script that allows you to vary the number of observations and select the policy by specifying its name as a string. Useful for running comparative experiments across different policies and observation settings.
# Simulation Model
A simulation model written in Python for the generation of spatial-temporal clinical data on infection outbreaks in hospitals. Source: https://www.nature.com/articles/s41598-023-47296-1

## Installation 
The source code is currently hosted on https://github.com/denissekim/Simulation-Model.

## Execution
To run the simulator, in the terminal go to the folder containing the simulator and run: `python main.py` this will open a prompt asking to choose 1 option: if you choose `1` the simulator is going to run a simulation that may or may not have an outbreak. If you choose `2` the simulator is going to run a simulation with an outbreak. The outcomes are saved in patients.csv and movements.csv.

## Outcomes
After running the program, the following outcomes are obtained:
- patients.csv: dataset containing the information of each patient that was hospitalized during the simulation.
- movements.csv: dataset containing for each step of time, the movements of the patients in the hospital and their state of health according to the SEIRD epidemiological model. 

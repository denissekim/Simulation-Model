# H-Outbreak
H-outbreak is a simulation model designed to generate a reliable spatial-temporal dataset detailing the activity of hospitalized patients and the progression of infection spread within hospitals due to significant bacterial factors. This model integrates a compartmental approach to illustrate the development of bacterial infections, an agent-based model to capture the dynamics and propagation of infections, along with individual actions, and spatial-temporal constraints imposed by the hospital's infrastructure. These constraints are determined by the representation of the hospital's layout, cleaning protocols, and staff schedules. Source: https://www.nature.com/articles/s41598-023-47296-1

## Hospital Structure
The hospital structure is organized into floors, each housing various services and wards. Each service, such as Radiology, Surgery, the ICU, and the ER, is equipped with a specific number of beds, while each ward comprises a set number of rooms. Notably, each room within a ward accommodates two beds.
The hospital that we implemented (but can be changed by the user) has an ER with 20 beds, 3 operating rooms, 5 radiology rooms, 4 wards with 14 rooms each, 3 wards with 10 rooms each, 1 ward with 5 rooms, and an ICU with 10 beds. Each room has 2 beds and there are 212 beds in total.

## Parameters
The input parameters that configure each simulation are:

**Population Parameters**
-	Patients rate: Daily occupancy rate
-	Arrival rate: Daily arrival rate
-	Arrival ER: Daily arrival rate at ER 
-	Occupancy ICU: Occupancy rate of the ICU  
-	Population: Hospital area of influence
-	Age: Patient's age distribution  
-	LOS: Patient's Length of Stay mean
  
**Epidemiological Model Parameters**
-	Arrival_S: Prob. of arrivals in Susceptible state 
-	Arrival_I: Prob. of arrivals in Infected state        
-	Arrival_NS: Prob. of arrivals in Non Susceptible state 
-	Arrival_C: Prob. of arrival in colonized state over the whole population                                 
-	P_pl: Prob. of patient infecting place. It follows a triangular distribution                           
-	P_lp: Prob. of place infecting patient. It follows a triangular distribution
-	P_pp: Prob. of patient infecting patient. It follows a triangular distribution
-	P_CI: Prob. of colonized patient becoming infected. It follows a triangular distribution
-	Incubation_time: Min. and max. incubation period (hours)   
-	P_qr: Prob. of quick recovery. It follows a triangular distribution                            
-	P_lr: Prob. of long recovery. It follows a triangular distribution
-	Treatment_days: Treatment duration. It follows a triangular distribution
-	P_death: Prob. of death
  
**Simulation configuration Parameters**
-	Step_time: Step duration (hours)                              
-	Max_patients_movements: Max. number of patients allowed by service per step                          
-	Max_time_infected: Max. infection duration of each place 
-	Steps_to_infect: Number of steps required to infect a place 

## Main Functions 
-	Patients generation: creation of new patients based on demographic input parameters (ie. Population parameters).
-	Patients movements: each step, the patients’ movements are constrained by a series of spatial-temporal rules to make them clinically realistic.
-	Patients contact: contacts between patients can happen when they share a room or a service. They can also interact with the environment. During these interactions, the infection can spread.
-	Cleaning control: contaminated places are cleaned after a predefined number of steps.

## Installation 
The source code is currently hosted on https://github.com/denissekim/Simulation-Model.

## Execution
To run the simulator, in the terminal go to the folder containing the simulator and run: `python main.py` this will open a prompt asking to choose 1 option: if you choose `1` the simulator is going to run a simulation that may or may not have an outbreak. If you choose `2` the simulator is going to run a simulation with an outbreak. The outcomes are saved in patients.csv and movements.csv.

## Outcomes
After running the program, the following outcomes are obtained:
- patients.csv: dataset containing the information of each patient that was hospitalized during the simulation.
- movements.csv: dataset containing for each step of time, the movements of the patients in the hospital and their state of health according to the SEIRD epidemiological model. 

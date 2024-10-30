""" INPUT PARAMETERS CONFIGURATION """

required_parameters = {
    'n_patients': 1,                  # hospital occupation rate
    'steps': 2190,                      # total duration of the simulation 1095
    'population': 680000,               # total population of the hospital area 170000
    'step_time': 8,                     # duration of each step (hours)
    'init_exposed': 1,                  # number of exposed patients at the beginning of the simulation
    'init_infected': 0,                 # number of infected patients at the beginning of the simulation
    'arrival_rate': 104.05,             # patients arrival rate per day 18.603
    'prob_arrival_ER': 0.7,             # proportion of arrivals through the ER
    'arrival_state_colonized': 0.076,   # proportion of colonized patients in the whole population of the hospital area
    'arrival_state_S': 0.9973429,       # proportion of arrivals in state S
    'arrival_state_I': 0.001563,        # proportion of arrivals in state I
    'arrival_state_NS': 0.0010941,      # proportion of arrivals in state NS
    'prob_p-env_min': 0.14,             # min probability of patient infecting environment
    'prob_p-env_max': 0.9,              # max probability of patient infecting environment
    'prob_p-env_mean': 0.52,            # mean (mode) probability of patient infecting environment
    'prob_env-p_min': 0.3262,           # min probability of environment infecting patient
    'prob_env-p_max': 0.5437,           # max probability of environment infecting patient
    'prob_env-p_mean': 0.435,           # mean probability of environment infecting patient
    'prob_pe_min': 0.18,                # min probability of patient exposure from interacting with infected patient
    'prob_pe_max': 0.3,                 # max probability of patient exposure from interacting with infected patient
    'prob_pe_mean': 0.24,               # mean probability of patient exposure from interacting with infected patient
    'prob_pc-i_min': 0,                 # min probability of colonized patient becoming infected
    'prob_pc-i_max': 0.0227,            # max probability of colonized patient becoming infected
    'prob_pc-i_mean': 0.0114,           # mean probability of colonized patient becoming infected
    'incubation_time_min': 48,          # min incubation time in hours
    'incubation_time_max': 72,          # max incubation time in hours
    'prob_quick_recov_min': 0.0,        # min probability of quick recovery
    'prob_quick_recov_max': 0.23,       # max probability of quick recovery
    'prob_quick_recov_mean': 0.115,     # mean probability of quick recovery
    'prob_long_recov_min': 0.5985,      # min probability of long recovery
    'prob_long_recov_max': 0.9975,      # max probability of long recovery
    'prob_long_recov_mean': 0.7981,     # mean probability of long recovery
    'treatment_days_min': 5,            # min duration of treatment in days
    'treatment_days_max': 15,           # max duration of treatment in days
    'treatment_days_mean': 10,          # mean duration of treatment in days
    'prob_death': 0.027,                # probability of death
    'max_patients_rx': 30,              # max number of patients that can go to radiology
    'max_patients_qx': 15,              # max number of patients that can go to surgery
    'min_steps_rx': 10,                 # number of steps during which a patient hasn't gone to radiology, for being allowed to go again (3-4 days)
    'min_steps_qx': 30,                 # number of steps during which a patient hasn't been to surgery, for being allowed to go again (10 days)
    'max_ward_movements': 2,            # max number of allowed movements inside a same ward
    'max_steps_er_icu': 3,              # number of steps that a patient has to remain in ER or ICU, for being changed to a room
    'max_movements_room': 5,            # max number of patients that can change to rooms
    'occupancy_icu': 0.46               # occupancy rate of the ICU
}

""" GLOBAL VARIABLES """

# MODEL VARIABLES
# patient unique counter
patient_id = 0
# patients' LOS and Ages for plotting (BORRAR? TODO)
patients_LOS = []
patients_ages = []
# dict to save patients' information
patients_log = {}

# HOSPITAL VARIABLES
# actual occupancy rate of the hospital
hosp_occupancy_rate = 0


def init():
    global hosp_occupancy_rate, patient_id, patients_LOS, patients_ages, patients_log
    hosp_occupancy_rate = 0
    patient_id = 0
    patients_LOS = []
    patients_ages = []
    patients_log = {}

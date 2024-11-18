""" SIMULATION MODEL FILE"""
from scipy.stats import bernoulli
from copy import copy
from random import normalvariate
import config
from hospital import *
from classes import *
import matplotlib.pyplot as plt
import plotly.express as px
import plotly.graph_objects as go
import pandas as pd

# Function that initializes the simulator (hospital, SEIRD-NS vector and patients)
def initialize():
    #L = initialize_hospital(20, 10, 8, [14, 10, 14, 10, 14, 10, 14, 5], 3, 5, 2, 4)
    L = initialize_hospital(20, 10, 7, [36, 7, 6, 6, 6, 6, 24], 3, 5, 2, 4) 
    seird_vector, non_susceptibles = initialize_seird()
    P, used_spaces, seird_vector = initialize_patients(L, seird_vector)
    P, seird_vector = expose_hospital(P, seird_vector)

    available_spaces = list(set(L) - set([p.l for p in P]))
    A = {}
    A[TypeLocalization.ER] = []
    A[TypeLocalization.Surgery] = []
    A[TypeLocalization.Radiology] = []
    A[TypeLocalization.Ward] = []
    A[TypeLocalization.Room] = []
    A[TypeLocalization.ICU] = []
    A[TypeLocalization.Bed] = []
    for a in available_spaces:
        A[a.name].append(a)

    # Dictionary with maximum steps that a type of localization will remain infected
    max_time_infected = {TypeLocalization.Surgery: 1, TypeLocalization.Radiology: 3, TypeLocalization.Bed: 1,
                         TypeLocalization.Room: 3, TypeLocalization.ICU: 3, TypeLocalization.Ward: 3,
                         TypeLocalization.ER: 6}

    # Steps that a patient has to be in a place to infect it
    steps_to_infect = {TypeLocalization.Surgery: 1, TypeLocalization.Radiology: 2, TypeLocalization.Bed: 1,
                       TypeLocalization.Room: 2, TypeLocalization.ICU: 2, TypeLocalization.Ward: 2,
                       TypeLocalization.ER: 1}

    # Steps that a patient has to be in a place to be contaminated by it
    steps_to_infect_p = {TypeLocalization.Surgery: 0, TypeLocalization.Radiology: 1, TypeLocalization.Bed: 0,
                         TypeLocalization.Room: 1, TypeLocalization.ICU: 1, TypeLocalization.Ward: 1,
                         TypeLocalization.ER: 1}

    return L, P, A, max_time_infected, steps_to_infect, used_spaces, seird_vector, non_susceptibles, steps_to_infect_p


# Function that initializes the SEIRD-NS vector
def initialize_seird():
    days = round(config.required_parameters['steps'] * config.required_parameters['step_time'] / 24)
    seird = [[0 for _ in range(5)] for _ in range(days)]
    non_susceptibles = [[] for _ in range(days)]
    return seird, non_susceptibles


# Function that initializes the patients
def initialize_patients(localizations, seird_vector):
    # we get all the available beds in the hospital
    bed_indexes = [i for i in range(len(localizations)) if
                   localizations[i].is_localization(TypeLocalization.Bed) and localizations[
                       i].located_in.name != TypeLocalization.Radiology]
    # we get the max amount of beds that have to be occupied at a given time
    config.hosp_occupancy_rate = int(len(bed_indexes) * config.required_parameters['n_patients'])
    # Dictionary with a localization that is being/has been used and number of steps that it remains infected.
    # If it isn't infected, value = -1
    used_spaces = {}

    P = []
    # total days necessary for the arrival of max_occupied_beds patients before the start of the simulation
    days = round(config.hosp_occupancy_rate / config.required_parameters['arrival_rate'])
    # while there isn't enough patients to cover the daily occupancy rate, we add more and we rest the remaining stay of the ones that are already in the hospital,
    # so that not everybody starts in the day 0 of the simulation with the same LOS (with this we prevent a massive discharge at the day 4-5 of the simulation)
    while len(P) < config.hosp_occupancy_rate:
        for day in range(days):
            for p in P:
                p.remaining_stay = p.remaining_stay - 1
            for i in range(round(config.required_parameters['arrival_rate'])):
                if (len(P) < config.hosp_occupancy_rate):
                    # we create a patient
                    P.append(Patient(config.patient_id, None, lognormal_distr_los(), normal_distr_age(), set_sex()))
                    config.patient_id = config.patient_id + 1
                else:
                    break

    # we assign a bed to each patient
    for p in P:
        # we get a random bed
        index = random.choice(bed_indexes)
        l = localizations[index]
        p.l = l
        # we add a susceptible patient to the day 0
        seird_vector[0][0] = seird_vector[0][0] + 1
        bed_indexes.remove(index)
        used_spaces[l] = -1
        # We add the patient p to the log on the day 0
        add_patient_log(p, 0)

    config.patients_ages = [x.age for x in P]
    config.patients_LOS = [x.LOS_initial for x in P]

    return P, used_spaces, seird_vector


# Function that adds exposed patients to the hospital at the beginning of the simulation
def expose_hospital(P, seird_vector):
    # we expose patients at the hospital
    i, j = 0, 0
    while i < required_parameters['init_exposed'] and j < len(P):
        if P[j].l.is_localization(TypeLocalization.Bed):
            P[j].state = 1
            i = i + 1
            j = j + 1
            seird_vector[0][1] = seird_vector[0][1] + 1
    seird_vector[0][0] = seird_vector[0][0] - seird_vector[0][1]

    # we infect patients at the hospital
    # Get init_infected patients in the ICU and infect them
    for i in range(0, required_parameters['init_infected']):
        p = next(x for x in P if x.l.is_located_in(TypeLocalization.ICU))
        p.state = 2
        seird_vector[0][2] = seird_vector[0][2] + 1
    seird_vector[0][0] = seird_vector[0][0] - seird_vector[0][2]

    return P, seird_vector


# For each step from 0 to max_steps, this function simulates the movement of patients in a hospital.
# P: list of patients to work with.
# A: available localizations where a patient can be moved to.
# used_spaces: dict that contains the nº of used and/or contaminated spaces in the hospital.
# max_time_infected: dict with time that a type of place remains infected
# steps_to_infect: min steps needed for a patient to infect a place
# seird_vector: SEIRD information
# non_susceptible: list with current non-susc patients
# steps_to_infect_p: min steps needed for a place to infect a patient
def simulator(P, A, used_spaces, max_time_infected, steps_to_infect, seird_vector, non_susceptibles,
              steps_to_infect_p):
    # information for the movements log (where they are in each step and in which health state)
    tuplas_patients = {}
    # total deceased patients during the simulation
    deceased_patients = []
    # saves the contaminated places to know when they get cleaned (for visualization purpose)
    spaces_to_log = {}
    # total population each day
    pobl_day = []
    # diccionario para guardar los pacientes que se infectaron por día
    infected = {}
    day = 0
    pobl_day.append(len(P))
    # step that falls in the afternoon (8-16hs)
    afternoon = 1
    P_d = []
    # Where do they get infected
    exposition = {'env': 0, 'people': 0}

    for step in range(config.required_parameters['steps']):
        spaces_to_log[step] = []
        if step > 0:
            # if step = multiple of 3, it's a new day
            # update SEIRD-NS vector information, save the population information and start a new day
            if step % 3 == 0:
                seird_vector[day], non_susceptibles[day] = update_seird(seird_vector[day], non_susceptibles[day],
                                                                        new_ns_id, P)
                pobl_day[day] = len(P)

                day = day + 1
                seird_vector[day] = copy(seird_vector[day - 1])
                # we reset the deaths count to show only the deaths by day
                seird_vector[day][4] = 0
                non_susceptibles[day] = copy(non_susceptibles[day - 1])
                pobl_day.append(len(P))

            # if we are in the Afternoon (8hs-16hs): admissions, discharges
            if step == afternoon:
                afternoon = afternoon + 24 / config.required_parameters['step_time']
                P = infect_colonized(P)
                P_a, A, used_spaces, new_ns_id = admit_patients(A, used_spaces, day)
                P.extend(P_a)
                P, A, P_d, actual_discharges, places_discharges = discharge_patients(P, P_d, A, day)

            # in each step, patients return from temporal places to their beds
            P, moved_patients, A, used_spaces = move_patients_temporal_places(P, A, used_spaces)

            # afternoon after moving patients from temporal places: indefinite places
            if step == (afternoon - 24 / config.required_parameters['step_time']):
                P, A, used_spaces = move_patients_beds(P, moved_patients, A, step, tuplas_patients, used_spaces)

        # in each step, patients can recover, die, being exposed and infected. Places can get contaminated and cleaned
        P, seird_vector[day] = recover_patients(P, seird_vector[day])

        discharged_patients, deceased_patients, seird_vector[day] = decease_patients(P, deceased_patients,
                                                                                     seird_vector[day])
        P, A, discharged_patients, actual_discharges, places_discharges = discharge_patients(P, discharged_patients, A, day)

        P, used_spaces, seird_vector[day], exposition = expose(P, step, tuplas_patients, used_spaces, steps_to_infect,
                                                               seird_vector[day], exposition, steps_to_infect_p)

        P, seird_vector[day], infected = infect(P, seird_vector[day], infected, day)

        used_spaces, spaces_to_log[step] = cleaning_control(used_spaces, max_time_infected, A, spaces_to_log[step])

        # we add the step to the time counters
        P, P_d = add_step_time_discharge(P, P_d)

        # we save the information of the patients in this step
        tuplas_patients[step] = []
        for p in P:
            # if the patient is non susceptible, we save state = 5
            if p.non_susceptible:
                t = (p.id, 5, copy(p.l))
            else:
                t = (p.id, p.state, copy(p.l))
            tuplas_patients[step].append(t)
        # if the patient has died in this step, we add it to the historic here, because they aren't in P anymore
        for ad in range(len(actual_discharges)):
            if actual_discharges[ad].state == 4:
                t = (actual_discharges[ad].id, actual_discharges[ad].state, places_discharges[ad])
                tuplas_patients[step].append(t)

    # save information in a CSV of patients movements and contaminated places
    print_csv(config.required_parameters['steps'], tuplas_patients, spaces_to_log)
    # save the patients' information
    print_patients_log()



    plot_simulation_runs(day, seird_vector, non_susceptibles, pobl_day)

    return tuplas_patients, seird_vector, non_susceptibles, pobl_day, infected, exposition


# Function that updates the SEIRD-NS vector
def update_seird(seird_vector, non_susc_vector, new_non_susceptibles, patients):
    # we remove those NS patients that were discharged and add the new admissions
    ids = [p.id for p in patients if p.state == 3 and p.non_susceptible]
    discharged_ns = []
    for ns in non_susc_vector:
        if ns not in ids:
            discharged_ns.append(ns)

    non_susc_vector = list(set(non_susc_vector) - set(discharged_ns))
    non_susc_vector.extend(new_non_susceptibles)

    # we update the seir vector without adding the NS
    s, e, i, r = 0, 0, 0, 0
    for p in patients:
        if p.state == 0:
            s = s + 1
        elif p.state == 1:
            e = e + 1
        elif p.state == 2:
            i = i + 1
        elif p.state == 3 and not p.non_susceptible:
            r = r + 1

    seird_vector[0] = s
    seird_vector[1] = e
    seird_vector[2] = i
    seird_vector[3] = r

    return seird_vector, non_susc_vector


# Function that admits new patients to the hospital
def admit_patients(available_spaces, used_spaces, day):
    # we get the beds from the ER and from Rooms
    available_beds = len(available_spaces[TypeLocalization.Bed])
    ER_beds_indexes = [i for i in range(available_beds) if
                       (available_spaces[TypeLocalization.Bed][i].is_located_in(TypeLocalization.ER))]
    room_beds_indexes = [i for i in range(available_beds) if
                         (available_spaces[TypeLocalization.Bed][i].is_located_in(TypeLocalization.Room))]

    # we get the number of patients that are going to arrive
    n = np.random.poisson(required_parameters['arrival_rate'])
    max_beds_available = len(ER_beds_indexes) + len(room_beds_indexes)
    n = min(n, max_beds_available)

    new_patients = []
    occupied_beds = []
    # new admissions by state -> S:0, I:0, NS:0
    adm_state = {0: 0, 1: 0, 2: 0, 3: 0}
    non_susceptibles_id = []
    for i in range(0, n):
        colonized, non_susceptible = False, False
        # we see if the patient arrives in state S, I or NS
        max_value = (required_parameters['arrival_state_NS'] + required_parameters['arrival_state_S'] +
                     required_parameters['arrival_state_I']) * 10000
        ran = random.randint(1, max_value)

        if ran <= required_parameters['arrival_state_NS'] * 10000:
            state = 3
            non_susceptible = True
        elif ran > required_parameters['arrival_state_NS'] * 10000 and ran <= (
                required_parameters['arrival_state_NS'] + required_parameters['arrival_state_S']) * 10000:
            state = 0
            # we check if the S patient is colonized
            prob_colonized = required_parameters['arrival_state_colonized'] / required_parameters['population']
            tn = truncnorm((0 - prob_colonized) / 0.0001, (1 - prob_colonized) / 0.0001, loc=prob_colonized, scale=0.0001)
            prob_c = tn.rvs()
            r = bernoulli.rvs(prob_c)
            if r == 1:
                colonized = True

        elif ran > (required_parameters['arrival_state_NS'] + required_parameters[
            'arrival_state_S']) * 10000 and ran <= max_value:
            state = 2

        # we see if the patients is admitted through ER or Ward
        r = bernoulli.rvs(required_parameters['prob_arrival_ER'])
        # arrives to ER
        if r == 1 and (len(ER_beds_indexes) > 0):
            index = random.choice(ER_beds_indexes)
            ER_beds_indexes.remove(index)
            p = Patient(config.patient_id, available_spaces[TypeLocalization.Bed][index], lognormal_distr_los(),
                        normal_distr_age(), set_sex(), state, colonized=colonized, non_susceptible=non_susceptible)
            new_patients.append(p)
            add_patient_log(p, day)
            occupied_beds.append(available_spaces[TypeLocalization.Bed][index])
            if state == 3:
                non_susceptibles_id.append(config.patient_id)
            config.patient_id = config.patient_id + 1
            adm_state[state] = adm_state[state] + 1
        # arrives to room
        else:
            if (len(room_beds_indexes) > 0):
                index = random.choice(room_beds_indexes)
                room_beds_indexes.remove(index)
                p = Patient(config.patient_id, available_spaces[TypeLocalization.Bed][index], lognormal_distr_los(),
                            normal_distr_age(), set_sex(), state, colonized=colonized, non_susceptible=non_susceptible)
                new_patients.append(p)
                add_patient_log(p, day)
                occupied_beds.append(available_spaces[TypeLocalization.Bed][index])
                if state == 3:
                    non_susceptibles_id.append(config.patient_id)
                config.patient_id = config.patient_id + 1
                adm_state[state] = adm_state[state] + 1

    # we remove the new occupied beds from A
    if occupied_beds:
        for bed in occupied_beds:
            available_spaces, used_spaces = occupy_localization(bed, available_spaces, used_spaces)

    config.patients_LOS.extend([x.LOS_initial for x in new_patients])
    config.patients_ages.extend([x.age for x in new_patients])

    return new_patients, available_spaces, used_spaces, non_susceptibles_id


# Function that adds a step to the patients' time attributes
def add_step_time_discharge(patients, discharges):
    # we rest the remaining LOS for each patient
    for p in patients:
        if p.state != 2:
            p.remaining_stay = p.remaining_stay - required_parameters['step_time'] / 24

        # if the LOS is reached and the patient isn't infected, they can be discharged
        if p.remaining_stay <= 0:
            if p.state != 2:
                if p not in discharges:
                    discharges.append(p)
            else:
                p.LOS_final = p.LOS_final + required_parameters['step_time'] / 24

        # if the patient is exposed, we add to the incubation time
        if p.state == 1:
            p.incubation_period = p.incubation_period + required_parameters['step_time']

        # if the patient is infected, we add a step counting the duration of the infection and we rest the treatment time
        if p.state == 2:
            p.duration_infection = p.duration_infection + 1
            if p.treatment_days > 0 and p.treatment_days < 1000:
                p.treatment_days = p.treatment_days - required_parameters['step_time'] / 24

        modify_patients_log(p, None)
    return patients, discharges


# Function that recovers infected patients
def recover_patients(patients, seird_vector):
    n_quick, n_long = 0, 0
    for p in patients:
        if p.state == 2:
            # if the patient hasn't been infected for 2-3 days, they can still have a quick recovery
            if p.duration_infection > 3 and p.duration_infection <= 9:
                p_quick = np.random.triangular(required_parameters['prob_quick_recov_min'],
                                               required_parameters['prob_quick_recov_mean'],
                                               required_parameters['prob_quick_recov_max'])
                r = bernoulli.rvs(p_quick)
                # the patient recovers
                if r == 1:
                    p.state = 3
                    n_quick = n_quick + 1
            # if the patient has been infected for more than 3 days, they can have a long recovery
            elif p.duration_infection > 9:
                # if we haven't set up a treatment, we do it now
                if p.treatment_days == 1000:
                    p.treatment_days = np.random.triangular(required_parameters['treatment_days_min'],
                                                            required_parameters['treatment_days_mean'],
                                                            required_parameters['treatment_days_max'])
                    p.max_treatment_days = p.treatment_days
                    # if there are more treatment days than the remaining stay, we increase the remaining stay
                    if p.remaining_stay - p.treatment_days < 0:
                        p.remaining_stay = p.remaining_stay + p.treatment_days
                        p.LOS_final = p.LOS_final + p.treatment_days
                    modify_patients_log(p, None)
                # we check if the patient is recovered
                else:
                    prob = np.random.triangular(required_parameters['prob_long_recov_min'],
                                                required_parameters['prob_long_recov_mean'],
                                                required_parameters['prob_long_recov_max'])
                    r = bernoulli.rvs(prob)
                    # the patient recovers
                    if r == 1:
                        p.state = 3
                        n_long = n_long + 1

    seird_vector[2] = seird_vector[2] - (n_long + n_quick)
    seird_vector[3] = seird_vector[3] + (n_long + n_quick)

    return patients, seird_vector


# Function that see if a patient has died
def decease_patients(patients, deceased_patients, seird_vector):
    # discharge_patients saves the patients that have died in this step
    # deceased_patients accumulates all the patients that have died in the simulation
    discharge_patients = []
    for p in patients:
        # if the patient has been infected for less than 30 days
        if p.state == 2 and p.duration_infection <= 90:
            r = bernoulli.rvs(required_parameters['prob_death'])
            # the patient dies
            if r == 1:
                p.state = 4
                discharge_patients.append(p)
                seird_vector[2] = seird_vector[2] - 1
                seird_vector[4] = seird_vector[4] + 1
    deceased_patients = deceased_patients + discharge_patients
    return discharge_patients, deceased_patients, seird_vector


# Function that discharges patients from the hospital
def discharge_patients(patients, patients_discharge, available_places, day):
    p_d = []
    l_d = []
    total_patients = len(patients)

    infected = []
    for p in patients_discharge:
        # if the patient was exposed at the end of their LOS, but then got infected, the patient remains in the hospital
        if p.state == 2:
            infected.append(p)
        else: #if p.state == 4 or (p.state != 4 and total_patients > normal_distr_occupancy_rate()):
            if p.b:
                available_places = free_localization(p.b, available_places)
                p.del_bed()
            if p.l:
                available_places = free_localization(p.l, available_places)
                l_d.append(p.l)
                p.l = None
            modify_patients_log(p, day)
            p_d.append(p)
            total_patients = total_patients - 1

    patients = list(set(patients) - set(p_d))
    patients_discharge = list(set(patients_discharge) - set(p_d))
    patients_discharge = list(set(patients_discharge) - set(infected))
    return patients, available_places, patients_discharge, p_d, l_d


# Function that moves the patients that are in temporary places
def move_patients_temporal_places(P, A, used_spaces):
    moved_patients = []
    for p in P:
        # 1-step localization movement
        if p.l.is_temporal_localization():
            # the current localization becomes an available space
            A = free_localization(p.l, A)
            if p.b is not None:
                # if patient's bed is in ER, we look for an available room bed
                if p.b.is_located_in(TypeLocalization.ER):
                    bed = set_bed(A)
                    if bed is not None:
                        p.l = bed
                        A, used_spaces = occupy_localization(p.l, A, used_spaces)
                        A = free_localization(p.b, A)
                        p.del_bed()
                    else:
                        p.l = p.b
                # if patient's bed wasn't in ER, he returns there
                else:
                    p.l = p.b

                moved_patients.append(p)

    # we return the patients that we moved now + the patients moved in move_patients_beds()
    return P, moved_patients, A, used_spaces


# Function that moves patients from indefinite places
def move_patients_beds(P, moved_patients, A, step, historic_patients, used_spaces):
    # we rest the patients that have already been moved during this step, so we don't move them twice
    patients_to_move = list(set(P) - set(moved_patients))
    # the remaining patients (those that hasn't been moved) are sorted by their length of stay in the same localization
    ranking = set_patients_ranking(patients_to_move, step, historic_patients)
    patients_2rx = 0
    patients_2qx = 0
    movements_in_service = 0
    movements_to_room = 0

    for p in ranking:
        new_bed_same_service = set_bed(A, p.l)
        new_bed = set_bed(A)
        # if the patient can move to Radiology
        if can_move_to(p, A, historic_patients, step, patients_2rx, config.required_parameters['max_patients_rx'],
                       config.required_parameters['min_steps_rx'], TypeLocalization.Radiology):
            # if the patient is in a bed, we save the bed for them to come back
            if p.l.is_localization(TypeLocalization.Bed) and p.l.located_in.name != TypeLocalization.Radiology:
                p.b = p.l
            # else the current localization becomes an available space
            else:
                A = free_localization(p.l, A)
            p.l = get_location(A, TypeLocalization.Radiology)
            A, used_spaces = occupy_localization(p.l, A, used_spaces)
            patients_2rx = patients_2rx + 1

        # if the patient can move to Surgery
        elif can_move_to(p, A, historic_patients, step, patients_2qx, config.required_parameters['max_patients_qx'],
                         config.required_parameters['min_steps_qx'], TypeLocalization.Surgery):
            # if the patient is in a bed, we save the bed for them to come back
            if p.l.is_localization(TypeLocalization.Bed) and p.l.located_in.name != TypeLocalization.Radiology:
                p.b = p.l
            # else the current localization becomes an available space
            else:
                A = free_localization(p.l, A)
            p.l = get_location(A, TypeLocalization.Surgery)
            A, used_spaces = occupy_localization(p.l, A, used_spaces)
            patients_2qx = patients_2qx + 1

        # if the patient can be moved from bed in the same service
        elif can_move_in_ward(p, movements_in_service,
                              config.required_parameters['max_ward_movements']) and new_bed_same_service:
            A, used_spaces = occupy_localization(new_bed_same_service, A, used_spaces)
            A = free_localization(p.l, A)
            p.l = new_bed_same_service
            if p.b is not None:
                A = free_localization(p.b, A)
                p.del_bed()
            movements_in_service = movements_in_service + 1

        # if the patient is in ER or UCI and can be transferred to a room
        elif can_move_to_room(p, historic_patients, step, config.required_parameters['max_steps_er_icu'],
                              config.required_parameters['max_movements_room'], movements_to_room) and new_bed:
            A, used_spaces = occupy_localization(new_bed, A, used_spaces)
            A = free_localization(p.l, A)
            p.l = new_bed
            if p.b is not None:
                A = free_localization(p.b, A)
                p.del_bed()
            movements_to_room = movements_to_room + 1

        # if the patient is in a room and hasn't changed of ward yet
        elif p.l.is_room_bed() and not p.changed_ward:
            # get random ward that's not the current one
            ward = get_ward(A, p.l.get_ward())
            if ward is not None:
                # get new bed in new service
                new_ward_bed = get_bed(A, ward)
                if new_ward_bed is not None:
                    A, used_spaces = occupy_localization(new_ward_bed, A, used_spaces)
                    A = free_localization(p.l, A)
                    p.l = new_ward_bed
                    p.changed_ward = True
                    if p.b is not None:
                        A = free_localization(p.b, A)
                        p.del_bed()

        # if the patient is in a room or in ER, they can go to ICU if the occupancy rate is lower than the ICU average occupancy rate
        elif get_icu_occupancy(A) < config.required_parameters['occupancy_icu']:
            if p.l.is_located_in(TypeLocalization.ER) or p.l.is_located_in(TypeLocalization.Room):
                bed = next((x for x in A[TypeLocalization.Bed] if x.is_located_in(TypeLocalization.ICU)), None)
                if bed is not None:
                    A = free_localization(p.l, A)
                    A, used_spaces = occupy_localization(bed, A, used_spaces)
                    p.l = bed
                    if p.b is not None:
                        A = free_localization(p.b, A)
                        p.del_bed()

    patients = moved_patients + ranking
    return patients, A, used_spaces


# Function that returns the current occupancy of the ICU
def get_icu_occupancy(available_places):
    # we get the available beds of the ICU
    icu = available_places[TypeLocalization.ICU][0]
    available_beds = [b for b in available_places[TypeLocalization.Bed] if b.located_in is icu]
    # we get the total occupied beds of the ICU
    oc_beds = len(icu.children) - len(available_beds)
    # we return the % of occupancy
    return oc_beds / len(icu.children)


# Function that returns an available bed
def set_bed(available_spaces, patient_localization=None):
    # if no patient's localization is indicated, we assign a free bed following a normal distribution
    if patient_localization is None:
        return normal_choice(available_spaces[TypeLocalization.Bed], True)
    # if the patient's bed is indicated, we search a bed in the same ward
    else:
        ward = patient_localization.get_ward()
        if ward is not None:
            return get_bed(available_spaces, ward)
    return None


# Function that chooses a bed in wards or a new ward following a normal distribution
def normal_choice(options, isbed, mean=None, stddev=None):
    # if the options are beds
    if isbed:
        lst = [x for x in options if x.is_room_bed()]
    # if the options are wards
    else:
        lst = options

    if len(lst) > 0:
        if mean is None:
            # if mean is not specified, use center of list
            mean = (len(lst) - 1) / 2

        if stddev is None:
            # if stddev is not specified, let list be -3 .. +3 standard deviations
            stddev = len(lst) / 6

        while True:
            index = int(normalvariate(mean, stddev) + 0.5)
            if 0 <= index < len(lst):
                return lst[index]
    else:
        return None


# Function that returns an available Bed in a Room that is in the ward
def get_bed(available_spaces, ward):
    for space in available_spaces[TypeLocalization.Bed]:
        if space.is_room_bed() and space.get_ward() is ward:
            return space
    return None


# Function that returns an available space of type type_location
def get_location(available_spaces, type_location):
    # if the localization is Radiology, we return a bed inside
    if type_location == TypeLocalization.Radiology:
        beds_radiology = [bed for bed in available_spaces[TypeLocalization.Bed] if
                          bed.located_in.name == TypeLocalization.Radiology]
        return beds_radiology[0]
    return available_spaces[type_location][0]


# Function that returns an available random ward following a normal distribution that isn't old_ward
def get_ward(available_spaces, old_ward):
    wards = [x for x in available_spaces[TypeLocalization.Ward] if x is not old_ward]
    return normal_choice(wards, False)


# Function that sorts the patients according to the nº of steps they've been in the same place
def set_patients_ranking(patients, step, historic_patients):
    ranking = {}
    for p in patients:
        ranking[p] = 0
        if step > 0:
            s = step - 1
            same_localiz = True

            while s >= 0 and same_localiz:
                # we search for the patient information at step s
                i = next((i for i, v in enumerate(historic_patients[s]) if v[0] == p.id), None)
                if i is not None:
                    # if they are in the same spot, we add a point
                    if p.l.id is historic_patients[s][i][2].id:
                        ranking[p] = ranking[p] + 1
                        s = s - 1
                    # else we stop looking back
                    else:
                        same_localiz = False
                # if we haven't found the patient at step s, we continue
                else:
                    same_localiz = False

    # we sort the ranking by the points accumulated for the patients
    ranking = dict(sorted(ranking.items(), key=lambda item: item[1], reverse=True))
    # we return the list of patients associated to that ranking
    return list(ranking.keys())


# Function that returns bool indicating whether a patient can be moved to a specific type of location or not
# patient = current patient
# historic_patients = list of previous steps to see if the patient had been to the location before
# localization_counts = dict of number of available spaces for each type of localization
# step = current step,
# patients_2location = number of patients moved to that location in the current step
# threshold_patients = max nº of patients that can be in that location at the same time
# threshold_steps = max nº of steps to look backwards
# location = type of location we are asking for
def can_move_to(patient, available_spaces, historic_patients, step, patients_2location, threshold_patients,
                threshold_steps, location):
    # if we haven't reach the top of patients moved to the new localization and the new localization has available spaces
    if patients_2location < threshold_patients and is_available(available_spaces, location):
        if step > 0:
            s = step - 1
            is_location = False
            while s >= 0 and threshold_steps > 0 and not is_location:
                # for each step we find the index of the patient in the history
                i = next((i for i, v in enumerate(historic_patients[s]) if v[0] == patient.id), None)

                # if the patient wasn't in the location at step s, we go on
                if i is not None and not historic_patients[s][i][2].is_localization(location):
                    s = s - 1
                    threshold_steps = threshold_steps - 1
                # else, if the patient was there, we continue
                else:
                    is_location = True

            if not is_location:
                return True
        else:
            return True

    return False


# Function that indicates if there are available places of a type or not
def is_available(A, location):
    if location == TypeLocalization.Radiology:
        beds_radiology = [bed for bed in A[TypeLocalization.Bed] if bed.located_in.name == TypeLocalization.Radiology]
        return len(beds_radiology) > 0
    return len(A[location]) > 0


# Function that indicates whether a patient can move to another bed of the same ward
def can_move_in_ward(patient, movements_in_ward, max_ward_movements):
    if patient.l.is_room_bed() and movements_in_ward < max_ward_movements:
        return True
    return False


# Function that indicates whether a patient can move from ER/ICU to a room
# patient = current patient
# historic_patients = list of previous steps to see if the patient had been to the location before
# step = current step
# threshold_steps = max nº of steps to look backwards
# threshold_move = max nº of allowed movements from ER/ICU to Room in this step
# movements_to_room = nº of movements made to rooms in this step
def can_move_to_room(patient, historic_patients, step, threshold_steps, threshold_move, movements_to_room):
    # if there are still movements available to rooms in this step and the patient is in ER or in ICU
    if movements_to_room < threshold_move and (
            patient.l.is_located_in(TypeLocalization.ER) or patient.l.is_located_in(TypeLocalization.ICU)):
        # we check if the patient has been in ER or ICU for the past threshold_steps steps:
        if step > 0:
            s = step - 1
            while s >= 0 and threshold_steps > 0:
                i = next((i for i, v in enumerate(historic_patients[s]) if v[0] == patient.id), None)
                if i is not None and patient.l.id is historic_patients[s][i][2].id:
                    s = s - 1
                    threshold_steps = threshold_steps - 1
                else:
                    return False

            if threshold_steps == 0:
                return True
    return False


# Function that exposes patients and places
def expose(patients, step, historic_patients, used_spaces, steps_to_infect, seird_vector, exposition,
           steps_to_infect_p):
    locations_checked = []
    for p in patients:
        ## Expose patient ##
        # if the patient isn't infected and their location is, the location infects the patient with probability p_env
        if p.state == 0:
            ## Expose patient by another patient ##
            # if the patient is in a bed, we search for the neighbor beds and see if there is an infected neighbour
            # if the neighbor is colonized, they can infect the other with the same probability
            if p.l.is_localization(TypeLocalization.Bed) and p.l.located_in.name != TypeLocalization.Radiology:
                beds = p.l.adjacent
                if beds is not None:
                    neighbours = [patient for i, patient in enumerate(patients) if (
                                patient.state == 2 or (patient.state == 0 and patient.colonized)) and patient.l in beds]
                    if neighbours is not None:
                        exposed = False
                        i = 0
                        while i < len(neighbours) and not exposed:
                            prob_pe = np.random.triangular(config.required_parameters['prob_pe_min'],
                                                           config.required_parameters['prob_pe_mean'],
                                                           config.required_parameters['prob_pe_max'])
                            r = bernoulli.rvs(prob_pe)
                            i = i + 1
                            # expose
                            if r == 1:
                                p.state = 1
                                seird_vector[0] = seird_vector[0] - 1
                                seird_vector[1] = seird_vector[1] + 1
                                # we save that the patient has been exposed by contact with other patient
                                exposition['people'] = exposition['people'] + 1
                                exposed = True

            # if the patient is still susceptible, we see if the place is infected
            if p.state == 0:
                if (p.l.infected and steps_in_same_localization(p, historic_patients, step) > steps_to_infect_p[
                    p.l.name]) or (p.l.is_localization(
                        TypeLocalization.Bed) and p.l.is_room_infected() and steps_in_same_localization(p,
                                                                                                        historic_patients,
                                                                                                        step) >
                                   steps_to_infect_p[p.l.located_in.name]):  # and riesgo de contagio
                    prob_e = np.random.triangular(config.required_parameters['prob_env-p_min'],
                                                  config.required_parameters['prob_env-p_mean'],
                                                  config.required_parameters['prob_env-p_max'])
                    r = bernoulli.rvs(prob_e)
                    # expose
                    if r == 1:
                        p.state = 1
                        seird_vector[0] = seird_vector[0] - 1
                        seird_vector[1] = seird_vector[1] + 1
                        # we save that the patient has been exposed by contact with environment
                        exposition['env'] = exposition['env'] + 1

        # if the patient is infected or colonized
        elif p.state == 2 or (p.state == 0 and p.colonized):
            ## Infect environment ##
            # if the patient is infected, and their location isn't, the patient infects the location with a bernoulli trial
            # if the location is a bed, it doesn't matter if the bed was infected or not, while it contains an infected patient, the bed will remain infected.
            # we assume that the surgeries don't get infected, as they are maintained sterile
            if steps_in_same_localization(p, historic_patients, step) >= steps_to_infect[p.l.name] and (
                    not p.l.infected or p.l.is_localization(TypeLocalization.Bed)) and not p.l.is_localization(
                    TypeLocalization.Surgery):
                if (p.treatment_days == 1000 or (p.max_treatment_days - p.treatment_days) < 3):
                    prob = np.random.triangular(0.6, 0.75, 0.9)
                    # prob = np.random.triangular(config.required_parameters['prob_p-env_mean'], np.mean([config.required_parameters['prob_p-env_mean'], config.required_parameters['prob_p-env_max']]), config.required_parameters['prob_p-env_max'])
                else:
                    prob = np.random.triangular(0.14, 0.37, 0.6)
                    # prob = np.random.triangular(config.required_parameters['prob_p-env_min'], 0.37, config.required_parameters['prob_p-env_mean'])

                r = bernoulli.rvs(prob)
                locations_checked.append(p.l.id)
                if p.l.is_localization(TypeLocalization.Bed):
                    locations_checked.append(p.l.located_in.id)
                # infects
                if r == 1:
                    p.l.infected = True
                    # starts counting the steps that the location remains infected
                    used_spaces[p.l] = 0

    return patients, used_spaces, seird_vector, exposition


# Function that infects patients
def infect(patients, seird_vector, infected, day):
    for p in patients:
        # if the patient has been exposed and the incubation time has passed, they are infected
        if p.state == 1 and p.incubation_period >= p.max_incubation_period:
            p.state = 2
            seird_vector[1] = seird_vector[1] - 1
            seird_vector[2] = seird_vector[2] + 1
            # we save the patient
            if day not in infected:
                infected[day] = []
            infected[day].append(p.id)

    return patients, seird_vector, infected


# Function that makes colonized patients develop infection
def infect_colonized(patients):
    for p in patients:
        if p.state == 0 and p.colonized:
            prob = np.random.triangular(config.required_parameters['prob_pc-i_min'],
                                        config.required_parameters['prob_pc-i_mean'],
                                        config.required_parameters['prob_pc-i_max'])
            r = bernoulli.rvs(prob)
            # if a colonized patient develops the disease, they become infected
            # to respect the seird transitions, we transform the patient to expose, with no incubation time, so that they inmediatelly become infected
            # infect
            if r == 1:
                p.state = 1
                p.incubation_period = p.max_incubation_period
    return patients


# Function that cleans the places
def cleaning_control(used_spaces, max_time_infected, available_spaces, spaces_to_log):
    # we compare the steps that a location has been infected with the max nº of steps that it has to be infected
    new_rooms_infected = []
    for space in used_spaces:
        # if the space is infected
        if used_spaces[space] != -1:
            # we add a new step
            used_spaces[space] = used_spaces[space] + 1
            # if it is a bed, it can infect the room/ER/ICU it's in
            if space.is_localization(TypeLocalization.Bed) and not space.is_room_infected():
                new_rooms_infected.append(space.located_in)

        # if the location has been infected longer than the max time determined for that type of space
        if used_spaces[space] > max_time_infected[space.name]:
            # the location is cleaned
            space.infected = False
            used_spaces[space] = -1
            if space in available_spaces[space.name]:
                spaces_to_log.append(space)

    # we add the rooms that have been infected in this step
    for room in new_rooms_infected:
        used_spaces[room] = 0
        room.infected = True

    return used_spaces, spaces_to_log


# Function that calculates how many steps a patient has been in the same place
def steps_in_same_localization(patient, historic_patients, step):
    # accum_steps starts with 1 for the present step
    accum_steps = 1
    if step > 0:
        s = step - 1
        same_location_as_current = True
        # for each step we find the index of the patient in the history
        while s >= 0 and same_location_as_current:
            i = next((i for i, v in enumerate(historic_patients[s]) if v[0] == patient.id), None)

            # if the patient was in the location at step s, we go on
            if i is not None and historic_patients[s][i][2] == patient.l:
                s = s - 1
                accum_steps = accum_steps + 1
            # else we continue
            else:
                same_location_as_current = False

    return accum_steps


# Function that frees a place
def free_localization(localization, available_spaces):
    if localization not in available_spaces[localization.name]:
        available_spaces[localization.name].append(localization)
    return available_spaces


# Function that occupies a place
def occupy_localization(localization, available_spaces, used_spaces):
    available_spaces[localization.name].remove(localization)
    if localization not in used_spaces.keys():
        used_spaces[localization] = -1
    return available_spaces, used_spaces


""" Print to logs """


# output: step | p_id | p_state | localization | p_id | p_state | localization | ...
def print_csv(max_steps, patients_record, spaces_to_log):
    fichero = open('movements.csv', 'w')
    for step in range(0, max_steps):
        s = str(step)
        new_list = sorted(patients_record[step])
        for t in new_list:
            if t[2] is not None:
                s = s + ';' + str(t[0]) + ';' + str(t[1]) + ';' + t[2].to_string()  # s = s + 'p'+ str(t[0]) + ';' + str(t[1]) + ';' + t[2].to_string() + ';'
            else:
                s = s + ';' + str(t[0]) + ';' + str(t[1]) + ';Null'
                
        if len(spaces_to_log[step]) > 0:
            s = s + ';' + "places" + ";" + str(spaces_to_log[step][0].id)
            
            for i in range(1, len(spaces_to_log[step])):
                s = s + ";" + str(spaces_to_log[step][i].id)
            
        s = s + '\n'
        fichero.write(s)


# Output: id | age | sex | LOS | incubation_time | infection_time | treatment_days | treatment | antibiotic | microorganism | adm_day | last_day | died | LOS_final | colonized | non_susceptible
def print_patients_log():
    fichero = open('patients.csv', 'w')
    for id, patient_data in config.patients_log.items():
        indice = 1
        s = str(id)
        values = "nextval('paciente_id_paciente_seq')"
        for v in patient_data.values():
            s = s + ';' + str(v)
            if v is not None:
                if indice == 3:
                    values = values + ",'" + str(v) + " days'"
                elif indice == 4:
                    values = values + ",'" + str(v) + " hours'"
                elif indice == 5:
                    values = values + ",'" + str(v * 8) + "'"
                else:
                    values = values + ",'" + str(v) + "'"
            indice = indice + 1
        s = s + '\n'
        fichero.write(s)


# Function that saves the patient information in the patients_log
def add_patient_log(patient, adm_day):
    config.patients_log[patient.id] = {'age': patient.age, 'sex': patient.sex, 'LOS': patient.LOS_initial,
                                       'incub_time': patient.max_incubation_period,
                                       'infect_time': patient.duration_infection,
                                       'treatment_days': patient.treatment_days, 'treatment': patient.treatment,
                                       'antibiotic': patient.antibiotic, 'microorg': patient.microorganism,
                                       'adm_day': adm_day,
                                       'last_day': -1, 'died': False, 'LOS_final': patient.LOS_final,
                                       'colonized': patient.colonized, 'non_susceptible': patient.non_susceptible}


# Function that modifies an existing entry of the patients_log
def modify_patients_log(patient, last_day):
    # if the patient has died or has been discharged
    if last_day is not None:
        # if the patient's LOS has increased due to the infection
        config.patients_log[patient.id]['LOS_final'] = last_day - config.patients_log[patient.id]['adm_day'] #patient.LOS_final
        config.patients_log[patient.id]['infect_time'] = patient.duration_infection
        config.patients_log[patient.id]['treatment'] = patient.treatment
        config.patients_log[patient.id]['antibiotic'] = patient.antibiotic
        config.patients_log[patient.id]['microorg'] = patient.microorganism
        config.patients_log[patient.id]['last_day'] = last_day
        config.patients_log[patient.id]['died'] = (patient.state == 4)
    # if no day was passed, we only modify the treatment days
    else:
        config.patients_log[patient.id]['treatment_days'] = patient.treatment_days
        config.patients_log[patient.id]['LOS_final'] = patient.LOS_final


""" Epidemiological indicators"""
# Function that returns the incidence density
def get_incidence_density(historic_patients):
    # saves the time each patient is in state S, E or I (is at risk and gives information)
    time_patients = {}
    infected = 0
    infected_week = []
    added_infected_patients = []
    day = 0
    week = 0
    for step in range(required_parameters['steps']):
        # we observe the patients' data only once per day
        day = day + 1
        # for each patient that day
        for p in historic_patients[step]:
            # we only consider those patients that are S, E or I
            if p[1] == 0 or p[1] == 1:
                # if there isn't any entry of time for the patient of id p[0]
                if p[0] not in time_patients:
                    time_patients[p[0]] = 0
                # we add 1 day to the time of the patient
                time_patients[p[0]] = time_patients[p[0]] + 1

            # we add the infected patient
            if p[1] == 2 and p[0] not in added_infected_patients:
                infected = infected + 1
                added_infected_patients.append(p[0])

    # we get the total time that patients are at risk of being infected in weeks
    week_patient = 0
    for tp in time_patients.values():
        week_patient = week_patient + tp  # this is in days

    week_patient = week_patient / 7  # this is in weeks

    ID = []
    for w in range(week):
        ID.append(infected_week[w] * 10000 / week_patient)

    # we get the ID for the total execution
    density = infected * 10000 * 3 / sum(time_patients.values())
    print('Incidence density (cases per 10000 days-person)', round(density, 2), ', cases:', infected, ', d-p:', round(sum(time_patients.values()) / 3, 2))
    #print(round(density, 2), ',', round(sum(time_patients.values()) / 3, 2), ',', infected)

    return density, infected


# Function that returns the incidence rate as is calculated in the literature
def get_incidence_rate(infected):
    days = 0
    for p_id, patient_data in config.patients_log.items():
        if patient_data['last_day'] != -1:
            last_day = patient_data['last_day']
        else:
            last_day = required_parameters['steps'] * required_parameters['step_time'] / 24
        days = days + (last_day - patient_data['adm_day'] + 1)

    rate = infected * 10000 / days

    print('Tasa de incidencia (casos por 10000 días)', round(rate, 2))
    return rate

def plot_simulation_runs(days, seird_vector, non_susceptibles, population_day):
    run = []
    s = [row[0] for row in seird_vector]
    e = [row[1] for row in seird_vector]
    i = [row[2] for row in seird_vector]
    r = [row[3] for row in seird_vector]
    d = [row[4] for row in seird_vector]
    ns = [len(row) for row in non_susceptibles]
    run = [s, e, i, r, d, ns, population_day]


    t = np.linspace(0, days+1, days+1)  # 160 puntos entre 0 y 160

    exec = 1
    
    # we separate the susceptibles and the rest to escalate them properly
    fig = plt.figure(figsize=(40, 10))
    ax = fig.add_subplot(222, axisbelow=True)
    simS = run[0]
    simPopulation = run[6]
    ax.plot(t, simS, lw=0.75, label='Susceptible', color='#3197f7', linestyle='dashed')
    ax.plot(t, simPopulation, lw=0.75, label="Total population", color='black', linestyle='dashed')
    ax.set_xlabel('Time (days)')
    ax.set_ylabel('Population')
    # Major ticks every 20, minor ticks every 5
    major_ticks = np.arange(0, days, 25)
    minor_ticks = np.arange(0, days, 5)
    ax.set_xticks(major_ticks)
    ax.set_xticks(minor_ticks, minor=True)
    ax.set_yticks(np.arange(90, 170, 20)) 
    ax.set_yticks(np.arange(90, 170, 5), minor=True)
    ax.grid(which='minor', alpha=0.2)
    ax.grid(which='major', alpha=0.5)
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.savefig(str(exec) + 'Simulator_S+TP.png', bbox_inches='tight')
    plt.tight_layout()

    fig = plt.figure(figsize=(40, 10))
    ax = fig.add_subplot(222, axisbelow=True)
    simE = run[1]
    simI = run[2]
    simR = run[3]
    simD = run[4]
    simNS = run[5]
    ax.plot(t, simE, lw=0.75, label='Exposed', color='orange', linestyle='dashed')
    ax.plot(t, simI, lw=0.75, label='Infected', color='#2fab48', linestyle='dashed')
    ax.plot(t, simR, lw=0.75, label='Recovered', color='#914fbc', linestyle='dashed')
    ax.plot(t, simD, lw=0.75, label='Deceased', color='red', linestyle='dashed')
    ax.plot(t, simNS, lw=0.75, label="Non-susceptible", color='#25249b', linestyle='dashed')
    ax.set_xlabel('Time (days)')
    ax.set_ylabel('Population')
    # Major ticks every 20, minor ticks every 5
    major_ticks = np.arange(0, days, 10)
    minor_ticks = np.arange(0, days, 1)
    ax.set_xticks(major_ticks)
    plt.xticks(fontsize=5)
    ax.set_xticks(minor_ticks, minor=True)
    ax.set_yticks(np.arange(0, 20, 5))
    ax.set_yticks(np.arange(0, 20, 1), minor=True)
    ax.grid(which='minor', alpha=0.2)
    ax.grid(which='major', alpha=0.5)
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    plt.savefig(str(exec) + 'Simulator_E+I+R+D+NS.png', bbox_inches='tight')
    plt.tight_layout()
    exec = exec +1

    data = {
        'Exposed': simE,
        'Infected': simI,
        'Recovered': simR,
        'Deceased': simD,
        'Non-susceptibles': simNS,
        'time': t
    }

    figure = go.Figure()
    figure.add_trace(go.Scatter(x=t, y=simE, name='Exposed', mode='lines', line_color="orange"))
    figure.add_trace(go.Scatter(x=t, y=simI, name='Infected', mode='lines', line_color='#2fab48'))
    figure.add_trace(go.Scatter(x=t, y=simR, name='Recovered', mode='lines', line_color='#914fbc'))
    figure.add_trace(go.Scatter(x=t, y=simD, name='Deceased', mode='lines', line_color='red'))
    figure.add_trace(go.Scatter(x=t, y=simNS, name='Non-susceptibles', mode='lines', line_color='#25249b'))
    figure.update_xaxes(title_text="Time (days)")
    figure.update_yaxes(title_text="Population")
    figure.layout.plot_bgcolor='white'
    figure.layout.yaxis.gridcolor='grey'
    figure.layout.xaxis.gridcolor='grey'
    #figure.update_layout(
    #    grid = {'rows': grid_rows, 'columns': grid_cols, 'pattern': "independent"}
    #)
    figure.show()
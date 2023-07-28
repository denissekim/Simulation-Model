# -*- coding: utf-8 -*-
from simulation import *



# Function that runs a simulation
# outbreak: bool indicating if there must be at least one outbreak or not
def run_simulation(with_outbreak):
    # reset the global parameters of the simulation
    config.init()

    if with_outbreak:
        while with_outbreak:
            L, P, A, max_time_infected, steps_to_infect, used_spaces, seird_vector, non_susceptibles, steps_to_infect_p = initialize()
            tuplas_patients, seird_vector, non_susceptibles, population_day, infected, exposition = simulator(P, A, used_spaces, max_time_infected, steps_to_infect, seird_vector, non_susceptibles, steps_to_infect_p)
            # if there's an outbreak, we stop
            outbreak, firstDay, outbreak_days, outbreaks = has_outbreak(infected)
            # reset the global parameters of the simulation
            config.init()
            if outbreak:
                with_outbreak = False

    else:
        L, P, A, max_time_infected, steps_to_infect, used_spaces, seird_vector, non_susceptibles, steps_to_infect_p = initialize()
        tuplas_patients, seird_vector, non_susceptibles, population_day, infected, exposition = simulator(P, A, used_spaces, max_time_infected, steps_to_infect, seird_vector, non_susceptibles, steps_to_infect_p)



# Function that indicates if there has been an outbreak in the simulation or not
def has_outbreak(infected_people):
    # returns first key (day) if the dict is ordered
    first_day = next(iter(infected_people), None)
    # if no outbreak
    if first_day is None:
        return False, first_day, [], 0

    outbreak_days = []
    first_outbreak = -1
    outbreaks = 0
    continue_outbreak = False
    for day, patients in infected_people.items():
        # if we are in the first day that appears an infected
        if first_day == day:
            # if on the first day there are 2 or more infected, it's an outbreak
            if len(patients) < 3:
                continue
            else:
                first_outbreak = day
                outbreak_days.append(day)
        # if we haven't found the first outbreak
        elif first_outbreak == -1:
            # if there hasn't been 7 days between 2 contagions
            # we save the day that the outbreak happens, not including the 1º day (because there wasn't any outbreak yet)
            if day - first_day <= 7:
                first_outbreak = day
                outbreak_days.append(day)
            elif len(patients) >= 3:
                first_outbreak = day
                outbreak_days.append(day)
                outbreaks = outbreaks + 1
            else:
                first_day = day
        # if we already have the first outbreak, we search for more outbreaks
        elif day - last_outbreak <= 7:
            outbreak_days.append(day)
            if last_outbreak in outbreak_days:
                continue_outbreak = True
            else:
                continue_outbreak = False

        else:
            continue_outbreak = False

            if len(patients) >= 3:
                outbreak_days.append(day)

        if not continue_outbreak and day in outbreak_days:
            outbreaks = outbreaks + 1

        last_outbreak = day

    if len(outbreak_days) > 0:
        return True, first_outbreak, outbreak_days, outbreaks

    return False, first_outbreak, outbreak_days, 0


if __name__ == '__main__':
    config.init()

    option = input("Choose one option: \n1 - Run a simulation \n2 - Run a simulation with at least one outbreak\n")
    with_outbreak = False if option == '1' else True
    print("Running...")

    run_simulation(with_outbreak)

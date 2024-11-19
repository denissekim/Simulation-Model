""" CLASSES PATIENT AND LOCALIZATION """

from enum import Enum
import random
from config import required_parameters
import numpy as np
from scipy.stats import truncnorm, bernoulli


class Patient:
    def __init__(self, id, localization, LOS, age, sex, state=0, record=None, bed=None, remaining_stay=None, colonized=False, non_susceptible=False):
        self.id = id
        self.l = localization
        self.r = record
        self.age = age
        self.sex = sex
        self.b = bed
        self.changed_ward = False
        # state = [0, 1, 2, 3, 4, 5] -> [S, E, I, R, D, NS]
        self.state = state
        self.LOS_initial = LOS
        self.LOS_final = LOS
        if remaining_stay:
            self.remaining_stay = remaining_stay
        else:
            self.remaining_stay = LOS
        self.incubation_period = 0
        self.max_incubation_period = random.randint(required_parameters['incubation_time_min'], required_parameters['incubation_time_max'])
        self.duration_infection = 0
        self.treatment_days = 1000
        self.max_treatment_days = 1000
        self.treatment = None
        self.antibiotic = None
        self.microorganism = None
        self.colonized = colonized
        self.non_susceptible = non_susceptible

    def del_bed(self):
        self.b = None

    def print_patient(self):
        print(self.to_string())

    def to_string(self):
        return 'patient:' + str(self.id) + ' ' + str(self.state) + '; age: ' + str(
            self.age) + '; sex: ' + self.sex + '; LOS: ' + str(self.LOS_initial) + '; Rem_time: ' + str(
            self.remaining_stay) + '; treatment_days: ' + str(self.treatment_days) + '; l: ' + self.l.to_string()


class TypeLocalization(Enum):
    ER = 1
    Surgery = 2
    Radiology = 3
    Room = 4
    ICU = 5
    Ward = 6
    Bed = 7


class Localization:
    def __init__(self, id, type_local, father_localization=None, infected=False, neighbor_localization=None,
                 children=None):
        self.id = id
        self.name = type_local
        self.infected = infected
        self.located_in = father_localization
        self.adjacent = neighbor_localization
        self.children = children

    def to_string(self):
        s = str(self.id) + ';' + str(self.infected)
        return s

    def print_localization(self):
        print(self.to_string())

    def print_to_csv(self):
        nombre = str(self.name.name) + "_" + str(self.id)
        # save the service name if it is a service
        if (self.name not in [TypeLocalization.Room, TypeLocalization.Bed]):
            service = self.name.name
        else:
            service = "NULL"

        # save the adjacents
        if (self.adjacent is not None):
            adjacents = "{" + str(self.adjacent[0].id) 
            for i in range(1, len(self.adjacent)):
                adjacents = adjacents + "," + str(self.adjacent[i].id)
            adjacents = adjacents + "}"
        else:
            adjacents = "NULL"

        # save the children
        if (self.children is not None):
            kids = "{" + str(self.children[0].id) 
            for i in range(1, len(self.children)):
                kids = kids + "," + str(self.children[i].id)
            kids = kids + "}"
        else:
            kids = "NULL"

        # save the parent
        if (self.located_in is not None):
            father = str(self.located_in.id)
        else:
            father = "NULL"


        s = str(self.id) + ';' + str(self.infected) + ';' + adjacents + ';' + kids + ';' + nombre + ';' + father + ';NULL;' + str(service)
        return s

    def is_temporal_localization(self):
        return self.name in [TypeLocalization.Surgery, TypeLocalization.Radiology] or self.located_in.name in [TypeLocalization.Surgery, TypeLocalization.Radiology]

    def is_room_bed(self):
        return self.name is TypeLocalization.Bed and self.located_in.name is TypeLocalization.Room

    def is_room_infected(self):
        return self.located_in.infected

    def is_located_in(self, location):
        return self.located_in is not None and self.located_in.name is location

    def is_localization(self, localization):
        if self.name is not None:
            if self.located_in is not None and localization == TypeLocalization.Radiology:
                return self.located_in.name is localization
            return self.name is localization
        return False

    def get_ward(self):
        if self.located_in and self.located_in.name is TypeLocalization.Room:
            return self.located_in.located_in
        else:
            return None


# Function that returns an age for an adult patient following a normal distribution
def normal_distr_age():
    ages = np.arange(18, 96, 1)
    mean = np.mean(ages) - 2.5
    sd = np.std(ages)
    tn = truncnorm((18 - mean) / sd, (95 - mean) / sd, loc=mean, scale=sd)
    return int(tn.rvs())


# Function that returns an initial LOS for a patient following a lognormal distribution
def lognormal_distr_los():
    mean = 4.254
    sd = 0.293
    normal_mean = np.log(mean) - sd ** 2 / 2
    dist = np.random.lognormal(normal_mean, sd, 1)
    return int(dist)


# Function that assigns a sex to a patient (M = Male, F = Female)
def set_sex():
    r = bernoulli.rvs(0.5)
    if r == 1:
        return 'M'
    else:
        return 'F'

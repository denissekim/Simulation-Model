""" HOSPITAL MODEL CREATION """
from scipy.stats import truncnorm
from classes import Localization, TypeLocalization
import config


# Function that initialize a new hospital
# er_nbeds = nº beds in the ER
# icu_nbeds = nº beds in the ICU
# nwards = nº of wards
# wards_nrooms = nº of rooms in each ward
# room_nbeds = nº of beds in each ward room
# sx_nrooms = nº of operation rooms
# rx_nrooms = nº of radiology rooms
# rx_nbeds = nº of beds in each radiology room
def initialize_hospital(er_nbeds, icu_nbeds, nwards, wards_nrooms, sx_nrooms, rx_nrooms, room_nbeds, rx_nbeds):
    # ER, ICU
    er = Localization(0, TypeLocalization.ER)
    icu = Localization(1, TypeLocalization.ICU)
    L = [er, icu]

    # WARDS
    index = 2
    wards = []
    for i in range(nwards):
        ward = Localization(index, TypeLocalization.Ward)
        wards.append(ward)
        L.append(ward)
        index = index + 1

    # OPERATING ROOMS
    for i in range(sx_nrooms):
        surgery = Localization(index, TypeLocalization.Surgery)
        L.append(surgery)
        index = index + 1

    # RADIOLOGY ROOMS
    rxs = []
    for i in range(rx_nrooms):
        rx = Localization(index, TypeLocalization.Radiology)
        rxs.append(rx)
        L.append(rx)
        index = index + 1

    # ER AND ICU BEDS
    beds_er, index = create_beds(er_nbeds, index, er, 0)
    L = L + beds_er
    beds_icu, index = create_beds(icu_nbeds, index, icu, 1)
    L = L + beds_icu

    # WARDS ROOMS AND BEDS
    for i in range(nwards):
        rooms, beds, index = create_rooms(wards_nrooms[i], room_nbeds, index, wards[i])
        L = L + rooms
        L = L + beds

    # RADIOLOGY ROOMS BEDS
    for i in range(rx_nrooms):
        beds, index = create_beds(rx_nbeds, index, rxs[i], 0)
        L = L + beds

    return L


# Function that returns new nbeds beds following the id index, located in located_in and in the floor floor
def create_beds(nbeds, index, located_in, floor):
    beds = []
    for i in range(nbeds):
        beds.append(Localization(index, TypeLocalization.Bed, located_in))
        index = index + 1
    for b in beds:
        b.adjacent = [bed for bed in beds if bed is not b]
    located_in.children = beds
    return beds, index


# Function that returns new nrooms rooms with room_nbeds beds following the id index, located in located_in
def create_rooms(nrooms, room_nbeds, index, located_in):
    rooms = []
    beds = []
    for i in range(nrooms):
        r = Localization(index, TypeLocalization.Room, located_in)
        rooms.append(r)
        index = index + 1
        bs, index = create_beds(room_nbeds, index, r, 1)
        beds.extend(bs)

    located_in.children = rooms
    return rooms, beds, index


# Function that returns the daily hospital occupancy rate following a normal distribution
def normal_distr_occupancy_rate():
    max_beds = round(config.hosp_occupancy_rate / config.required_parameters['n_patients'])
    mean = config.hosp_occupancy_rate
    sd = 30
    tn = truncnorm((0 - mean) / sd, (max_beds - mean) / sd, loc=mean, scale=sd)
    return int(tn.rvs())

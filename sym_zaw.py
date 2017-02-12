#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  sym_zaw.py

from read_file import Measurment
from matplotlib import pyplot as plt
from hydrostatic_drive import *
import json



# Utworzenie obiektu, który jest instancją klasy Measurment
# Plik z pomiarami zaworu proporcjonalnego
ZawProp = Measurment('85_2500.txt')

# Wczytanie parametrów zaworu z pliku json
filename = 'zawór_ciśnieniowy.json'

with open(filename) as f_obj:
	valve_parameters = json.load(f_obj)

print(valve_parameters)

# utworzenie obiektu zaworu ciśnieniowego
zawSym = PressureReliefValveComplex(**valve_parameters)

# utworzenie układu hydrostatycznego
naped_hydr = HydrostaticDrive(a=1.0, c=0.1, prv=zawSym)

# pętla symulacji
wyniki_sym = naped_hydr.simulation(0.09, 0., 4800, ZawProp.values[4])

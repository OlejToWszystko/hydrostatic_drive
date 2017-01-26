#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  symulacja_zaworu.py
#  
#  Copyright 2017 Maciej Olejnik <maciej@debian>

from read_file import Measurment
from matplotlib import pyplot as plt
from hydrostatic_drive import *


# Utworzenie obiektu, który jest instancją klasy Measurment
# Plik z pomiarami zaworu proporcjonalnego
ZawProp = Measurment('85_2500.txt')

# Utworzenie obiektów modeli symulacyjnych
ZawSym = PressureReliefValve(Tz=0.0085, hz=1.85, pz=85)
naped_hydr = HydrostaticDrive(0.8, 0.146, prv=ZawSym)
wyniki_sym = naped_hydr.simulation(0.09, 0., 4800, ZawProp.values[4])

# Utworzenie obiektów klasy Figure (1 obiekt) oraz Axes (dwa obiekty)
fig1, (sp1, sp2) = plt.subplots(1, 2)

# utworzenie serii danych p(t) na obiekcie sp1
sp1.plot(ZawProp.values[0], ZawProp.values[5])

# utworzenie nowego obiektu Axes, ze wspólną osią x z sp1
sp1_other = sp1.twinx()

# dodanie serii danych do obiektu sp1_other
sp1_other.plot(ZawProp.values[0], ZawProp.values[4], c='r')
sp1_other.plot(ZawProp.values[0], ZawProp.values[3], c='green')
sp1_other.plot(ZawProp.values[0], ZawProp.values[2], c='pink')
sp1_other.plot(wyniki_sym['t'], wyniki_sym['Q_act'], c='yellow')

# utworzenie serii danych p(t) na obiekcie sp2
sp2.plot(ZawProp.values[0], ZawProp.values[5])

# utworzenie nowego obiektu Axes, ze wspólną osią x z sp1 oraz
# utworzenie serii danych n(t) na tym obiekcie
sp2.twinx().plot(ZawProp.values[0], ZawProp.values[1], c='green')

# ustawienie układu na obiekcie fig1
fig1.tight_layout()

# Utworzenie obiektów klasy Figure (1 obiekt) oraz Axes (1 obiekt)
fig2, sp3 = plt.subplots()

# utworzenie serii danych p(t) na obiekcie sp3
sp3.plot(ZawProp.values[0], ZawProp.values[5])
sp3.plot(wyniki_sym['t'], wyniki_sym['p'], c='red')
sp3_other = sp3.twinx()
sp3_other.plot(wyniki_sym['t'], wyniki_sym['Q_act'], c='yellow')
sp3_other.plot(wyniki_sym['t'], wyniki_sym['Qz'], c='green')


# ustawienie układu na obiekcie fig2
fig2.tight_layout()

# wyświetlenie wykresów
plt.show()





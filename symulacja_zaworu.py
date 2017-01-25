#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  symulacja_zaworu.py
#  
#  Copyright 2017 Maciej Olejnik <maciej@debian>

from read_file import Measurment
from matplotlib import pyplot as plt
# from array import array

ZawProp = Measurment('85_2500.txt')
# ZawProp.print_channels()
# ZawProp.print_units()
'''
for value in ZawProp.values:
	print(value)
'''


# print(ZawProp.values[5][:100])
# print("Ilość elementów na liście : " + str(len(ZawProp.values)))
# ZawProp.measurment_to_dict()
# print(ZawProp.dictMeasurment)

settings_pressure = {
				'x' : ZawProp.values[0],
				'y' : ZawProp.values[5],
				's' : 1,
				'c' : 'red',
				'edgecolor' : 'None',
				}
				
#plt.scatter(**settings_pressure)
#plt.csd(array('f', ZawProp.values[0]), array('f', ZawProp.values[5]))
'''
pressure_chart = plt.subplot(121)
pressure_chart.plot(ZawProp.values[0], ZawProp.values[5])
others_chart = pressure_chart.twinx()
others_chart.plot(ZawProp.values[0], ZawProp.values[4], c='r')
others_chart.plot(ZawProp.values[0], ZawProp.values[3], c='green')
others_chart.plot(ZawProp.values[0], ZawProp.values[2], c='pink')
p_v = plt.subplot(122)
p_v.plot(ZawProp.values[0], ZawProp.values[5])
p_v.twinx().plot(ZawProp.values[0], ZawProp.values[1], c='green')
'''
fig1, (sp1, sp2) = plt.subplots(1, 2)
sp1.plot(ZawProp.values[0], ZawProp.values[5])
sp1_other = sp1.twinx()
sp1_other.plot(ZawProp.values[0], ZawProp.values[4], c='r')
sp1_other.plot(ZawProp.values[0], ZawProp.values[3], c='green')
sp1_other.plot(ZawProp.values[0], ZawProp.values[2], c='pink')
sp2.plot(ZawProp.values[0], ZawProp.values[5])
sp2.twinx().plot(ZawProp.values[0], ZawProp.values[1], c='green')
fig1.tight_layout()
fig2, sp3 = plt.subplots()
sp3.plot(ZawProp.values[0], ZawProp.values[5])
#plt.plot(ZawProp.values[0], ZawProp.values[5])
#plt.plot(ZawProp.values[0], ZawProp.values[4])
#plt.plot(ZawProp.values[0], ZawProp.values[3])
#plt.plot(ZawProp.values[0], ZawProp.values[2])
#plt.plot(ZawProp.values[0], ZawProp.values[1])
plt.show()





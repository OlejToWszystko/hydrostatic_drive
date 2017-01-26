#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  hydrostatic_drive.py

from math import pi

class PressureReliefValve():
	'''
	Obiekt reprezentujący zawór przelewowy.
	'''
	
	def __init__(self, Tz=0.0085, hz=1.85, pz=85):
		'''
		Inicjalizacja parametrów zaworu.
		
		Tz - stała czasowa zaworu [s]
		hz - współczynnik wzmocnienia zaworu
		pz - ciśnienie otwarcia zaworu [Pa]
		''' 
		self.Tz = Tz				
		self.hz = hz * 10**(-10)	
		self.pz = pz * 10**5		
		
		
	def diff_flow_rate(self, p=0, Qz=0):
		'''
		Metoda obliczająca pochodną po czasie natężenia przepływu przez 
		zawór ciśnieniowy.
		'''
		if p >= self.pz:
			# dQz/dt = (hz/Tz)*p - (1/Tz)*Qz - (hz/Tz)*pz
			dQz = (self.hz/self.Tz) * (p-self.pz) - (Qz/self.Tz)
		else:
			dQz = -(Qz/self.Tz)
			
		return dQz
		
		
class LinearReceiver():
	'''
	Klasa reprezentująca liniowy odbiornik, np. siłownik nurnikowy.
	'''
	
	def __init__(self, dt, f, mq):
		''' 
		Inicjalizacja parametrów odbiornika.
		
		dt - średnica tłoka [m]
		f - współczynnik oporów wiskotycznych [N*s/m]
		mq - masa ładunku [kg]
		'''
		self.dt = dt
		self.f = f
		self.mq = mq
		self.At = (pi*(self.dt)**2)/4
		
		
	def acceleration(self, p, v):
		'''
		Metoda obliczająca przyspieszenie odbiornika zgodnie z Drugą
		Zasadą Dynamiki Newton'a.
		'''
		dv = (self.At/self.mq)*p - (self.f/self.mq)*v - 9.81
		
		return dv
		
		
class HydrostaticDrive():
	'''
	Klasa reprezentująca układ hydrostatyczny.
	'''
	
	def __init__(
		self, a, c, if_prv=True, if_receiver=False, prv_kwargs=None,
		receiver_kwargs=None,
		):
		'''
		Inicjalizacja układu hydrostatycznego.
		
		a - współczynnik przecieków
		c - kapacytancja układau
		if_prv - czy obecny jest zawór ciśnieniowy
		if_receiver - czy obecny jest odbiornik
		'''
		
		self.a = a
		self.c = c
		self.if_prv = if_prv
		self.if_receiver = if_receiver
		self.prv_kwargs = prv_kwargs
		self.receiver_kwargs = receiver_kwargs
		self.components()
		
		
	def components(self):
		'''
		Tworzy opcjonalne obiekty układu hydrostatycznego.
		'''
		if self.if_prv:
			if self.prv_kwargs:
				prv_parameters_list = ['Tz', 'pz', 'hz']
				prv_parameters = {}
				for prv_parameter in prv_parameters_list:
					if prv_parameter in self.prv_kwargs.keys():
						prv_parameters[prv_parameter] = ( 
							self.prv_kwargs[prv_parameter]
							)
				self.prv = PressureReliefValve(**prv_parameters)
			else:
				self.prv = PressureReliefValve()
		
		if self.if_receiver:
			if self.receiver_kwargs:
				receiver_parameters_list = ['dt', 'f', 'mq']
				receiver_parameters = {}
				for receiver_parameter in receiver_parameters_list:
					if (receiver_parameter in 
						self.receiver_kwargs.keys()
						):
						receiver_parameters[receiver_parameter] = ( 
							self.receiver_kwargs[receiver_parameter]
							)
				self.receiver = LinearReceiver(**receiver_parameters)
			else:
				self.receiver = LinearReceiver()		 
			
		
		


def main(args):
    return 0

if __name__ == '__main__':
    import sys
    
    Zawor = PressureReliefValve(Tz=0.01, hz=2, pz=150)
    print('Tz = ' + str(Zawor.Tz))
    print('hz = ' + str(Zawor.hz))
    print('pz = ' + str(Zawor.pz))
    
    valveParam = {'pz' : 250, 'hz' : 4, 'Tz' : 0.07, 'Tz2' : 1}
    receiverParam = {'f' : 2500, 'dt' : 0.069, 'mq' : 900}
    hydr_sys = HydrostaticDrive(
		1, 1, prv_kwargs=valveParam, if_receiver=True, 
		receiver_kwargs=receiverParam,
		)
		
    print('Tz = ' + str(hydr_sys.prv.Tz))
    print('hz = ' + str(hydr_sys.prv.hz))
    print('pz = ' + str(hydr_sys.prv.pz))
    
    
    print('f = ' + str(hydr_sys.receiver.f))
    print('dt = ' + str(hydr_sys.receiver.dt))
    print('mq = ' + str(hydr_sys.receiver.mq))
    sys.exit(main(sys.argv))

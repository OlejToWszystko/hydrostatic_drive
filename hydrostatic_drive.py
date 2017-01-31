#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  hydrostatic_drive.py

from math import pi, sin, cos, radians


def euler(x, dx, dt):
	x_new = x + dx * dt
	return x_new


class PressureReliefValveComplex():
	'''
	Szczegółowy model zaworu przelewowego.
	'''
	
	def __init__(
		self, dg, lg, alfa, Cg, ksg, kmi, mg, msg, adj_force, dtl, Cst,
		xt0, kst, dd, Cd, kmi_t, mt, mst, ro,
		):
		'''Inicjalizacja parametrów zaworu.'''
		self.dg = dg * 0.001
		self.lg = lg * 0.001
		self.alfa = radians(alfa)
		self.tetha = self.alfa + 0.5*(pi/2 - self.alfa)
		self.Cg = Cg
		self.ksg = ksg
		self.kmi = kmi
		self.mg = mg
		self.msg = msg
		self.mgz = self.mg + 1/3 * self.msg
		self.poppetArea = self.area(self.dg)
		self.adj_force = adj_force
		self.fgz = self.equivalent_area()
		
		self.dtl = dtl * 0.001
		self.Cst = Cst
		self.xt0 = xt0
		self.kst = kst
		self.dd = dd * 0.001
		self.Cd = Cd
		self.fd = self.area(self.dd)
		self.kmi_t = kmi_t
		self.mt = mt
		self.mst = mst
		self.mtz = self.mt + 1/3 * self.mst 
		self.spoolArea = self.area(self.dtl)
		
		self.ro = ro
		
		
	def area(self, diameter):
		'''Metoda obliczająca pole koła'''
		A = pi * (diameter ** 2) / 4
		
		return A
		
	def equivalent_area(self):
		'''Metoda obliczająca powierzchnię zastępczą grzybka'''
		slotArea = (
			pi * ((self.dg / 2 + self.lg * sin(self.alfa)) ** 2
			- (self.dg / 2) ** 2)
			)
		equivalentArea = self.poppetArea + 0.5 * slotArea
		
		return equivalentArea
		
		
	def hydrostatic_force(self, area, pressure_drop):
		'''Metoda obliczająca siłę hydrostatyczną'''
		hydrostaticForce = pressure_drop * area
		
		return hydrostaticForce
		
		
	def spring_force(self, x0, x, k):
		'''Metoda obliczająca siłę sprężyny'''
		Ps = k * (x0 + x)
		
		return Ps
		
		
	def poppet_hydrodynamic_force(self, pressure_drop, 
		poppet_displacement):
		'''
		Metoda obliczająca siłę hydrodynamiczną działającą na 
		grzybek zaworu.
		'''
		l2z = (
			pi * (sin(self.alfa) / cos(self.tetha - self.alfa)) 
			* (self.dg + self.lg * sin(self.alfa))
			)
		xg = poppet_displacement
		pt = pressure_drop
		
		P_hdg = (
			2 * (self.Cg ** 2) * (l2z ** 2) * 
			((xg * cos(self.tetha)) / l2z 
			- (xg ** 2) / self.poppetArea) * pt
			)
			
		return P_hdg 
		
		
	def spool_hydrodynamic_force(self, p, xt):
		''' 
		Metoda obliczająca siłę hydrodynamiczną działającą na tłoczek
		zaworu ciśnieniowego
		'''
		
		P_hdt = 0.72 * self.Cst * pi * self.dtl * xt * p
		
		return P_hdt 
		
		
	def spool_damping_force(self, vt):
		'''
		Metoda obliczająca siłę tłumiącą, wynikającą z przepływu cieczy
		przez dławik, działającą na tłoczek głównego stopnia zaworu.
		'''
		P_damp = (
			(self.ro * self.spoolArea ** 3) / (2 * (self.Cd ** 2) * 
			(self.fd ** 2)) * vt ** 2
			)
			
		return P_damp
		
		
	def friction_force(self, coefficient, velocity):
		'''
		Metoda obliczająca siłę tarcia lepkiego.
		'''
		return (coefficient * velocity)
		
		
	def diff_equations(self, p, pt, xg, vg, xt, vt):
		'''
		Równania różniczkowe zaworu przelewowego
		'''
		dvg = (
			self.hydrostatic_force(self.fgz, pt) - 
			self.poppet_hydrodynamic_force(pt, xg) - self.adj_force - 
			self.friction_force(self.kmi, vg)
			) * 1 / self.mgz
		dvt = (
			self.hydrostatic_force(self.spoolArea, (p-pt)) - 
			self.spool_hydrodynamic_force(p, xt) - 
			self.spring_force(self.xt0, xt, self.kst) -
			self.spool_damping_force(vt) - 
			self.friction_force(self.kmi_t, vt)
			) * 1 / self.mtz 
		 	
			

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
		
		
	def flow_rate(self, p, Qz, dt):
		'''
		Metoda obliczająca natężenie przepływu przez zawór.
		'''
		dQz = self.diff_flow_rate(p, Qz)
		flowRate = euler(Qz, dQz, dt)
		
		
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
	
	def __init__(self, a, c, prv=PressureReliefValve(), receiver=None):
		'''
		Inicjalizacja układu hydrostatycznego.
		
		a - współczynnik przecieków
		c - kapacytancja układau
		'''
		
		self.a = a * 10 ** (-12)
		self.c = c * 10 ** (-12)
		self.prv = prv
		self.receiver = receiver
	
		
	def diff_pressure(self, Qp, p, vt=0, Qz=0):
		''' 
		Metoda obliczająca pochodną ciśnienia po czasie w układzie, na 
		podstawie równania bilansu przepływu w układzie
		'''
		if self.receiver:
			dp = (
					-(self.a/self.c)*p - (self.receiver.At/self.c)*vt
					- (Qz/self.c) + (Qp/self.c)
					)
		else:
			dp = -(self.a/self.c)*p - (Qz/self.c) + (Qp/self.c)
			
		return dp
	
	
	def actuation(self, t, t_start, tr, Qp_meas):
		'''
		Funkcja generująca sygnał wymuszający w postaci natężenia
		przepływu.
		'''
		if t < t_start:
			Q_act = 0
		elif t < (t_start + tr):
			Q_act = (t - t_start)/(tr) * Qp_meas/60000
		else:
			Q_act = Qp_meas/60000
			
		return Q_act
		
		
	def simulation(self, t_start, tr, freq, Qp_meas):
		'''
		Główna pętla symulacji
		'''
		
		i_max = len(Qp_meas)
		t = 0
		dt = 1.0/freq
		if self.receiver and self.prv:
			p, Qz, vt = 0, 0, 0
			# utworzenie pustego słownika z wynikami symulacji
			sim_res = {
				't' : [],
				'Q_act' : [],
				'p' : [],
				'Qz' : [],
				'vt' : [],
				}
			for i in range(i_max):
				# wymuszenie
				Q_act = self.actuation(t, t_start, tr, Qp_meas[i])
				# pochodna ciśnienia
				dp = self.diff_pressure(Q_act, p, vt, Qz)
				# pochodna prędkości
				dvt = self.receiver.acceleration(p, vt)
				# pochodna przepływu przez zawór
				dQz = self.prv.diff_flow_rate(p, Qz)
				# zapis wartości do słownika
				sim_res['t'].append(t)
				sim_res['Q_act'].append(Q_act * 60000)
				sim_res['p'].append(p / 10**5)
				sim_res['Qz'].append(Qz * 60000)
				sim_res['vt'].append(vt)
				# metoda Eulera
				p = euler(p, dp, dt)
				vt = euler(vt, dvt, dt)
				Qz = euler(Qz, dQz, dt)
				if Qz >= (Qp_meas[i]/60000):
					Qz = Qp_meas[i]/60000
				t += dt
				
		elif self.receiver:
			p, vt = 0, 0
			# utworzenie pustego słownika z wynikami symulacji
			sim_res = {
				't' : [],
				'Q_act' : [],
				'p' : [],
				'vt' : [],
				}
			for i in range(i_max):
				# wymuszenie
				Q_act = self.actuation(t, t_start, tr, Qp_meas[i])
				# pochodna ciśnienia
				dp = self.diff_pressure(Q_act, p, vt)
				# pochodna prędkości
				dvt = self.receiver.acceleration(p, vt)
				# zapis wartości do słownika
				sim_res['t'].append(t)
				sim_res['Q_act'].append(Q_act * 60000)
				sim_res['p'].append(p / 10**5)
				sim_res['vt'].append(vt)
				# metoda Eulera
				p = euler(p, dp, dt)
				vt = euler(vt, dvt, dt)
				t += dt
				
		elif self.prv:
			p, Qz = 0, 0
			# utworzenie pustego słownika z wynikami symulacji
			sim_res = {
				't' : [],
				'Q_act' : [],
				'p' : [],
				'Qz' : [],
				}
			for i in range(i_max):
				# wymuszenie
				Q_act = self.actuation(t, t_start, tr, Qp_meas[i])
				# pochodna ciśnienia
				dp = self.diff_pressure(Q_act, p, 0, Qz)
				# pochodna przepływu przez zawór
				dQz = self.prv.diff_flow_rate(p, Qz)
				# zapis wartości do słownika
				sim_res['t'].append(t)
				sim_res['Q_act'].append(Q_act * 60000)
				sim_res['p'].append(p / 10**5)
				sim_res['Qz'].append(Qz * 60000)
				# metoda Eulera
				p = euler(p, dp, dt)
				Qz = euler(Qz, dQz, dt)
				if Qz >= (Qp_meas[i]/60000):
					Qz = Qp_meas[i]/60000
				t += dt
		
		return sim_res
		
			
			
			


def main(args):
    return 0

if __name__ == '__main__':
    import sys
    
    Zawor = PressureReliefValve(Tz=0.01, hz=2, pz=150)
    print('Tz = ' + str(Zawor.Tz))
    print('hz = ' + str(Zawor.hz))
    print('pz = ' + str(Zawor.pz))
    
    valveParam = {'pz' : 250, 'hz' : 4, 'Tz' : 0.07}
    receiverParam = {'f' : 2500, 'dt' : 0.069, 'mq' : 900}
    Zawor2 = PressureReliefValve(**valveParam)
    silownik = LinearReceiver(**receiverParam)
    hydr_sys = HydrostaticDrive(1, 1, prv=Zawor2, receiver=silownik)
		
    print('Tz = ' + str(hydr_sys.prv.Tz))
    print('hz = ' + str(hydr_sys.prv.hz))
    print('pz = ' + str(hydr_sys.prv.pz))
    
    
    print('f = ' + str(hydr_sys.receiver.f))
    print('dt = ' + str(hydr_sys.receiver.dt))
    print('mq = ' + str(hydr_sys.receiver.mq))
    sys.exit(main(sys.argv))

#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  hydrostatic_drive.py

from math import pi, sin, cos, radians, sqrt


def euler(x, dx, dt):
	x_new = x + dx * dt
	return x_new


class PressureReliefValveComplex():
	'''
	Szczegółowy model zaworu przelewowego.
	'''
	
	def __init__(
		self, dg, lg, alfa, Cg, ksg, kmi, mg, msg, adj_force, dtl, Cst,
		xt0, kst, dd, Cd, kmi_t, mt, mst, ckt, ro,
		):
		'''Inicjalizacja parametrów zaworu.'''
		self.dg = dg * 0.001
		self.lg = lg * 0.001
		self.alfa = radians(alfa)
		self.tetha = self.alfa + 0.5*(pi/2 - self.alfa)
		self.Cg = Cg
		self.ksg = ksg
		self.kmi = kmi
		self.mg = mg * 0.001
		self.msg = msg * 0.001
		self.mgz = self.mg + 1/3 * self.msg
		self.poppetArea = self.area(self.dg)
		self.adj_force = adj_force
		self.fgz = self.equivalent_area()
		self.l2z = self.equivalent_lenght_l2z()
		
		self.dtl = dtl * 0.001
		self.Cst = Cst
		self.xt0 = xt0 * 0.001
		self.kst = kst
		self.dd = dd * 0.001
		self.Cd = Cd
		self.fd = self.area(self.dd)
		self.kmi_t = kmi_t
		self.mt = mt * 0.001
		self.mst = mst * 0.001
		self.mtz = self.mt + 1/3 * self.mst 
		self.ckt = ckt * 10 ** (-12)
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
		
		
	def equivalent_lenght_l2z(self):
		''' Metoda obliczająca długość zastępczą l2z. '''
		l2z = (
			pi * (sin(self.alfa) / cos(self.tetha - self.alfa)) 
			* (self.dg + self.lg * sin(self.alfa))
			)
			
		return l2z
		
		
	def hydrostatic_force(self, area, pressure_drop):
		'''Metoda obliczająca siłę hydrostatyczną'''
		hydrostaticForce = pressure_drop * area
		
		return hydrostaticForce
		
		
	def spring_force(self, x0, x, k):
		'''Metoda obliczająca siłę sprężyny'''
		Ps = k * (x0 + x)
		
		return Ps
		
		
	def poppet_hydrodynamic_force(self, pt, xg):
		'''
		Metoda obliczająca siłę hydrodynamiczną działającą na 
		grzybek zaworu.
		'''
		P_hdg = (
			2 * (self.Cg ** 2) * (self.l2z ** 2) * 
			((xg * cos(self.tetha)) / self.l2z 
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
			(self.ro * self.spoolArea ** 3) / (2 * (self.Cd ** 2) * (self.fd ** 2)) * vt ** 2
			)
			
		return P_damp
		
		
	def friction_force(self, coefficient, velocity):
		'''
		Metoda obliczająca siłę tarcia lepkiego.
		'''
		return (coefficient * velocity)
	
	
	def throttle_flow_rate(self, p, pt):
		'''
		Metoda obliczająca natężenie przepływu w dławiku tłoczka.
		'''
		if pt > p:
			pt = p
		
		Qd = self.Cd * self.fd * sqrt(2 * (p - pt) / self.ro)
		
		return Qd
		
		
	def spool_flow_rate(self, xt, p):
		'''
		Metoda obliczająca natężenie przepływu przez szczelinę 
		utworzoną przez tłoczek.
		'''
		if p < 0:
			p = 0
			
		Qt = self.Cst * pi * self.dtl * xt * sqrt(2 * p / self.ro)
		
		return Qt
		
		
	def poppet_flow_rate(self, xg, pt):
		'''
		Metoda obliczająca natężenie przepływu przez grzybek zaworu.
		'''
		if pt < 0:
			pt = 0
			
		Qg = self.Cg * xg * self.l2z * sqrt(2 * pt / self.ro)
		
		return Qg
		
		
	def valve_flow_rate(self, p, pt, xt, xg):
		'''
		Metoda obliczająca natężenie przepływu przez zawór.
		'''
		Qz = self.poppet_flow_rate(xg, pt) + self.spool_flow_rate(xt, p)
		
		return Qz
		
		
	def diff_equations(self, p, pt, xg, vg, xt, vt):
		'''
		Równania różniczkowe zaworu przelewowego
		'''
		dvg = (
			self.hydrostatic_force(self.fgz, pt) - 
			self.poppet_hydrodynamic_force(pt, xg) - self.adj_force - 
			self.friction_force(self.kmi, vg)
			) / self.mgz
			
		dvt = (
			self.hydrostatic_force(self.spoolArea, (p-pt)) - 
			self.spool_hydrodynamic_force(p, xt) - 
			self.spring_force(self.xt0, xt, self.kst) -
			self.spool_damping_force(vt) - 
			self.friction_force(self.kmi_t, vt)
			) / self.mtz
		 
		dxg = vg
		 
		dxt = vt
		
		dpt = (
			self.throttle_flow_rate(p, pt) + self.spoolArea * vt -
			self.poppet_flow_rate(xg, pt) - self.poppetArea * vg
			) / self.ckt
			
		return dpt, dxg, dvg, dxt, dvt
				

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
		#flowRate = euler(Qz, dQz, dt)
		
		
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
			if type(self.prv) == PressureReliefValve:
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
					#if Qz >= (Qp_meas[i]/60000):
					#	Qz = Qp_meas[i]/60000
					t += dt
			
			elif type(self.prv) == PressureReliefValveComplex:
				p, pt, xg, vg, xt, vt, Qz = 0, 0, 0, 0, 0, 0, 0
				# utworzenie pustego słownika z wynikami symulacji
				sim_res = {
					't' : [],
					'Q_act' : [],
					'p' : [],
					'Qz' : [],
					'xg' : [],
					'vg' : [],
					'xt' : [],
					'vt' : [],
					}
				for i in range(i_max):
					# wymuszenie
					Q_act = self.actuation(t, t_start, tr, Qp_meas[i])
					# pochodna ciśnienia
					dp = self.diff_pressure(Q_act, p, 0, Qz)
					# pochodne od zaworu
					dpt, dxg, dvg, dxt, dvt = ( 
						self.prv.diff_equations(p, pt, xg, vg, xt, vt)
						)
					# zapis wartości do słownika
					sim_res['t'].append(t)
					sim_res['Q_act'].append(Q_act * 60000)
					sim_res['p'].append(p / 10**5)
					sim_res['Qz'].append(Qz * 60000)
					sim_res['xg'].append(xg * 1000)
					sim_res['vg'].append(vg)
					sim_res['xt'].append(xt)
					sim_res['vt'].append(vt)
					# metoda Eulera
					p = euler(p, dp, dt)
					pt = euler(pt, dpt, dt)
					xg = euler(xg, vg, dt)
					vg = euler(vg, dvg, dt)
					xt = euler(xt, vt, dt)
					vt = euler(vt, dvt, dt)
					Qz = self.prv.valve_flow_rate(p, pt, xt, xg)
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

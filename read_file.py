#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  read_file.py
#  
#  Copyright 2017 Maciej Olejnik <maciej@debian>

from collections import OrderedDict

class Measurment():
	'''Tworzy obiekt pomiaru z pliku tekstowego'''
	
	def __init__(self, file_name):
		'''Inicjalizacja atrybutów obiektu'''
		self.file_name = file_name
		self.units = []
		self.channels = []
		self.values = []
		self.dictMeasurment = OrderedDict()
		
		# Wczytanie pliku i utworzenie odpowiednich list
		self.read_measurment()
		# usunięcie białych znaków po lewej i prawej stronie w values
		self.clear_lines(self.values)
		# zamiana ',' na '.'
		self.replace_char(self.values, ',', '.')
		# zapakowanie każdej serii pomiarowej do osobnej listy i 
		# umieszczenie każdej tej listy do listy self.values
		self.series_to_lists()
		

	def read_measurment(self):
		'''
		Wczytuje plik z pomiarami, gdzie w pierwszej linijce są opisane 
		jednostki, a w drugiej linii opisane są kanały 
		'''
		with open(self.file_name) as f:
			self.units = f.readline().strip().split('\t')
			self.channels = f.readline().strip().split('\t')
			self.values = f.readlines()
		
			
	def clear_lines(self, list_of_lines):		
		'''
		Usuwa białe znaki po lewej i prawej stronie każdej linii 
		tekstu, jaka znajduje się na liście
		'''
		for i in range(len(list_of_lines)):
			list_of_lines[i] = list_of_lines[i].strip()	
			
	
	def replace_char(self, list_of_lines, char_1, char_2):
		''' 
		Zamienia znak char_1 na char_2 w każdym wierszu listy 
		list_of_lines
		'''
		for i in range(len(list_of_lines)):
			list_of_lines[i] = list_of_lines[i].replace(char_1, char_2)
		
	
	def list_list(self, list_of_values):
		'''
		Rozdziela, przechowywany w elemencie listy, wiersz tekstu, 
		tworząc z niego listę. Znakiem podziału jest biały znak '\t'.
		'''
		for i in range(len(list_of_values)):
			list_of_values[i] = list_of_values[i].split('\t')
	
	
	def series_to_lists(self):
		'''
		Zmierzone wartości z każdego kanału zapisuje w osobnej liście
		'''	
		# utworzenie pustej, pomocniczej listy	
		aux_list = []
		# dodaje tyle pustych list do listy pomocniczej, ile jest 
		# kanałów pomiarowych
		for i in range(len(self.channels)):
			aux_list.append([])
		# utworzenie listy list z values
		self.list_list(self.values)
		
		# dla każdej listy, odpowiadającej chwili czasowej, z listy 
		# values ...
		for channels_values in self.values:
			i = 0 # iterator
			# dla każdej wartości z listy, odpowiadającej chwili 
			# czasowej
			for channel_value in channels_values:
				# dołącz do i-tej listy w aux_list
				aux_list[i].append(float(channel_value))
				i += 1
				
		# po wypełnieniu aux_list, nadpisz listę values
		# najpierw wyczyszczona zostaje lista self.values
		del self.values[:]
		# skopiowanie aux_list do pustej listy self.values
		self.values[:] = aux_list[:] 
					
			
	def print_units(self):
		'''Drukuje listę jednostek'''
		units_string = ''
		for unit in self.units:
			units_string += unit + '\t'
		# usunięcie znaku '\t' z ostatniej kolumny
		units_string = units_string.rstrip()
		print(units_string)
		
		
	def print_channels(self):
		'''drukuje listę kanałów'''
		channels_string = ''
		for channel in self.channels:
			channels_string += channel + '\t'
		# usunięcie znaku '\t' z ostatniej kolumny
		channels_string = channels_string.rstrip()
		print(channels_string)
		
		

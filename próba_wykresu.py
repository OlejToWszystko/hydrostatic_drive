#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  próba_wykresu.py
#  
#  Copyright 2017 Maciej Olejnik <maciej@debian>

from math import cos
import pygal
xy_chart = pygal.XY()
xy_chart.title = 'XY Cosinus'
xy_chart.add('x = cos(y)', [(cos(x / 10.), x / 10.) for x in range(-50, 50, 5)])
xy_chart.add('y = cos(x)', [(x / 10., cos(x / 10.)) for x in range(-50, 50, 5)])
xy_chart.add('x = 1',  [(1, -5), (1, 5)])
xy_chart.add('x = -1', [(-1, -5), (-1, 6)])
xy_chart.add('y = 1',  [(-5, 1), (5, 1)])
xy_chart.add('y = -1', [(-5, -1), (5, -1)])
lista_x = [x for x in range(1, 11)]
lista_y = [2*x for x in lista_x]
lista_xy = []
for i in range(len(lista_x)):
	lista_xy.append((lista_x[i], lista_y[i]))
xy_chart.add('y = 2*x', lista_xy, secondary=True)

#plik = xy_chart.render()
#print(plik)
xy_chart.render_to_file('próba.svg')
#xy_chart.render_to_png('próba.png')

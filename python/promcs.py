#!/usr/bin/env python3                                                                                                                           
# -*- coding: utf-8 -*- 
import numpy as np
import matplotlib.pyplot as plt

# Leer el archivo
ruta = '/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos/datos/promcs.dat'
data = np.loadtxt(ruta)

# Separar columnas
T_grid = data[:,0]
avg_xs = data[:,1]

# Graficar
plt.figure(figsize=(8,5))
plt.plot(T_grid, avg_xs, marker='o', linestyle='-', color='blue', label='avg_xs vs T_grid')
plt.yscale('log')  # opcional si los valores son muy pequeños
plt.xlabel('T_grid')
plt.ylabel('avg_xs')
plt.title('Gráfica de avg_xs vs T_grid')
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.legend()
plt.tight_layout()
plt.show()

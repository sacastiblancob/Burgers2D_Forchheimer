""""
Programa para graficar la malla que se utilizara para el dominio 2D de la ecua-
cion de Burgers. 
Los parametros se toman del archivo geometry.f95 que genera un archivo .dat
El archivo se lee y este programa grafica el dominio con puntos en un plano xy

Antonio Preziosi Ribero
December 2016
El archivo se incluye en la carpeta principal del programa 
"""

# Obteniendo las librerias necesarias
import numpy as np
import os
import matplotlib.pyplot as plt
import sys

# Imprimir parametros para ver donde esta trabajando cambiar con cada carpeta
# os.chdir('/home/toni/Documents/IAHR_2017/170118_Full2D')

# Definir salto o parte del archivo Geometry que no se va a leer
salto = 302
salto += 2 # Por las entradas de BC que tiene el archivo 

# Leyendo el archivo .dat de geometria (solo espacio con dos columnas)
data = np.genfromtxt('PruebaGeom.dat', skip_footer=salto)

# Haciendo el scatterplot del dominio computacional que se va a modelar
f = plt.figure()
plt.scatter(data[:,0], data[:,1], s = 1.5, color='blue')
plt.xlim(0,60)
plt.ylim(-28, 0)
plt.xlabel(r'Streamwise direction $(cm)$')
plt.ylabel(r'Vertical direction $(cm)$')
plt.show()
f.savefig("grid2D.eps")

print(sys.maxsize)

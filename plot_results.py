import matplotlib.pyplot as plt 
import sys

'''
Mejia Juarez Nancy Arlette
Tarea 03 Larga


Script que permite graficar los resultados obtenidos con una ejecucion
para resolver las funciones WFG2, WFG3 y WFG8.
'''


# Se grafican los datos del optimo de pareto
function = sys.argv[1]
scale_factor = [0.1, 0.005, 0.05, 0.01]

plt.figure()
file = open("results/" + function.lower() + ".dat")
x_vals = []
y_vals = []
for line in file:
	line = [float(x) for x in line.split()]
	x_vals.append(line[0])
	y_vals.append(line[1]) 
file.close()
plt.plot(x_vals, y_vals, "*", label="Optimo de pareto")

# Se grafican los resultados obtenidos con los diferentes factores de escala
for scale in scale_factor:
	file = None
	if scale == 0.1:
		file = open("results/" + function + "_solution0_scale_" + str(scale) + "00000.txt")
	elif scale == 0.005:
		file = open("results/" + function + "_solution0_scale_" + str(scale) + "000.txt")

	elif scale == 0.05:
		file = open("results/" + function + "_solution0_scale_" + str(scale) + "0000.txt")

	elif scale == 0.01:
		file = open("results/" + function + "_solution0_scale_" + str(scale) + "0000.txt")
	
	x_vals = []
	y_vals = []
	for line in file:
		line = [float(x) for x in line.split()]
		x_vals.append(line[0])
		y_vals.append(line[1]) 

	plt.plot( x_vals, y_vals, "*" , label="factor: " + str(scale))

	file.close()
plt.legend()
plt.savefig("results/frentes_"  + function + ".png")
plt.show()
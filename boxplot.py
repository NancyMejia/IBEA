import matplotlib.pyplot as plt 
import sys

'''
Mejia Juarez Nancy Arlette
Tarea 03 Larga

Script para crear boxplot de los datos, los datos
se tomaron como los datos de la ejecuciones realizadas
Si se realizan nuevas ejecuciones habria que modificar los datos
del calculo del hipervolumen
'''

WFG2_scale1 = [11.0870, 10.3372, 10.5979, 10.5518, 10.6018  ]
WFG2_scale2 = [10.8286, 11.3852, 10.5478, 10.6013, 10.7324 ]
WFG2_scale3 = [10.6179, 10.2895, 10.6034, 10.6190,  10.3326] 
WFG2_scale4 = [10.5588, 10.5860, 10.5938, 10.2099, 10.5797] 


WFG3_scale1 = [10.0301 ,9.8968 ,9.7992 ,9.9373 ,9.9846]
WFG3_scale2 = [9.7964 ,9.8439 ,9.8837 ,9.9025 ,9.8956]
WFG3_scale3 = [9.8766 ,9.8869 ,9.8454 ,9.8986 ,9.8681]
WFG3_scale4 = [9.8950 ,9.9105 ,9.9909 ,9.7581 ,9.9891]

WFG8_scale1 = [6.7207 , 6.6988 , 6.8138 , 7.2172 , 6.7042]
WFG8_scale2 = [6.7050 , 6.8070 , 7.0193 , 6.7073 , 6.7249] 
WFG8_scale3 = [6.5987 , 6.6416 , 6.6572 , 6.6151 , 6.5335]
WFG8_scale4 = [6.4715 , 6.5144 , 6.5510 , 6.4745 , 6.5504]          

plt.figure()
plt.boxplot([WFG2_scale1,WFG2_scale2, WFG2_scale3, WFG2_scale4])
plt.xticks([1, 2, 3, 4], ['0.1', '0.05', '0.01', '0.005'])
plt.title("Boxplot ejecuciones WFG2")
plt.xlabel("Factor de escala")
plt.ylabel("Hipervolumen")
plt.savefig("results/boxplot_wfg2.png")
plt.show()

plt.figure()
plt.boxplot([WFG3_scale1,WFG3_scale2, WFG3_scale3, WFG3_scale4])
plt.xticks([1, 2, 3, 4], ['0.1', '0.05', '0.01', '0.005'])
plt.title("Boxplot ejecuciones WFG3")
plt.xlabel("Factor de escala")
plt.ylabel("Hipervolumen")
plt.savefig("results/boxplot_wfg3.png")
plt.show()

plt.figure()
plt.boxplot([WFG8_scale1,WFG8_scale2, WFG8_scale3, WFG8_scale4])
plt.xticks([1, 2, 3, 4], ['0.1', '0.05', '0.01', '0.005'])
plt.title("Boxplot ejecuciones WFG8")
plt.xlabel("Factor de escala")
plt.ylabel("Hipervolumen")
plt.savefig("results/boxplot_wfg8.png")
plt.show()

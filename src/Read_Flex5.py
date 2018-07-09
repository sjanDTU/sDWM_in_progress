
import csv, re
import math as m
import numpy as np
import matplotlib.pyplot as plt

def read_Pitch_Flex5():
    fichier = 'Flex5_data/Wind_minPitch_pow.txt'
    f = open(fichier, "rb")
    fichier_csv = csv.reader(f, delimiter=" ")
    tab = list(fichier_csv)
    f.close()

    #print tab

    DATA=[]

    for i in range(len(tab)):
        L = []
        for c in tab[i]:
            x = re.findall(r"[-+]?\d*\.\d+|\d+", c)
            # print x
            if x != []:
                for j in range(len(x)):
                    x[j] = float(x[j])
                #L.append(float(x))
                L = L + x

        DATA = DATA + [L]
    #print DATA

    WS = []
    PITCH = []
    POWER = []

    for i in range(len(DATA)):
        WS.append(DATA[i][0])                # m/s
        PITCH.append(DATA[i][1])             # degree
        POWER.append(DATA[i][2])             # kW



    """
    plt.figure(1)
    plt.subplot(121), plt.plot(WS,PITCH,label='Pitch(U)'), plt.legend()
    plt.subplot(122), plt.plot(WS,POWER,label='Power(U)'), plt.legend()
    plt.show()
    #"""
    return [WS+[20.],PITCH + [PITCH[-1]], POWER + [POWER[-1]]]
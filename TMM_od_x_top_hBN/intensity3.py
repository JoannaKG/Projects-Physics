"""intensity2 zwraca R w funkcji wavelength uwzgledniajac zaleznosc od x w warstwie aktywnej
"""
from __future__ import division, print_function, absolute_import

from numpy import arange, array, pi, linspace, inf, zeros, exp

import numpy as np
from scipy.integrate import quad

import matplotlib.pyplot as plt

from scipy.interpolate import interp1d
import csv
degree = pi/180

def import_n(layer_csv):
    wavelength = []
    n = []
    k = []
    
    with open(layer_csv) as op:
        reader = csv.reader(op, delimiter=";")
        wavelength_list = list(zip(*reader))[0]
    for line in wavelength_list:
        wavelength.append(float(line)*1000)
        
    with open(layer_csv) as op:
        reader = csv.reader(op, delimiter=";")
        n_list = list(zip(*reader))[1]
    for line in n_list:
        n.append(float(line))
        
    with open(layer_csv) as op:
        reader = csv.reader(op, delimiter=";")
        k_list = list(zip(*reader))[2]
    for line in k_list:
        k.append(float(line)*1j)
    #print(wavelength, n, k)
        
    wavelength_array = array(wavelength)
    n_array = array(n)
    k_array = array(k)
        
    layer_n_list = []
    for i in range (0, len(wavelength)-1):
        layer_n_list.append([wavelength_array[i],n_array[i]-k_array[i]])

    layer_n_data = array(layer_n_list)
    layer_n_fn = interp1d(layer_n_data[:,0], layer_n_data[:,1], kind='linear')
    return layer_n_fn


def intensity(TMD, d_SiO2=300, d_hBN=20, lam_0=514, lambda_min = 360, lambda_max = 1000, number_of_layers = 1, th_0=0, dtype=complex): 
       #funkcja liczaca kontrast w funkcji dlugosci fali dla zadanej grubosci SiO2
    if TMD == 'MoS2':
        TMD_n_fn = lambda wavelenth : 4.024-1j*1.032
        TMD_thickness = number_of_layers*0.65
    if TMD == 'MoS2-ML':
        TMD_n_fn = import_n('MoS2-ML.csv')
        TMD_thickness = number_of_layers*0.65
    if TMD == 'MoSe2':
        TMD_n_fn = import_n('MoSe2.csv')
        TMD_thickness = number_of_layers*0.7
    if TMD == 'MoSe2-ML':
        TMD_n_fn = import_n('MoSe2-ML.csv')
        TMD_thickness = number_of_layers*0.7
        
    if TMD == 'WS2':
        TMD_n_fn = import_n('WS2.csv')
        TMD_thickness = number_of_layers*0.8
    if TMD == 'WS2-ML':
        TMD_n_fn = import_n('WS2-ML.csv')
        TMD_thickness = number_of_layers*0.8
    if TMD == 'WSe2':
        TMD_n_fn = import_n('WSe2.csv')
        TMD_thickness = number_of_layers*0.7
    if TMD == 'WSe2-ML':
        TMD_n_fn = import_n('WSe2-ML.csv')
        TMD_thickness = number_of_layers*0.7
        
    def beta_j(d,layer,lam, dtype=complex):
        
        if layer == 'air':
            n_layer = lambda wavelength : 1
        if layer == 'MoS2':
            n_layer = lambda wavelenth : 4.024-1j*1.032
        if layer == 'hBN':
            n_layer = lambda wavelength : 2.2
        if layer == 'SiO2':
            n_layer = import_n('SiO2.csv')
        if layer == 'Si':
            n_layer = lambda wavelenth : 4.149+1j*0.0426
        
        beta_j=2*pi*n_layer(lam)*d/lam
       # print("n1=" +str(n_layer(lam)))
        return beta_j
        
    def beta_x(x,lam, dtype=complex):
        
        n_layer = lambda wavelenth : 4.024-1j*1.032
        
        beta_x=2*pi*n_layer(lam)*x/lam
        return beta_x
        
    def r(layer_1, layer_2, lam, dtype=complex):
        if layer_1 == 'air':
            n_layer_1 = lambda wavelength : 1
        if layer_1 == 'MoS2':
            n_layer_1 = lambda wavelenth : 4.024-1j*1.032
        if layer_1 == 'hBN':
            n_layer_1 = lambda wavelength : 2.2
        if layer_1 == 'SiO2':
            n_layer_1 = import_n('SiO2.csv')
        if layer_1 == 'Si':
            n_layer_1 = lambda wavelenth : 4.149+1j*0.0426
            
        if layer_2 == 'air':
            n_layer_2 = lambda wavelength : 1
        if layer_2 == 'MoS2':
            n_layer_2 = lambda wavelenth : 4.024-1j*1.032
        if layer_2 == 'hBN':
            n_layer_2 = lambda wavelength : 2.2
        if layer_2 == 'SiO2':
            n_layer_2 = import_n('SiO2.csv')
        if layer_2 == 'Si':
            n_layer_2 = lambda wavelenth : 4.149+1j*0.0426
            
            
        r=(n_layer_1(lam)-n_layer_2(lam))/(n_layer_1(lam)+n_layer_2(lam))
        return r
    
    def t(layer_1, layer_2, lam, dtype=complex):
        if layer_1 == 'air':
            n_layer_1 = lambda wavelength : 1
        if layer_1 == 'MoS2':
            n_layer_1 = lambda wavelenth : 4.024-1j*1.032
        if layer_1 == 'hBN':
            n_layer_1 = lambda wavelength : 2.2
        if layer_1 == 'SiO2':
            n_layer_1 = import_n('SiO2.csv')
        if layer_1 == 'Si':
            n_layer_1 = lambda wavelenth : 4.149+1j*0.0426
            
        if layer_2 == 'air':
            n_layer_2 = lambda wavelength : 1
        if layer_2 == 'MoS2':
            n_layer_2 = lambda wavelenth : 4.024-1j*1.032
        if layer_2 == 'hBN':
            n_layer_2 = lambda wavelength : 2.2
        if layer_2 == 'SiO2':
            n_layer_2 = import_n('SiO2.csv')
        if layer_2 == 'Si':
            n_layer_2 = lambda wavelenth : 4.149+1j*0.0426
            
        t=(2*n_layer_1(lam))/(n_layer_1(lam)+n_layer_2(lam))
        return t
    
    def r_prim(lam, dtype=complex):
        r_prim=(r("hBN", "SiO2",lam)+r("SiO2", "Si",lam)*exp(-2j*beta_j(d_SiO2,"SiO2",lam)))/(1+r("hBN", "SiO2",lam)*r("SiO2", "Si",lam)*exp(-2j*beta_j(d_SiO2,"SiO2",lam)))
      #  print("r_prim=" + str(r_prim)) 
        return r_prim
        
    def r_bis(lam, dtype=complex):
        r_bis=(r("MoS2","hBN",lam)+r_prim(lam)*exp(-2j*beta_j(d_hBN,"hBN",lam)))/(1+r("MoS2", "hBN",lam)*r_prim(lam)*exp(-2j*beta_j(d_hBN,"hBN",lam)))
      #  print("r_bis=" + str(r_bis)) 
        return r_bis
        
    def r_tri(lam, dtype=complex):
        r_tri=(r("air","MoS2",lam)+r_bis(lam)*exp(-2j*beta_j(TMD_thickness,"MoS2",lam)))/(1+r("air", "MoS2",lam)*r_bis(lam)*exp(-2j*beta_j(TMD_thickness,"MoS2",lam)))
      #  print("r_tri=" + str(r_tri)) 
      #  print("r_bis=" + str(r_bis)) 
      #  print("r_01=" + str(r("air","MoS2",lam)))
      #  print("beta_j1=" + str(beta_j(TMD_thickness,"MoS2",lam))) 
        return r_tri
    
    def F_abs(x,lam):
        #print("r(air,MoS2,lam)="+ str(r("air","MoS2",lam)))
        F_abs=(t("air","MoS2",lam)*((exp(-1j*beta_x(x,lam))+r_bis(lam)*exp(-1j*(2*beta_j(TMD_thickness,"MoS2",lam)-beta_x(x,lam))))/(1+r_bis(lam)*r("air","MoS2",lam)*exp(-2j*beta_j(TMD_thickness,"MoS2",lam)))))
        return F_abs
    
    def F_sc(x,lam):
        F_sc=(t("MoS2","air",lam)*((exp(-1j*beta_x(x,lam))+r_bis(lam)*exp(-1j*(2*beta_j(TMD_thickness,"MoS2",lam)-beta_x(x,lam))))/(1+r_bis(lam)*r("air","MoS2",lam)*exp(-2j*beta_j(TMD_thickness,"MoS2",lam)))))
        return F_sc

    num=lambda_max-lambda_min+1    
    lambdas = linspace(lambda_min, lambda_max, num)
    
    intensity=[]

    file = open('enchancement_factor-literature-od_x.txt', 'w')
    for lam in lambdas:
        def integrand(x):
           # return (abs(F_abs(x,lam_0)*F_sc(x,lam)))**2
            return (abs(F_abs(x,lam_0)*F_sc(x,lam)))**2
            
        ans, err = quad(integrand,0,TMD_thickness)
        intensity.append(ans)
        print("integrand="+str(ans))
        #intensity.append(abs(r_tri(lam))**2)
        file.write(str(lam))
        file.write(';')
        file.write(str(ans))
        file.write('\n')
        
    plt.figure()
    plt.plot(lambdas, intensity, 'blue')
    plt.xlabel('Wavelength (nm)')
    plt.xlim(lambda_min, lambda_max)
    plt.ylabel('Intensity')
    plt.title('Intensity from air/TMD/hBN/SiO2/Si')
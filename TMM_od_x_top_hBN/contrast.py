from __future__ import division, print_function, absolute_import

from numpy import pi, linspace, inf, array, zeros
from scipy.interpolate import interp1d
import matplotlib
from . import color
import numpy as np
from matplotlib import pyplot as plt
import csv
degree = pi/180
from .tmm_core import coh_tmm


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
        layer_n_list.append([wavelength_array[i],n_array[i]+k_array[i]])

    layer_n_data = array(layer_n_list)
    layer_n_fn = interp1d(layer_n_data[:,0], layer_n_data[:,1], kind='linear')
    return layer_n_fn

def intensity_TMD_hBN_SiO2_d(TMD, number_of_layers, hBN_bottom_d = 300, SiO2_d=300, lambda_min = 360, lambda_max = 1000, lam0=514, nhBN=2.2): 
    #funkcja liczaca kontrast w funkcji dlugosci fali dla zadanej grubosci hBN
    Si_n_fn = import_n('Si.csv')
    SiO2_n_fn = import_n('SiO2_2.csv')
    hBN_n_fn = lambda wavelength : nhBN#zmienilam
    # air refractive index
    air_n_fn = lambda wavelength : 1
    
    if TMD == 'MoS2':
        TMD_n_fn = import_n('MoS2.csv')
        layer_thickness = number_of_layers*0.65
    if TMD == 'MoS2-ML':
        TMD_n_fn = import_n('MoS2-ML.csv')
        layer_thickness = number_of_layers*0.65
    if TMD == 'MoSe2':
        TMD_n_fn = import_n('MoSe2.csv')
        layer_thickness = number_of_layers*0.7
    if TMD == 'MoSe2-ML':
        TMD_n_fn = import_n('MoSe2-ML.csv')
        layer_thickness = number_of_layers*0.7
        
    if TMD == 'WS2':
        TMD_n_fn = import_n('WS2.csv')
        layer_thickness = number_of_layers*0.8
    if TMD == 'WS2-ML':
        TMD_n_fn = import_n('WS2-ML.csv')
        layer_thickness = number_of_layers*0.8
    
    if TMD == 'WSe2':
        TMD_n_fn = import_n('WSe2.csv')
        layer_thickness = number_of_layers*0.7
    if TMD == 'WSe2-ML':
        TMD_n_fn = import_n('WSe2-ML.csv')
        layer_thickness = number_of_layers*0.7
    
    n_fn_list = [air_n_fn, hBN_n_fn, TMD_n_fn, hBN_n_fn, SiO2_n_fn, Si_n_fn]
    th_0 = 0
     
    num=lambda_max-lambda_min+1
    
    lambdas = linspace(lambda_min, lambda_max, num)
    
    d_list = [inf, 0,layer_thickness, hBN_bottom_d, SiO2_d, inf]
    #d_list_1 = [inf, hBN_top_d, hBN_bottom_d, SiO2_d, inf]
    d_list_1 = [inf,0, hBN_bottom_d, SiO2_d, inf]
    #reflectances_BG = color.calc_reflectances(n_fn_list_1, d_list_1, th_0)
    reflectances = color.calc_reflectances(TMD, number_of_layers, n_fn_list, d_list, th_0,'s', lambda_min, lambda_max, lam0)
    #contrast=[]
    file = open('contrast_TMD_hBN_SiO2_d-enhancement_factor.txt', 'w')
    R=[]
    #R_BG=[]
    
    for i,lam in enumerate(reflectances[:,0]):
        file.write(str(lam))
        file.write(';')
        file.write(str((reflectances[i,1])))
       # print("Ench="+str(reflectances[i,1]))
        file.write('\n')
        #contrast.append((reflectances[i,1]-reflectances_BG[i,1])/(reflectances[i,1]+reflectances_BG[i,1]))
        R.append(reflectances[i,1])
        #R_BG.append(reflectances_BG[i,1])
    file.close()

    font = {'size'   : 6}

    matplotlib.rc('font', **font)
    
    plt.figure(num=None, figsize=(2.953, 3.15), dpi=120, facecolor='w', edgecolor='k')
   # plt.plot(lambdas, R, 'blue', lambdas, R_BG,'purple', lambdas, contrast, 'red')
    plt.plot(lambdas, R)
    plt.xlabel('Wavelength (nm)', fontsize=6)
    plt.xlim(lambda_min, lambda_max)
    plt.ylabel('Enhancement function',fontsize=6)
    plt.title('Enhancement function for air/'+ str(number_of_layers)+'L-MoS2/'+ str(hBN_bottom_d) +'nm-hBN/'+ str(SiO2_d) +'nm-SiO2/Si \n(R-blue, R_BG-purple, Contrast-red)', fontsize=6)
    plt.savefig("hBN_thickness=" + str(hBN_bottom_d) + "nm.png")

def intensity_TMD_hBN_wavelength(TMD, nhBN=2.2, wavelength=550, lam0=514, min_hBN_thickness = 0, max_hBN_thickness = 400, SiO2_d=200, number_of_layers = 1):
    #funkcja zwracająca kontrast w funkcji grubosci hBN dla zadanej dlugosci fali
    Si_n_fn = import_n('Si.csv')
    # SiO2 refractive index (approximate): 1.46 regardless of wavelength
    SiO2_n_fn = import_n('SiO2_2.csv')
    # air refractive index
    if TMD == 'MoS2':
        TMD_n_fn = import_n('MoS2.csv')
        layer_thickness = number_of_layers*0.65
    if TMD == 'MoS2-ML':
        TMD_n_fn = import_n('MoS2-ML.csv')
        layer_thickness = number_of_layers*0.65
    if TMD == 'MoSe2':
        TMD_n_fn = import_n('MoSe2.csv')
        layer_thickness = number_of_layers*0.7
    if TMD == 'MoSe2-ML':
        TMD_n_fn = import_n('MoSe2-ML.csv')
        layer_thickness = number_of_layers*0.7
        
    if TMD == 'WS2':
        TMD_n_fn = import_n('WS2.csv')
        layer_thickness = number_of_layers*0.8
    if TMD == 'WS2-ML':
        TMD_n_fn = import_n('WS2-ML.csv')
        layer_thickness = number_of_layers*0.8
    
    if TMD == 'WSe2':
        TMD_n_fn = import_n('WSe2.csv')
        layer_thickness = number_of_layers*0.7
    if TMD == 'WSe2-ML':
        TMD_n_fn = import_n('WSe2-ML.csv')
        layer_thickness = number_of_layers*0.7
    
    air_n_fn = lambda wavelength : 1
    hBN_n_fn = lambda wavelength : nhBN #zmienilam

    n_fn_list = [air_n_fn, hBN_n_fn,TMD_n_fn, hBN_n_fn, SiO2_n_fn, Si_n_fn]
    th_0 = 0
    
    num = max_hBN_thickness
    hBN_thickness_list = linspace(min_hBN_thickness, max_hBN_thickness, num)
    
    intensity = []
    num_layers = len(n_fn_list)
    
    file = open('intensity_TMD_'+str(TMD)+'_nhBN_'+str(nhBN)+'_SiO2_d_'+str(SiO2_d)+'_lam0_'+str(lam0)+'_lam_'+str(wavelength)+'.txt', 'w')
    
    for hBN_d in hBN_thickness_list:
        d_list = [inf, 0, layer_thickness, hBN_d, SiO2_d, inf]
        n_list = [n_fn_list[i](wavelength) for i in range(num_layers)]
        n_list_abs = [n_fn_list[i](lam0) for i in range(num_layers)]
       
        I=coh_tmm(TMD, number_of_layers,'s', n_list, d_list, th_0, wavelength, n_list_abs, lam0)['Enh']
        intensity.append(I)
       # R_BG = coh_tmm('s', n_list_1, d_list_1, th_0, wavelength)['R']
        
        file.write(str(hBN_d))
        file.write(';')
        file.write(str(I))
        #file.write(str(-(reflectances[i,1]-reflectances_BG[i,1])/(reflectances[i,1]+reflectances_BG[i,1])))
        file.write('\n')
       # contrast.append(-(reflectances[i,1]-reflectances_BG[i,1])/(reflectances[i,1]+reflectances_BG[i,1]))
        
        
    file.close()    
         
    plt.figure()
    plt.plot(hBN_thickness_list, intensity)
    """   plt.xlabel('hBN thickness (nm)')
    plt.xlim(min_hBN_thickness, max_hBN_thickness)
    plt.ylabel('Enhancement function')
    plt.title('Enhancement function for air/TMD/hBN/'+ str(SiO2_d)+ 'nm-SiO2/Si')
    """
def intensity_TMD_hBN_SiO2(TMD, hBN_d=300, SiO2_d=300, lambda_min = 360, lambda_max = 1000, number_of_layers = 1, nhBN=2.1, lam0=532): 
    #funkcja liczaca intensywnosc w funkcji dlugosci fali dla zadanej grubosci hBN
    
    Si_n_fn = import_n('Si.csv')
    SiO2_n_fn = import_n('SiO2.csv')
    # SiO2 refractive index (approximate): 1.46 regardless of wavelength
    hBN_n_fn = lambda wavelength : nhBN#1.65
    # air refractive index
    air_n_fn = lambda wavelength : 1
    
    if TMD == 'MoS2':
        TMD_n_fn = import_n('MoS2.csv')
        layer_thickness = number_of_layers*0.65
    if TMD == 'MoS2-ML':
        TMD_n_fn = import_n('MoS2-ML.csv')
        layer_thickness = number_of_layers*0.65
    
    if TMD == 'MoSe2':
        TMD_n_fn = import_n('MoSe2.csv')
        layer_thickness = number_of_layers*0.7
    if TMD == 'MoSe2-ML':
        TMD_n_fn = import_n('MoSe2-ML.csv')
        layer_thickness = number_of_layers*0.7
        
    if TMD == 'WS2':
        TMD_n_fn = import_n('WS2.csv')
        layer_thickness = number_of_layers*0.8
    if TMD == 'WS2-ML':
        TMD_n_fn = import_n('WS2-ML.csv')
        layer_thickness = number_of_layers*0.8
    
    if TMD == 'WSe2':
        TMD_n_fn = import_n('WSe2.csv')
        layer_thickness = number_of_layers*0.7
    if TMD == 'WSe2-ML':
        TMD_n_fn = import_n('WSe2-ML.csv')
        layer_thickness = number_of_layers*0.7
    
    n_fn_list = [air_n_fn, hBN_n_fn, TMD_n_fn, hBN_n_fn, SiO2_n_fn, Si_n_fn]
    th_0 = 0
   
    num=lambda_max-lambda_min+1   
    
    lambdas = linspace(lambda_min, lambda_max, num)
    
    d_list = [inf, 0, layer_thickness, hBN_d, SiO2_d, inf]
    reflectances = color.calc_reflectances(TMD, number_of_layers, n_fn_list, d_list, th_0,'s', lambda_min, lambda_max, lam0)
    contrast=[]
    file = open('Intensity from air_'+str(TMD)+'_'+str(hBN_d)+'nm-hBN_'+str(SiO2_d)+'nm-SiO2_Si.txt', 'w')
    
    for i,lam in enumerate(reflectances[:,0]):
        file.write(str(lam))
        file.write(';')
        file.write(str((reflectances[i,1])))
        #file.write(str((reflectances[i,1]-reflectances_BG[i,1])/(reflectances[i,1]+reflectances_BG[i,1])))
        file.write('\n')
        contrast.append((reflectances[i,1]))
        #contrast.append((reflectances[i,1]-reflectances_BG[i,1])/(reflectances[i,1]+reflectances_BG[i,1]))
    
    file.close()
     
    plt.figure()
    plt.plot(lambdas, contrast, 'blue')
    plt.xlabel('Wavelength (nm)')
    #plt.xlim(1.9,2.07)
    plt.xlim(lambda_min, lambda_max)
    plt.ylabel('Intensity')
    plt.title('Intensity from air/'+str(TMD)+'/'+str(hBN_d)+'nm-hBN/'+str(SiO2_d)+'nm-SiO2/Si')  
    plt.savefig('Intensity from air_'+str(TMD)+'_'+str(hBN_d)+'nm-hBN_'+str(SiO2_d)+'nm-SiO2_Si.png')

def intensity_TMD_SiO2_wavelength(TMD, nhBN=2.1, wavelength=532, lam0=532, min_SiO2_thickness = 0, max_SiO2_thickness = 400, number_of_layers = 1):
    #funkcja zwracająca kontrast w funkcji grubosci SiO2 dla zadanej dlugosci fali
    Si_n_fn = import_n('Si.csv')
    # SiO2 refractive index (approximate): 1.46 regardless of wavelength
    SiO2_n_fn = import_n('SiO2.csv')
    # air refractive index
    if TMD == 'MoS2':
        TMD_n_fn = import_n('MoS2.csv')
        layer_thickness = number_of_layers*0.65
    if TMD == 'MoS2-ML':
        TMD_n_fn = import_n('MoS2-ML.csv')
        layer_thickness = number_of_layers*0.65
    if TMD == 'MoSe2':
        TMD_n_fn = import_n('MoSe2.csv')
        layer_thickness = number_of_layers*0.7
    if TMD == 'MoSe2-ML':
        TMD_n_fn = import_n('MoSe2-ML.csv')
        layer_thickness = number_of_layers*0.7
        
    if TMD == 'WS2':
        TMD_n_fn = import_n('WS2.csv')
        layer_thickness = number_of_layers*0.8
    if TMD == 'WS2-ML':
        TMD_n_fn = import_n('WS2-ML.csv')
        layer_thickness = number_of_layers*0.8
    
    if TMD == 'WSe2':
        TMD_n_fn = import_n('WSe2.csv')
        layer_thickness = number_of_layers*0.7
    if TMD == 'WSe2-ML':
        TMD_n_fn = import_n('WSe2-ML.csv')
        layer_thickness = number_of_layers*0.7
    
    air_n_fn = lambda wavelength : 1
    hBN_n_fn = lambda wavelength : nhBN #zmienilam

    n_fn_list = [air_n_fn, hBN_n_fn,TMD_n_fn, hBN_n_fn, SiO2_n_fn, Si_n_fn]
    th_0 = 0
    
    num = max_SiO2_thickness
    SiO2_thickness_list = linspace(min_SiO2_thickness, max_SiO2_thickness, num)
    
    intensity = []
    num_layers = len(n_fn_list)
    hBN_top_d=0
    hBN_bottom_d=0
    file = open('intensity_abs_TMD_'+str(TMD)+'_lam0_'+str(lam0)+'_lam_'+str(wavelength)+'.txt', 'w')
    
    for SiO2_d in SiO2_thickness_list:
        d_list = [inf, hBN_top_d, layer_thickness, hBN_bottom_d, SiO2_d, inf]
        n_list = [n_fn_list[i](wavelength) for i in range(num_layers)]
        n_list_abs = [n_fn_list[i](lam0) for i in range(num_layers)]
       
        I=coh_tmm(TMD, number_of_layers,'s', n_list, d_list, th_0, wavelength, n_list_abs, lam0)['Enh']
        intensity.append(I)
       # R_BG = coh_tmm('s', n_list_1, d_list_1, th_0, wavelength)['R']
        
        file.write(str(SiO2_d))
        file.write(';')
        file.write(str(I))
        #file.write(str(-(reflectances[i,1]-reflectances_BG[i,1])/(reflectances[i,1]+reflectances_BG[i,1])))
        file.write('\n')
       # contrast.append(-(reflectances[i,1]-reflectances_BG[i,1])/(reflectances[i,1]+reflectances_BG[i,1]))
        
        
    file.close()    
         
    plt.figure()
    plt.plot(SiO2_thickness_list, intensity)
    """   plt.xlabel('hBN thickness (nm)')
    plt.xlim(min_hBN_thickness, max_hBN_thickness)
    plt.ylabel('Enhancement function')
    plt.title('Enhancement function for air/TMD/hBN/'+ str(SiO2_d)+ 'nm-SiO2/Si')
    """

def intensity_map_hBN_SiO2(TMD,number_of_layers=1, SiO2_d=300, nhBN=2.1, max_hBN_thickness = 400, lam0=532, lambda_min=450, lambda_max=800): 
    #funkcja rysujaca mape intensywnosci w funkcji grubosci hBN dla roznych dlugosci fali na SiO2 o podanej grubosci
    
    Si_n_fn = import_n('Si.csv')
    SiO2_n_fn = import_n('SiO2_2.csv')
    # SiO2 refractive index (approximate): 1.46 regardless of wavelength
    hBN_n_fn = lambda wavelength : nhBN
    # air refractive index
       
    air_n_fn = lambda wavelength : 1
    
    if TMD == 'MoS2':
        TMD_n_fn = import_n('MoS2.csv')
        layer_thickness = number_of_layers*0.65
    if TMD == 'MoS2-ML':
        TMD_n_fn = import_n('MoS2-ML.csv')
        layer_thickness = number_of_layers*0.65
    if TMD == 'MoSe2':
        TMD_n_fn = import_n('MoSe2.csv')
        layer_thickness = number_of_layers*0.7
    if TMD == 'MoSe2-ML':
        TMD_n_fn = import_n('MoSe2-ML.csv')
        layer_thickness = number_of_layers*0.7
        
    if TMD == 'WS2':
        TMD_n_fn = import_n('WS2.csv')
        layer_thickness = number_of_layers*0.8
    if TMD == 'WS2-ML':
        TMD_n_fn = import_n('WS2-ML.csv')
        layer_thickness = number_of_layers*0.8
    
    if TMD == 'WSe2':
        TMD_n_fn = import_n('WSe2.csv')
        layer_thickness = number_of_layers*0.7
    if TMD == 'WSe2-ML':
        TMD_n_fn = import_n('WSe2-ML.csv')
        layer_thickness = number_of_layers*0.7
    
    n_fn_list = [air_n_fn, hBN_n_fn, TMD_n_fn, hBN_n_fn, SiO2_n_fn, Si_n_fn]
    th_0 = 0   
    
    num = max_hBN_thickness
    hBN_thickness_list = linspace(0, max_hBN_thickness, num)
    
    nrows = lambda_max-lambda_min+1
    ncols = len(hBN_thickness_list)

    grid = np.zeros((nrows, ncols), dtype=np.float)
    
    j=0
    file = open('intensity_map_'+str(TMD)+'_lam0_'+str(lam0)+'_SiO2_d_'+str(SiO2_d)+'.txt', 'w')
    for hBN_d in hBN_thickness_list:
        
        d_list = [inf, 0, layer_thickness, hBN_d, SiO2_d, inf]
        reflectances = color.calc_reflectances(TMD, number_of_layers, n_fn_list, d_list, th_0,'s', lambda_min, lambda_max, lam0)
        
        for i,lam in enumerate(reflectances[:,0]):
            grid[i,j] = reflectances[i,1]
            
            file.write(str(j))
            file.write(';')
            file.write(str(i+lambda_min))
            file.write(';')
            file.write(str(grid[i,j]))
            #file.write(str((reflectances[i,1]-reflectances_BG[i,1])/(reflectances[i,1]+reflectances_BG[i,1])))
            file.write('\n')
        j+=1    
    file.close()  
    
    plt.figure()
    plt.title('SiO2 thickness = '+ str(SiO2_d)+ ' nm, lambda0=' + str(lam0) + ' nm')
    plt.imshow(grid, interpolation='nearest', 
            extent=(0, max_hBN_thickness, lambda_max, lambda_min)) #cmap=cm.gist_rainbow
    plt.gca().invert_yaxis()
    plt.xlabel('hBN thickness (nm)')
    plt.ylabel('Wavelength (nm)')
    plt.colorbar()
    plt.show()
    
def intensity_map_hBN_top_SiO2(TMD,number_of_layers=1, SiO2_d=300, nhBN=2.1, max_hBN_thickness = 400, lam0=532, lambda_min=450, lambda_max=800): 
    #funkcja rysujaca mape intensywnosci w funkcji grubosci górnego hBN dla roznych dlugosci fali na SiO2 o podanej grubosci
    
    Si_n_fn = import_n('Si.csv')
    SiO2_n_fn = import_n('SiO2_2.csv')
    # SiO2 refractive index (approximate): 1.46 regardless of wavelength
    hBN_n_fn = lambda wavelength : nhBN
    # air refractive index
       
    air_n_fn = lambda wavelength : 1
    
    if TMD == 'MoS2':
        TMD_n_fn = import_n('MoS2.csv')
        layer_thickness = number_of_layers*0.65
    if TMD == 'MoS2-ML':
        TMD_n_fn = import_n('MoS2-ML.csv')
        layer_thickness = number_of_layers*0.65
    if TMD == 'MoSe2':
        TMD_n_fn = import_n('MoSe2.csv')
        layer_thickness = number_of_layers*0.7
    if TMD == 'MoSe2-ML':
        TMD_n_fn = import_n('MoSe2-ML.csv')
        layer_thickness = number_of_layers*0.7
        
    if TMD == 'WS2':
        TMD_n_fn = import_n('WS2.csv')
        layer_thickness = number_of_layers*0.8
    if TMD == 'WS2-ML':
        TMD_n_fn = import_n('WS2-ML.csv')
        layer_thickness = number_of_layers*0.8
    
    if TMD == 'WSe2':
        TMD_n_fn = import_n('WSe2.csv')
        layer_thickness = number_of_layers*0.7
    if TMD == 'WSe2-ML':
        TMD_n_fn = import_n('WSe2-ML.csv')
        layer_thickness = number_of_layers*0.7
    
    n_fn_list = [air_n_fn, hBN_n_fn, TMD_n_fn, SiO2_n_fn, Si_n_fn]
    th_0 = 0   
    
    num = max_hBN_thickness
    hBN_thickness_list = linspace(0, max_hBN_thickness, num)
    
    nrows = lambda_max-lambda_min+1
    ncols = len(hBN_thickness_list)

    grid = np.zeros((nrows, ncols), dtype=np.float)
    
    j=0
    
    file = open('intensity_map_hBN_top_'+str(TMD)+'_lam0_'+str(lam0)+'_SiO2_d_'+str(SiO2_d)+'.txt', 'w')
    for hBN_d in hBN_thickness_list:
        
        d_list = [inf, hBN_d, layer_thickness, SiO2_d, inf]
        reflectances = color.calc_reflectances(TMD, number_of_layers, n_fn_list, d_list, th_0,'s', lambda_min, lambda_max, lam0)
        
        for i,lam in enumerate(reflectances[:,0]):
            grid[i,j] = reflectances[i,1]
            
            file.write(str(j))
            file.write(';')
            file.write(str(i+lambda_min))
            file.write(';')
            file.write(str(grid[i,j]))
            #file.write(str((reflectances[i,1]-reflectances_BG[i,1])/(reflectances[i,1]+reflectances_BG[i,1])))
            file.write('\n')
        j+=1    
    file.close()
    
    plt.figure()
    plt.title('SiO2 thickness = '+ str(SiO2_d)+ ' nm, lambda0=' + str(lam0) + ' nm')
    plt.imshow(grid, interpolation='nearest', 
            extent=(0, max_hBN_thickness, lambda_max, lambda_min)) #cmap=cm.gist_rainbow
    plt.gca().invert_yaxis()
    plt.xlabel('hBN thickness (nm)')
    plt.ylabel('Wavelength (nm)')
    plt.colorbar()
    plt.show()    
 #---------------------------NOWE WERSJE PROGRAMU--------------------------------   
def intensity_hBN_TMD_hBN_SiO2(TMD, hBN_top_d=30, hBN_bottom_d=300, SiO2_d=300, lambda_min = 360, lambda_max = 1000, number_of_layers = 1, nhBN=2.1, lam0=532): 
    #funkcja liczaca intensywnosc w funkcji dlugosci fali dla zadanej grubosci hBN
    
    Si_n_fn = import_n('Si.csv')
    SiO2_n_fn = import_n('SiO2.csv')
    # SiO2 refractive index (approximate): 1.46 regardless of wavelength
    hBN_n_fn = lambda wavelength : nhBN#1.65
    # air refractive index
    air_n_fn = lambda wavelength : 1
    
    if TMD == 'MoS2':
        TMD_n_fn = import_n('MoS2.csv')
        layer_thickness = number_of_layers*0.65
    if TMD == 'MoS2-ML':
        TMD_n_fn = import_n('MoS2-ML.csv')
        layer_thickness = number_of_layers*0.65
    
    if TMD == 'MoSe2':
        TMD_n_fn = import_n('MoSe2.csv')
        layer_thickness = number_of_layers*0.7
    if TMD == 'MoSe2-ML':
        TMD_n_fn = import_n('MoSe2-ML.csv')
        layer_thickness = number_of_layers*0.7
        
    if TMD == 'WS2':
        TMD_n_fn = import_n('WS2.csv')
        layer_thickness = number_of_layers*0.8
    if TMD == 'WS2-ML':
        TMD_n_fn = import_n('WS2-ML.csv')
        layer_thickness = number_of_layers*0.8
    
    if TMD == 'WSe2':
        TMD_n_fn = import_n('WSe2.csv')
        layer_thickness = number_of_layers*0.7
    if TMD == 'WSe2-ML':
        TMD_n_fn = import_n('WSe2-ML.csv')
        layer_thickness = number_of_layers*0.7
    
    n_fn_list = [air_n_fn, hBN_n_fn, TMD_n_fn, hBN_n_fn, SiO2_n_fn, Si_n_fn]
    th_0 = 0
   
    num=lambda_max-lambda_min+1   
    
    lambdas = linspace(lambda_min, lambda_max, num)
    
    d_list = [inf, hBN_top_d, layer_thickness, hBN_bottom_d, SiO2_d, inf]
    reflectances = color.calc_reflectances(TMD, number_of_layers, n_fn_list, d_list, th_0,'s', lambda_min, lambda_max, lam0)
    contrast=[]
    file = open('Intensity from air_'+str(hBN_top_d)+'nm-hBN_'+str(TMD)+'_'+str(hBN_bottom_d)+'nm-hBN_'+str(SiO2_d)+'nm-SiO2_Si.txt', 'w')
    
    for i,lam in enumerate(reflectances[:,0]):
        file.write(str(lam))
        file.write(';')
        file.write(str((reflectances[i,1])))
        #file.write(str((reflectances[i,1]-reflectances_BG[i,1])/(reflectances[i,1]+reflectances_BG[i,1])))
        file.write('\n')
        contrast.append((reflectances[i,1]))
        #contrast.append((reflectances[i,1]-reflectances_BG[i,1])/(reflectances[i,1]+reflectances_BG[i,1]))
    
    file.close()
     
    plt.figure()
    plt.plot(lambdas, contrast, 'blue')
    plt.xlabel('Wavelength (nm)')
    #plt.xlim(1.9,2.07)
    plt.xlim(lambda_min, lambda_max)
    plt.ylabel('Intensity')
    plt.title('Intensity from air/'+str(hBN_top_d)+'nm-hBN/'+str(TMD)+'/'+str(hBN_bottom_d)+'nm-hBN/'+str(SiO2_d)+'nm-SiO2/Si')  
    plt.savefig('Intensity from air_'+str(hBN_top_d)+'nm-hBN_'+str(TMD)+'_'+str(hBN_bottom_d)+'nm-hBN_'+str(SiO2_d)+'nm-SiO2_Si.png')

def intensity_map_TMD_SiO2(TMD,number_of_layers=1, max_SiO2_thickness = 400, lam0=532, lambda_min=450, lambda_max=800): 
    #funkcja rysujaca mape intensywnosci w funkcji grubosci hBN dla roznych dlugosci fali na SiO2 o podanej grubosci
    
    Si_n_fn = import_n('Si.csv')
    SiO2_n_fn = import_n('SiO2.csv')
    hBN_n_fn = lambda wavelength : 2.1
    # air refractive index       
    air_n_fn = lambda wavelength : 1
    
    if TMD == 'MoS2':
        TMD_n_fn = import_n('MoS2.csv')
        layer_thickness = number_of_layers*0.65
    if TMD == 'MoS2-ML':
        TMD_n_fn = import_n('MoS2-ML.csv')
        layer_thickness = number_of_layers*0.65
    if TMD == 'MoSe2':
        TMD_n_fn = import_n('MoSe2.csv')
        layer_thickness = number_of_layers*0.7
    if TMD == 'MoSe2-ML':
        TMD_n_fn = import_n('MoSe2-ML.csv')
        layer_thickness = number_of_layers*0.7
        
    if TMD == 'WS2':
        TMD_n_fn = import_n('WS2.csv')
        layer_thickness = number_of_layers*0.8
    if TMD == 'WS2-ML':
        TMD_n_fn = import_n('WS2-ML.csv')
        layer_thickness = number_of_layers*0.8
    
    if TMD == 'WSe2':
        TMD_n_fn = import_n('WSe2.csv')
        layer_thickness = number_of_layers*0.7
    if TMD == 'WSe2-ML':
        TMD_n_fn = import_n('WSe2-ML.csv')
        layer_thickness = number_of_layers*0.7
    
    n_fn_list = [air_n_fn, hBN_n_fn, TMD_n_fn, hBN_n_fn, SiO2_n_fn, Si_n_fn]
    th_0 = 0   
    
    num = max_SiO2_thickness
    SiO2_thickness_list = linspace(0, max_SiO2_thickness, num)
    
    nrows = lambda_max-lambda_min+1
    ncols = len(SiO2_thickness_list)

    grid = np.zeros((nrows, ncols), dtype=np.float)
    
    j=0
    
    hBN_top_d=0
    hBN_bottom_d=0
    
    file = open('Enhancement map WS2 SiO2 lam0=532nm.txt', 'w')

    for SiO2_d in SiO2_thickness_list:
        
        d_list = [inf, hBN_top_d, layer_thickness, hBN_bottom_d, SiO2_d, inf]
        reflectances = color.calc_reflectances(TMD, number_of_layers, n_fn_list, d_list, th_0,'s', lambda_min, lambda_max, lam0)
        
        for i,lam in enumerate(reflectances[:,0]):
            grid[i,j] = reflectances[i,1]
            file.write(str(j))
            file.write(';')
            file.write(str(i+lambda_min))
            file.write(';')
            file.write(str(grid[i,j]))
            #file.write(str((reflectances[i,1]-reflectances_BG[i,1])/(reflectances[i,1]+reflectances_BG[i,1])))
            file.write('\n')
        j+=1    
    file.close()
    
    plt.figure()
    plt.title('lambda0=' + str(lam0) + ' nm')
    plt.imshow(grid, interpolation='nearest', 
            extent=(0, max_SiO2_thickness, lambda_max, lambda_min)) #cmap=cm.gist_rainbow
    plt.gca().invert_yaxis()
    plt.xlabel('SiO2 thickness (nm)')
    plt.ylabel('Wavelength (nm)')
    plt.colorbar()
    plt.show()
    
def intensity_map_TMD_Al_SiO2(TMD,number_of_layers=1, max_Al_thickness = 400, lam0=532, lambda_min=450, lambda_max=800, SiO2_d=300): 
    #funkcja rysujaca mape intensywnosci w funkcji grubosci hBN dla roznych dlugosci fali na SiO2 o podanej grubosci
    
    Si_n_fn = import_n('Si.csv')
    SiO2_n_fn = import_n('SiO2.csv')
    Al_n_fn = import_n('Al.csv')
    hBN_n_fn = lambda wavelength : 2.1
    # air refractive index       
    air_n_fn = lambda wavelength : 1
    
    if TMD == 'MoS2':
        TMD_n_fn = import_n('MoS2.csv')
        layer_thickness = number_of_layers*0.65
    if TMD == 'MoS2-ML':
        TMD_n_fn = import_n('MoS2-ML.csv')
        layer_thickness = number_of_layers*0.65
    if TMD == 'MoSe2':
        TMD_n_fn = import_n('MoSe2.csv')
        layer_thickness = number_of_layers*0.7
    if TMD == 'MoSe2-ML':
        TMD_n_fn = import_n('MoSe2-ML.csv')
        layer_thickness = number_of_layers*0.7
        
    if TMD == 'WS2':
        TMD_n_fn = import_n('WS2.csv')
        layer_thickness = number_of_layers*0.8
    if TMD == 'WS2-ML':
        TMD_n_fn = import_n('WS2-ML.csv')
        layer_thickness = number_of_layers*0.8
    
    if TMD == 'WSe2':
        TMD_n_fn = import_n('WSe2.csv')
        layer_thickness = number_of_layers*0.7
    if TMD == 'WSe2-ML':
        TMD_n_fn = import_n('WSe2-ML.csv')
        layer_thickness = number_of_layers*0.7
    
    n_fn_list = [air_n_fn, hBN_n_fn, TMD_n_fn, Al_n_fn, SiO2_n_fn, Si_n_fn]
    th_0 = 0   
    
    num = max_Al_thickness
    Al_thickness_list = linspace(0, max_Al_thickness, num)
    
    nrows = lambda_max-lambda_min+1
    ncols = len(Al_thickness_list)

    grid = np.zeros((nrows, ncols), dtype=np.float)
    
    j=0
    
    hBN_top_d=0
    
    file = open('Enhancement map MoS2_Al_SiO2 lam0=532nm.txt', 'w')

    for Al_d in Al_thickness_list:
        
        d_list = [inf, hBN_top_d, layer_thickness, Al_d, SiO2_d, inf]
        reflectances = color.calc_reflectances(TMD, number_of_layers, n_fn_list, d_list, th_0,'s', lambda_min, lambda_max, lam0)
        
        for i,lam in enumerate(reflectances[:,0]):
            grid[i,j] = reflectances[i,1]
            file.write(str(j))
            file.write(';')
            file.write(str(i+lambda_min))
            file.write(';')
            file.write(str(grid[i,j]))
            #file.write(str((reflectances[i,1]-reflectances_BG[i,1])/(reflectances[i,1]+reflectances_BG[i,1])))
            file.write('\n')
        j+=1    
    file.close()
    
    plt.figure()
    plt.title('lambda0=' + str(lam0) + ' nm')
    plt.imshow(grid, interpolation='nearest', 
            extent=(0, max_Al_thickness, lambda_max, lambda_min)) #cmap=cm.gist_rainbow
    plt.gca().invert_yaxis()
    plt.xlabel('Al thickness (nm)')
    plt.ylabel('Wavelength (nm)')
    plt.colorbar()
    plt.show()
    
def intensity_map_TMD_air_Si(TMD,number_of_layers=1, max_SiO2_thickness = 400, lam0=532, lambda_min=450, lambda_max=800): 
    #funkcja rysujaca mape intensywnosci w funkcji grubosci hBN dla roznych dlugosci fali na SiO2 o podanej grubosci
    
    Si_n_fn = import_n('Si.csv')
    SiO2_n_fn = import_n('SiO2.csv')
    hBN_n_fn = lambda wavelength : 2.1
    # air refractive index       
    air_n_fn = lambda wavelength : 1
    
    if TMD == 'MoS2':
        TMD_n_fn = import_n('MoS2.csv')
        layer_thickness = number_of_layers*0.65
    if TMD == 'MoS2-ML':
        TMD_n_fn = import_n('MoS2-ML.csv')
        layer_thickness = number_of_layers*0.65
    if TMD == 'MoSe2':
        TMD_n_fn = import_n('MoSe2.csv')
        layer_thickness = number_of_layers*0.7
    if TMD == 'MoSe2-ML':
        TMD_n_fn = import_n('MoSe2-ML.csv')
        layer_thickness = number_of_layers*0.7
        
    if TMD == 'WS2':
        TMD_n_fn = import_n('WS2.csv')
        layer_thickness = number_of_layers*0.8
    if TMD == 'WS2-ML':
        TMD_n_fn = import_n('WS2-ML.csv')
        layer_thickness = number_of_layers*0.8
    
    if TMD == 'WSe2':
        TMD_n_fn = import_n('WSe2.csv')
        layer_thickness = number_of_layers*0.7
    if TMD == 'WSe2-ML':
        TMD_n_fn = import_n('WSe2-ML.csv')
        layer_thickness = number_of_layers*0.7
    
    n_fn_list = [air_n_fn, hBN_n_fn, TMD_n_fn, hBN_n_fn, air_n_fn, Si_n_fn]
    th_0 = 0   
    
    num = max_SiO2_thickness
    SiO2_thickness_list = linspace(0, max_SiO2_thickness, num)
    
    nrows = lambda_max-lambda_min+1
    ncols = len(SiO2_thickness_list)

    grid = np.zeros((nrows, ncols), dtype=np.float)
    
    j=0
    
    hBN_top_d=0
    hBN_bottom_d=0
    
    file = open('Enhancement map MoS2 air_Si lam0=532nm.txt', 'w')

    for SiO2_d in SiO2_thickness_list:
        
        d_list = [inf, hBN_top_d, layer_thickness, hBN_bottom_d, SiO2_d, inf]
        reflectances = color.calc_reflectances(TMD, number_of_layers, n_fn_list, d_list, th_0,'s', lambda_min, lambda_max, lam0)
        
        for i,lam in enumerate(reflectances[:,0]):
            grid[i,j] = reflectances[i,1]
            file.write(str(j))
            file.write(';')
            file.write(str(i+lambda_min))
            file.write(';')
            file.write(str(grid[i,j]))
            #file.write(str((reflectances[i,1]-reflectances_BG[i,1])/(reflectances[i,1]+reflectances_BG[i,1])))
            file.write('\n')
        j+=1    
    file.close()
    
    plt.figure()
    plt.title('lambda0=' + str(lam0) + ' nm')
    plt.imshow(grid, interpolation='nearest', 
            extent=(0, max_SiO2_thickness, lambda_max, lambda_min)) #cmap=cm.gist_rainbow
    plt.gca().invert_yaxis()
    plt.xlabel('air thickness (nm)')
    plt.ylabel('Wavelength (nm)')
    plt.colorbar()
    plt.show()

    
   
def NL_map_TMD_SiO2(TMD,max_number_of_layers=20, max_SiO2_thickness = 100, lam0=532, lam_vac=532): 
    #funkcja rysujaca mape intensywnosci w funkcji grubosci SiO2 dla roznych ilosc wasrstw TMD
    
    Si_n_fn = import_n('Si.csv')
    SiO2_n_fn = import_n('SiO2.csv')
    hBN_n_fn = lambda wavelength : 2.1
    # air refractive index       
    air_n_fn = lambda wavelength : 1
    
    if TMD == 'MoS2':
        TMD_n_fn = import_n('MoS2.csv')
        layer_thickness = 0.65
    if TMD == 'MoS2-ML':
        TMD_n_fn = import_n('MoS2-ML.csv')
        layer_thickness = 0.65
    if TMD == 'MoSe2':
        TMD_n_fn = import_n('MoSe2.csv')
        layer_thickness = 0.7
    if TMD == 'MoSe2-ML':
        TMD_n_fn = import_n('MoSe2-ML.csv')
        layer_thickness = 0.7
        
    if TMD == 'WS2':
        TMD_n_fn = import_n('WS2.csv')
        layer_thickness = 0.8
    if TMD == 'WS2-ML':
        TMD_n_fn = import_n('WS2-ML.csv')
        layer_thickness = 0.8
    
    if TMD == 'WSe2':
        TMD_n_fn = import_n('WSe2.csv')
        layer_thickness = 0.7
    if TMD == 'WSe2-ML':
        TMD_n_fn = import_n('WSe2-ML.csv')
        layer_thickness = 0.7
    
    n_fn_list = [air_n_fn, hBN_n_fn, TMD_n_fn, hBN_n_fn, SiO2_n_fn, Si_n_fn]
    th_0 = 0   
    

    num = max_SiO2_thickness+1
    num_d = max_number_of_layers+1
    SiO2_thickness_list = linspace(0, max_SiO2_thickness, num)
    number_of_layers_list = linspace(0,max_number_of_layers,num_d)
    
    nrows = len(number_of_layers_list)
    ncols = len(SiO2_thickness_list)

    grid = np.zeros((nrows, ncols), dtype=np.float)
    
    i=0
    j=0
    
    hBN_top_d=0
    hBN_bottom_d=0
    num_layers = len(n_fn_list)
    
    lam_vac_list=[]
    lam_vac_list.append(lam_vac)
    for lam_vac in lam_vac_list:
        #lam_vac=514
        n_list = [n_fn_list[i](lam_vac) for i in range(num_layers)]
        #print(n_list)
        n_list_abs = [n_fn_list[i](lam0) for i in range(num_layers)]
        #print(n_list_abs)
    pol='s'
    
    file = open('Intensity_map_NL_SiO2_d_NL_max='+str(max_number_of_layers)+'SiO2_d_max'+str(max_SiO2_thickness)+'lam0'+str(lam0)+'lam_vac'+str(lam_vac)+'.txt', 'w')
     
    for number_of_layers in number_of_layers_list:
        j=0
        for SiO2_d in SiO2_thickness_list:
            #print(number_of_layers)
            d_list = [inf, hBN_top_d, number_of_layers*layer_thickness, hBN_bottom_d, SiO2_d, inf]
            
            I=coh_tmm(TMD, number_of_layers, pol, n_list, d_list, th_0, lam_vac, n_list_abs, lam0)['Enh']
            grid[i,j] = I
            file.write(str(j))
            file.write(';')
            file.write(str(i))
            file.write(';')
            file.write(str(I))
            #file.write(str((reflectances[i,1]-reflectances_BG[i,1])/(reflectances[i,1]+reflectances_BG[i,1])))
            file.write('\n')
           # print(I)
            j+=1    
        i+=1
        
    file.close()
    
    plt.figure()
    plt.title('excitation wavelength=' + str(lam0) + ' nm, emission wavelength=' + str(lam_vac) + ' nm')
    plt.imshow(grid, interpolation='nearest', 
            extent=( 0,max_SiO2_thickness, max_number_of_layers, 0)) #cmap=cm.gist_rainbow
    plt.gca().invert_yaxis()
    plt.xlabel('SiO2 thickness (nm)')
    plt.ylabel('Number of layers')    
    plt.colorbar()
    plt.show()
    

def intensity_NL_TMD_SiO2(TMD,max_number_of_layers=200, SiO2_d = 300, lam0=532, lam_vac=532): 
    #funkcja rysujaca zależnoć intensywnosci w funkcji iloci warstw TMD dla zadanych lam0 i lam
    
    Si_n_fn = import_n('Si.csv')
    SiO2_n_fn = import_n('SiO2.csv')
    hBN_n_fn = lambda wavelength : 2.1
    # air refractive index       
    air_n_fn = lambda wavelength : 1
    
    if TMD == 'MoS2':
        TMD_n_fn = import_n('MoS2.csv')
        layer_thickness = 0.65
    if TMD == 'MoS2-ML':
        TMD_n_fn = import_n('MoS2-ML.csv')
        layer_thickness = 0.65
    if TMD == 'MoSe2':
        TMD_n_fn = import_n('MoSe2.csv')
        layer_thickness = 0.7
    if TMD == 'MoSe2-ML':
        TMD_n_fn = import_n('MoSe2-ML.csv')
        layer_thickness = 0.7
        
    if TMD == 'WS2':
        TMD_n_fn = import_n('WS2.csv')
        layer_thickness = 0.8
    if TMD == 'WS2-ML':
        TMD_n_fn = import_n('WS2-ML.csv')
        layer_thickness = 0.8
    
    if TMD == 'WSe2':
        TMD_n_fn = import_n('WSe2.csv')
        layer_thickness = 0.7
    if TMD == 'WSe2-ML':
        TMD_n_fn = import_n('WSe2-ML.csv')
        layer_thickness = 0.7
    
    n_fn_list = [air_n_fn, hBN_n_fn, TMD_n_fn, hBN_n_fn, SiO2_n_fn, Si_n_fn]
    th_0 = 0   
    
    num_d = max_number_of_layers+1
    number_of_layers_list = linspace(0,max_number_of_layers,num_d)
    
    hBN_top_d=0
    hBN_bottom_d=0
    num_layers = len(n_fn_list)
    
    lam_vac_list=[]
    lam_vac_list.append(lam_vac)
    for lam_vac in lam_vac_list:
        #lam_vac=514
        n_list = [n_fn_list[i](lam_vac) for i in range(num_layers)]
        #print(n_list)
        n_list_abs = [n_fn_list[i](lam0) for i in range(num_layers)]
        #print(n_list_abs)
    pol='s'
    intensity=[]
    file = open('Intensity_NL_SiO2_d_NL_max='+str(max_number_of_layers)+'SiO2_d'+str(SiO2_d)+'.txt', 'w')
     
    for number_of_layers in number_of_layers_list:

            #print(number_of_layers)
        d_list = [inf, hBN_top_d, number_of_layers*layer_thickness, hBN_bottom_d, SiO2_d, inf]
        #print(number_of_layers*layer_thickness)
        I=coh_tmm(TMD, number_of_layers, pol, n_list, d_list, th_0, lam_vac, n_list_abs, lam0)['Enh']
        intensity.append(I)
        file.write(str(number_of_layers))
        file.write(';')
        file.write(str(I))
        #file.write(str((reflectances[i,1]-reflectances_BG[i,1])/(reflectances[i,1]+reflectances_BG[i,1])))
        file.write('\n')

    file.close()
    
    plt.figure()
    plt.title('excitation wavelength=' + str(lam0) + ' nm, emission wavelength=' + str(lam_vac) + ' nm')
    plt.plot(number_of_layers_list,intensity)
    plt.xlabel('Number of layers')
    plt.ylabel('Intensity')    
    plt.show()
    
def intensity_NL_TMD_Al_SiO2(TMD,max_number_of_layers=200, Al_d=150, SiO2_d = 300, lam0=532, lam_vac=532): 
    #funkcja rysujaca zależnoć intensywnosci w funkcji iloci warstw TMD dla zadanych lam0 i lam
    
    Si_n_fn = import_n('Si.csv')
    SiO2_n_fn = import_n('SiO2.csv')
    Al_n_fn=import_n('Al.csv')
    
    hBN_n_fn = lambda wavelength : 2.1
    # air refractive index       
    air_n_fn = lambda wavelength : 1
    
    if TMD == 'MoS2':
        TMD_n_fn = import_n('MoS2.csv')
        layer_thickness = 0.65
    if TMD == 'MoS2-ML':
        TMD_n_fn = import_n('MoS2-ML.csv')
        layer_thickness = 0.65
    if TMD == 'MoSe2':
        TMD_n_fn = import_n('MoSe2.csv')
        layer_thickness = 0.7
    if TMD == 'MoSe2-ML':
        TMD_n_fn = import_n('MoSe2-ML.csv')
        layer_thickness = 0.7
        
    if TMD == 'WS2':
        TMD_n_fn = import_n('WS2.csv')
        layer_thickness = 0.8
    if TMD == 'WS2-ML':
        TMD_n_fn = import_n('WS2-ML.csv')
        layer_thickness = 0.8
    
    if TMD == 'WSe2':
        TMD_n_fn = import_n('WSe2.csv')
        layer_thickness = 0.7
    if TMD == 'WSe2-ML':
        TMD_n_fn = import_n('WSe2-ML.csv')
        layer_thickness = 0.7
    
    n_fn_list = [air_n_fn, hBN_n_fn, TMD_n_fn, Al_n_fn, SiO2_n_fn, Si_n_fn]
    th_0 = 0   
    
    num_d = max_number_of_layers+1
    number_of_layers_list = linspace(0,max_number_of_layers,num_d)
    
    hBN_top_d=0
    num_layers = len(n_fn_list)
    
    lam_vac_list=[]
    lam_vac_list.append(lam_vac)
    for lam_vac in lam_vac_list:
        #lam_vac=514
        n_list = [n_fn_list[i](lam_vac) for i in range(num_layers)]
        #print(n_list)
        n_list_abs = [n_fn_list[i](lam0) for i in range(num_layers)]
        #print(n_list_abs)
    pol='s'
    intensity=[]
    file = open('Intensity_NL_Al_d_SiO2_NL_max='+str(max_number_of_layers)+'SiO2_d'+str(SiO2_d)+'.txt', 'w')
     
    for number_of_layers in number_of_layers_list:

            #print(number_of_layers)
        d_list = [inf, hBN_top_d, number_of_layers*layer_thickness, Al_d, SiO2_d, inf]
        #print(number_of_layers*layer_thickness)
        I=coh_tmm(TMD, number_of_layers, pol, n_list, d_list, th_0, lam_vac, n_list_abs, lam0)['Enh']
        intensity.append(I)
        file.write(str(number_of_layers))
        file.write(';')
        file.write(str(I))
        #file.write(str((reflectances[i,1]-reflectances_BG[i,1])/(reflectances[i,1]+reflectances_BG[i,1])))
        file.write('\n')

    file.close()
    
    plt.figure()
    plt.title('excitation wavelength=' + str(lam0) + ' nm, emission wavelength=' + str(lam_vac) + ' nm')
    plt.plot(number_of_layers_list,intensity)
    plt.xlabel('Number of layers')
    plt.ylabel('Intensity')    
    plt.show()

def intensity_NL_TMD_air_Si(TMD,max_number_of_layers=200, air_d=1500, lam0=532, lam_vac=532): 
    #funkcja rysujaca zależnoć intensywnosci w funkcji iloci warstw TMD dla zadanych lam0 i lam
    
    Si_n_fn = import_n('Si.csv')
    SiO2_n_fn = import_n('SiO2.csv')
    hBN_n_fn = lambda wavelength : 2.1
    # air refractive index       
    air_n_fn = lambda wavelength : 1
    
    if TMD == 'MoS2':
        TMD_n_fn = import_n('MoS2.csv')
        layer_thickness = 0.65
    if TMD == 'MoS2-ML':
        TMD_n_fn = import_n('MoS2-ML.csv')
        layer_thickness = 0.65
    if TMD == 'MoSe2':
        TMD_n_fn = import_n('MoSe2.csv')
        layer_thickness = 0.7
    if TMD == 'MoSe2-ML':
        TMD_n_fn = import_n('MoSe2-ML.csv')
        layer_thickness = 0.7
        
    if TMD == 'WS2':
        TMD_n_fn = import_n('WS2.csv')
        layer_thickness = 0.8
    if TMD == 'WS2-ML':
        TMD_n_fn = import_n('WS2-ML.csv')
        layer_thickness = 0.8
    
    if TMD == 'WSe2':
        TMD_n_fn = import_n('WSe2.csv')
        layer_thickness = 0.7
    if TMD == 'WSe2-ML':
        TMD_n_fn = import_n('WSe2-ML.csv')
        layer_thickness = 0.7
    
    n_fn_list = [air_n_fn, hBN_n_fn, TMD_n_fn, air_n_fn, SiO2_n_fn, Si_n_fn]
    th_0 = 0   
    
    num_d = max_number_of_layers+1
    number_of_layers_list = linspace(0,max_number_of_layers,num_d)
    
    hBN_top_d=0
    SiO2_d=0
    num_layers = len(n_fn_list)
    
    lam_vac_list=[]
    lam_vac_list.append(lam_vac)
    for lam_vac in lam_vac_list:
        #lam_vac=514
        n_list = [n_fn_list[i](lam_vac) for i in range(num_layers)]
        #print(n_list)
        n_list_abs = [n_fn_list[i](lam0) for i in range(num_layers)]
        #print(n_list_abs)
    pol='s'
    intensity=[]
    file = open('Intensity_NL_air_Si_NL_max='+str(max_number_of_layers)+'.txt', 'w')
     
    for number_of_layers in number_of_layers_list:

            #print(number_of_layers)
        d_list = [inf, hBN_top_d, number_of_layers*layer_thickness, air_d, SiO2_d, inf]
        #print(number_of_layers*layer_thickness)
        I=coh_tmm(TMD, number_of_layers, pol, n_list, d_list, th_0, lam_vac, n_list_abs, lam0)['Enh']
        intensity.append(I)
        file.write(str(number_of_layers))
        file.write(';')
        file.write(str(I))
        #file.write(str((reflectances[i,1]-reflectances_BG[i,1])/(reflectances[i,1]+reflectances_BG[i,1])))
        file.write('\n')

    file.close()
    
    plt.figure()
    plt.title('excitation wavelength=' + str(lam0) + ' nm, emission wavelength=' + str(lam_vac) + ' nm')
    plt.plot(number_of_layers_list,intensity)
    plt.xlabel('Number of layers')
    plt.ylabel('Intensity')    
    plt.show()
    
   
def intensity_NL_TMD_suspended(TMD,max_number_of_layers=200, lam0=532, lam_vac=532): 
    #funkcja rysujaca zależnoć intensywnosci w funkcji iloci warstw TMD dla zadanych lam0 i lam
    
    Si_n_fn = import_n('Si.csv')
    SiO2_n_fn = import_n('SiO2.csv')
    hBN_n_fn = lambda wavelength : 2.1
    # air refractive index       
    air_n_fn = lambda wavelength : 1
    
    if TMD == 'MoS2':
        TMD_n_fn = import_n('MoS2.csv')
        layer_thickness = 0.65
    if TMD == 'MoS2-ML':
        TMD_n_fn = import_n('MoS2-ML.csv')
        layer_thickness = 0.65
    if TMD == 'MoSe2':
        TMD_n_fn = import_n('MoSe2.csv')
        layer_thickness = 0.7
    if TMD == 'MoSe2-ML':
        TMD_n_fn = import_n('MoSe2-ML.csv')
        layer_thickness = 0.7
        
    if TMD == 'WS2':
        TMD_n_fn = import_n('WS2.csv')
        layer_thickness = 0.8
    if TMD == 'WS2-ML':
        TMD_n_fn = import_n('WS2-ML.csv')
        layer_thickness = 0.8
    
    if TMD == 'WSe2':
        TMD_n_fn = import_n('WSe2.csv')
        layer_thickness = 0.7
    if TMD == 'WSe2-ML':
        TMD_n_fn = import_n('WSe2-ML.csv')
        layer_thickness = 0.7
    
    n_fn_list = [air_n_fn, hBN_n_fn, TMD_n_fn, air_n_fn]
    th_0 = 0   
    
    num_d = max_number_of_layers+1
    number_of_layers_list = linspace(0,max_number_of_layers,num_d)
    
    hBN_top_d=0
    num_layers = len(n_fn_list)
    
    lam_vac_list=[]
    lam_vac_list.append(lam_vac)
    for lam_vac in lam_vac_list:
        #lam_vac=514
        n_list = [n_fn_list[i](lam_vac) for i in range(num_layers)]
        #print(n_list)
        n_list_abs = [n_fn_list[i](lam0) for i in range(num_layers)]
        #print(n_list_abs)
    pol='s'
    intensity=[]
    file = open('Intensity_NL_air_Si_NL_max='+str(max_number_of_layers)+'.txt', 'w')
     
    for number_of_layers in number_of_layers_list:

            #print(number_of_layers)
        d_list = [inf, hBN_top_d, number_of_layers*layer_thickness, inf]
        #print(number_of_layers*layer_thickness)
        I=coh_tmm(TMD, number_of_layers, pol, n_list, d_list, th_0, lam_vac, n_list_abs, lam0)['Enh']
        intensity.append(I)
        file.write(str(number_of_layers))
        file.write(';')
        file.write(str(I))
        #file.write(str((reflectances[i,1]-reflectances_BG[i,1])/(reflectances[i,1]+reflectances_BG[i,1])))
        file.write('\n')

    file.close()
    
    plt.figure()
    plt.title('excitation wavelength=' + str(lam0) + ' nm, emission wavelength=' + str(lam_vac) + ' nm')
    plt.plot(number_of_layers_list,intensity)
    plt.xlabel('Number of layers')
    plt.ylabel('Intensity')    
    plt.show()
  

def map_hBN_TMD_hBN_SiO2(TMD,max_hBN_top_d=20, max_hBN_bottom_d=200, SiO2_d = 300, lam0=532, lam_vac=543, number_of_layers=1, nhBN=2.2): 
    #funkcja rysujaca mape intensywnosci w funkcji grubosci gornej i dolnej warstwy hBN
    
    Si_n_fn = import_n('Si.csv')
    SiO2_n_fn = import_n('SiO2.csv')
    hBN_n_fn = lambda wavelength : nhBN
    # air refractive index       
    air_n_fn = lambda wavelength : 1
    
    if TMD == 'MoS2':
        TMD_n_fn = import_n('MoS2.csv')
        layer_thickness = 0.65
    if TMD == 'MoS2-ML':
        TMD_n_fn = import_n('MoS2-ML.csv')
        layer_thickness = 0.65
    if TMD == 'MoSe2':
        TMD_n_fn = import_n('MoSe2.csv')
        layer_thickness = 0.7
    if TMD == 'MoSe2-ML':
        TMD_n_fn = import_n('MoSe2-ML.csv')
        layer_thickness = 0.7
        
    if TMD == 'WS2':
        TMD_n_fn = import_n('WS2.csv')
        layer_thickness = 0.8
    if TMD == 'WS2-ML':
        TMD_n_fn = import_n('WS2-ML.csv')
        layer_thickness = 0.8
    
    if TMD == 'WSe2':
        TMD_n_fn = import_n('WSe2.csv')
        layer_thickness = 0.7
    if TMD == 'WSe2-ML':
        TMD_n_fn = import_n('WSe2-ML.csv')
        layer_thickness = 0.7
    
    n_fn_list = [air_n_fn, hBN_n_fn, TMD_n_fn, hBN_n_fn, SiO2_n_fn, Si_n_fn]
    th_0 = 0   
    

    num = max_hBN_bottom_d+1
    num_d = max_hBN_top_d+1
    hBN_bottom_d_list = linspace(0, max_hBN_bottom_d, num)
    hBN_top_d_list = linspace(0,max_hBN_top_d,num_d)
    
    nrows = len(hBN_top_d_list)
    ncols = len(hBN_bottom_d_list)

    grid = np.zeros((nrows, ncols), dtype=np.float)
    
    i=0
    j=0

    num_layers = len(n_fn_list)
    
    lam_vac_list=[]
    lam_vac_list.append(lam_vac)
    for lam_vac in lam_vac_list:
        #lam_vac=514
        n_list = [n_fn_list[i](lam_vac) for i in range(num_layers)]
        n_list_abs = [n_fn_list[i](lam0) for i in range(num_layers)]
    pol='s'
    
    file = open('Intensity_map_hBN_top_hBN_bottom_hBN__SiO2_d='+str(SiO2_d)+'lam0'+str(lam0)+'lam'+str(lam_vac)+'.txt', 'w')
    for hBN_top_d in hBN_top_d_list:
        j=0
        for hBN_bottom_d in hBN_bottom_d_list:
            #print(number_of_layers)
            d_list = [inf, hBN_top_d, number_of_layers*layer_thickness, hBN_bottom_d, SiO2_d, inf]
            
            I=coh_tmm(TMD, number_of_layers, pol, n_list, d_list, th_0, lam_vac, n_list_abs, lam0)['Enh']
            grid[i,j] = I
            file.write(str(j))
            file.write(';')
            file.write(str(i))
            file.write(';')
            file.write(str(I))
            #file.write(str((reflectances[i,1]-reflectances_BG[i,1])/(reflectances[i,1]+reflectances_BG[i,1])))
            file.write('\n')
           # print(I)
            j+=1    
        i+=1
    file.close()
    
    plt.figure()
    plt.title('excitation wavelength=' + str(lam0) + ' nm, emission wavelength=' + str(lam_vac) + ' nm')
    plt.imshow(grid, interpolation='nearest', 
            extent=( 0,max_hBN_bottom_d,max_hBN_top_d, 0)) #cmap=cm.gist_rainbow
    plt.gca().invert_yaxis()
    plt.xlabel('bottom hBN thickness (nm)')
    plt.ylabel('top hBN thickness (nm)')    
    plt.colorbar()
    plt.show()
   
def map_SiO2_TMD_SiO2(TMD,max_SiO2_top_d=200, max_SiO2_bottom_d=200, lam0=532, lam_vac=532, number_of_layers=1): 
    #funkcja rysujaca mape intensywnosci w funkcji grubosci gornej i dolnej warstwy SiO2
    
    Si_n_fn = import_n('Si.csv')
    SiO2_n_fn = import_n('SiO2.csv')
    hBN_n_fn = lambda wavelength : 2.12
    # air refractive index       
    air_n_fn = lambda wavelength : 1
    
    if TMD == 'MoS2':
        TMD_n_fn = import_n('MoS2.csv')
        layer_thickness = 0.65
    if TMD == 'MoS2-ML':
        TMD_n_fn = import_n('MoS2-ML.csv')
        layer_thickness = 0.65
    if TMD == 'MoSe2':
        TMD_n_fn = import_n('MoSe2.csv')
        layer_thickness = 0.7
    if TMD == 'MoSe2-ML':
        TMD_n_fn = import_n('MoSe2-ML.csv')
        layer_thickness = 0.7
        
    if TMD == 'WS2':
        TMD_n_fn = import_n('WS2.csv')
        layer_thickness = 0.8
    if TMD == 'WS2-ML':
        TMD_n_fn = import_n('WS2-ML.csv')
        layer_thickness = 0.8
    
    if TMD == 'WSe2':
        TMD_n_fn = import_n('WSe2.csv')
        layer_thickness = 0.7
    if TMD == 'WSe2-ML':
        TMD_n_fn = import_n('WSe2-ML.csv')
        layer_thickness = 0.7
    
    n_fn_list = [air_n_fn, SiO2_n_fn, TMD_n_fn, SiO2_n_fn, SiO2_n_fn, Si_n_fn]
    th_0 = 0   
    

    num = max_SiO2_bottom_d+1
    num_d = max_SiO2_top_d+1
    SiO2_bottom_d_list = linspace(0, max_SiO2_bottom_d, num)
    SiO2_top_d_list = linspace(0,max_SiO2_top_d,num_d)
    
    nrows = len(SiO2_top_d_list)
    ncols = len(SiO2_bottom_d_list)

    grid = np.zeros((nrows, ncols), dtype=np.float)
    
    i=0
    j=0

    num_layers = len(n_fn_list)
    
    lam_vac_list=[]
    lam_vac_list.append(lam_vac)
    for lam_vac in lam_vac_list:
        #lam_vac=514
        n_list = [n_fn_list[i](lam_vac) for i in range(num_layers)]
        n_list_abs = [n_fn_list[i](lam0) for i in range(num_layers)]
    pol='s'
    SiO2_d=0
    
    file = open('Intensity_map_SiO2_top_SiO2_bottom_SiO2_top_max='+str(max_SiO2_top_d)+'SiO2_bottom_max'+str(max_SiO2_bottom_d)+'lam0'+str(lam0)+'.txt', 'w')
      
    for SiO2_top_d in SiO2_top_d_list:
        j=0
        for SiO2_bottom_d in SiO2_bottom_d_list:
            #print(number_of_layers)
            d_list = [inf, SiO2_top_d, number_of_layers*layer_thickness, SiO2_bottom_d, SiO2_d, inf]
            
            I=coh_tmm(TMD, number_of_layers, pol, n_list, d_list, th_0, lam_vac, n_list_abs, lam0)['Enh']
            grid[i,j] = I
            file.write(str(j))
            file.write(';')
            file.write(str(i))
            file.write(';')
            file.write(str(I))
            #file.write(str((reflectances[i,1]-reflectances_BG[i,1])/(reflectances[i,1]+reflectances_BG[i,1])))
            file.write('\n')
           # print(I)
            j+=1    
        i+=1
    file.close()
    plt.figure()
    plt.title('excitation wavelength=' + str(lam0) + ' nm, emission wavelength=' + str(lam_vac) + ' nm')
    plt.imshow(grid, interpolation='nearest', 
            extent=( 0,max_SiO2_bottom_d, max_SiO2_top_d, 0)) #cmap=cm.gist_rainbow
    plt.gca().invert_yaxis()
    plt.xlabel('bottom SiO2 thickness (nm)')
    plt.ylabel('top SiO2 thickness (nm)')    
    plt.colorbar()
    plt.show()
    
def intensity_SiO2_TMD_SiO2(TMD,max_SiO2_top_d=200, SiO2_bottom_d=260, lam0=532, lam_vac=532, number_of_layers=1): 
    #funkcja rysujaca mape intensywnosci w funkcji grubosci gornej i dolnej warstwy SiO2
    
    Si_n_fn = import_n('Si.csv')
    SiO2_n_fn = import_n('SiO2.csv')
    hBN_n_fn = lambda wavelength : 2.12
    # air refractive index       
    air_n_fn = lambda wavelength : 1
    
    if TMD == 'MoS2':
        TMD_n_fn = import_n('MoS2.csv')
        layer_thickness = 0.65
    if TMD == 'MoS2-ML':
        TMD_n_fn = import_n('MoS2-ML.csv')
        layer_thickness = 0.65
    if TMD == 'MoSe2':
        TMD_n_fn = import_n('MoSe2.csv')
        layer_thickness = 0.7
    if TMD == 'MoSe2-ML':
        TMD_n_fn = import_n('MoSe2-ML.csv')
        layer_thickness = 0.7
        
    if TMD == 'WS2':
        TMD_n_fn = import_n('WS2.csv')
        layer_thickness = 0.8
    if TMD == 'WS2-ML':
        TMD_n_fn = import_n('WS2-ML.csv')
        layer_thickness = 0.8
    
    if TMD == 'WSe2':
        TMD_n_fn = import_n('WSe2.csv')
        layer_thickness = 0.7
    if TMD == 'WSe2-ML':
        TMD_n_fn = import_n('WSe2-ML.csv')
        layer_thickness = 0.7
    
    n_fn_list = [air_n_fn, SiO2_n_fn, TMD_n_fn, SiO2_n_fn, SiO2_n_fn, Si_n_fn]
    th_0 = 0   
    
    num_d = max_SiO2_top_d+1
    SiO2_top_d_list = linspace(0,max_SiO2_top_d,num_d)

    num_layers = len(n_fn_list)
    
    lam_vac_list=[]
    lam_vac_list.append(lam_vac)
    for lam_vac in lam_vac_list:
        #lam_vac=514
        n_list = [n_fn_list[i](lam_vac) for i in range(num_layers)]
        n_list_abs = [n_fn_list[i](lam0) for i in range(num_layers)]
    pol='s'
    SiO2_d=0
    
    intensity=[]
    file = open('Intensity_SiO2_top_SiO2_bottom_SiO2_top_max='+str(max_SiO2_top_d)+'SiO2_bottom_max'+str(SiO2_bottom_d)+'lam0'+str(lam0)+'.txt'+'lam'+str(lam_vac), 'w')
      
    for SiO2_top_d in SiO2_top_d_list:
        d_list = [inf, SiO2_top_d, number_of_layers*layer_thickness, SiO2_bottom_d, SiO2_d, inf]
        
        I=coh_tmm(TMD, number_of_layers, pol, n_list, d_list, th_0, lam_vac, n_list_abs, lam0)['Enh']
        intensity.append(I)
        file.write(str(SiO2_top_d))
        file.write(';')
        file.write(str(I))
        #file.write(str((reflectances[i,1]-reflectances_BG[i,1])/(reflectances[i,1]+reflectances_BG[i,1])))
        file.write('\n')

    file.close()
    plt.figure()
    plt.plot(SiO2_top_d_list, intensity)
    plt.title('excitation wavelength=' + str(lam0) + ' nm, emission wavelength=' + str(lam_vac) + ' nm')
    plt.xlabel('Capping SiO2 thickness (nm)')
    plt.ylabel('Intensity')    
    plt.show()
    
    
def intensity_TMD_SiO2(TMD,max_SiO2_bottom_d=200, lam0=532, lam_vac=752, number_of_layers=1): 
    #funkcja rysujaca mape intensywnosci w funkcji grubosci dolnej warstwy SiO2
    
    Si_n_fn = import_n('Si.csv')
    SiO2_n_fn = import_n('SiO2.csv')
    hBN_n_fn = lambda wavelength : 2.12
    # air refractive index       
    air_n_fn = lambda wavelength : 1
    
    if TMD == 'MoS2':
        TMD_n_fn = import_n('MoS2.csv')
        layer_thickness = 0.65
    if TMD == 'MoS2-ML':
        TMD_n_fn = import_n('MoS2-ML.csv')
        layer_thickness = 0.65
    if TMD == 'MoSe2':
        TMD_n_fn = import_n('MoSe2.csv')
        layer_thickness = 0.7
    if TMD == 'MoSe2-ML':
        TMD_n_fn = import_n('MoSe2-ML.csv')
        layer_thickness = 0.7
        
    if TMD == 'WS2':
        TMD_n_fn = import_n('WS2.csv')
        layer_thickness = 0.8
    if TMD == 'WS2-ML':
        TMD_n_fn = import_n('WS2-ML.csv')
        layer_thickness = 0.8
    
    if TMD == 'WSe2':
        TMD_n_fn = import_n('WSe2.csv')
        layer_thickness = 0.7
    if TMD == 'WSe2-ML':
        TMD_n_fn = import_n('WSe2-ML.csv')
        layer_thickness = 0.7
    
    n_fn_list = [air_n_fn, SiO2_n_fn, TMD_n_fn, SiO2_n_fn, SiO2_n_fn, Si_n_fn]
    th_0 = 0   
    
    num_d = max_SiO2_bottom_d+1
    SiO2_bottom_d_list = linspace(0,max_SiO2_bottom_d,num_d)

    num_layers = len(n_fn_list)
    
    lam_vac_list=[]
    lam_vac_list.append(lam_vac)
    for lam_vac in lam_vac_list:
        #lam_vac=514
        n_list = [n_fn_list[i](lam_vac) for i in range(num_layers)]
        n_list_abs = [n_fn_list[i](lam0) for i in range(num_layers)]
    pol='s'
    SiO2_d=0
    SiO2_top_d=0
    
    intensity=[]
    file = open('Intensity_SiO2_bottom_SiO2_bottom_max'+str(max_SiO2_bottom_d)+'lam0'+str(lam0)+'lam'+str(lam_vac)+'.txt', 'w')
      
    for SiO2_bottom_d in SiO2_bottom_d_list:
        d_list = [inf, SiO2_top_d, number_of_layers*layer_thickness, SiO2_bottom_d, SiO2_d, inf]
        
        I=coh_tmm(TMD, number_of_layers, pol, n_list, d_list, th_0, lam_vac, n_list_abs, lam0)['Enh']
        intensity.append(I)
        file.write(str(SiO2_top_d))
        file.write(';')
        file.write(str(I))
        #file.write(str((reflectances[i,1]-reflectances_BG[i,1])/(reflectances[i,1]+reflectances_BG[i,1])))
        file.write('\n')

    file.close()
    plt.figure()
    plt.plot(SiO2_bottom_d_list, intensity)
    plt.title('excitation wavelength=' + str(lam0) + ' nm, emission wavelength=' + str(lam_vac) + ' nm')
    plt.xlabel('Substrate SiO2 thickness (nm)')
    plt.ylabel('Intensity')    
    plt.show()
"""
     
def intensity_TMD_SiO2_wavelength(TMD, nhBN=2.2, wavelength=550, lam0=532, min_hBN_thickness = 0, max_hBN_thickness = 400, SiO2_d=200, number_of_layers = 1):
    #funkcja zwracająca kontrast w funkcji grubosci SiO2 dla zadanej dlugosci fali
    Si_n_fn = import_n('Si.csv')
    # SiO2 refractive index (approximate): 1.46 regardless of wavelength
    SiO2_n_fn = import_n('SiO2.csv')
    # air refractive index
    if TMD == 'MoS2':
        TMD_n_fn = import_n('MoS2.csv')
        layer_thickness = number_of_layers*0.65
    if TMD == 'MoS2-ML':
        TMD_n_fn = import_n('MoS2-ML.csv')
        layer_thickness = number_of_layers*0.65
    if TMD == 'MoSe2':
        TMD_n_fn = import_n('MoSe2.csv')
        layer_thickness = number_of_layers*0.7
    if TMD == 'MoSe2-ML':
        TMD_n_fn = import_n('MoSe2-ML.csv')
        layer_thickness = number_of_layers*0.7
        
    if TMD == 'WS2':
        TMD_n_fn = import_n('WS2.csv')
        layer_thickness = number_of_layers*0.8
    if TMD == 'WS2-ML':
        TMD_n_fn = import_n('WS2-ML.csv')
        layer_thickness = number_of_layers*0.8
    
    if TMD == 'WSe2':
        TMD_n_fn = import_n('WSe2.csv')
        layer_thickness = number_of_layers*0.7
    if TMD == 'WSe2-ML':
        TMD_n_fn = import_n('WSe2-ML.csv')
        layer_thickness = number_of_layers*0.7
    
    air_n_fn = lambda wavelength : 1
    hBN_n_fn = lambda wavelength : nhBN #zmienilam

    n_fn_list = [air_n_fn, TMD_n_fn, hBN_n_fn, SiO2_n_fn, Si_n_fn]
    th_0 = 0
    
    num = max_hBN_thickness
    hBN_thickness_list = linspace(min_hBN_thickness, max_hBN_thickness, num)
    
    intensity = []
    num_layers = len(n_fn_list)
    
    file = open('intensity_TMD_'+str(TMD)+'_nhBN_'+str(nhBN)+'_SiO2_d_'+str(SiO2_d)+'_lam0_'+str(lam0)+'_lam_'+str(wavelength)+'.txt', 'w')
    
    for hBN_d in hBN_thickness_list:
        d_list = [inf, layer_thickness, hBN_d, SiO2_d, inf]
        n_list = [n_fn_list[i](wavelength) for i in range(num_layers)]
        n_list_abs = [n_fn_list[i](lam0) for i in range(num_layers)]
       
        I=coh_tmm(TMD, number_of_layers,'s', n_list, d_list, th_0, wavelength, n_list_abs, lam0)['Enh']
        intensity.append(I)
       # R_BG = coh_tmm('s', n_list_1, d_list_1, th_0, wavelength)['R']
        
        file.write(str(hBN_d))
        file.write(';')
        file.write(str(I))
        #file.write(str(-(reflectances[i,1]-reflectances_BG[i,1])/(reflectances[i,1]+reflectances_BG[i,1])))
        file.write('\n')
       # contrast.append(-(reflectances[i,1]-reflectances_BG[i,1])/(reflectances[i,1]+reflectances_BG[i,1]))
  
    file.close()    
         
    plt.figure()
    plt.plot(hBN_thickness_list, intensity)
"""
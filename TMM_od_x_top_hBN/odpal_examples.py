# -*- coding: utf-8 -*-
"""
Created on Thu Oct 12 09:40:57 2017

@author: Joanna
"""
from numpy import arange
import numpy as np

import TMM_od_x_top_hBN
#import tmm.contrast
#import tmm.examples
#import tmm.tmm2
#import tmm.tmm2.contrast #tmm2 umozliwia liczenie koloru, zakres 360-830nm bo taki zakres ma iluminant D65
import TMM_od_x_top_hBN.contrast #tmm3 umozliwia obliczenia dla WS2, WSe2 i ML TMD ale dla wezszego zakresu: 415-805
#import TMM_od_x_top_hBN.intensity
import TMM_od_x_top_hBN.intensity2
"""
tmm.examples.sample1()
tmm.examples.sample2()
tmm.examples.sample3()
tmm.examples.sample4()
#tmm.examples.sample5()
tmm.examples.sample6()
"""
#tmm.examples.sample8()
#tmm.contrast.reflectance_angle('MoS2')
#tmm.contrast.contrast_hBN_TMD_hBN('MoS2',10,1,140, 620, 700)
#tmm.contrast.contrast_TMD_hBN('MoS2',140)
#tmm.contrast.contrast_hBN_TMD_hBN_SiO2_d('MoS2',14,1,10,285, 620, 700)
#tmm.contrast.contrast_TMD_SiO2_d('MoS2', 2030, 460, 720, 1,20)
#tmm.contrast.transmitance_TMD('MoS2', 2030, 450, 700, 1)
#tmm.contrast.contrast_TMD_function_of_wavelength('MoS2',315)
#tmm.contrast.contrast_TMD_SiO2_d('MoS2', 320, 450, 900, 1)

#tmm.tmm2.contrast.contrast('MoS2')

#tmm.contrast.contrast('MoS2')
#tmm.tmm2.contrast.calc_color('MoS2',0)

#tmm.tmm2.contrast.calc_color_hBN_SiO2_Si(200)

#tmm.contrast.calc_color('MoS2',0)

#tmm.tmm3.contrast.contrast_hBN_SiO2('WS2-ML', 200)
#tmm.tmm3.contrast.reflectance_map_hBN_SiO2(200)

#tmm.tmm3.contrast.contrast_TMD_hBN('WS2',200, 500, 550, 1)
#tmm.tmm3.contrast.contrast_hBN_TMD_hBN_SiO2_d('WS2-ML',0, 1, 0, 320, 1.95, 2.1)

#tmm.tmm3.contrast.contrast_TMD_hBN_wavelength('MoS2',514,550,600,285,1)


#TMM_od_x.contrast.intensity_TMD_hBN_SiO2_d('WSe2-ML', number_of_layers=1, hBN_bottom_d = 63, SiO2_d=300, lambda_min = 415, lambda_max = 825, lam0=532, nhBN=2.12)
#TMM_od_x_top_hBN.intensity2.intensity('MoS2', d_SiO2=270, d_hBN=20, lam_0=514, lambda_min = 500, lambda_max = 501, number_of_layers = 10, th_0=0)
#TMM_od_x_top_hBN.contrast.intensity_TMD_hBN_SiO2('MoS2', hBN_d=20, SiO2_d=270, lambda_min = 500, lambda_max = 501, number_of_layers = 10, nhBN=2.2, lam0=514)
#TMM_od_x.intensity.intensity('MoS2',285, 280, 514, 550, 555)

#TMM_od_x.contrast.intensity_TMD_hBN_wavelength('WSe2',542,532,0,300,300,1,2.1) #nhBN=2.13 lub 2.17
#TMM_od_x.contrast.intensity_TMD_hBN_wavelength('WSe2-ML',644,633,0,300,300,1,2.1)
"""
max_nhBN=2.3
min_nhBN=2.2
nhBN_list=np.arange(min_nhBN,max_nhBN+0.001, 0.02)
for nhBN in nhBN_list:

    TMM_od_x.contrast.intensity_TMD_hBN_wavelength('WSe2', nhBN, wavelength=539, lam0=532, min_hBN_thickness = 0, max_hBN_thickness = 550, SiO2_d=300, number_of_layers = 1)
"""
#hereTMM_od_x_top_hBN.contrast.intensity_TMD_hBN_wavelength('MoSe2-ML', nhBN=2.2, wavelength=756, lam0=712, min_hBN_thickness = 50, max_hBN_thickness = 400, SiO2_d=83, number_of_layers = 1)
#TMM_od_x_top_hBN.intensity.intensity('MoS2', d_SiO2=300, d_hBN=0, lam_0=532, lambda_min = 700, lambda_max = 710, number_of_layers = 1, th_0=0, dtype=complex)
#TMM_od_x_top_hBN.contrast.intensity_hBN_TMD_hBN_SiO2('MoS2', hBN_top_d=100, hBN_bottom_d=300, SiO2_d=300, lambda_min = 400, lambda_max = 750, number_of_layers = 1, nhBN=2.1, lam0=532)

#TMM_od_x_top_hBN.contrast.intensity_TMD_SiO2_wavelength('WSe2',nhBN=2.1, wavelength=532, lam0=732, min_SiO2_thickness = 0, max_SiO2_thickness = 400, number_of_layers = 1)
#TMM_od_x_top_hBN.contrast.intensity_map_TMD_SiO2('WS2-ML',number_of_layers=1, max_SiO2_thickness = 400, lam0=532, lambda_min=430, lambda_max=820)
#TMM_od_x_top_hBN.contrast.intensity_NL_TMD_air_Si('MoS2-ML',max_number_of_layers=120, air_d=1500, lam0=532, lam_vac=542)
#TMM_od_x_top_hBN.contrast.intensity_NL_TMD_suspended('MoS2-ML',max_number_of_layers=120, lam0=532, lam_vac=542)

"""
TMM_od_x_top_hBN.contrast.intensity_NL_TMD_SiO2('MoS2-ML',max_number_of_layers=10, SiO2_d = 310, lam0=532, lam_vac=542)
TMM_od_x_top_hBN.contrast.intensity_NL_TMD_air_Si('MoS2-ML',max_number_of_layers=10, air_d=1500, lam0=532, lam_vac=542)
TMM_od_x_top_hBN.contrast.intensity_NL_TMD_Al_SiO2('MoS2-ML',max_number_of_layers=10, Al_d=150, SiO2_d = 300, lam0=532, lam_vac=542)

TMM_od_x_top_hBN.contrast.intensity_NL_TMD_SiO2('MoS2-ML',max_number_of_layers=5, SiO2_d = 310, lam0=532, lam_vac=660)
TMM_od_x_top_hBN.contrast.intensity_NL_TMD_air_Si('MoS2-ML',max_number_of_layers=5, air_d=1500, lam0=532, lam_vac=660)
TMM_od_x_top_hBN.contrast.intensity_NL_TMD_Al_SiO2('MoS2-ML',max_number_of_layers=5, Al_d=150, SiO2_d = 300, lam0=532, lam_vac=660)
#TMM_od_x_top_hBN.contrast.intensity_map_TMD_SiO2('MoS2-ML',number_of_layers=1, max_SiO2_thickness = 400, lam0=532, lambda_min=430, lambda_max=820)
#TMM_od_x_top_hBN.contrast.intensity_map_TMD_air_Si('MoS2',number_of_layers=1, max_SiO2_thickness = 2000, lam0=532, lambda_min=430, lambda_max=820)
#TMM_od_x_top_hBN.contrast.intensity_map_TMD_Al_SiO2('MoS2',number_of_layers=1, max_Al_thickness = 200, lam0=532, lambda_min=430, lambda_max=800, SiO2_d=820)
"""

#TMM_od_x_top_hBN.contrast.intensity_SiO2_TMD_SiO2('WSe2',max_SiO2_top_d=300, SiO2_bottom_d=260, lam0=532, lam_vac=752, number_of_layers=1)
#TMM_od_x_top_hBN.contrast.intensity_TMD_SiO2('WSe2',max_SiO2_bottom_d=300, lam0=532, lam_vac=752, number_of_layers=1)


#TMM_od_x_top_hBN.contrast.NL_map_TMD_SiO2('WSe2',max_number_of_layers=120, max_SiO2_thickness = 560, lam0=633, lam_vac=752) 
#TMM_od_x_top_hBN.contrast.intensity_NL_TMD_SiO2('WSe2',max_number_of_layers=300, SiO2_d = 270, lam0=532, lam_vac=532)
TMM_od_x_top_hBN.contrast.map_hBN_TMD_hBN_SiO2('MoS2-ML',max_hBN_top_d=600, max_hBN_bottom_d=600, SiO2_d = 300, lam0=532, lam_vac=653, number_of_layers=1, nhBN=2.2)
#TMM_od_x_top_hBN.contrast.map_hBN_TMD_hBN_SiO2('WS2-ML',max_hBN_top_d=600, max_hBN_bottom_d=600, SiO2_d = 300, lam0=532, lam_vac=620, number_of_layers=1, nhBN=2.12)
#TMM_od_x_top_hBN.contrast.map_hBN_TMD_hBN_SiO2('MoSe2-ML',max_hBN_top_d=600, max_hBN_bottom_d=600, SiO2_d = 300, lam0=532, lam_vac=756, number_of_layers=1, nhBN=2.12)
#TMM_od_x_top_hBN.contrast.map_hBN_TMD_hBN_SiO2('WSe2-ML',max_hBN_top_d=600, max_hBN_bottom_d=600, SiO2_d = 300, lam0=532, lam_vac=752, number_of_layers=1, nhBN=2.12)
#TMM_od_x_top_hBN.contrast.map_SiO2_TMD_SiO2('WSe2',max_SiO2_top_d=600, max_SiO2_bottom_d=600, lam0=532, lam_vac=532, number_of_layers=1)
#TMM_od_x_top_hBN.contrast.map_hBN_TMD_hBN_SiO2('WSe2',max_hBN_top_d=600, max_hBN_bottom_d=600, SiO2_d = 300, lam0=532, lam_vac=532, number_of_layers=1)
#TMM_od_x_top_hBN.contrast.map_hBN_TMD_hBN_SiO2('WSe2',max_hBN_top_d=600, max_hBN_bottom_d=600, SiO2_d = 100, lam0=532, lam_vac=752, number_of_layers=1)
"""
max_hBN_thickness=300
num = max_hBN_thickness/10+1
hBN_thickness_list = linspace(0, max_hBN_thickness, num)
for hBN_d in hBN_thickness_list:
    TMM_od_x.contrast.intensity_TMD_hBN_SiO2('WS2-ML', hBN_d, SiO2_d=300, lambda_min = 450, lambda_max = 800, number_of_layers = 1, nhBN=2.1, lam0=532)
"""

#TMM_od_x.contrast.intensity_TMD_hBN_SiO2('MoS2', hBN_d=290, SiO2_d=285, lambda_min = 500, lambda_max = 750, number_of_layers = 1, nhBN=2.1, lam0=514)


#TMM_od_x_top_hBN.contrast.intensity_map_hBN_SiO2('MoSe2-ML',number_of_layers=1, SiO2_d=80, nhBN=2.12, max_hBN_thickness = 600, lam0=532, lambda_min=430, lambda_max=820) 
#hereTMM_od_x_top_hBN.contrast.intensity_map_hBN_SiO2('MoSe2-ML',number_of_layers=1, SiO2_d=80, nhBN=2.2, max_hBN_thickness = 400, lam0=712, lambda_min=500, lambda_max=800) 
#TMD,number_of_layers=1, SiO2_d=300, nhBN=2.1, max_hBN_thickness = 400, lam0=532

#intensity(TMD, d_SiO2=300, d_hBN=20, lam_0=532, lambda_min = 360, lambda_max = 1000, number_of_layers = 1, th_0=0
#tmm.tmm3.contrast.contrast_hBN_TMD_hBN_SiO2_d('MoS2',0,1,280,285, 415, 800)
"""
tmm.examples.sample8()
tmm.examples.sample9()
tmm.examples.sample10()
"""

#tmm.examples.sample15()


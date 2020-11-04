from __future__ import division, print_function, absolute_import

from numpy import arange, array

from .tmm_core import coh_tmm


inf = float('inf')

def calc_reflectances(TMD, number_of_layers, n_fn_list, d_list, th_0, pol='s', lambda_min = 360, lambda_max = 1000, lam0=514):
    lam_vac_list = arange(lambda_min, lambda_max+1)

    num_layers = len(n_fn_list)

    Enh=[]

    for lam_vac in lam_vac_list:
        #lam_vac=514
        n_list = [n_fn_list[i](lam_vac) for i in range(num_layers)]
        n_list_abs = [n_fn_list[i](lam0) for i in range(num_layers)]
        #R = coh_tmm(pol, n_list, d_list, th_0,lam_vac)['R']
        Enhancement=coh_tmm(TMD, number_of_layers, pol, n_list, d_list, th_0, lam_vac, n_list_abs, lam0)['Enh']
        #print('Ench='+str(Enchancement))
        #final_answer.append([lam_vac,R])
        Enh.append([lam_vac,Enhancement])
    #final_answer = array(final_answer)
    Enh=array(Enh)

    return Enh
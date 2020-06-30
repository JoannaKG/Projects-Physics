import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd
from math import pi as PI
import traceback

def lorentzian(x, xc, A, w):
    return (2*A/PI)*(w/(4*(x-xc)**2 + w**2))

def func(x, *params):
    """
    Parameters:
        f - Single Lorentzian function to be repeated several times 
            depending on the length of *params
        x - M-length sequence, independent variables (horizontal axis)
        *params - array, parameters of curve in the form of the list with arguments in order:
                  [center_0, amplitude_0, width_0, center_1 ...] 
    Returns: 
        Multi-peak Lorentzian function 
    """
    y = np.zeros_like(x)
    for i in range(1, len(params), 3):
        xc = params[i]
        A = params[i+1]
        w = params[i+2]
        y = y + lorentzian(x, xc, A, w)
    y = y + params[0]
    return y

## Uses scipy curve_fit to optimise the lorentzian fitting
def fit_function(func, x, y, guess):
    """
    Parameters: 
        guess - array, parameters used for first iteration in scipy curve_fit 
                with arguments in order:
                [y_0, center_0, amplitude_0, width_0, center_1 ...]
        func - the model function to be fitted
        x - M-length sequence, independent variables (horizontal axis)
        y - M-length sequence, dependent variables (vertical axis)
        *params - array, parameters of curve with arguments in order:
                  [center_0, amplitude_0, width_0, center_1 ...] 
    Returns: 
        popt - optimal values 
    """
    popt, pcov = curve_fit(func, x, y, p0=guess, maxfev=4000, bounds=(lower_bounds(),upper_bounds()))
    fit = func(x, *popt)
    return (popt, fit)

def predict_fitting_parameters(filepath, n_peaks_to_find = 9, delta_E = 0.001):
    """
        Parameters:
            filepath - filepath to input file with data. Format of input data: x, y0, y1, y2 [...]
            n_peaks_to_find - integer, number of peaks 
            delta_E - float, value of energy distance between incident laser energies used to collect data;
                      change of energy affecting guess of new energy of Raman modes and new laser energy
        Returns:
            out_df - pandas DataFrame, with values of fitted parameters for each column
            Report.txt - file with reports: successed/not succesed fitting for each column
            Params_extracted.csv - file with values of fitted parameters for each successfully fitted column
    """
    df = pd.read_csv(filepath, sep=',', header = None)
    number_of_columns = df.shape[1]
    xs = df[0]
    #create labels for output DataFrame
    rows=[]
    rows.append('y0')
    for i in range(0,n_peaks_to_find): 
        rows.append('xc_'+str(i+1))
        rows.append('A_'+str(i+1))
        rows.append('w_'+str(i+1))
    out_df=pd.DataFrame(rows)
    
    folder = 'Figures bounds' #folder for output reports
    file = open(folder + '/Report.txt', 'w') #file with reports: successed/not succesed fitting for each column
    cols = ['parameters']
    for col in range(1,number_of_columns):
        ys = df[col]
        try:
            if col == 1:
                params, fit = fit_function(func, xs, ys, first_guess() )
            else:
                params[22]=params[22] + delta_E #change guess - Raman mode energy
                params[25]=params[25] + delta_E #change guess - excitation laser energy
                params, fit = fit_function(func, xs, ys, params)
            out_df[str(params[25])]=params
            cols.append(col)
            print('Column No. '+ str(col)+' fitted successfully! \n')
            file.write('Column No. '+ str(col)+' fitted successfully! \n')           
        except RuntimeError:
            print('Failure in fitting column No. '+ str(col)+'\n')
            file.write('Failure in fitting column No. '+ str(col)+'\n')
            print(traceback.format_exc())
    file.close()
    out_df.columns = cols
    writer = ('C:\\Users\\Joanna\Desktop\\PLE\\' + folder +'\\Params_extracted.csv')
    out_df.to_csv(writer,header=False,index=False)
    return out_df

def plot_fitted_curves(filepath, df_fitted_parameters):
    """
        Parameters:
            filepath - filepath to input file with data. Format of input data: x, y0, y1, y2 [...]
            df_fitted_parameters - DataFrame with fitting parameters 
                                   returned by function predict_fitting_parameters
    """
    df = pd.read_csv(filepath, sep=',', header = None)
    number_of_columns = df.shape[1]
    number_of_rows = df_fitted_parameters.shape[0]
    xs = df[0]
    folder = 'Figures bounds'

    for col in range(1,number_of_columns):
        if col == df_fitted_parameters.columns[col]:
            ys = df[col]
            plt.figure(figsize=(8,6))
            plt.title("Raman Scattering")
            
            for row in range(1, number_of_rows, 3):
                xc = df_fitted_parameters.loc[row,col] 
                A = df_fitted_parameters.loc[row+1,col] 
                w = df_fitted_parameters.loc[row+2,col] 
                plt.plot(xs, lorentzian(xs, xc, A, w)+df_fitted_parameters.loc[0,col], ls='-',label=row//3+1)
            plt.plot(xs,ys, lw=1, label='data', c='black')
            #plt.plot(xs, fit, 'r-', label='fit', c='red', lw=2, ls='--')
            plt.legend()
            plt.ylabel("Intensity (arb. units)")
            plt.xlabel("Raman Shift ($cm^{-1}$)")
            plt.title('E_exc = ' + str(round(df_fitted_parameters.loc[25,col],4)) + ' eV')
            plt.savefig(folder + '/E_exc = ' + str(round(df_fitted_parameters.loc[25,col],4)) + ' eV.png')
            plt.show()
        

def first_guess():
    """
        Returns:
            2-tuple of array-like, parameters used for first iteration in scipy curve_fit 
            in the form of the list with arguments in order:
            [y_0, center_0, amplitude_0, width_0, center_1 ...]. 
            
            For best result of curve_fit those parameters should be taken 
            from normal fitting in other program, eg. from OriginPro
    """
    guess = []
    #y0 - baseline
    guess.append(592)
    #peak1
    guess.append(1.67) #xc
    guess.append(3.5) #A
    guess.append(0.011) #w 
    #peak2
    guess.append(1.68) #xc
    guess.append(0.131) #A
    guess.append(0.0021) #w 
    #peak3
    guess.append(1.684) #xc
    guess.append(1.0) #A
    guess.append(0.005) #w 
    #peak4
    guess.append(1.691) #xc
    guess.append(1.0) #A
    guess.append(0.0056) #w 
    #peak5
    guess.append(1.6997) #xc
    guess.append(1.0) #A
    guess.append(0.007) #w 
    #peak6
    guess.append(1.703) #xc
    guess.append(0.75) #A
    guess.append(0.00395) #w 
    #peak7
    guess.append(1.720) #xc
    guess.append(1.75) #A
    guess.append(0.0062) #w 
    #peak8
    guess.append(1.7281) #xc
    guess.append(1.75) #A
    guess.append(0.0055) #w 
    #peak9 - laser
    guess.append(1.8271) #xc
    guess.append(0.58) #A
    guess.append(0.005) #w 
    
    return guess

def lower_bounds():
    """
        Returns:
            Array, parameters used for lower bounds during scipy curve_fit 
            in the form of the list with arguments in order:
            [y_0, center_0, amplitude_0, width_0, center_1 ...]. 
            
            For best result of curve_fit those parameters should be close to 
            expected values of fitting parameters
    """
    bound = []
    #y0 - baseline
    bound.append(-np.inf)
    #peak1
    bound.append(1.668) #xc
    bound.append(0) #A
    bound.append(0.01) #w 
    #peak2
    bound.append(1.679) #xc
    bound.append(0.13) #A
    bound.append(0.002) #w 
    #peak3
    bound.append(1.6824) #xc
    bound.append(0.8) #A
    bound.append(0.0046) #w 
    #peak4
    bound.append(1.6895) #xc
    bound.append(0.98) #A
    bound.append(0.0055) #w 
    #peak5
    bound.append(1.6996) #xc
    bound.append(0) #A
    bound.append(0.0001) #w 
    #peak6
    bound.append(1.698) #xc
    bound.append(0) #A
    bound.append(0.0038) #w 
    #peak7
    bound.append(1.719) #xc
    bound.append(0.1) #A
    bound.append(0.0061) #w 
    #peak8
    bound.append(1.728) #xc
    bound.append(0.2) #A
    bound.append(0.0017) #w 
    #peak9 - laser
    bound.append(1.827) #xc
    bound.append(0.47) #A
    bound.append(0.0011) #w

    return bound

def upper_bounds():
    """
        Returns:
            Array, parameters used for upper bounds during scipy curve_fit 
            in the form of the list with arguments in order:
            [y_0, center_0, amplitude_0, width_0, center_1 ...]. 
            
            For best result of curve_fit those parameters should be close to 
            expected values of fitting parameters
    """
    bound = []
    #y0 - baseline
    bound.append(np.inf)
    #peak1
    bound.append(1.672) #xc
    bound.append(130) #A
    bound.append(0.016) #w 
    #peak2
    bound.append(1.681) #xc
    bound.append(2) #A
    bound.append(0.003) #w 
    #peak3
    bound.append(1.685) #xc
    bound.append(2.0) #A
    bound.append(0.0062) #w 
    #peak4
    bound.append(1.6915) #xc
    bound.append(25) #A
    bound.append(0.008) #w 
    #peak5
    bound.append(1.703) #xc
    bound.append(100) #A
    bound.append(0.012) #w 
    #peak6
    bound.append(1.722) #xc
    bound.append(22) #A
    bound.append(0.007) #w 
    #peak7 
    bound.append(1.738) #xc
    bound.append(100) #A
    bound.append(0.0065) #w 
    #peak8
    bound.append(1.772) #xc
    bound.append(6) #A
    bound.append(0.006) #w 
    #peak9 - laser
    bound.append(1.868) #xc
    bound.append(1.34) #A
    bound.append(0.0051) #w 

    return bound 

name = 'PLE_all' #format of input data: x, y0, y1, y2 [...]
filepath = 'C:\\Users\\Joanna\Desktop\\PLE\\'+name+'.dat'
#if __name__ == "__main__":
    #out = predict_fitting_parameters(filepath)

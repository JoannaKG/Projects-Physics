import pandas as pd
import numpy as np

df = pd.read_csv('01_ml_circ_polarizer_1power_6sx10_linescan.txt', sep='\t', skiprows=[0], header=None)

a = 0
size_of_acquisitions = 1015
number_of_rows = df.shape[0]
number_of_acquisitions = int(number_of_rows / size_of_acquisitions)
print(number_of_acquisitions)

def cutter(a, i):
    df2 = df[[2, 3]].loc[a:a + 1014]
    writer = ('01_ml_circ_polarizer_1power_6sx10_linescan_' + str(i + 1) + '_X=' + str(np.round(df.at[a, 0],2)) + '_Y=' + str(np.round(df.at[a, 1],2)) + '.txt')
    df2.to_csv(writer, index=False, header=[str(np.round(df.at[a, 0],2)),str(np.round(df.at[a, 1],2))])


for i in range(number_of_acquisitions):
    cutter(a, i)
    a += size_of_acquisitions

import pickle

import pandas as pd

from plotting_functions import *
from df_functions import *

### LOADING VARIABLES ############################3

""""
# Getting back the objects:
with open('objs.pkl') as f:  # Python 3: open(..., 'rb')
    obj0, obj1, obj2 = pickle.load(f)
"""

# with open("CD59FL/ER_time_array") as fp:
#    CD59FL_ER_time_array= pickle.load(fp)

# ER

with open('CD59FL/ER_MFI', 'rb') as f:
    CD59FL_ER_array = pickle.load(f)

with open('GPIFRA/ER_MFI', 'rb') as f:
    GPIFRA_ER_array = pickle.load(f)

with open('GPICD59/ER_MFI', 'rb') as f:
    GPICD59_ER_array = pickle.load(f)

with open('PCX/ER_MFI', 'rb') as f:
    PCX_ER_array = pickle.load(f)

with open('Ecad/ER_MFI', 'rb') as f:
    Ecad_ER_array = pickle.load(f)

# time

with open('CD59FL/ER_time_array', 'rb') as f:
    CD59FL_ER_time = pickle.load(f)

with open('GPIFRA/ER_time_array', 'rb') as f:
    GPIFRA_ER_time = pickle.load(f)

with open('GPICD59/ER_time_array', 'rb') as f:
    GPICD59_ER_time = pickle.load(f)

with open('PCX/ER_time_array', 'rb') as f:
    PCX_ER_time = pickle.load(f)

with open('Ecad/ER_time_array', 'rb') as f:
    Ecad_ER_time = pickle.load(f)

# Golgi
with open('CD59FL/Golgi_MFI', 'rb') as f:
    CD59FL_Golgi_array = pickle.load(f)

with open('GPIFRA/Golgi_MFI', 'rb') as f:
    GPIFRA_Golgi_array = pickle.load(f)

with open('GPICD59/Golgi_MFI', 'rb') as f:
    GPICD59_Golgi_array = pickle.load(f)

with open('PCX/Golgi_MFI', 'rb') as f:
    PCX_Golgi_array = pickle.load(f)

with open('Ecad/Golgi_MFI', 'rb') as f:
    Ecad_Golgi_array = pickle.load(f)

# Golgi
with open('CD59FL/Golgi_time_array', 'rb') as f:
    CD59FL_Golgi_time = pickle.load(f)

with open('GPIFRA/Golgi_time_array', 'rb') as f:
    GPIFRA_Golgi_time = pickle.load(f)

with open('GPICD59/Golgi_time_array', 'rb') as f:
    GPICD59_Golgi_time = pickle.load(f)

with open('PCX/Golgi_time_array', 'rb') as f:
    PCX_Golgi_time = pickle.load(f)

with open('Ecad/Golgi_time_array', 'rb') as f:
    Ecad_Golgi_time = pickle.load(f)

#####

df_CD59FL_golgi=pd.read_csv('CD59FL/CD59FL Golgi Exit.csv')
df_GPIFRA_golgi=pd.read_csv('GPIFRA/GPIFRA Golgi Exit.csv')
df_GPICD59_golgi=pd.read_csv('GPICD59/GPICD59 Golgi Exit.csv')
df_PCX_golgi=pd.read_csv('PCX/PCX Golgi Exit.csv')
df_Ecad_golgi=pd.read_csv('Ecad/Ecad Golgi Exit.csv')


df_CD59FL_ER=pd.read_csv('CD59FL/CD59FL RE Exit.csv')
df_GPIFRA_ER=pd.read_csv('GPIFRA/GPIFRA RE Exit.csv')
df_GPICD59_ER=pd.read_csv('GPICD59/GPICD59 RE Exit.csv')
df_PCX_ER=pd.read_csv('PCX/PCX RE Exit.csv')
df_Ecad_ER=pd.read_csv('Ecad/Ecad RE Exit.csv')

#### CALCULATING MEAN AND STD #############

CD59FL_ER_array_mean=mean(CD59FL_ER_array)
GPIFRA_ER_array_mean=mean(GPIFRA_ER_array)
GPICD59_ER_array_mean=mean(GPICD59_ER_array)
PCX_ER_array_mean=mean(PCX_ER_array)
Ecad_ER_array_mean=mean(Ecad_ER_array)

CD59FL_ER_array_std=deviation(CD59FL_ER_array)
GPIFRA_ER_array_std=deviation(GPIFRA_ER_array)
GPICD59_ER_array_std=deviation(GPICD59_ER_array)
PCX_ER_array_std=deviation(PCX_ER_array)
Ecad_ER_array_std=deviation(Ecad_ER_array)


CD59FL_Golgi_array_mean=mean(CD59FL_Golgi_array)
GPIFRA_Golgi_array_mean=mean(GPIFRA_Golgi_array)
GPICD59_Golgi_array_mean=mean(GPICD59_Golgi_array)
PCX_Golgi_array_mean=mean(PCX_Golgi_array)
Ecad_Golgi_array_mean=mean(Ecad_Golgi_array)

CD59FL_Golgi_array_std=deviation(CD59FL_Golgi_array)
GPIFRA_Golgi_array_std=deviation(GPIFRA_Golgi_array)
GPICD59_Golgi_array_std=deviation(GPICD59_Golgi_array)
PCX_Golgi_array_std=deviation(PCX_Golgi_array)
Ecad_Golgi_array_std=deviation(Ecad_Golgi_array)



###### CALCULATING MAXIMUN PEAK ######

print('CD59FL_ER_peaks')
CD59FL_ER_peaks = peak_estimador(CD59FL_ER_array,CD59FL_ER_time)
print('GPIFRA_ER_peaks')
GPIFRA_ER_peaks = peak_estimador(GPIFRA_ER_array,GPIFRA_ER_time)
print('GPICD59_ER_peaks')
GPICD59_ER_peaks = peak_estimador(GPICD59_ER_array,GPICD59_ER_time)
print('PCX_ER_peaks')
PCX_ER_peaks = peak_estimador(PCX_ER_array,PCX_ER_time)
print('Ecad_ER_peaks')
Ecad_ER_peaks = peak_estimador(Ecad_ER_array,Ecad_ER_time)

peaks_frames_ER=[CD59FL_ER_peaks,GPIFRA_ER_peaks,GPICD59_ER_peaks,PCX_ER_peaks,Ecad_ER_peaks]
df_peaks_frames_ER=pd.concat(peaks_frames_ER)
df_peaks_frames_ER.to_csv('Results/ER_peaks_time.csv')

print('CD59FL_Golgi_peaks')
CD59FL_Golgi_peaks = peak_estimador(CD59FL_Golgi_array,CD59FL_Golgi_time)
print('GPIFRA_Golgi_peaks')
GPIFRA_Golgi_peaks = peak_estimador(GPIFRA_Golgi_array,GPIFRA_Golgi_time)
print('GPICD59_Golgi_peaks')
GPICD59_Golgi_peaks = peak_estimador(GPICD59_Golgi_array,GPICD59_Golgi_time)
print('PCX_Golgi_peaks')
PCX_Golgi_peaks = peak_estimador(PCX_Golgi_array,PCX_Golgi_time)
print('Ecad_Golgi_peaks')
Ecad_Golgi_peaks = peak_estimador(Ecad_Golgi_array,Ecad_Golgi_time)

peaks_frames_golgi=[CD59FL_Golgi_peaks,GPIFRA_Golgi_peaks,GPICD59_Golgi_peaks,
                    PCX_Golgi_peaks,Ecad_Golgi_peaks]
df_peaks_frames_golgi=pd.concat(peaks_frames_golgi)
df_peaks_frames_golgi.to_csv('Results/Golgi_peaks_time.csv')

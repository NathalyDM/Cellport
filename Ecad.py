import pickle
from golgi_functions import *
from plotting_functions import *
from df_functions import *
import warnings
warnings.filterwarnings("ignore")
warnings.simplefilter(action='ignore', category=FutureWarning)
sheet_name_golgi='ECad Golgi Exit'
sheet_name_ER='ECad ER Exit'

plt.close('all')
###############  IMPORTING DATAFRAME   ###############
df_golgi = pd.read_excel(r'Z:/Paris/Quantification/Prim All Cargoes from Time point 8 included.xlsx', sheet_name=sheet_name_golgi)
df_golgi.dropna(axis=0,inplace=True)

df_ER = pd.read_excel(r'Z:/Paris/Quantification/Prim All Cargoes from Time point 8 included.xlsx', sheet_name=sheet_name_ER)
df_ER.dropna(axis=0,inplace=True)


###############  SPLITTING MFI AND TIME   ###############
step=0.4375
golgi_time=df_golgi.iloc[:, 0]*step
golgi_MFI=df_golgi.iloc[:, 1:len(df_golgi.columns)-1]

ER_time=df_ER.iloc[:, 0]*step
ER_MFI=df_ER.iloc[:, 1:len(df_golgi.columns)-1]



###############  FILTERING DATA   ###############
ER_time_array,ER_array,ER_xlim_time=filtering_ER(ER_time,ER_MFI,0,int((21-8*step)/step),0.58)
Golgi_time_array,Golgi_array,Golgi_xlim_time=filtering_Golgi(golgi_time,golgi_MFI,int((10-8*step)/step),int((55-8*step)/step),0.5)

print('ER before filtering',len(ER_MFI.columns))
print('ER after filtering',len(ER_array))

print('GoLgi before filtering',len(golgi_MFI.columns))
print('ER after filtering',len(Golgi_array))


###############  PLOTING VALUES   ###############
plot_ind_protein(ER_time,ER_MFI, sheet_name_golgi)
plot_ind_protein(golgi_time,golgi_MFI, sheet_name_golgi)

###############  PLOTING FILTERED VALUES   ###############
plot_filtered_protein_ER(ER_time_array,ER_array,ER_xlim_time,step,sheet_name_ER,ER_xlim_time)
plot_filtered_protein_ER(Golgi_time_array,Golgi_array,Golgi_xlim_time,step,sheet_name_golgi,Golgi_xlim_time)


###############  ANALYZING DATAFRAME   ###############
df_max = peak_estimador(ER_array,ER_time_array)
df_max = peak_estimador(Golgi_array,Golgi_time_array)

###############  FITTING MODELLING GOLGI  ###############
df_solutions_RE= pd.DataFrame(columns=['System','Organelle','Cell','D(min-1)','V(min-1)','r(min-1)','kL(min-1)',
                                         'kR(min-1)', 'dt(min)', 'Dt(min)','Error'])
cell=1
for protein_mfi in ER_array:
    #print(column)
    desired_output = protein_mfi
    best_specimen, best_specimen_error = start(1000, desired_output, 10, 0.3)
    # [D,V,r,kL, kR, dt, Dt]
    df2 = pd.DataFrame([{'System':'Ecad','Organelle':'Golgi','Cell':cell,'D(min-1)':best_specimen[0],'V(min-1)':best_specimen[1],
                             'r(min-1)':best_specimen[2],'kL(min-1)':best_specimen[3],
                             'kR(min-1)':best_specimen[4],'dt(min)':best_specimen[5],
                             'Dt(min)':best_specimen[6], 'Error':best_specimen_error}])
    df_solutions_RE = pd.concat([df_solutions_RE, df2], axis=0, ignore_index=True)
    cell=cell+1

###############  FITTING MODELLING GOLGI  ###############
df_solutions_Golgi = pd.DataFrame(columns=['System','Organelle','Cell','D(min-1)','V(min-1)','r(min-1)','kL(min-1)',
                                         'kR(min-1)', 'dt(min)', 'Dt(min)','Error'])

cell=1
for protein_mfi in Golgi_array:
    #print(column)
    desired_output = protein_mfi
    best_specimen, best_specimen_error = start(1000, desired_output, 10, 0.3)
    # [D,V,r,kL, kR, dt, Dt]
    df2 = pd.DataFrame([{'System':'Ecad','Organelle':'Golgi','Cell':cell,'D(min-1)':best_specimen[0],'V(min-1)':best_specimen[1],
                             'r(min-1)':best_specimen[2],'kL(min-1)':best_specimen[3],
                             'kR(min-1)':best_specimen[4],'dt(min)':best_specimen[5],
                             'Dt(min)':best_specimen[6], 'Error':best_specimen_error}])
    df_solutions_Golgi = pd.concat([df_solutions_Golgi, df2], axis=0, ignore_index=True)
    cell=cell+1

###############  EXPORTING RESULTS   ###############

with open("Ecad/ER_time_array", "wb") as fp:
    pickle.dump(ER_time_array, fp)

with open("Ecad/ER_MFI", "wb") as fp:
    pickle.dump(ER_array, fp)

with open("Ecad/Golgi_time_array", "wb") as fp:
    pickle.dump(Golgi_time_array, fp)

with open("Ecad/Golgi_MFI", "wb") as fp:
    pickle.dump(Golgi_array, fp)

df_solutions_Golgi.to_csv('Ecad/Ecad Golgi Exit.csv')
df_solutions_RE.to_csv('Ecad/Ecad RE Exit.csv')


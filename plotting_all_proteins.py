from plotting_functions import *
from variables import *


##### PLOTTING
means=[CD59FL_ER_array_mean,GPIFRA_ER_array_mean,GPICD59_ER_array_mean,PCX_ER_array_mean,Ecad_ER_array_mean]
deviations=[CD59FL_ER_array_std,GPIFRA_ER_array_std,GPICD59_ER_array_std,PCX_ER_array_std,Ecad_ER_array_std]

means_g=[CD59FL_Golgi_array_mean,GPIFRA_Golgi_array_mean,GPICD59_Golgi_array_mean,PCX_Golgi_array_mean,Ecad_Golgi_array_mean]
deviations_g=[CD59FL_Golgi_array_std,GPIFRA_Golgi_array_std,GPICD59_Golgi_array_std,PCX_Golgi_array_std,Ecad_Golgi_array_std]


# MAximun length
max_length=FindMaxLength(means)
max_length_g=FindMaxLength(means_g)

step=0.4375
plot_all_proteins_final('ER exit', means, max_length,deviations,step)
plot_all_proteins_final('Golgi exit', means_g, max_length_g,deviations_g,step)





import pandas as pd
import numpy as np

def peak_estimador(MFI_arrays,time_arrays):
    df_max = pd.DataFrame(columns=['colname', 'max_value','min_value','mean_protein_level',
                                   'time_max_value', 'time_min_value', 'mean_time_max_min'])
    cell=1

    for protein_mfi in MFI_arrays:
        if len(protein_mfi)>0:
            #print(protein_mfi)
            #protein_mfi=np.array(protein_mfi)
            max_index = np.argmax(protein_mfi)
            min_index = np.argmin(protein_mfi)
            peak_value=protein_mfi[max_index]
            low_value=protein_mfi[min_index]
            time_max_value = time_arrays[cell-1][max_index]
            time_min_value = time_arrays[cell - 1][min_index]
            mean_time_max_min=time_arrays[cell-1][int((max_index+min_index)/2)]
            mean_protein_level=protein_mfi[int((max_index+min_index)/2)]

            df1 = pd.DataFrame([{'colname': cell, 'max_value': peak_value, 'min_value': low_value, 'mean_protein_level': mean_protein_level,
                                 'time_max_value': time_max_value, 'time_min_value': time_min_value, 'mean_time_max_min': mean_time_max_min}])
            df_max = pd.concat([df_max, df1], axis=0, ignore_index=True)

        cell=cell+1

    df_max['std_max_value'] = df_max['max_value'].std()
    df_max['std_mean_protein_level'] = df_max['mean_protein_level'].std()
    df_max['std_mean_time_max_min'] = df_max['mean_time_max_min'].std()

    df_max['std_min_value']=df_max['min_value'].std()
    df_max['std_time_max_value'] = df_max['time_max_value'].std()
    df_max['std_time_min_value'] = df_max['time_min_value'].std()

    print('The max mean value of the protein is',round( df_max['max_value'].mean(),2), ' +/-', round(df_max['max_value'].std(),2),
          ' at time:', round(df_max['time_max_value'].mean(),2), ' +/-', round(df_max['time_max_value'].std(),2), 's.')
    #print('\n')
    print('The middle mean value of the protein is', round(df_max['mean_protein_level'].mean(),2), ' +/-',
          round(df_max['mean_protein_level'].std(),2),
          ' at time:', round(df_max['mean_time_max_min'].mean(),2), ' +/-', round(df_max['mean_time_max_min'].std(),2), 's.')

    print('The min mean value of the protein is', round(df_max['min_value'].mean(),2), ' +/-', round(df_max['min_value'].std(),2),
          ' at time:', round(df_max['time_min_value'].mean(),2), ' +/-', round(df_max['time_min_value'].std(),2), 's.')
   # print('\n')

    print('\n')
    return df_max

def filtering_ER(ER_time,ER_MFI,min,max,trh):
    empty_ER_array=[]

    empty_ERtime_array = []
    xlim_list = []
    cnt=0
    for column in ER_MFI:
        sample=ER_MFI[column].to_numpy()
        ER_time_array = ER_time.to_numpy()
        #print(ER_time_array[-1])
        minpos = np.argmin(sample[min:max])
        maxpos = np.argmax(sample[min:max])


        if sample[minpos]<trh:
            #print(sample[minpos])
            #print(column)
            #print(minpos)
            empty_ER_array.append(sample[maxpos:minpos + 1])
            empty_ERtime_array.append(ER_time_array[maxpos:minpos + 1])
            xlim_list.append(ER_time_array[minpos + 1])
            #print(ER_time_array[minpos ])

        #else:
            #print(sample[minpos])
            #print(column)
            #print(minpos)
    xlim_max = np.amax(xlim_list)
    return empty_ERtime_array,empty_ER_array,xlim_max


def filtering_Golgi(ER_time,ER_MFI,min,max,trh):
    empty_ER_array=[]

    empty_ERtime_array = []
    xlim_list=[]
    cnt=0

    for column in ER_MFI:
        sample=ER_MFI[column].to_numpy()
        ER_time_array = ER_time.to_numpy()
        #print(ER_time_array[-1])
        minpos = np.argmin(sample[min:max])
        maxpos = np.argmax(sample[min:max])

        print(sample[maxpos:minpos])
        if len(sample[maxpos:minpos])>0:
            if sample[minpos] < trh and np.amax(sample[:minpos][-20:]) < np.amax(sample[min:max]):
                # print(sample[minpos])
                # print(column)
                # print(minpos)
                empty_ER_array.append(sample[0:minpos + 1])
                empty_ERtime_array.append(ER_time_array[0:minpos + 1])
                xlim_list.append(ER_time_array[minpos + 1])
                # print(ER_time_array[max])

        #else:
            #print(sample[minpos])
            #print(column)
            #print(minpos)

    xlim_max=np.amax(xlim_list)
    return empty_ERtime_array,empty_ER_array,xlim_max


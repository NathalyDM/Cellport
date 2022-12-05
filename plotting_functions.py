import matplotlib.pyplot as plt
import numpy as np

def FindMaxLength(lst):
    maxLength = max(len(x) for x in lst)

    return maxLength




def plot_all_proteins_final(title, means, max_length,deviations,step):
    fig = plt.figure()
    ax = fig.gca()
    ax.set_xticks(np.arange(0,48, 3))
    ax.set_yticks(np.arange(0, 1.1, 0.1))
    plt.grid()

    for idx in range(len(means)):
        mean_value = np.array(means[idx])
        std_value = np.array(deviations[idx])

        if len(mean_value) < max_length:
            mean_value = np.append(mean_value, np.repeat(mean_value[-1], max_length - len(mean_value)))
            std_value = np.append(std_value, np.repeat(std_value[-1], max_length - len(std_value)))

        X_minus_sigma = mean_value - std_value
        X_plus_sigma = mean_value + std_value

        #colors = ['blue', 'black', 'red', 'green']
        labels = ['CD59FL', 'GPIFRA', 'GPICD59', 'PCX', 'Ecad']

        time = np.arange(0, len(mean_value) * step, step)
        # print(len(mean_value))
        # print(len(time)) , color=colors[idx] , color=colors[idx]

        plt.plot(time, mean_value, alpha=0.8, label=labels[idx])
        plt.fill_between(time, X_minus_sigma, X_plus_sigma, alpha=0.1)

    # Add title and axis names
    plt.title(title)
    plt.xlabel('Time (min)')
    plt.ylabel('Mean Fluorescence Intensity (MFI)')
    plt.legend(loc="upper right")
    plt.show()

def plot_ind_protein(GPIFRA_golgi_time,GPIFRA_golgi_MFI, title):
    mean_GPIFRA_golgi_MFI = GPIFRA_golgi_MFI.mean(axis=1)
    std_GPIFRA_golgi_MFI = GPIFRA_golgi_MFI.std(axis=1)
    X_plus_sigma = mean_GPIFRA_golgi_MFI + std_GPIFRA_golgi_MFI
    X_minus_sigma = mean_GPIFRA_golgi_MFI - std_GPIFRA_golgi_MFI

    ###############  PLOTING VALUES   ###############
    fig = plt.figure()
    ax = fig.gca()
    ax.set_xticks(np.arange(0, GPIFRA_golgi_time.iloc[-1], 5))
    ax.set_yticks(np.arange(0, 1.1, 0.1))
    plt.grid()

    for column in GPIFRA_golgi_MFI.columns:
        plt.scatter(GPIFRA_golgi_time, GPIFRA_golgi_MFI[column], alpha=0.05, color='blue')

    plt.plot(GPIFRA_golgi_time, mean_GPIFRA_golgi_MFI, color='black')
    # plt.fill_between(GPIFRA_golgi_time, X_minus_sigma, X_plus_sigma ,alpha=0.3, color = 'blue')

    # Add title and axis names
    plt.title(title )
    plt.xlabel('Time (min)')
    plt.ylabel('Mean Fluorescence Intensity (MFI)')
    plt.ylim(0, 1)
    # Showing plot adding grid
    plt.show()

def Average(lst):
    return sum(lst) / len(lst)

def mean(arrays):
    lengths=[]
    for array in arrays:
        lengths.append(len(array))
    max_length=max(lengths)
    #print(max_length)

    mean_array=[]
    for idx in range(max_length):
        column_mean = []
        for array in arrays:
            #print(array)
            if idx>=len(array):
                continue
            column_mean.append(array[idx])

        mean_array.append(Average(column_mean))
    #print(mean_array)
    return mean_array

def deviation(arrays):
    lengths = []
    for array in arrays:
        lengths.append(len(array))
    max_length = max(lengths)
    #print(max_length)

    std_array = []
    for idx in range(max_length):
        column_std = []
        for array in arrays:
            # print(array)
            if idx >= len(array):
                continue
            column_std.append(array[idx])

        std_array.append(np.std(column_std))
   # print(std_array)
    return std_array


################## Ploting ##################

def plot_filtered_protein_ER(ER_time_array,ER_array,xlim_time,step,title, xmaxtick):
    flat_time = [item for sublist in ER_time_array for item in sublist]
    mean_line=mean(ER_array)

    print(min(flat_time))
    print(max(flat_time))
    mean_time=np.arange(min(flat_time), max(flat_time)+step,step)
    print(mean_line)
    print(mean_time)

    fig = plt.figure()
    ax = fig.gca()
    ax.set_xticks(np.arange(0, xlim_time, 1))
    ax.set_yticks(np.arange(0, 1.1, 0.1))
    plt.grid()

    for array in range(len(ER_array)):
        plt.scatter(ER_time_array[array], ER_array[array], alpha=0.1, color='blue')

    plt.plot(mean_time, mean_line, color='black')

    # Add title and axis names
    plt.title(title)
    plt.xlabel('Time (min)')
    plt.ylabel('Mean Fluorescence Intensity (MFI)')
    plt.ylim(0, 1.05)
    plt.xlim(8*step, xmaxtick)
    # Showing plot adding grid
    plt.show()






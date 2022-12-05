#!/usr/bin/env python
# coding: utf-8
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import random

# In[10]:


def generate_chromosome(): 
    
    D_range=np.arange(0.01,0.5,0.1)
    V_range=np.arange(0.01,0.5,0.01)
    r_range=np.arange(0,0.5,0.1)
    kL_range=np.arange(0.0,0.6,0.01)
    kR_range=np.arange(0.0,0.6,0.01)
    dt=[0.4375]
    Dt=0.4375*np.arange(1,121,2)


    # [D,V,r,kL, kR, dt, Dt]
    chromosome=[]
    chromosome.append(random.choice(D_range))
    chromosome.append(random.choice(V_range))
    chromosome.append(random.choice(r_range))
    chromosome.append(random.choice(kL_range))
    chromosome.append(random.choice(kR_range))
    chromosome.append(random.choice(dt))
    chromosome.append(random.choice(Dt))
    
    return chromosome


# In[1]:


def golgi_intra_transport(solution, A, flux):

    #Fitting
    D = solution[0]          # [m^2/sec] diffusion constant 
    V = solution[1]           # [m/sec] velocity 
    r=solution[2]            # rate of exit rn from cisterna n.
    kL=solution[3]             # rates of intercisternal exchange backward <-
    kR=solution[4]             # rates of intercisternal exchange forward ->
    dt=solution[5]         # Step time of difussion
    Dt=solution[6]     # Step time of convection

    Tf=len(A)*dt  # Final time
    t=0
    
    Lx=len(A)
    new_flux= np.zeros_like(A)
    new_flux[0:len(flux)]= flux
    
    # Parameters are linked by the relations 
    dx=Dt*V            # spatial spacing
    n=D*Dt/dx**2
    N=Dt/dt
    cnt=1
    
    while t<Tf :
        t=t+Dt
        Jon= new_flux #0.3  
        A[1:Lx-1]=A[0:Lx-2]
        A[0]=Jon[cnt-1] #%Jon;
        k=0
        cnt=cnt+1
        for k  in range(1,int(N)):
            #print(k)
            AL=A[1]
            AR=A[Lx-2]
            A[1:Lx-2]= A[1:Lx-2]*(1-2/n-r*dt)+A[0:Lx-3]/n+A[2:Lx-1]/n
            A[0]=A[1]*(1-1 / n-kL/n)+AL/n+kL*Jon[cnt-1]*dt
            A[Lx-1]=A[Lx-1]*(1 -1/n-kR/n)+AR/n
    
    return A 


# In[3]:


def f_error(solution, desired_output, flux):
    
    #A = pd.read_csv('data_GOLGI.csv')
    #A =np.asarray(A['Normalized to the Background_C1'][1:68])
    A=desired_output.copy()
    output = golgi_intra_transport(solution, A, flux)
    f_error = sum(abs(output - desired_output))
    return f_error


# In[4]:


def generate_population(population_size): 
    population_size = int(population_size)
    population = []
    
    for i in range(population_size):
        while True:
            chromosome = generate_chromosome()
            if chromosome not in population:
                break
        population.append(chromosome)
    
    return population


# In[5]:


def start_Golgi(population_size, desired_output, generations_size, top_percentage, flux, plot_TF=True):
    population = generate_population(population_size)
    pop_error = []

    for chromosome in population:
        pop_error.append(f_error(chromosome, desired_output, flux))

    generations_size = int(generations_size)
    for gen in range(generations_size - 1):
        new_population = []
        best_parents = fitness(population, pop_error, top_percentage)
        pop_error.clear()
        for chromosome in population:
            [parent_1, parent_2] = _parents_selection(best_parents)
            child = _crossover(parent_1, parent_2, False, 0.5)
            pop_error.append(f_error(child, desired_output, flux))
            new_population.append(child)
        population = new_population

    (best_specimen, best_specimen_error)=_show_best_specimen(population, pop_error, desired_output, plot_TF, flux)
    return (best_specimen, best_specimen_error)

# In[6]:


def _show_best_specimen(population, pop_error, desired_output,plot_TF, flux):
    best_specimen = fitness(population, pop_error)[0]
    best_specimen_error = f_error(best_specimen, desired_output, flux)

    print('The best solution is: ', best_specimen)
    print('Its error is: ', best_specimen_error)

    A = desired_output.copy()
    output=golgi_intra_transport(best_specimen, A, flux)

    if plot_TF == True: 
        # Plotting
        plt.plot(output, 'r')
        plt.plot(desired_output, 'b')
        plt.show()

    return (best_specimen, best_specimen_error)

# In[7]:


def fitness(population,pop_error,top_percentage=1):
    dataset = pd.DataFrame({'chr': population, 'cost': pop_error}, columns=['chr', 'cost'])
    dataset=dataset.sort_values(by='cost', ascending=True)
    dataset=dataset.head(int(len(dataset)*top_percentage))
    #print(dataset['chr'].to_numpy())
    return(dataset['chr'].to_numpy())


# In[8]:


def _parents_selection(best_parents): 
    parent_1 = random.choice(best_parents)
    parent_2 = random.choice(best_parents)
    while parent_2 == parent_1:
        parent_2 = random.choice(best_parents)
    return([parent_1,parent_2])


# In[9]:


import math
def _crossover(parent1, parent2,crossover_split_random,crossover_split_size):
    if crossover_split_random:
        split_size = random.randint(0, len(parent1))

    else:
        fraction = float(crossover_split_size)
        split_size = math.floor(fraction * len(parent1))
    
    heritage_p1=random.sample(range(0, len(parent1)-1), int(split_size))
    child=[]
    for gene in range(len(parent1)):
        if gene in heritage_p1: 
            child.append(parent1[gene])
        else: 
            child.append(parent2[gene])
    return child


# In[ ]:





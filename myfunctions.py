#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Functions for MobiVirus Simulator
"""
#Packages
import numpy as np
import random
import matplotlib.pyplot as plt
import pandas as pd
import statistics
import imageio
import os
import sys
from configparser import ConfigParser

"""
=================================
READ THE PARAMETERS FROM .INI FILE
=================================
"""

# Specify the directory and file name that contains the parameters
directory = './'
file_name = 'parameters.ini'
file_path = os.path.join(directory, file_name)

# Check if the file exists
if not os.path.exists(file_path):
    print(f"Error: {file_name} must be in the current directory.")
    sys.exit("Exiting the program.")

# File exists, proceed with parsing
config = ConfigParser()
config.read(file_path)


"""
-----------------------------------------------
Read directory from the Initial_Parameters section
-----------------------------------------------
"""
directory = config.get('Directories', 'directory').strip('"')                 # Directory where the scripts exist


"""
===================================
BUILDING THE DATASET OF INDIVIDUALS 
===================================
"""

''' 
--------------------
Coordinates Function 
--------------------
'''

def coords_func(n, bound_l, bound_h):
    
    ## Creates a stacked column array of the (x,y) 2D coordinates of each individual in the dataset ##
    ## The coordinates are assingned randomly using a uniform distribution ##
    
    # n = total number of individuals in the sample
    # bound_l =  lower bound of the spacial axes
    # bound_h = higher bound of the spacial axes
    # x = random coordinate in the x axis 
    # y = random coordinate in the y axis
    
    # Default values:  The defined space is a rectangle with the same bounds in the x and y axis 
    
    
    x = np.random.uniform(low=bound_l, high=bound_h, size=n)
    y = np.random.uniform(low=bound_l, high=bound_h, size=n)
    
    return np.column_stack((x, y))

'''
------------------------
Infection Label Function
------------------------
'''

def label(n, i):
    
    ## Creates a stacked column array where every cell corresponds to an individual, with a label depending on if they are infected with the virus or not ##
    ## In the initial stage of the simulation, there is a sample of infected individuals among the healthy ones ##
    
    # n = total number of individuals in the sample
    # i = number of infectd individuals
    # label 0 = healthy individual
    # label 1 = infected individual
    
    
    z = np.array([1]*i)
    o = np.array([0]*(n-i))
    s = np.concatenate((z, o))
    random.shuffle(s)
    l = np.column_stack(s)
    return l.T

'''
---------------
Genome Function
---------------
'''

def genome(n, l, r_m, n_i):
    
    ## Creates an empty data frame where the rows correspond to the genomes of different individuals and the columns to the genome positions ##
    ## Returns the dataframe together with the total mutation rate of a genome
    
    # n = total number of individuals in the sample
    # l = total length of the genome
    # r_m = mutation rate of a position in the genome
    # n_i = number of important genome positions that form a specific strain (they will be in the begining of the genome)
    
    # Default result: If the genome of an individual has a specific value in at least one of the important positions, the individual is infected by the Super Strain
    
    g =  pd.DataFrame(index=range(n), columns=range(l))
    r_tot = r_m * l #Total rate of mutation of genome
    
    #probi = sum(g[:n_i]) #Rate of infection of individual
    
    return g, r_tot

'''
----------------------
Movement Rate Function
----------------------
'''

def probs(coords, rm_i, rm_h):
    
    ## Creates a stacked column array where every cell corresponds to an individual with a value for the rate of movement for that individual ##
    ## The rate of movement depends on the infection label (see parameters.ini) ##
    

    # coords = the full coordinate table
    # rm_i = rate of movement for infected individual
    # rm_h = rate of movement for healthy individual
    
    probm = np.zeros(len(coords)) #random rate of movement for each individual

    mask = coords[:, 2] == 1.0 #for those who are infected
    probm[mask] = rm_i         #they get rm_i in the corresponding positions in probm
    #probi[mask] = 0.5

    mask = coords[:, 2] != 1.0 #for those that are healthy
    probm[mask] = rm_h          #they get rm_h in the corresponding positions in probm
    
    return np.column_stack(probm).T #, np.column_stack(probi).T


'''
-----------------------
Infection Rate Function
-----------------------
'''

def infectivity(in_probi, g, n_i, ri_s):
    
    ## Gives the specific value of the infection rate of an individual, depending on the strain of the virus ##
    ## The strain is decided by the n_i important positions in the individual's genome ##
    
    # g = genome of individual
    # n_i = important positions in the genome
    # in_probi = infection rate of the infector
    # ri_s = infectivity rate of super strain
    
    
    if np.any(np.array(g[:n_i])==1):
        probi = ri_s
    else:
        probi = in_probi
        
    #probi = np.where(g[:n_i].any() == 1 , 1.5, probi)
    
    return probi



'''
--------------------------------
Initial Distance Matrix Function
--------------------------------
'''

def ini_distm(coords):
    
    ## Creates a dataframe containing the distance between all the individuals in the sample ##
    ## The rows and the columns of the dataframe both correspond to the individuals ##
    
    # coords = dataframne of the individual's spacial coordinates
    
    n = len(coords)
    df = np.zeros((n, n))
    for i in range(n):
        for j in range(i+1, n):
            x1, y1 = coords[i, :2]
            #print(x1, y1)
            x2, y2 = coords[j, :2]
            #print(x2, y2)
            df[j,i] = np.sqrt((x2-x1)**2 + (y2-y1)**2)
            #print(df[j, i])
    # return pd.DataFrame(df, index=range(n), columns=range(n))
    return df




"""
================================================
SIMULATING THE SPACIAL MOVEMENT OF AN INDIVIDUAL 
================================================
"""

'''
-----------------
Movement Function
-----------------
'''

def movement(coords, bound_l, bound_h, c, std_x, std_y, rim): 

    ## Calculates the new spacial coordinates of the individual that will move ##
    ## Returns the full data table but with new coordinates ##
    ## The new coordinates are randomly selected from two different Gaussian distributions (one for each coordinate) ##
    ## The mean value of the distribution is the initial coordinate (x or y accordingly) and the standard deviation is taken from the standard deviation of the distribution of all of the coordinates from our sample (x or y accordingly) ##
    ## If the individual is infected, the standard deviation of the Gaussian distribution is 3 times smaller than the one for a healthy individual ##
    ## In that way we add the condition that an infected individual is more likely to stay in the same place or move closer to its initial position ##
    
    
    # coords = the full coordinate table in order to get the standard deviations
    # bound_l =  lower bound of the spacial axes
    # bound_h = higher bound of the spacial axes
    # c = the specific individual that will move
    # std_x = standard deviation of x coordinate
    # std_y = standard deviation of y coordinate
    # If these standard deviations are defined, we use them and not the standard deviation from all of the coordinates in the sample
    # rim = relative infected mobility. The parameter that differentiates between the mobility of a healthy individual and an infected one
    ## Individuals are not allowed to move outside the box (or space) ##
    
    #n = len(coords[c, :])  
    x_f = 0
    y_f = 0
    label_f = coords[c, 2]
    probm_f = coords[c, 3]
    probi_f = coords[c, 4]
    mut_f = coords[c, 5]
    sus_f = coords[c, 6]

    r_x = statistics.stdev(coords[:, 0])
    r_y = statistics.stdev(coords[:, 1])

    if coords[c, 2] == 0:

        a = np.random.normal(coords[c,0], r_x )
        while a > bound_h or a< bound_l:
            a = np.random.normal(coords[c,0], r_x )
                
        b = np.random.normal(coords[c,1], r_y )
        while b > bound_h or b< bound_l:
            b = np.random.normal(coords[c,1], r_y )
                 
        x_f = a
        y_f = b
    else:
            
        a = np.random.normal(coords[c,0], r_x*rim )
        while a > bound_h or a< bound_l:
            a = np.random.normal(coords[c,0], r_x*rim )
                
        b = np.random.normal(coords[c,1], r_y*rim )
        while b > bound_h or b< bound_l:
            b = np.random.normal(coords[c,1], r_y*rim )
                 
        x_f = a
        y_f = b

    return np.column_stack((x_f, y_f, label_f, probm_f, probi_f, mut_f, sus_f))


'''
---------------------------------------
Distance Matrix After Movement Function
---------------------------------------
'''

## Changes only the column of the initial distance matrix that coresponds to the individual that moved ##
## For that column onlly, it calculates the new distances between the individual c and the rest in the sample ##


def new_dist(coords, coords_f, initial_dist, c): #where c in the index of the coords that corresponds the individual that moved
    n = len(coords)
    df = initial_dist.copy()
    # dff = ini_distm(coords_f) 
    for i in range(n):
        for j in range(i+1, n):
            x1, y1 = coords_f[0,:2]
            x2, y2 = coords[j, :2]
            df[j,i] = np.where((i==c or j==c), np.sqrt((x2-x1)**2 + (y2-y1)**2), df[j,i])  #Change only the distance matrix cells of the individul that moved
    return df


def ind_probi(df, c):
    
    #Gives an array of length as many as the individuals and shape (n,1), only for the one individual that infects in that moment
    #In the cells, 1 is for those who are in the "infection distance" and 0 for the others
    #1 is also the probability of infection
    
    n = len(df) #the number of individuals
    ipi = np.zeros((n))
    for i in range(n):
        ipi[i] = np.where((df[i,c]<0.4)&(df[i,c]!=0), 1, 0) #The first condition defines the infection distance
    return ipi


    
def mutation(g_i, r_tot): #Mutation of an individuals'  genome
    
    #g_i = genome of an individual 
    #r_tot = Total rate of mutation of genome
    
    p = np.random.poisson(r_tot) #Random Poisson number with lambda = r_tot to give me how many mutations will happen
    N = np.random.randint(0, len(g_i), size = p) #Positions in the genome that mutations will happen (as many as the poisson)
    for i in range(len(N)):
        g_i[N[i]] = np.where(g_i[N[i]]==0, 1, 0)
    
    return g_i


def plot_map(plots_directory, coords, tt, ts, mv, ss, un):
    
    #Create plot with all of the regions

    #xc1, yc1, r1 = 100, 800, 200
    #xc2, yc2, r2 = 200, 200, 200
    #xc3, yc3, r3 = 800, 400, 200
    #circle1 = plt.Circle((xc1, yc1), r1, color='k', fill=False)
    #circle2 = plt.Circle((xc2, yc2), r2, color='k', fill=False)
    #circle3 = plt.Circle((xc3, yc3), r3, color='k', clip_on=False, fill=False)
    fig, ax = plt.subplots() # note we must use plt.subplots, not plt.subplot
    # (or if you have an existing figure)
    # fig = plt.gcf()
    # ax = fig.gca()
    #ax.add_patch(circle1)
    #ax.add_patch(circle2)
    #ax.add_patch(circle3)
    #ax.annotate('P1', (200,200), c='r')
    #ax.annotate('P2', (800,400), c='r')
    #ax.annotate('P3', (100,800), c='r')
    #ax.annotate('P4', (800,800), c='r')
    ax.scatter(200,200, c='k', s=0)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_xlim(-2, 8)
    ax.set_ylim(-0.5, 10)
    ax.set_title("Simulation time =""{0:.4f}".format(ts))
    #im=[ax.scatter(list(coords.iloc[0]), list(coords.iloc[1]), s=3)]
    #ax.text(-0.9, 5.9, "Infected people ="+str(len(coords[coords[:, 2] == 1][:,0])-ss))
    #ax.text(-0.9, 5.6, "Cumulative Movements ="+str(mv))
    #ax.text(-0.9, 5.3, "Super Spreaders ="+str(ss))
    ax.scatter(list(coords[(coords[:, 2] == 1)&(coords[:, 5] == 1)][:, 0]), list(coords[(coords[:, 2] == 1)&(coords[:, 5] == 1)][:, 1]), s=3, c='indianred', label='Normal infections = %i' %(len(coords[coords[:, 5] == 1][:,0]))) #normal strain
    ax.scatter(list(coords[(coords[:, 2] == 1)&(coords[:, 5] == 2)][:, 0]), list(coords[(coords[:, 2] == 1)&(coords[:, 5] == 2)][:, 1]), s=3, c='b', label='Super spreaders = %i' %(len(coords[coords[:, 5] == 2][:,0]))) #super strain
    #ax.scatter(list(coords[(coords_t[:, 3] == 1)&(coords_t[:, 6] == 2)][:, 0]), list(coords[(coords_t[:, 3] == 1)&(coords_t[:, 6] == 2)][:, 1]), s=3, c='lightcoral')
    ax.scatter(list(coords[coords[:, 2] == 0][:, 0]), list(coords[coords[:, 2] == 0][:, 1]), s=3, c='gold')
    ax.scatter([], [], label="Total Un-infections = %i" %un, s=0)
    #box = ax.get_position()
    #plt.subplots_adjust(right=0.6)
    #ax.set_position([box.x0, box.y0, box.width * 0.6, box.height])
    plt.legend(title = "Cumulative Movements = %i" %mv, loc='upper center')#, bbox_to_anchor=(1.04, 0.5),  borderaxespad=0)
    plt.savefig(plots_directory+'/'+str(tt) +'.jpg', dpi=1200)
    # ax.scatter(list(coords[i].iloc[0]), list(coords[i].iloc[1]), s=3, c='r')
    

def create_gif(directory):
    images = []
    names = []
    for filename in os.listdir(directory):
        if int(filename[:-4])==0 or int(filename[:-4])%100==0:
            names.append(os.path.join(directory, filename))
        
    names.sort()
    names.sort(key=len)
    
    for g in names:
        images.append(imageio.imread(g, format='jpg'))
    
    imageio.mimsave(directory+'figs.gif', images, duration=0.8)
    
def sample_data(samples_directory, genomes_directory, g, tt, coords_t, all_inf, sample_times):

    if tt==0 or tt%sample_times==0 or sum(coords_t[:,5]==1)==0:
        
        all_inf = pd.DataFrame(data=all_inf, columns=['Total infected', 'Super spreaders', 'Normal spreaders', 'Time'])
        all_inf.to_csv(samples_directory+'/all_inf_'+str(tt)+'.csv', header=True, index=False)
        
        g.to_csv(genomes_directory+'/genomes_'+str(tt)+'.csv',header=False, index=False)
        
        coords_t = pd.DataFrame(data=coords_t, columns=["x","y", "label", "rate of movement", "rate of infection", "mutation", "susceptibility"])
        coords_t.to_csv(samples_directory+'/coords_'+str(tt)+'.csv', header=True, index=False)
    
def do_plots(plots_directory, coords_t, tt, t_s, mv, ss, un):
    
    #Make plots according to some conditions
    #Here the condition is to make a plot if its the 0th generation of the simulation,
    #if the number of generation can be devided by 1000 or if the "Super Strain" disappears.
    
    if tt==0 or tt%1000==0 or sum(coords_t[:,5]==1)==0: 
        plot_map(plots_directory, coords_t, tt, t_s, mv, ss, un)
    
    
def save_data(samples_directory, genomes_directory, coords_2, coords_t, tt, g, all_inf, unin, hah, ss, ns, mv, t_un):
    
    all_inf = pd.DataFrame(data=all_inf, columns=['Total infected', 'Super spreaders', 'Normal spreaders', 'Time'])
    all_inf.to_csv(samples_directory+'/all_inf_'+'final'+'.csv', header=True, index=False)
    
    g.to_csv(genomes_directory+'/genomes_'+'final'+'.csv',header=False, index=False)
    
    unin = np.concatenate([np.column_stack(np.array((unin, t_un), dtype=float))], axis=0)
    unin = pd.DataFrame(data=unin, columns=['Uninfected individual', 'Time of uninfection'])
    unin.to_csv(samples_directory+'/uninfections.csv', header=True, index=False)
    #np.savetxt(directory+'uninfections.txt', pd.DataFrame(data = unin), fmt=['%3d','%.3f'])
    
    hah = pd.DataFrame(data=hah, columns=['Infecting', 'Infected', 'Time', 'Mutation'])
    hah.to_csv(samples_directory+'/infections.csv', header=True, index=False)
    #np.savetxt(directory+'infections.txt', pd.DataFrame(data = hah, columns=['Infects', 'Infected', 'Infection time', 'Infection Probability']), fmt=['%3d','%3d','%.3f', '%.1f'])
    
    extra_data = np.column_stack(np.array((ss, ns, mv), dtype=int))
    extra_data = pd.DataFrame(data=extra_data, columns=['Total Super spreaders', 'Total Normal spreaders', 'Total movements'])
    extra_data.to_csv(samples_directory+'/extra_data.csv', header=True, index=False)
    #np.savetxt(directory+'extra_data.txt', pd.DataFrame(data=np.column_stack(np.array((ss, ns, mv),dtype=int)), columns=['Total Super spreaders', 'Total Normal spreaders', 'Total movements']), fmt=['%1d', '%1d', '%1d'], header='Total Super infections, Total Normal infections , Total movements', comments='')

    coords_t = pd.DataFrame(data=coords_t, columns=["x","y", "label", "rate of movement", "rate of infection", "mutation", "susceptibility"])
    coords_t.to_csv(samples_directory+'/final_coords.csv', header=True, index=False)
    
    coords_2 = pd.DataFrame(data=coords_2, columns=["x","y", "label", "rate of movement", "rate of infection", "mutation", "susceptibility"])
    coords_2.to_csv(samples_directory+'/initial_coords.csv', header=True, index=False)
    

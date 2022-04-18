#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 9 10:34:25 2022

@author: laurabarbara
"""

import numpy as np
from numpy import random as rd
from matplotlib import pyplot as plt
from scipy import optimize as opt


import warnings  
warnings.filterwarnings("ignore")  # supress runtime warnings for tidyness
    

import time
startTime = time.time()  # keep track of how long each simulation takes to run


walks = 1  # number of simulations


J=20  # length of well
walkers = 1000  # number of walkers

J_values = np.array([10,20,35,50,90,150])  # length of wells to compare
walker_values = np.array([1000, 10000, 100000])  # number of walkers to compare

bessel_zeros = np.array([[2.4048, 5.5201, 8.6537, 11.7915, 14.9309], \
                         [3.8317, 7.0156, 10.1735, 13.3237, 16.4706], \
                         [5.1356, 8.4172, 11.6198, 14.7960, 17.9598], \
                         [6.3802, 9.7610, 13.0152, 16.2235, 19.4094], \
                         [7.5883, 11.0647, 14.3725, 17.6160, 20.8269], \
                         [8.7715, 12.3386, 15.7002, 18.9801, 22.2178]])

    
def exp_fit(s, q_j, rate):
    '''

    Parameters
    ----------
    s : STEPS.
    q_j : SPATIAL
        DISTRIBUTION.
    rate : RATE
        OF RREST.

    Returns
    -------
    None.

    '''
    
    p_js = q_j * np.exp(-1 * rate * s)
    
    return p_js


def bar_plot(x_data, y_data, x_label, y_label):
    '''

    Parameters
    ----------
    x_data : X COORDINATES
        OF THE BARS.
    y_data : HEIGHTS
        OF THE BARS.
    x_label : LABEL
        OF X-AXIS.
    y_label : LABEL
        OF Y-AXIS.

    Returns
    -------
    None.

    '''

    plt.bar(x_data, y_data, width=2, align='center', color='r', alpha=0.3)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.tight_layout()
    
    return


def log_plot(x_data, y_data, x_label, y_label):
    '''

    Parameters
    ----------
    x_data : X COORDINATES
        OF THE DATA POINTS.
    y_data : Y COORDINATES
        OF THE DATA POINTS.
    x_label : LABEL
        OF X-AXIS.
    y_label : LABEL
        OF Y-AXIS.

    Returns
    -------
    None.

    '''
    
    plt.yscale("log")
    plt.plot(x_data, y_data, 'r. ', markersize = 4, alpha=0.3)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.tight_layout()
    
    return


def walk(J, J_values, walkers, walker_values, bessel_zeros, n=1, nx=1, ny=1, m=1, nr=1, \
             well_1D = False, square_well_2D = False, circular_well_2D = False, \
             J_compare = False, walker_compare = False): 
    '''    

    Parameters
    ----------
    J : WELL
        LENGTH.
    J_values : WELL
        LENGTHS TO COMPARE.
    walkers : NUMBER
        OF WALKERS.
    walker_values : NUMBER
        OF WALKERS TO COMPARE.
    bessel_zeros : ZEROS 
        OF BESSEL FUNCTION OF THE FIRST KIND.
    n : QUANTUM
        NUMBER 1D WELL. Optional. The default is 1.
    nx : QUANTUM
        NUMBER 2D SQUARE WELL. Optional. The default is 1.
    ny : QUANTUM
        NUMBER 2D SQUARE WELL. Optional. The default is 1.
    m : QUANTUM
        NUMBER 2D CIRCULAR WELL. Optional. The default is 1.
    nr : QUANTUM
        NUMBER 2D CIRCULAR WELL. Optional. The default is 1.
    well_1D : 1D
        WELL. Optional. The default is False.
    square_well_2D : 2D
        SQUARE WELL. Optional. The default is False.
    circular_well_2D : 2D
        CIRCULAR WELL. Optional. The default is False.
    J_compare : COMPARE
        WELL LENGTHS. Optional. The default is False.
    walker_compare : COMPARE
        NUMBER OF WALKERS. Optional. The default is False.

    Returns
    -------
    None.

    ''' 
    
    if well_1D == True:
        
        results = np.zeros(walks)  # save result of each simulation run
        
        for walk in range(walks):
        
            steps = []  # save number of steps each walker takes before arrest
                                   
            for walker in range(walkers):  # compute walk for each walker
            
                j = J * (1 - 1/nx)  # current position
                s = 0  # steps taken
        
                caught = 'False'  # not caught yet
        
                while caught == 'False':  # as long as the walker is inside the well
            
                    s += 1  # increase step count
                    j += rd.choice([-1,1])  # take step with equal probability of stepping left or right     
                
                    if j >= J or j <= J * (1 - 2/nx):  # check if the walker leaves the well
                    
                        caught = 'True'  # arrest walker            
                        steps.append(s)  # note how many steps the walker took
            
            results[walk], r2d = analysis_wells(np.array(steps),walkers)  # data analysis function
        
            print('r squared = ', r2d)
        
        average = np.average(results)  # find mean result of all simulations
        std_dev = np.std(results)  # find standard deviation of simulations
        
        print('lambda J^2 =', average * J**2, '±', std_dev*J**2)
        print('expected =', (np.pi**2) * (nx**2)/8)
        
                
    elif square_well_2D == True:
        
        results = np.zeros(walks)  # save result of each simulation run
        
        for walk in range(walks):
        
            steps = []  # save number of steps each walker takes before arrest
                                   
            for walker in range(walkers):  # compute walk for each walker

                j1 = J * (1 - 1/nx)  # current x position
                j2 = J * (1 - 1/ny)  # current y position
                s = 0  # steps taken
    
                caught = 'False'  # not caught yet

                while caught == 'False':  # as long as the walker is inside the well
        
                    s += 1  # increase step count
                    direction = rd.choice([1,2])  # decide direction of step with equal probability                                 
                    
                    if direction == 1:  # x direction
                        j1 += rd.choice([-1,1])

                    elif direction == 2:  # y direction
                        j2 += rd.choice([-1,1])
            
                    if j1 >= J or j1 <= J * (1 - 2/nx) or j2 >= J or j2 <= J * (1 - 2/ny):  # check if the walker leaves the well
                
                        caught = 'True'  # arrest walker            
                        steps.append(s)  # note how many steps the walker took
    
            results[walk], r2d = analysis_wells(steps)  # data analysis function

            print('r squared = ', r2d)
        
        average = np.average(results)  # find mean result of all simulations
        std_dev = np.std(results)  # find standard deviation of simulations
        
        print('2 lambda J^2 =', 2 * average * J**2, '±', 2 * std_dev * J**2)
        print('expected =', (np.pi**2) * (nx**2 + ny**2)/8)


    elif circular_well_2D == True:

        results = np.zeros(walks)  # save result of each simulation run
        
        for walk in range(walks):
        
            steps = []  # save number of steps each walker takes before arrest
                                   
            for walker in range(walkers):  # compute walk for each walker

                j1 = 0  # current x position
                j2 = 0  # current y position
                s = 0  # steps taken
    
                caught = 'False'  # not caught yet

                while caught == 'False':  # as long as the walker is inside the well
        
                    s += 1  # increase step count
                    direction = rd.choice([1,2])  # decide direction of step with equal probability                                 
                    
                    if direction == 1:  # x direction
                        j1 += rd.choice([-1,1])

                    elif direction == 2:  # y direction
                        j2 += rd.choice([-1,1])
            
                    if j1**2 + j2**2 >= (J/nr)**2:  # check if the walker leaves the well
                
                        caught = 'True'  # arrest walker            
                        steps.append(s)  # note how many steps the walker took
    
            results[walk], r2d = analysis_wells(steps)  # data analysis function
            
            print('r squared = ', r2d)
        
        average = np.average(results)  # find mean result of all simulations
        std_dev = np.std(results)  # find standard deviation of simulations
        
        print('2 lambda J^2 =', 2 * average * J**2, '±', 2 * std_dev * J**2)
        print('expected =', (bessel_zeros[m-1][nr-1]**2)/2)
         
           
    elif J_compare == True:
        
        lambdas = np.zeros(len(J_values))  # save arrest rate for each well length
        
        for i in range(len(J_values)):  
            
            steps = []  # save number of steps each walker takes before arrest
        
            for walker in range(walkers):  # compute walk for each walker

                j = 0  # current position
                s = 0  # steps taken
        
                caught = 'False'  # not caught yet
        
                while caught == 'False':  # as long as the walker is inside the well
            
                    s += 1  # increase step count
                    j += rd.choice([-1,1])  # take a step with equal probability of stepping left or right     
                    
                    if j >= J_values[i] or j <= -J_values[i]:  # check if the walker leaves the well
                        
                        caught = 'True'  # arrest walker            
                        steps.append(s)  # note how many steps the walker took
                    
            lambdas[i] = analysis_compare(steps)  # data analysis function
    
        plot_J_compare(J_values, lambdas)
           

    elif walker_compare == True:

        lambdas = np.zeros(len(walker_values))  # save arrest rate for each well length
        
        for i in range(len(walker_values)):  

            steps = []  # save number of steps each walker takes before arrest
        
            for walker in range(walker_values[i]):  # compute walk for each walker

                j = 0  # current position
                s = 0  # steps taken
        
                caught = 'False'  # not caught yet
        
                while caught == 'False':  # as long as the walker is inside the well
            
                    s += 1  # increase step count
                    j += rd.choice([-1,1])  # take a step with equal probability of stepping left or right     

                    if j >= J or j <= -J:  # check if the walker leaves the well
                    
                        caught = 'True'  # arrest walker            
                        steps.append(s)  # note how many steps the walker took
                    
            lambdas[i] = analysis_compare(steps)  # data analysis function
            
        plot_walker_compare(walker_values, lambdas, J)

    return


def analysis_wells(steps, walkers):
    '''

    Parameters
    ----------
    steps : TOTAL
        STEPS TAKEN BY WALKERS.
    walkers : NUMBER
        OF WALKERS.

    Returns
    -------
    popt[1]: RATE
        OF ARREST.
    r2d : R-SQUARED
        VALUE.

    '''
    
    plt.figure(1)     
    # collect data into histogram to get number of arrests per number of steps
    hist_data = plt.hist(steps, density = False, color = 'r', alpha=0.5, \
        bins = int((max(steps) - min(steps) + 1)/2), align = 'mid')  
    plt.close(1)
    
    s = hist_data[1][:-1]  # array with number of steps
    arrests = hist_data[0][:]  # array with number of arrests at each s
    
    cutoff = np.argmax(arrests)  # define peak of arrests as low cutoff
    
    survivals = walkers-np.cumsum(arrests)  # number of surviving walkers at each s
    prob_survivals = survivals/max(survivals)  # probability of a walker surviving until each s

    s = s[cutoff:]  # remove low n data points
    prob_survivals = prob_survivals[cutoff:]  # remove low n data points
    
    # fitting an exp curve to data
    popt, pcov = opt.curve_fit(exp_fit, s, prob_survivals, p0=[1,0])  
    fitted = exp_fit(s, popt[0], popt[1])  # get data from fit
    
    # calculate r-squared
    residuals = prob_survivals - fitted
    sumsquares_residuals = np.sum(residuals**2)  
    sumsquares_total = np.sum((prob_survivals-np.mean(prob_survivals))**2)   
    r2d = 1 - (sumsquares_residuals/sumsquares_total)
    
    plot_wells(s, prob_survivals, fitted, popt, pcov)    
    
    return popt[1], r2d


def plot_wells(s, prob_survivals, fitted, popt, pcov):
    '''

    Parameters
    ----------
    s : STEP
        AMOUNTS.
    prob_survivals : PROBABILITY
        OF SURVIVING UP TO S STEPS.
    fitted : DATA
        FROM FIT.
    popt : FIT
        PARAMETERS.
    pcov : COVARIANCE
        MATRIX.

    Returns
    -------
    None.

    '''
    
    plt.figure(2)    
    # define x and y labels for the plots
    x_label = 'Number of Steps Taken (s)'
    y_label = 'Probability of Survival' + r'$(p(j_1,j_2,s))$'                
    # bar plot of probability of survival
    bar_plot(s, prob_survivals, x_label, y_label)
    plt.plot(s, fitted, 'k-', label=r'$p(j_1,j_2,s)=($' + str(np.round(popt[0],5)) \
             + r'$)\exp(-$' + str(np.round(popt[1],5)) + r'$ s)$')
    plt.legend()
    #plt.savefig('bar plot.png', dpi=1000, bbox_inches='tight')
    plt.show()        
    plt.close(2)
    
    plt.figure(3)
    # scatter log plot of probability of survival    
    log_plot(s, prob_survivals, x_label, y_label)
    plt.plot(s, fitted, 'k-', label=r'$p(j_1,j_2,s)=($' + str(np.round(popt[0],5)) \
             + r'$)\exp(-$' + str(np.round(popt[1],5)) + r'$ s)$')
    plt.legend()
    #plt.savefig('log plot.png', dpi=1000, bbox_inches='tight')
    plt.show()        
    plt.close(3)

    return


def analysis_compare(steps):

    '''

    Parameters
    ----------
    steps : TOTAL
        STEPS TAKEN BY WALKERS.

    Returns
    -------
    popt[1]: RATE
        OF ARREST.

    '''

    plt.figure(1)     
    # collect data into histogram to get number of arrests per number of steps
    hist_data = plt.hist(steps, density = False, color = 'r', alpha=0.5, \
        bins = int((max(steps) - min(steps) + 1)/2), align = 'mid')  
    plt.close(1)
    
    s = hist_data[1][:-1]  # array with number of steps
    arrests = hist_data[0][:]  # array with number of arrests at each s
    
    cutoff = np.argmax(arrests)  # define peak of arrests as low cutoff
    
    survivals = walkers-np.cumsum(arrests)  # number of surviving walkers at each s
    prob_survivals = survivals/max(survivals)  # probability of a walker surviving until each s

    s = s[cutoff:]  # remove low n data points
    prob_survivals = prob_survivals[cutoff:]  # remove low n data points
    
    # fitting an exp curve to data
    popt, pcov = opt.curve_fit(exp_fit, s, prob_survivals, p0=[1,0])  

    return popt[1]


def plot_J_compare(J_values, lambdas):
    '''

    Parameters
    ----------
    J_values : WELL
        LENGTHS TO COMPARE.
    lambdas : ARREST
        RATE FOR EACH WELL LENGTH.

    Returns
    -------
    None.

    '''

    plt.figure(7)
    # define x and y labels for the plots
    plt.xlabel('Well width (J)')
    plt.ylabel(r'$\lambda J^2$')
    # scatter plot of calculated energies for each well length 
    plt.plot(J_values, lambdas*J_values**2, '.r ', markersize = 13)
    plt.axhline((np.pi**2)/8, color='k', label=r'$\pi^2/8$')
    plt.legend()
    plt.tight_layout()
    #plt.savefig('J comparison.png', dpi=1000, bbox_inches='tight')
    plt.show()

    return


def plot_walker_compare(walker_values, lambdas, J):
    '''

    Parameters
    ----------
    walker_values : NUMBER
        OF WALKERS TO COMPARE.
    lambdas : ARREST
        RATE FOR EACH WELL LENGTH.
    J : WELL
        LENGTH.

    Returns
    -------
    None.

    '''

    plt.figure(7)  # clear initialized plot
    # define x and y labels for the plots
    plt.xlabel('Number of Drunkards')
    plt.ylabel(r'$\lambda J^2$')
    # scatter plot of calculated energies for each number of walkers
    plt.plot(walker_values, lambdas*J**2, '.r ', markersize = 13)
    plt.axhline((np.pi**2)/8, color='k', label=r'$\frac{\pi^2}{2}$')
    plt.tight_layout()
    #plt.savefig('walker comparison.png', dpi=1000, bbox_inches='tight')
    plt.show()

    return


walk(J, J_values, walkers, walker_values, bessel_zeros)

executionTime = (time.time() - startTime)
print('Execution time in seconds: ' + str(executionTime))
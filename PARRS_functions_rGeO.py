# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 10:49:49 2024

@author: Moritz Meißner (@PDI)

This script contains all the functions necessary for the fit. You must change the Raman Tensor function (or angular functions) for each material and measurement.
The naming given by me contains the mode and the scattering geometry.
For example the B1g mode in the (100) plane with perpendicular polarized backscattered light is called: B1g_A_par
You can choose the naming of your functions as you like, but they must be selectable when going through the for loop of the modes and geometries.
That means if your mode is labeled 'A1g_4' and you geometry (combination of plane and polarization) is (010)_perp, both should be in the function name.

The general Fit function is at the end of this script.
"""
import numpy as np
from numpy import pi, sin, cos, sqrt
import matplotlib.pyplot as plt

import inspect

def B1g_A_par(phi,c ,min_ampl_B1g_ap,phi_offset):
    phi=phi+phi_offset
    I = c**2*sin(phi)**4+min_ampl_B1g_ap
    return I
def B1g_A_perp(phi,c ,min_ampl_B1g_as,phi_offset):
    phi=phi+phi_offset
    I = c**2*cos(phi)**2 *sin(phi)**2+min_ampl_B1g_as
    return I
def B1g_C_par(phi,c ,min_ampl_B1g_cp,phi_offset):
    phi=phi+phi_offset
    I = c**2*cos(2*phi)**2 + min_ampl_B1g_cp
    return I

def B1g_C_perp(phi,c, min_ampl_B1g_cs,phi_offset ):
    phi=phi+phi_offset
    I = c**2*sin(2*phi)**2 + min_ampl_B1g_cs
    return I
def B1g_110_par(phi,min_ampl_B1g_dp,phi_offset):
    phi=phi+phi_offset
    return min_ampl_B1g_dp
def B1g_110_perp(phi,min_ampl_B1g_ds,phi_offset):
    phi=phi+phi_offset
    return min_ampl_B1g_ds
#### Eg1 ###
def E1g_A_par(phi,e , min_ampl_Eg_ap,phi_offset):
    phi=phi+phi_offset
    I = e**2*sin(2*phi)**2 + min_ampl_Eg_ap
    return I

def E1g_A_perp(phi,e, min_ampl_Eg_as ,phi_offset):
    phi=phi+phi_offset
    I = e**2*cos(2*phi)**2+ min_ampl_Eg_as
    return I
def E1g_110_par(phi,e , min_ampl_Eg_dp,phi_offset):
    phi=phi+phi_offset
    I = 2*e**2*cos(phi)**2*sin(phi)**2+ min_ampl_Eg_dp
    return I
def E1g_110_perp(phi,e, min_ampl_Eg_ds,phi_offset ):
    phi=phi+phi_offset
    I = (1/2)*e**2*cos(2*phi)**2+ min_ampl_Eg_ds
    return I
def E1g_C_par(phi,min_ampl_Eg_cp,phi_offset):
    phi=phi+phi_offset
    return min_ampl_Eg_cp
def E1g_C_perp(phi,min_ampl_Eg_cs,phi_offset):
    phi=phi+phi_offset
    return min_ampl_Eg_cs

### A1g ###
def A1g_A_par(phi,a,b, min_ampl_A1g_ap,phi_offset ):
    phi=phi+phi_offset
    I =(b*cos(phi)**2 + a*sin(phi)**2)**2 + min_ampl_A1g_ap
    return I
def A1g_A_perp(phi,a,b , min_ampl_A1g_as,phi_offset):
    phi=phi+phi_offset
    I =(a - b)**2*cos(phi)**2*sin(phi)**2+ min_ampl_A1g_as
    return I

def A1g_C_par(phi,a,min_ampl_A1g_cp,phi_offset ):
    phi=phi+phi_offset
    I =a**2 + min_ampl_A1g_cp
    return I
def A1g_C_perp(phi,min_ampl_A1g_cs,phi_offset):
    phi=phi+phi_offset
    return min_ampl_A1g_cs

def A1g_110_par(phi,a,b , min_ampl_A1g_dp,phi_offset):
    phi=phi+phi_offset
    I =(b*cos(phi)**2 + a*sin(phi)**2)**2+ min_ampl_A1g_dp
    return I
def A1g_110_perp(phi,a,b,min_ampl_A1g_ds,phi_offset ):
    phi=phi+phi_offset
    I =(a - b)**2*cos(phi)**2*sin(phi)**2+ min_ampl_A1g_ds
    return I

### B2g ###
def B2g_A_par(phi,min_ampl_B2g_ap,phi_offset):
    phi=phi+phi_offset
    return min_ampl_B2g_ap
def B2g_A_perp(phi,min_ampl_B2g_as,phi_offset):
    phi=phi+phi_offset
    return min_ampl_B2g_as
def B2g_C_par(phi,d,min_ampl_B2g_cp,phi_offset ):
    phi=phi+phi_offset
    I =d**2*sin(2*phi)**2+min_ampl_B2g_cp
    return I
def B2g_C_perp(phi,d ,min_ampl_B2g_cs,phi_offset):
    phi=phi+phi_offset
    I =d**2*cos(2*phi)**2+min_ampl_B2g_cs
    return I
def B2g_110_par(phi,d ,min_ampl_B2g_dp,phi_offset):
    phi=phi+phi_offset
    phi = phi -1.55046023
    I =d**2*sin(phi)**4+min_ampl_B2g_dp
    return I
def B2g_110_perp(phi,d ,min_ampl_B2g_ds,phi_offset):
    phi=phi+phi_offset
    phi = phi -1.55046023
    I =d**2*cos(phi)**2*sin(phi)**2+min_ampl_B2g_ds
    return I

# These functions are called in the inner for loop of the general function to calculate the intensity along the energy axis.
def gauss(x, A, x0, sigma):
    return A/(np.sqrt(2*np.pi*sigma**2))* np.exp(-(x - x0)**2 / (2 * sigma**2)) #Gauß

def exp(x,a,t0):
    return a*np.exp(-x*t0)

def lorentzian(x, amplitude, center, sigma):
    return amplitude / (1 + ((x - center) / sigma) ** 2)

######################## General Fit Function ##########################

#This function gets called by the lmfit.model.fit function in the main script.
def Lorentzian_2D_general(energy:list, angle:list, mode_dict:dict,**params):# params is a lmfit Parameters class
    #get the angles and convert them into radians from degree
    phi = np.radians(angle)
    # get the modes that will be fitted
    mode_list = mode_dict['modes']
    I_dict={} #define the dictionary in which the calculated intesities for each mode in each scattering geometry are stored (local)
    geometries = mode_dict['geometries']
    suffix = mode_dict['suffix']     #suffis lists
    #get all functions which are contained in this file
    current_module = inspect.getmodule(inspect.currentframe())
    functions = inspect.getmembers(current_module, inspect.isfunction)
    functions = [func for func in functions]
    
    # Collect all variables in the current scope, they are stored inb a dict
    available_variables = params
    available_variables['energy'] = energy
    available_variables['phi'] = phi
    available_variables['phi_offset'] = 0

    #make troubbleshooting easier with some exceptions:
    if len(geometries) == 0:
        raise Exception('ERROR: 0 scattering geometries were passed.')
    if len(geometries) != len(suffix):
        raise Exception(f'ERROR: The length of geometries ({len(geometries)}) must be the same as suffix({len(suffix)})')
    if len(available_variables.keys()) <= 3:
        raise Exception(f'Something went wrong with the available_variables, there length is {len(available_variables.keys())}')
    if len(params) ==0:
        raise Exception('ERROR: The length of the Parameters is 0!')
        
    # Go through each mode and each scattering intensity and calculate their 2d intensity.
    for mode in mode_list:
        counter=0 #used to go through the suffix lists
        for geometry in  geometries:
            for func in functions: # Go through each availablefunction and check if its the correct one, i.e. right mode and geometry.
                if mode in func[0] and geometry in func[0]:  
                    func=func[1]
                    if not callable(func):
                        print(f"Skipping non-callable object: {func}")
                        continue
                    #Now the by the current function required parameters are calles from "available_variables".
                    phi_offset = [available_variables[i]  for i in available_variables if suffix[counter][0] in i and 'phi_offset' in i ][0]
                    available_variables['phi_offset'] = phi_offset
                    # Get the signature of the function
                    sig = inspect.signature(func)
                    # Extract the required parameters
                    required_params = sig.parameters.keys()
                    # Filter available variables to match required parameters
                    arguments = {param: available_variables[param] for param in required_params if param in available_variables}
                    # Call the function with the filtered arguments
                    amplitude = func(**arguments)
                    # get the center, sigma,offset parameters from the parameters
                    center =  [available_variables[i]  for i in available_variables if mode in i and 'center' in i][0]
                    sigma = [available_variables[i] for i in available_variables if mode in i and 'sigma' in i][0]
                    offset = [available_variables[i] for i in available_variables if suffix[counter] in i and 'y' in i][0]
  
                    ### run the calculation  ###                
                    Intensity = lorentzian(energy, amplitude, center, sigma) + offset
                    I_dict[f'{mode}_{geometry}'] = Intensity #save the calculation
                    
            counter += 1
    if len(I_dict) == 0:
        raise Exception('ERROR: Intensities could not be calculated, the dictionary is empty.')
        
    # Intenstiy_stacked = []
    results=[]
    for geometry in geometries:
        # Collect lists where 'A_par' is in the key
        selected_lists = [value for key, value in I_dict.items() if geometry in key]
        
        # Convert the lists to a numpy array and sum element-wise along axis 0
        results.append(np.sum(np.array(selected_lists), axis=0) )

    Intensity_stacked = np.vstack(results) 
      
    if len(geometries) != len(Intensity_stacked[:,0]):
        raise Exception(f'ERROR: Intensities are not stacked correctly. There are {len(geometries)} but the the vstack has only dim{len(Intensity_stacked[:,0])}')

    return Intensity_stacked






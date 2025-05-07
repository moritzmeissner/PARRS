# -*- coding: utf-8 -*-
"""
Created on Mon Nov  4 16:04:54 2024

@author: Moritz Meißner (@PDI Berlin)

This is the main script. 
Attention: You MUST change the file names of the input and functions files to the ones the are required for the current analysis!
"""

import numpy as np
from numpy import pi, sin, cos, sqrt
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import lmfit
from lmfit import Parameters, Minimizer, Model
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os
from os.path import join
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from matplotlib import cbook, cm
from matplotlib.colors import LightSource
from lmfit.models import LorentzianModel
from datetime import datetime



############ !!! change .py file name accordingly!!! #########
from PARRS_Plot_functions import  Plot_arrays,make_single_plot, plot_data_array,  Plot_3D_arrays, Plot_ang_dependency, plot_single_spectra,PDF_merge
from PARRS_Input_rGeO import mode, normalize,plot_2D_Fit,plot_1D_ang_dependency,plot_3D_Fit,save_plots,save_fit_report,merge_plots,mode_list
from PARRS_Input_rGeO import Energy_start, Energy_end, material, data_path
from PARRS_Input_rGeO import start_values

import PARRS_functions_rGeO as fct
from PARRS_functions_rGeO import Lorentzian_2D_general


class read_data:
    def __init__(self, data_path: str,material:str,  limits: list, normalize: bool):
        print(f'Reading data from {data_path}, the following files are read:')
        file_list, data1D_list = [],[]
        data_dic={}
      
        
        energy_raw = np.genfromtxt(join(data_path,f'{material}_PARRS_energy.txt'))
        lower_limit_idx = np.where(energy_raw >= limits[0])[0][0]
        upper_limit_idx = np.where(energy_raw >= limits[1])[0][0]
        energy_list = energy_raw[lower_limit_idx:upper_limit_idx]
  
        angle_list = np.genfromtxt(join(data_path,f'{material}_PARRS_angle.txt'))
        
        ###1. make energies and angles 2D data, 2. flatten them into long 1d arrays ###
        energy_2d, angle_2d = np.meshgrid(energy_list, angle_list) # 2D x & y
        size = energy_2d.shape
        # make 1D data from 2D for fitting
        energy_flattend = energy_2d.reshape((1, np.prod(size)))
        angle_flattend = angle_2d.reshape((1, np.prod(size)))
        # make it into a vstack for the fit and save it in the class
        energy_stacked = np.vstack(energy_flattend)
        angle_stacked = np.vstack(angle_flattend)
        
        self.energy_for_fit =  energy_stacked
        self.angle_for_fit = angle_stacked
        self.array_size = size
        self.energy2D = energy_2d
        self.angle2D = angle_2d        
   
        
        for file in os.listdir(data_path):
            if not '.txt' in file:
                continue
            if 'energy' in file or 'angle' in file:
                continue
            file_name=file[:-4]
            print(file)
            data = np.genfromtxt(join(data_path,file), delimiter=',')
            # cut the data according to limits
            data = data[:,lower_limit_idx : upper_limit_idx]
            
            file_list.append(file_name)
            data_dic[file_name] = data
            data_flatten = np.ravel(data)
            data1D_list.append(data_flatten)

                
                
        self.energy = energy_list
        self.angle = angle_list
        self.data1D = np.vstack(data1D_list)
        self.files = np.vstack(file_list)
        self.data2D = data_dic
        
        print('Finished data import.')


limits = [Energy_start, Energy_end] #in cm^-1      
normalize = False
save_fit_report = False
     
print(f'Material: {material}, data path: {data_path}')   
PARRS = read_data(data_path, material, limits, normalize)  

Data = PARRS.data1D    
energy_for_fit  = PARRS.energy_for_fit
angle_for_fit  = PARRS.angle_for_fit

geometries_list = [i[0] for i in PARRS.files]
       
    
def save_result_txt(mode,result,normalize,geometries_list):
    # Save fit report to a text file
    if normalize == False:
        output_file = f"fit_report_{mode}.txt"
    elif normalize == True:
        output_file = f"fit_report_{mode}_normalized.txt"

    with open(output_file, "w") as file:
        file.write(f'Mode: {mode} \n')
        file.write(f'normalized = {normalize} \n')
        file.write(f'planes: {geometries_list}\n')
        file.write(result.fit_report())
    
    print(f"Fit report saved to {output_file}")

def write_files_txt(data_dic):
    # Save fit report to a text file
    for key in data_dic:
        output_file = os.path.join('r-GeO-korrigierte_Spektren',f'{key}_korrigiert.txt')
        energy = data_dic[key][0]
        intensity = data_dic[key][1]
        
        with open(output_file, "w") as file:
            for i in range(len(intensity)):
                file.write(f'{energy[i]} {intensity[i]}\n ')
        
    print(f"Fit report saved to {output_file}")    
# write_files_txt(data_dic)

def adjust_vary_params(params, modes):
    """
    Adjust the 'vary' attribute of parameters in an lmfit Parameters object.

    Args:
        params (Parameters): The Parameters object from model.make_params().
        modes (list): A list of strings. If a parameter name contains any of these strings,
                      the 'vary' attribute remains unchanged.
    
    Returns:
        Parameters: The modified Parameters object.
    """
    for param_name, param in params.items():
        # print(modes)
        # print(param_name,param)
        # Check if "min_ampl" is in the name
        if "min_ampl" in param_name :
            # Check if none of the modes are in the parameter name
            if not any(mode in param_name for mode in modes['modes']):
                print('no vary:',param_name)
                param.vary = False  # Set 'vary' to False

    return params

#################################    Fit  ###################################################
model = Model(Lorentzian_2D_general, independent_vars=['energy','angle', 'mode_dict'])

# Adjust parameters
adjusted_params = adjust_vary_params(start_values, mode_list)

#### make the Fit ####
print('Fit start')

startTime = datetime.now()
mode_list['geometries'] = geometries_list

result = model.fit(data=Data, energy  = energy_for_fit, angle  = angle_for_fit, 
                   params = start_values, method='powell', mode_dict = mode_list)

print(f'Fit finished in {datetime.now() - startTime}')
# least_squared, , differential_evolution, cobyla(ggod),’tnc’: Truncated Newton(not good),'bfgs' good,powell fast,karismos mobark bad,cg bad,lbfgsb good,
z_fit = result.best_fit #result is 1d

# Save parameters containing "y" as variables
variables = {}
y_list, min_ampl_list=[],[]
for param_name, param_value in result.params.items():
    if "y" in param_name:
        var_name = param_name
        variables[var_name] = param_value.value  # Save value
        # Dynamically create variables
        exec(f"{var_name} = {param_value.value}")
        y_list+=[param_value.value]
    if "min_ampl" in param_name:
        var_name = param_name
        variables[var_name] = param_value.value  # Save value
        # Dynamically create variables
        exec(f"{var_name} = {param_value.value}")
        min_ampl_list+=[param_value.value]


print(result.fit_report())

if save_fit_report == True:
    save_result_txt(mode,result,normalize, geometries_list)


peak_energies = [start_values[f'center{mode}'] for mode in mode_list['modes']]

#### Plot Results of 2D Fit for each Raman mode ####
pdfs=[]
if plot_2D_Fit == True:
    for i in range(len(geometries_list)):
      plane = geometries_list[i]
      intensity_2d = Data[i].reshape(PARRS.array_size) #- y_list[i] - min_ampl_list[i]
      Z_fit = z_fit[i].reshape(PARRS.array_size)# - y_list[i]- min_ampl_list[i]
      pdfs = Plot_arrays(PARRS.energy2D, PARRS.angle2D, intensity_2d, Z_fit, plane,mode,pdfs,normalize,save_plot=save_plots) #2D array plot

if plot_1D_ang_dependency == True:
    for energy,mode in  zip(peak_energies,mode_list['modes']):
        index1 = min(range(len(PARRS.energy)), key=lambda i: abs(PARRS.energy[i] - energy))
        for i in range(len(geometries_list)):
          plane = geometries_list[i]
          intensity_2d = Data[i].reshape(PARRS.array_size)# - y_list[i]- min_ampl_list[i]
          Z_fit = z_fit[i].reshape(PARRS.array_size)# - y_list[i]- min_ampl_list[i]
          pdfs = Plot_ang_dependency(PARRS.energy2D,PARRS.angle2D,intensity_2d,Z_fit,index1,plane,mode,pdfs,normalize,save_plot=save_plots) #ang. dependency plot

if plot_3D_Fit == True:
    for i in range(len(geometries_list)):
      plane = geometries_list[i]
      intensity_2d = Data[i].reshape(PARRS.array_size)# - y_list[i]- min_ampl_list[i]
      Z_fit = z_fit[i].reshape(PARRS.array_size) #- y_list[i]- min_ampl_list[i]
      pdfs = Plot_3D_arrays(PARRS.energy2D, PARRS.angle2D, intensity_2d, Z_fit, plane,mode,pdfs,normalize,save_plot=save_plots) # 3D array plot

if merge_plots == True:
    if normalize == False:
        PDF_merge(pdfs,name= os.path.join('Plots',f'{mode}_merged.pdf'))
    elif normalize == True:
        PDF_merge(pdfs,name= os.path.join('Plots',f'{mode}_merged_normalized.pdf'))




# -*- coding: utf-8 -*-
"""
Created on Dec 2024

@author: Moritz Meißner (@PDI)

This is the input file. This has to be changed for every project.
Here you define what modes should be fitted, what plots should be shown and you create the initual values for the fit.

"""
data_path = 'rGeO' # The directory in which the data is stored, in a relative to this script path. 
material = 'rGeO'  # for the naming of angle and energy list.

#for normalizing, this works only for rGeo at the moment. Therfore set normalie to False
mode = 'A1g'
normalize = False

### Plotting Options ###
plot_2D_Fit = True
plot_1D_ang_dependency = False
plot_3D_Fit = False

save_plots = False #save the plots
save_fit_report = False
merge_plots = False

# You can choose the mode names as needed, make sure they match the chosen names for the Raman tensor functions.
mode_list = {'modes':['A1g']}#'B1g','Eg','A1g','B2g'
# the suffix is needed in the general fit function to select the right parameter for each Raman function. "ap" stands for a plane - perpendicular.
# This is used to shorten the variable names. The suffix can be chosen freely, make sure they match the parameter names in the 'start_values'.
# There should be as many suffixes as there are scattering geometries, i.e. 6 if you have a (100),(001)&(110) plane in parallel and perpendicular backscattering geometry.
mode_list['suffix'] =  ['ap','cp','dp','as','cs','ds']

# Set up the energy range in which you want to fit. Start with only one mode and select the energy range accordingly. Increase the amount of modes needed.
#upper limit should be 900 cm^-1 for rGeO
Energy_start, Energy_end = 600, 900 #in cm^-1

from lmfit import Parameters
start_values = Parameters()

# The Parameter class is set up like this : # add with tuples: (NAME VALUE VARY MIN  MAX  EXPR  BRUTE_STEP)
# The parameter names are given by the raman tensor functions, which are defined in PARRS_functions.py. You can change them to whatever you want.
# The 'center ' and 'sigma' names are given by the lorentzian function.  
start_values.add_many(
    ('centerA1g',696.72,True, 650,730 ),('sigmaA1g', 4.64, True, 1, 7),  
    ('a', -82),('b', 154), ('c',117 ), ('d', 67), ('e', 1), #Raman Tensor Elements
    ('phi_offset_ap',0.01134185),('phi_offset_as',0),('phi_offset_cp',-0.05694673),('phi_offset_cs',0),('phi_offset_dp',0.02303618),('phi_offset_ds',0),
    ('y_ap',206,True,0),('y_as',90,True,0),('y_cp',1066,True,0),('y_cs',846,True,0),('y_dp',237,True,0),('y_ds',123,True,0),
    ('min_ampl_A1g_ap',10298,True,10000), ('min_ampl_A1g_as',1832,True,1000),('min_ampl_A1g_cp',24834,True,20000),('min_ampl_A1g_cs',3458,True,3200),('min_ampl_A1g_dp',11272,True,10000),('min_ampl_A1g_ds',2482,True,2000),

                    )

   # centerB1g=dict(value=PE1,min=PE1-20,max=PE1+20), centerEg=dict(value=PE2,min=PE2-20,max=PE2+20), centerA1g=dict(value=PE3,min=PE3-20,max=PE3+20),centerB2g=dict(value=PE4,min=PE4-20,max=PE4+20),
   # sigmaB1g=dict(value=4.64,min = 1, max = 7),sigmaEg=dict(value=4.64,min = 1, max = 7),sigmaA1g=dict(value=4.64,min = 1, max = 7), sigmaB2g=dict(value=4.64,min = 1, max = 7),
   # a= -82, b=154, c=dict(value=117,min=90,max=130),d = dict(value=67,min =50,max=80), e=dict(value=1,max=5,min=0),
   
   # min_ampl_B1g_ap=dict(value=302,min=300,max=450), min_ampl_B1g_as=dict(value=67),min_ampl_B1g_cp=dict(value=6460),min_ampl_B1g_cs=dict(value=4000,min=0),min_ampl_B1g_dp=dict(value=1955,min=0),min_ampl_B1g_ds=dict(value=191,min=0),
   # min_ampl_Eg_ap=dict(value=170,min=0), min_ampl_Eg_as=dict(value=40,min=0),min_ampl_Eg_cp=dict(value=2285,min=0),min_ampl_Eg_cs=dict(value=1740,min=0),min_ampl_Eg_dp=dict(value=263,min=0),min_ampl_Eg_ds=dict(value=111,min=0),
   # min_ampl_A1g_ap=dict(value=10298,min=10000), min_ampl_A1g_as=dict(value=1832,min=1000),min_ampl_A1g_cp=dict(value=24834,min=20000),min_ampl_A1g_cs=dict(value=3458,min=3200),min_ampl_A1g_dp=dict(value=11272,min=10000),min_ampl_A1g_ds=dict(value=2482,min=2000),
   # min_ampl_B2g_ap=dict(value=75,min=0), min_ampl_B2g_as=dict(value=60,min=0),min_ampl_B2g_cp=dict(value=1,min=0),min_ampl_B2g_cs=dict(value=337,min=0),min_ampl_B2g_dp=dict(value=208,min=0),min_ampl_B2g_ds=dict(value=340,min=0),
   
   # phi_offset_ap=0.01134185,phi_offset_as=0,phi_offset_cp=-0.05694673,phi_offset_cs=0,phi_offset_dp=0.02303618,phi_offset_ds=0,

  
   # y_ap=dict(value=206,min=0),y_as=dict(value=90,min=0),y_cp=dict(value=1066,min=0),y_cs=dict(value=846,min=0),y_dp=dict(value=237,min=0),y_ds=dict(value=123,min=0)
   # )

"""
# Energies for slicing, in cm^-1
# target1,target2= 750,920 #B2g
# target1,target2= 650,730 #A1g
# target1,target2= 500,550 #E1g
# target1,target2= 140,300 #B1g

# Peak_energy_B1g = 166
# Peak_energy_B2g = 869.54
# Peak_energy_A1g =696.72
# Peak_energy_E1g =521
"""



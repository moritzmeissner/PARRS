# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 16:39:35 2024

@author: meissner
"""

import numpy as np
from numpy import pi, sin, cos, sqrt
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import lmfit
from lmfit import Parameters, Minimizer, Model
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

from matplotlib import cbook, cm
from matplotlib.colors import LightSource

def make_colorbar(fig,im,ax):
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(im, cax=cax, orientation='vertical')

def Plot_arrays(X, Y, intensity_2d,Z_fit,plane,mode,pdfs,normalize,save_plot = False):
    limits = [X[0,0],X[0,-1],Y[0,0],Y[-1,0]]
    fig, (ax1, ax2,ax3) = plt.subplots(1,3, figsize=(12, 6))
    ax1.set_title("Data",fontsize=12)
    im1 = ax1.pcolormesh(X, Y, intensity_2d)
    ax1.axis(limits)
    make_colorbar(fig,im1,ax1)
      
    ax2.set_title("Fit",fontsize=12)
    im2 = ax2.pcolormesh(X, Y, Z_fit)
    ax2.axis(limits)
    make_colorbar(fig,im2,ax2)
    
    ax3.set_title("Residual",fontsize=12)
    im3 = ax3.pcolormesh(X, Y, intensity_2d- Z_fit)
    ax3.axis(limits)
    make_colorbar(fig,im3,ax3)

    ax1.set_xlabel(r'Raman Shift ($cm^{-1}$)',fontsize=12),ax2.set_xlabel(r'Raman Shift ($cm^{-1}$)',fontsize=12), ax3.set_xlabel(r'Raman Shift ($cm^{-1}$)',fontsize=12)
    ax1.set_ylabel(r'Polarization Angle (°)',fontsize=12)
    plt.suptitle(f'r-GeO {plane}, Mode {mode} - 2D Fit -Result ',fontsize=12)
    fig.tight_layout(h_pad=0.0)#, plt.subplots_adjust(wspace=0.01, hspace=0.01) # make smaller / remove space between plots
    if normalize == False:
        path = os.path.join('Plots',f'rGeO{plane}_{mode}_2DFit.pdf')
    else:
        path = os.path.join('Plots',f'rGeO{plane}_{mode}_2DFit_normalized.pdf')

    if save_plot != False:
        plt.savefig(path, bbox_inches='tight')
        print('plot saved')
    plt.show()
    pdfs += [path]
    return pdfs

def make_single_plot(energy,intensity,measurement,save_plot = False): 
    fig, (ax1) = plt.subplots(1,figsize=(8,5),dpi=800)
    # ax2 = ax1.twinx()
    
    ax1.plot(energy,intensity,label=measurement[:4])
    ax1.set_xlabel(r'Energy ($cm^{-1}$)'),  ax1.set_ylabel('Intensity (arb. units)')#, ax1.set_ylim(np.min(intensity_2d[:,Index ])-50)
    ax1.tick_params(which='both',axis='both', direction = 'in',labelsize=8),
    ax1.legend(fontsize=8,loc='upper left')#,ax2.legend(fontsize=8,loc='upper right'),
    plt.title(f'{measurement}')
    ax1.xaxis.set_major_locator(MultipleLocator(50)), ax1.xaxis.set_major_formatter('{x:0.0f}'),ax1.xaxis.set_minor_locator(MultipleLocator(10))

    # ax2.set_ylabel('Residual',fontsize=8),ax2.patch.set_visible(False) 
    # ax1.patch.set_visible(False)#,    
    
    if save_plot != False:
        plt.savefig(f'test.pdf', bbox_inches='tight')
        print('plot saved')
    plt.show() 
    
def plot_data_array(energy_list, angle_list,data_array,plane,mode, save=False) :
    X_2d, Y_2d = np.meshgrid(energy_list, angle_list) #now we have an array
     
    limits = [energy_list[0],energy_list[-1],angle_list[0],angle_list[-1]]

    fig, (ax1) = plt.subplots(1,1, figsize=(12, 6))
    ax1.set_title("Data")
    im1 = ax1.pcolormesh(X_2d, Y_2d, data_array)
    ax1.axis(limits)
    make_colorbar(fig,im1,ax1)

    ax1.set_xlabel(r'Raman Shift ($cm^{-1}$)',fontsize=12),ax1.set_ylabel(r'Polarization Angle (°)',fontsize=12)
    plt.suptitle(f'r- GeO {plane} plane - Mode {mode}  ',fontsize=12)
    fig.tight_layout(h_pad=0.0)#, plt.subplots_adjust(wspace=0.01, hspace=0.01) # make smaller / remove space between plots
    path = os.path.join('bGaO Fit Plots',f'bGaO_{plane}_{mode}_data.pdf')
    if save != False:
        plt.savefig(path, bbox_inches='tight')
        print(f'saved {plane} 2D Plot')
    plt.show()
    
def Plot_3D_arrays(X, Y, intensity_2d,Z_fit,plane,mode,pdfs,normalize,save_plot = False):
    limits = [X[0,0],X[0,-1],Y[0,0],Y[-1,0]]
    fig, (ax1, ax2,ax3) = plt.subplots(1,3, figsize=(18, 7),subplot_kw=dict(projection='3d'),dpi=600)

    ls = LightSource(270, 45)
    rgb = ls.shade(intensity_2d, cmap=cm.gist_earth, vert_exag=0.1, blend_mode='soft')
   
    ax1.set_title("Data",fontsize=12)
    # To use a custom hillshading mode, override the built-in shading and pass
    # in the rgb colors of the shaded surface calculated from "shade".
    surf1 = ax1.plot_surface(X, Y, intensity_2d, rstride=1, cstride=1, facecolors=rgb,
                           linewidth=0, antialiased=False, shade=False, edgecolor='none')
  
      
    ax2.set_title("Fit",fontsize=12)
    surf2 = ax2.plot_surface(X, Y, Z_fit, rstride=1, cstride=1, facecolors=rgb,
                           linewidth=0, antialiased=False, shade=False)
    
    ax3.set_title("Residual",fontsize=12)
    surf3 = ax3.plot_surface(X, Y, intensity_2d- Z_fit, rstride=1, cstride=1, facecolors=rgb,
                           linewidth=0, antialiased=True, shade=False)

    for ax in [ax1, ax2,ax3]:
        ax.set_xlabel(r'Raman Shift ($cm^{-1}$)',fontsize=12),  ax.set_zlabel(r'Intensity',fontsize=12),  ax.set_ylabel(r'Polarization Angle (°)',fontsize=12)
        ax.set_box_aspect(None, zoom=0.9)
    plt.suptitle(f'r-GeO {plane} - Mode {mode} -  3D Fit - Result ',fontsize=15)
    fig.tight_layout(), plt.subplots_adjust(wspace=0.01, hspace=0.01) # make smaller / remove space between plots
    if normalize == False:
        path = os.path.join('Plots',f'rGeO{plane}_{mode}_3DFit.pdf')
    else:
        path = os.path.join('Plots',f'rGeO{plane}_{mode}_3DFit_normalized.pdf')

    if save_plot != False:
        plt.savefig(path, bbox_inches='tight')
        print('plot saved')
    plt.show()
    pdfs += [path]
    return pdfs

def Plot_ang_dependency(X,Y,intensity_2d,Z_fit,Index,plane,mode,pdfs,normalize,save_plot=False):
    # Index = 585
    #### Plot angular dependency along one energy ####
    angles = Y[:,Index]
    fig, (ax1) = plt.subplots(1,figsize=(8,5),dpi=800)
    ax2 = ax1.twinx()
    ax2.plot(angles,intensity_2d[:,Index ]-Z_fit[:,Index ],marker='h',label='Residual', 
             markerfacecolor='none',ls='-',color='dodgerblue',zorder=1,alpha=0.3)
    ax1.fill_between(angles,Z_fit[:,Index],label='Fit',color='darkseagreen', alpha =0.6,zorder=1,linewidth=2)
    # intensity_2d = [np.sum(intensity_2d[i,Index-3:Index+3]) for i in range(len(intensity_2d[:,Index]))]

    ax1.plot(angles,intensity_2d[:,Index], color='firebrick', label=f'Intensity \nalong {round(X[0,Index],0)}cm^-1',
             ls='-',marker='o', markerfacecolor='none',zorder=2)
    
    plt.xlabel(r'Angle (°)'),  ax1.set_ylabel('Intensity (arb. units)'), 
    # ax1.set_ylim(np.min(intensity_2d[:,Index ])-50)
    ax1.tick_params(which='both',axis='both', direction = 'in',labelsize=8),
    ax1.legend(fontsize=8,loc='upper left'),ax2.legend(fontsize=8,loc='upper right'), plt.title(f'r-GeO {plane} - Mode {mode} - 2D Fit: {round(X[0,Index],0)}cm^-1')
    ax2.set_ylabel('Residual',fontsize=8), ax1.set_xlabel(r'Angle (°)',fontsize=8)
    ax1.set_zorder(1)  ,    ax2.set_zorder(2)  
    ax1.patch.set_visible(False),    ax2.patch.set_visible(False)
    if normalize == False:
        path = os.path.join('Plots',f'rGeO{plane}_{mode}_ang_dependency.pdf')
    else:
        path = os.path.join('Plots',f'rGeO{plane}_{mode}_ang_dependency_normalized.pdf')
    if save_plot != False:
        plt.savefig(path, bbox_inches='tight')
        print('plot saved')
    plt.show()   
    pdfs += [path]
    return pdfs

def plot_single_spectra(X_2d, plane,target_angle,lower_limit,upper_limit,label,Z_fit=False):
    #### Plot single Raman spectra along an angle ###
    fig, (ax1) = plt.subplots(1,figsize=(8,5),dpi=800)
    # plt.title(f'bGaO (010) par: simulated Raman spectra at {angle}°',fontsize=12)
    plt.plot(X_2d[0,:], plane[target_angle,lower_limit:upper_limit],label=f'data {label} ')
    if not Z_fit == False:
        plt.plot(X_2d[0,:], Z_fit[target_angle,lower_limit:upper_limit],label=f'Fit {label} ')

    plt.legend(fontsize=12),plt.tick_params(which='both',axis='both', direction = 'in',labelsize=8)
    plt.xlabel(r'Energy ($cm^{-1}$)'),  plt.ylabel('Intensity (arb. units)')
    ax1.xaxis.set_major_locator(MultipleLocator(50)), ax1.xaxis.set_major_formatter('{x:0.0f}'),ax1.xaxis.set_minor_locator(MultipleLocator(10))
    plt.show()

def PDF_merge(PDFs,name):
    from PyPDF2 import PdfMerger
    merger = PdfMerger()
    for pdf in PDFs:
        merger.append(f'{pdf}')
    merger.write(f"{name}")
    merger.close()
    print('PDFs merged')


# -*- coding: utf-8 -*-
"""
Created on Fri May  3 14:54:08 2019

@author: PAC | pac94@pitt.edu
"""
#based upon: https://rochefort-lab.github.io/fissa/examples/cNMF%20example.html


#%% INPUT PARAMS

######################   ONLY EDIT THIS (after 'r')  #########################
#animalDataPath = r'D:\Data\AA0340'
animalDataPath = '/media/user/sutter2P_data/Data/AA0343'
##############################################################################




#%% IMPORT DEPENDENCIES
import os
from scipy.io import loadmat
import numpy as np
import fissa
import re
import glob


#%% LOAD ROI for each sequence of tifs for experiment save FISSA formatted ROI in .npy file

#load ROI | Must be formatted as ordered points along a curve
#obtained via a matlab script that finds ROI outline points (mask2polygonCoord.m) then orders them along a curve (orderEllipsePtOnCurve.m)
os.chdir(animalDataPath)
anml = re.search(r'[A-Z]{2}\d{4}',animalDataPath)[0]

tiff_folder = os.path.join(animalDataPath,'NoRMCorred')
if not os.path.isdir(tiff_folder):
    os.mkdir(tiff_folder)

output_folder = os.path.join(tiff_folder,'FISSAoutput')
if not os.path.isdir(output_folder):
    os.mkdir(output_folder)

ROIcoor, nTif, rois_FISSA = [[] for _ in range(3)]
#[::-1] reverses list so pre comes before post
for fileN,ROIfile in enumerate(glob.glob(os.path.join(animalDataPath,'*_moCorrROI_*.mat'))[::-1]):
    ROIcoor.append(loadmat(ROIfile)['moCorROI']['ROIcurveOrderedXY'].transpose())
    nTif.append(loadmat(ROIfile)['nTifs'][0][0])
    rois_FISSA.extend([[[ROIcoor[fileN][i,0][0],\
                        ROIcoor[fileN][i,0][1]] for i in range(ROIcoor[fileN].shape[0])]]*nTif[fileN])

#This checks whether ROI files all have the same number of ROI:
#this works for however long ROIcorr is:
if not all([coor.shape[0] is ROIcoor[0].shape[0] for coor in ROIcoor]):
    raise NameError('Different number of ROI in some set')

#save ROIs in python format | can write code to split using nTif variable
#rois_FISSA is a list of ROIs of length number of tifs
np.save(os.path.join(tiff_folder,'FISSA_ROIs.npy'),rois_FISSA)


# %% RUN FISSA

#initiate experiment
exp = fissa.Experiment(tiff_folder,rois_FISSA,output_folder)

#run FISSA
exp.separate(redo_prep=True)

#save FISSA output to matlab format
os.chdir(output_folder)
exp.save_to_matlab()
print('DONE!')
#in 'result', for a given cell and trial there is a n x numTraceFrames double; row 1 is ROI trace, rows 2->n are traces from neuropil regions around ROI
#in 'ROIs', for a given cell and trial there is a n x 1 cell of doubles, the doubles are nPoints x 2, col1=Y,col2=X; 1st cell is ROI,  2->n are neuropil regions around ROI
# %%

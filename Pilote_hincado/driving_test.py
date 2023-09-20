#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  1 14:26:23 2021
@author: SebastianG
"""

from cProfile import label
import os
from pyclbr import Function
from turtle import color
from xml.etree.ElementTree import PI
import math 
import matplotlib.pyplot as plt
from matplotlib.legend_handler import HandlerLine2D
from matplotlib import pyplot
import numpy as np

parameters['allow_extrapolation'] = True

def create_mesh(mesh, cell_type, prune_z=False):
        cells = mesh.get_cells_type(cell_type)
        cell_data = mesh.get_cell_data("gmsh:physical", cell_type)
        out_mesh = meshio.Mesh(points=mesh.points, cells={cell_type: cells}, cell_data={"name_to_read":[cell_data]})
        if prune_z:
            out_mesh.prune_z_0()
        return out_mesh
print('creando malla..')

nombre = 'pilote_refinado_estructurado_axi'


   


#%%

steps =18000

Ti=900#tiempo total
delta= Ti/steps
dt=delta
t=0 # tiempo inicial

I=0
v=0

data = np.array([[I,v,t]])

for pot in range(steps):
    I=(v)*(delta)+I

    if I<4.7:
        v=0.54/60 
    elif I>=4.7 and I<8:
        v=1.65/60
    else:
        v=0.4
    
    if I<9:
        
    else:
        I=9
    np.append()
    
    
    t=t+delta
#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Miscellaneous utility functions

@author: C. J. Stoneking - cjstoneking@gmail.com
"""

import numpy as np

# get_grid_coords
#
# compute the subdivision of a given volume into cells of a given size
# return coordinates of corners of those cells
# 
# input:
#   volume_size - size of volume to subdivide
#   cell_size   - size of cells to divide it in (can be scalar, in which case replicated)
#   cell_offset - offset of cells (in absolute units)
#   center      - if true: have grid of cells centered in volume
#                           the margin which we don't cover will be equal on either side
#               - if false: have no margin on either side , achieve this by appending cells at end
# output:
#   coordinates of corners of the cells 
#   in array where zeroth dimension is volume dimension and first dimension is grid cell

def get_grid_coords(volume_size, cell_size, cell_offset, center):

    cell_size = int(cell_size)
    cell_offset = int(cell_offset)
    
    #if cell size is a scalar and volume size is not,
    #repeat cell size across all dimensions
    if(len(cell_size)==1 and len(volume_size) > 1):
        cell_size = np.tile(cell_size, len(volume_size))
        

    coords = []

    for i in range(len(volume_size)):
        c_i = np.arange(1, volume_size[i], cell_offset)
        c_i = c_i[c_i + cell_size[i] - 1 < volume_size[i]]
        if(center):
            if(len(c_i)>0):
                c_i.append(int((volume_size[i] - (c_i[-1] + cell_size[i] - 1))/2))
        else:
            if((c_i[-1] + cell_size[i]) < volume_size[i]):
                c_i.append( volume_size[i] - cell_size[i] + 1)
        coords.append(c_i)
        
    return np.array(coords)

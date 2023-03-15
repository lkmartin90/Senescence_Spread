import numpy as np
import pandas as pd
import scipy
import scipy.special as sp
from mpl_toolkits import mplot3d
import plotly.express as px
import plotly.graph_objects as go
from scipy.stats import norm
from scipy import spatial
import math as math
import datetime
from circles import *
from diffusion_class import *
from diffusion_functions import *
import multiprocessing as mp
import tqdm 
import parmap
import csv


def one_run(sim_number, r_cell, D_L, density, media_height, size, R_total, k_on):
    """" Tracks one ligand through the media until it binds to a cell
        
    Parameters
    ----------
    sim_numer: float
        number of simulation out of total when many are run
    r_cell: float
        radius of a cell [m]
    D_L: float
        x coordinate of second point
    density: float
        Density of cells in the 2D plane, between 0 and 1
    media_height:
        Height of the layer of medium [m]
    size: int
        Size of the plane of cells under consideration, in number of cell radii
    R_total: float
        Total number of receptors on one cell
    k_on: float
        Ligand receptor binding rate constant, [M^-1 min^-1], where M =  mol/L
            
    Returns
    ------
    list, bool
        list containing the time at which the ligand bound, the distance from (0,0) at which binding happened, and 
    the distance from the point the ligand was emitted that it binds at. For ligands which reach a distance greater
    the specified size, this list is empty.
    The boolean flags whether the trajectory was paracrine or not. If it is a 1 it was paracrine, if it was 0 it 
    was autocrine.
    """
    # set seed
    np.random.seed(sim_number)
    
    # Initialise simulation
    signal_sim = Particle_track(r_cell, D_L, density, media_height, size, R_total, k_on)
    
    count = 0 
    
    while count == 0:

        # Update ligand location
        status = signal_sim.update()

        # if status is True particle has been trapped
        if status == True:

            # distance of absorbtion of ligand from (0,0)
            trap_distance = ((signal_sim.x_loc[-1]**2 + signal_sim.y_loc[-1]**2)**0.5)/r_cell
            
            # distance of absorbition of ligand from that point at which it was released
            trap_distance_point = (((signal_sim.x_loc[-1] - signal_sim.x_loc[0])**2 \
                                    + (signal_sim.y_loc[-1] - signal_sim.y_loc[0])**2)**0.5)/r_cell

            
            count = 1
            if trap_distance < 1:
                # autocrine trajectory
                para_flag = 0

            else:
                # paracrine trajectory
                para_flag = 1
            
            return [signal_sim.time[-1], trap_distance, trap_distance_point], para_flag
            
        elif signal_sim.x_loc[-1]**2 + signal_sim.y_loc[-1]**2 > (size * r_cell/2.)**2:
            # particle has left the region in the simulation
            count = 1
                
            para_flag = 0
            
            trap_distance = ((signal_sim.x_loc[-1]**2 + signal_sim.y_loc[-1]**2)**0.5)/r_cell
            #print(trap_distance)
                
            return [], para_flag
                
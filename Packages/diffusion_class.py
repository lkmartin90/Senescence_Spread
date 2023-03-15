import numpy as np
import scipy
import scipy.special as sp
from scipy.stats import norm
from scipy import spatial
import math as math
from timeit import default_timer as timer
import datetime
from circles import *

class Particle_track():
    """This class contains the functions to update the position of a ligand emitted by a cell as it moves through a medium
    
    ...
    
    Attributes 
    ----------
    x_loc: list
        The list of x positions occupied by the ligand over time
    y_loc: list
        The list of y positions occupied by the ligand over time
    z_loc: list
        The list of z positions occupied by the ligand over time
    timestamp: list
        List containing the time taken for each jump in ligand location
    time: list
        List contaning the time of each jump in ligand location assuming the ligand was emitted at time 0
    r_cell: float
        Radius of the cell [m]
    density: float
        Density of cells in the 2D plane, between 0 and 1
    D_L: float
        Diffusion constant for the ligand, [m^2/s]
    h: float
        Height of the layer of medium
    boundary: float
        Height of the boundary layer in which a different propogation algorithm is used
    size: int
        Size of the plane of cells under consideration, in number of cell radii (one side of the square)
    cells: list
        Cells in the 2D plane
    R_total: float
        Total number of receptors on one cell
    k_on: float
        Ligand receptor binding rate constant, [M^-1 min^-1], where M =  mol/L
    kappa: float
        Rate constant
        
        
    Methods
    -------
    
    euclidean_distance
        
    populate_plane
    
    radius_of_sphere
        
    update_loc_boundary
    
    update_loc_med
    
    update_time_boundary
        
    update_time_boundary
    
    update_time_med
        
    lateral_advance_boundary
    
    boundary_test
        
    trap_test
        
    nearest_trap
        
    update
    
    """
    
    def __init__(self, r_cell, D_L, density, media_height, size, R_total, k_on):
        """
        Parameters
        ----------
        r_cell: float
            Radius of the cell [m]
        density: float
            Density of cells in the 2D plane, between 0 and 1
        D_L: float
            Diffusion constant for the ligand, [m^2/s]
        media_height: float
            Height of the layer of medium
        size: int
            Size of the plane of cells under consideration, in number of cell radii
        R_total: float
            Total number of receptors on one cell
        k_on: float
            Ligand receptor binding rate constant, [M^-1 min^-1], where M =  mol/L
        """

        # want the random source of the ligand to be located somewhere on the first cell with a uniform distribution, do this by generating a random r and theta
        r = np.random.rand()*r_cell
        theta = np.random.rand()* np.pi * 2
        self.x_loc = [r * np.sin(theta)]
        self.y_loc = [r * np.cos(theta)]
        self.z_loc = [0.0]
        self.timestamp = []
        self.time = [0.0]

        self.r_cell = r_cell
        self.density = density
        self.D_L = D_L
        self.h = media_height
        self.boundary = 0.001*self.r_cell # Optimal distance for algorithm speed
        self.d = ((np.pi*self.r_cell**2)/self.density)**0.5
        self.size = size
        self.cells = []
        self.R_total = R_total
        self.k_on = k_on
        self.kappa = (self.k_on*self.R_total)/(np.pi*self.r_cell**2*6.0221409e+23)
        
        self.populate_plane()
        #print("populated with cells")
    
    def euclidean_distance(self, x1, y1, x2, y2):
        """ Returns the Euclidean distance between two points
        
        Parameters
        ----------
        x1: float
            x coordinate of first point
        y1: float
            y coordinate of first point
        x2: float
            x coordinate of second point
        y1: float
            y coordinate of second point
            
        Returns
        -------
        float
            Euclidean distance between two points
        """

        return math.hypot((x1 - x2), (y1 - y2))
    
        
    def populate_plane(self):
        """Populates a 2D plane with non overlapping cells with the same radius.
        
        Used the density of cells to create the correct number of cells in the 2D plane (z = 0). Uses the Circles class, imported from circles.py. The coordinates of the centre of these cells are added to the list self.cells and a histogram of cell locations is assigned to self.boxes.
        """
        # create instance of the Cirles class
        circle = Circles()
        # fill a plane of size self.size and density self.density
        circle.fill(self.size, self.density)
        # fill self.cells with the coordinates of each cell
        self.cells = np.array(list(circle))*self.r_cell
        # store the dictionary containing the binned cells for later use
        self.boxes = circle.boxes
                
        
    def radius_of_sphere(self):
        """ Calculates the radius of the sphere on which the particle moves to when it is outside of the boundary layer
        
        Takes the last known z position of the ligand (the final value in the list self.z_loc) and calculates the radius that this particle can move to
        
        Returns
        -------
        float
            The distance in meters that the ligand can move
        """
        R = np.min([self.z_loc[-1], self.h - self.z_loc[-1]])
        if R > self.boundary:
            return R
        else:
            return self.boundary
        
        
    def update_loc_boundary(self):
        """ Updates the location of a ligand when it is in the boundary layer
        
        Checks whether the ligand is over a cell or an empty region of dish. If empty region of dish then the ligand is reflected from the surface, if cell, the ligand is absorbed with some probability
        
        Returns
        -------
        Boolean
            True when particle is trapped and False when it is not
        """
        
        t = self.timestamp[-1]
        z_0 = self.z_loc[-1]

        # Update x and y positions
        self.lateral_advance_boundary()
        # Update vertical (depends on whether the particle is above a trap or not - check if this is the case)
        trap_stat = self.trap_check()
        
        sd = (2*self.D_L*t)**0.5

                
        if trap_stat:
            value_check = 0
            # This method involes z being updated directly, instead of adding the value to the previous z
            while value_check == 0:
                
                # Sample a new zvalue from the distribution that you expect to see at a reflecting boundary
                z = np.random.normal(z_0, sd) + np.random.normal(-1*z_0, sd)
                
                # Only want positive z values, as the boundary is in the way...
                if z > 0.0:
                    # Calculate the probability that for this new z value the ligand is absorbed and not reflected.
                
                    p_a = norm.pdf(z, z_0, sd) + norm.pdf(z, -1*z_0, sd) \
                    - (self.kappa/self.D_L)*np.exp(-((z+z_0)**2)/(4*self.D_L*t)) \
                    * sp.erfcx((z+z_0+2*self.kappa*t)/((4 * self.D_L *t)**0.5))

                    reflect = (norm.pdf(z, z_0, sd) + norm.pdf(z, -1*z_0, sd))
                    
                    test = np.random.rand() * reflect

                    if test > p_a:
                        # The ligand has been absorbed by the trap
                        #print("absorbed")
                        value_check = 2
                        self.z_loc.append(0.0)
                    else:
                        # Particle is reflected and arrives at the new z position
                        value_check = 1
                        self.z_loc.append(z)
                

            if value_check == 2:
                # ligand absorbed
                return True
            else:
                # ligand not absorbed
                return False
                    
        else:
            # As this is the reflecting boundary, we want only values of z greater than 0
            value_check = 0
            while value_check == 0:
                z = np.random.normal(z_0, sd) + np.random.normal(-1*z_0, sd)
                if z > 0:
                    self.z_loc.append(z)
                    value_check = 1
            # ligand not absorbed
            return False
        
    def update_loc_med(self):
        """Updates the the location when the ligand is not in the boundary layer
        
        This is done by picking coordinates on a sphere with a uniform probabilty, with the radius of this sphere calculated by the 'radius_of_sphere' method, depending on the ligands location. The lists self.x_loc, self.y_loc, and self.z_loc are updated, with the new coordinates added.
        """
        
        r = self.radius_of_sphere()
        theta_1 = np.random.rand() * np.pi * 2
        theta_2 = np.random.rand() * np.pi * 2
        
        self.x_loc.append(self.x_loc[-1] + r*np.cos(theta_1)*np.sin(theta_2))
        self.y_loc.append(self.y_loc[-1] + r*np.sin(theta_1)*np.sin(theta_2))
        trial_z = self.z_loc[-1] + r*np.cos(theta_2)
        if trial_z < self.h:
            self.z_loc.append(trial_z)
        else:
            #print("reflected ", self.h - (trial_z - self.h))
            self.z_loc.append(self.h - (trial_z - self.h))
            
    
        
    def update_time_boundary(self):
        """ Updatse the timestamp when the ligand is in the boundary layer.
        
        The updated position depends on whether over trap or not. The lists self.timestamp and self.time are updated with the time taken for the particle to move to the next location, and the current total simulation time, respectively
        """
        trap_stat = self.trap_check()
        if trap_stat:
            # Over trap (0.01*r_cell)**2/ 2*D_L
            t = (0.01*self.r_cell)**2/(2*self.D_L)
        else:
            # Over reflective part of surface (0.1*d_nearest)**2/ 2*D_l
            d_near = self.nearest_trap()
            t = (0.1*d_near)**2/(2*self.D_L)
            
        self.timestamp.append(t)
        self.time.append(self.time[-1] + t)
        
    def update_time_med(self):
        """Updates the timestamp when the location is not in the boundary.
        
        The lists self.timestamp and self.time are updated with the time taken for the particle to move to the next location, and the current total simulation time, respectively.
        """
        t = self.radius_of_sphere()**2/(6*self.D_L)
        self.timestamp.append(t)
        self.time.append(self.time[-1] + t)
    
    def lateral_advance_boundary(self):
        """ Upadates the ligands x and y coordinates when it is in the boundary layer.
       
        Updates the lists self.x_loc and self.y_loc
        """
        # (4 pi DL delta_t)^-0.5 * exp(-delta_x**2/ 4 DL delta_t)
        t = self.timestamp[-1]
        x = np.random.normal(0, (2*self.D_L*t)**0.5)
        y = np.random.normal(0, (2*self.D_L*t)**0.5)
       
        self.x_loc.append(self.x_loc[-1] + x)
        self.y_loc.append(self.y_loc[-1] + y)
        
    def boundary_test(self):
        """ Check if particle is in the boundary or not
        
        Returns
        -------
        Boolean
            True is the particle is in the boundary layer, False if it is not
        """
        if len(self.z_loc) > 0:
            if self.z_loc[-1] < self.boundary:
                return True
            else:
                return False
        else:
            return True
        
    def trap_check(self):
        """Check if particle is over a cell or not, using the function nearest_trap
        
        Returns
        -------
        Boolean
            True is the particle is over a cell, False if it is not
        """
        dist_to_closest = self.nearest_trap()

        if dist_to_closest > self.r_cell:
            return False
        else:
            return True
        
    def nearest_trap(self):
        """Finds the nearest cell to a ligand for update_time_boundary

        Uses the binned locations of particles stored in self.boxes, only bins around current ligand location must be searched.
        
        Returns
        -------
        float
            Distance to the centre of closest cell
        """
        x, y = (self.x_loc[-1]/ self.r_cell, self.y_loc[-1]/self.r_cell)
        i, j = int(x), int(y)
        closest_array = []
        
        for di in steps:
            for dj in steps:
                c = self.boxes.get((i+di, j+dj))
                if c:
                    px, py = c
                    closest_array.append(self.euclidean_distance(x, y, px, py))
                   
        # backup method if no cells are in the imediate vicinity of the ligand
        if len(closest_array) ==0:
            cell_list = self.cells
            tree = spatial.KDTree(cell_list)
            closest = tree.query([(self.x_loc[-1], self.y_loc[-1])])
            return closest[0][0]

        return min(closest_array)* self.r_cell
            
        
    def update(self):
        """ Update the ligand coordinates and time
        
        The timestep for the move is first calculated and reccorded, and this value is then used to updte the ligand coordinates
        
        Returns
        -------
        Boolean
            True if ligand has bound to a cell, False if not
        """
        
        boundary_stat = self.boundary_test()
        if boundary_stat:
            self.update_time_boundary()
            trapped = self.update_loc_boundary()
        else:
            self.update_time_med()
            self.update_loc_med()
            trapped = False
            
        return trapped
            

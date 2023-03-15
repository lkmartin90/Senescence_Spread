import numpy as np
from itertools import *

# The number of points either side of the point of interest to search for overlapping cells
steps = (-2, -1, 0, 1, 2)

class Circles:
    """ A class to populate a plane with non overlapping circles of the same radius.
    
    ...
    
    Attributes
    ----------
    boxes: dictionary
        dictionary, where the key describes the integer bin in which the coordinate falls and the value is the
        coordinate of the point
    count: int
        keeps track of the number of cells created by the code
    
    
    Methods
    -------
    add:
        Adds an individual cell to the population of cells
    fill:
        Fills the plane with the correct number of cells given the density
    
    
    """
    
    def __init__(self):
        """
        Constructs parameters for the circles object
        """
        self.boxes = {}
        self.count = 0
        
    def fill(self, size, density):
        """
        Fills a 2D plane with a number of non overlapping circles to a given density.
        
        Parameters
        ----------
        size: int
            an integer specifying the length of one size of the square plane to be filled with cells
            (in units of cell radii)
        density: float
            the density of the plane which cells should occupy
        """
        # calculate the number of cells required
        n = int(size*size*density/np.pi)
        # always require central cell at (0,0) to make simulation easier later on
        self.add((0, 0))
        # while the number of cells is smaller than the required number, generate a new random coordinate on the plane
        while len(self) < n:
            x = np.random.rand() * size - size/2.
            y = np.random.rand() * size - size/2.
            self.add((x, y))

    def fill_lattice(self, size, density):
        """
         Fills a 2D plane with a number of non overlapping circles to a given density, on a lattice.

         Parameters
         ----------
         size: int
             an integer specifying the length of one size of the square plane to be filled with cells
             (in units of cell radii)
         density: float
             the density of the plane which cells should occupy
         """
        # calculate the required lattice spacing
        spacing = (np.pi/density)**0.5
        size_array = np.unique(np.concatenate((-1*np.arange(0, np.ceil(size/2), spacing), np.arange(0, np.ceil(size/2), spacing))), 0)
        for coord in permutations(size_array, 2):
            x = coord[0]
            y = coord[1]
            self.add((x, y))
        for i in size_array:
            self.add((i, i))
        
    def add(self, point):
        """
        Adds an individual point to the population of cells on the plane stored in the boxes dictionary if it does not
        overlap with other, pre-existing cells.
        
        Parameters
        ----------
        point: tuple
            coordinates of the point to be added to the plane
        """
        x, y = point
        # return the integer above which the coordinate values lies above
        i, j = int(x), int(y)
        # steps through nearest bins looking for overlap
        for di in steps:
            for dj in steps:
                c = self.boxes.get((i+di, j+dj))
                # if there is a cell in the bin then check the distance from this cell
                if c:
                    px, py = c
                    # if the newly generated coordinate is too close, do not add the cell
                    if (px-x)**2 + (py-y)**2 < 4:
                        return
        # otherwise add the cell
        self.boxes[(i, j)] = point
        self.count += 1
            
    def __len__(self):
        # allows length of class to be called, returning the number of cells
        return self.count
    
    def __repr__(self):
        # special method to represent class' objects as a string, so that calling the class returns some
        # useful information about it
        return f'Circles(count={len(self)})'
    
    def __iter__(self):
        # returns iterable version of the coordinates of each cell
        return iter(self.boxes.values())

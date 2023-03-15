from functools import lru_cache
import scipy.stats as stats
from circles import *
import random

class Senescence_grid():
    """This class contains the location of senescent cells and the strength of the SASP
    
    ...
    Attributes
    ----------
    x_range: int
        x axis ranges from 0 to this value
    y_range: int
        y axis ranges from 0 to this value
    r_cell: float
        cell radius in meters
    R_tot: int
        total number of receptors on each cell
    k_on: float
        forward binding constant,[m^3]/[mol][s]
    density: float
        density of cells in the 2D plane
    k_e: float
        rate of endocytosis min^-1
    k_off: float
        rate of dissociation min^-1
    D_L: float
        Diffuction constant m^2/s
    N_E: int
        number of ligands emitted per unit time
    N_D: int
        number of ligands absorbed per unit time to become senesent
    k: float
        kappa, rate constant
    Da: float
        Damkohler number, dimensionless    
    nu: float
        probability of internalisation, dimensionless
    keff: float
        effective binding rate [m]/[s]
    grid: list
        list of cell locations
    N_A: float
        avagadros number [mol^-1]
    scells: dict
        dictionary containing senescent cells, key is tuple of cell coordinate, value is the time at which the cell
        became senescent
    jux_scells: dict
        dictionary containing juxtacrine secondary senescent cells. key is a tuple of cell coordinates and value is
        the time at which the cell became senescent
    paracrine: dict
        dictionary containing paracrine secondary senescent cells. key is a tuple of cell coordinates and value is
        the time at which the cell became senescent
    sim_list: list
        a list containing the cells on the grid which have yet to become senescent
    probs: list
        a list containing a probability for each cell in sim_list to become senescent
    Tstruct: dict
        Structure to store the cells which are in the process of becoming senescent. Using a dict as in this version
        of python it will be ordered. The key is the delay until the cell becomes fully senescent. The value returned
        is a tuple, the first item in the tuple is a list of coordinates of cells which become senescent at this time.
        The second entry is the type of cell, 0 seeded senescence, 1 for paracrine senescent and 2 for Juxtacrine.
    cells_which_have_waited: list
        Used in the delay simulations to keep track of cells which have started the proces of becoming senescent.
        
    Methods
    -------
    cumulative_func_int
    
    cumulative func
    
    seed_senesence
    
    update paracrine
    
    update juxt
    
    update_senescent_list
    
    assign_prob
    
    prop_cacl
    
    F_calc
    
    time_calc
    
    cell_calc
    
    update_juxt_cells
    
    update_juxt_cells_delay
    
    one_run
    
    """
    
    def __init__(self, size, r_cell, R_tot, k_on, density, k_e, k_off, D_L, N_E, N_D,
                 delay, j_frac, j_prob, p_frac, lattice):
        """ 
        Paramters
        ---------
        size: int
            length of the side of the grid of cells in units of cell radii
        r_cell: float
            cell radius in meters
        R_tot: int
            total number of receptors on each cell
        k_on: float
            forward binding constant,[m^3]/[mol][s]
        density: float
            density of cells in the 2D plane
        k_e: float
            rate of endocytosis min^-1
        k_off: float
            rate of dissociation min^-1
        D_L: float
            Diffusion constant m^2/s
        N_E: int
            number of ligands emitted per unit time
        N_D: int
            number of ligands absorbed per unit time to become senescent
        delay: float
            number of secondss taken for a cell to become fully senescent
        j_frac: float
            fraction of the SASP produced by the original senescent cells produced by the juxtacrine senescent cells
        p_frac: float
            fraction of the SASP produced by the original senescent cells produced by the paracrine senescent cells
        j_prob: float
            probability oc a cell inducing juxtacrine senescence in its neighbours
        """

        self.size = size
        self.r_cell = r_cell
        self.N_A = 6.0221409*10**23
        self.R_tot = R_tot
        self.k_on = k_on
        self.sigma = density
        self.k_e = k_e
        self.k_off = k_off
        self.D_L = D_L
        self.N_E = N_E
        self.N_D = N_D
        self.delay = delay
        self.sim_time = 0
        self.j_prob = j_prob
        self.lattice = lattice
        self.cells_that_have_waited = []
        
        self.juxt_SASP_fraction = j_frac
        self.para_SASP_fraction = p_frac
        
        k = (self.k_on*self.R_tot)/(np.pi * self.r_cell**2 * self.N_A)
        Da =  k*self.r_cell/self.D_L
        
        # if values are 0 then there is no dissociation or unbinding of any sort
        try:
            nu = self.k_e/(self.k_off + self.k_e)
        except:
            nu = 1
        
        self.k = k
        self.Da = Da
        self.nu = nu
        
        k_eff = self.k*self.sigma/(1 + (np.pi*self.Da/4))
        
        self.k_eff = k_eff
        
        # create cells on 2D plane
        self.populate_plane()
        
        # create structure to store the cells which are in the process of becoming senescent.
        self.Tstruct = {}
        
        return

    def populate_plane(self):
        """Populates a 2D plane with non overlapping cells with the same radius.
        
        Used the density of cells to create the correct number of cells in the 2D plane (z = 0). Uses
        the Circles class, imported from circles.py. The coordinates of the centre of these cells are
        added to the list self.grid and a histogram of cell locations is assigned to self.boxes.
        """
        # create instance of the Cirles class
        circle = Circles()
        # fill a plane of size self.size and density self.density
        if self.lattice:
            circle.fill_lattice(self.size, self.sigma)
        else:
            circle.fill(self.size, self.sigma)
        # fill self.cells with the coordinates of each cell
        self.grid = [(x*self.r_cell, y*self.r_cell) for x, y in list(circle)]
        # store the dictionary containing the binned cells for later use
        self.boxes = circle.boxes
        
        return

    # If running not as a Jupyter notebook want to use @lru_cache() around the following 2 functions
    @lru_cache()
    def cumulative_func(self, r):
        """cumulative probability distribution for the probability of binding a ligand at distance r 
        
        Parameters
        ----------
        r: float
            distance from cell which emmits ligand [m]
        
        Returns
        -------
        out: float
            Cumulative probability for a ligand binding before this distance
        """
        # subtract r_cell from distance as the equation is true for distance outside the emitting cell
        r = r - self.r_cell
        # probability of autocrine binding
        P_au = self.nu*self.Da/(self.nu*self.Da + 4/np.pi)
        out = r/(r + (1.1 * self.D_L / (self.nu * self.k_eff)))
        return (1 - P_au)*out

    def seed_senescence(self, seed_r):
        """Seeds the senescent cells at the start of the simulation, input is an array of 2d coordinates. 
        Senescent cells are entered into the dict self.scells and any surrounding juxtacrine cells are entered 
        into self.juc_cells
        
        Parameters
        ----------
        seed_r: float
            radius which cells are senescent inside
        """
        
        self.scells = {}
        self.jux_scells = {}
        self.paracrine = {}
        for coord in self.grid:
            x = coord[0]
            y = coord[1]
            if x**2 + y**2 < seed_r**2:
                # create a dictionary containing senescent cells
                # no matter the delay this will be the first thing to happen so don't have
                # to put these values into Tstruct
                self.scells.update({(x, y): self.delay})

        # create a dictionary with juxtocrine senescent cells
        if self.delay == 0:
            for coord in self.scells:
                self.update_juxt_cells(coord[0], coord[1], 0)
        else:
            self.update_juxt_cells_delay(list(self.scells.keys()), self.delay)
            
        self.sim_time = self.delay
        return

    def seed_senescence_random(self, seed_dens):
        """Seeds the senescent cells at the start of the simulation, input is an array of 2d coordinates.
        Senescent cells are entered into the dict self.scells and any surrounding juxtacrine cells are entered
        into self.juc_cells

        Parameters
        ----------
        seed_dens: float
            density of randomly seeded senescent cells, where a density of 1 would mean that all cells are senescent
        """

        self.scells = {}
        self.jux_scells = {}
        self.paracrine = {}
        no_of_sen_cells = int(seed_dens*len(self.grid))
        for coord in random.sample(self.grid, no_of_sen_cells):
            x = coord[0]
            y = coord[1]

            self.scells.update({(x, y): self.delay})

        # create a dictionary with juxtocrine senescent cells
        if self.delay == 0:
            for coord in self.scells:
                self.update_juxt_cells(coord[0], coord[1], 0)
        else:
            self.update_juxt_cells_delay(list(self.scells.keys()), self.delay)

        self.sim_time = self.delay
        return
            
    def update_paracrine(self, new_cell, t):
        """Update the list of paracrine cells with a newly paracrine senescent cell
        
        Parameters
        ----------
        new_cell: tuple
            tuple of the coordinates of the new paracrine senescent cell
        t: float
            time that the cell becomes senescent
        """

        if new_cell not in self.jux_scells and new_cell not in self.paracrine and new_cell is not None:
            self.paracrine.update({(new_cell[0], new_cell[1]): [t]})
            
        return


    def update_senescent_list(self):
        """subtract paracrine and juxtocrine cells from the grid of cells that can become senescent
        """

        # find the coordinates of the cells that are in the process of becoming senescent
        struct_cells = []
        for key in self.Tstruct.keys():
            struct_cells = struct_cells + self.Tstruct[key][0]
            
        diff1 = list(set(self.grid) - set(self.scells))
        diff2 = list(set(diff1) - set(self.jux_scells))
        diff3 = list(set(diff2) - set(self.paracrine))
        diff4 = list(set(diff3) - set(struct_cells))
        self.sim_list = sorted(diff4)
        
        return

    def assign_prob(self):
        """assign each grid point a probability of becoming senescent
        """

        # create array containing 0s for each coordinate in the grid - this will be filled by the rest of the function
        #self.probs = [0 for ix, iy in it.product(range(self.x_range),range(self.y_range))]
        self.probs = [0 for i in self.sim_list]

        # cycle through each cell which could become senescent
        for i, point in enumerate(self.sim_list):
            # array that will contain the n*p from the binomial, from each already senescent cell
            # n = number of ligands emitted by senescent cell, p = probability of binding
            # np < 5 for Poisson approximation to hold.
            prob_array = []
            for cell in self.scells:
                # calculate the distance from the source
                x_dist = cell[0] - point[0]
                y_dist = cell[1] - point[1]
                r_dist = (x_dist**2 + y_dist**2)**0.5
                
                # to get to the probability want to integrate the prob dist at a given distance with the limits
                # dist+r and dist-r. Then need to multiply by pi*r_cell**2/(pi*((r_dist+r_cell)**2 - (r_dist-r_cell)**2))
                c = r_dist-self.r_cell
                b = r_dist+self.r_cell

                probability = (self.cumulative_func(b)-self.cumulative_func(c))*(self.r_cell**2/((r_dist+self.r_cell)**2 - (r_dist-self.r_cell)**2))
  
                # N_E = number of ligands emitted in a time period 
                # N_D = the number of ligands that must be detected in a time period to get a response 
                # binom_prob = the probability of a cell absorbing N_D or more things in a given time frame.
                
                # Expectation of the number of ligands bound from this source cell
                # binomial n*p

                prob_array.append(probability*self.N_E)
                
            if self.juxt_SASP_fraction != 0:
                for cell in self.jux_scells:
                    # calculate the distance from the source
                    x_dist = cell[0] - point[0]
                    y_dist = cell[1] - point[1]
                    r_dist = (x_dist**2 + y_dist**2)**0.5

                    # to get to the probability want to integrate the prob dist at a given distance with the limits
                    # dist+r and dist-r. Then need to multiply by pi*r_cell**2/(pi*((r_dist+r_cell)**2 - (r_dist-r_cell)**2))
                    c = r_dist-self.r_cell
                    b = r_dist+self.r_cell

                    probability = (self.cumulative_func(b)-self.cumulative_func(c))*(self.r_cell**2/((r_dist+self.r_cell)**2 - (r_dist-self.r_cell)**2))

                    # N_E = number of things emitted in a time period 
                    # N_D = the number of things that must be detected in a time period to get a response 
                    # binom_prob = the probability of a cell absorbing N_D or more things in a given time frame.

                    # Expectation of the number of ligands bound from this source cell
                    prob_array.append(probability*self.N_E*self.juxt_SASP_fraction)
                    
            if self.para_SASP_fraction != 0:
                for cell in self.paracrine:
                    # calculate the distance from the source
                    x_dist = cell[0] - point[0]
                    y_dist = cell[1] - point[1]
                    r_dist = (x_dist**2 + y_dist**2)**0.5

                    # to get to the probability want to integrate the prob dist at a given distance with the limits
                    # dist+r and dist-r. Then need to multiply by
                    # pi*r_cell**2/(pi*((r_dist+r_cell)**2 - (r_dist-r_cell)**2))
                    c = r_dist-self.r_cell
                    b = r_dist+self.r_cell

                    probability = (self.cumulative_func(b)-self.cumulative_func(c))*(self.r_cell**2/((r_dist+self.r_cell)**2 - (r_dist-self.r_cell)**2))

                    # N_E = number of things emitted in a time period 
                    # N_D = the number of things that must be detected in a time period to get a response 
                    # binom_prob = the probability of a cell absorbing N_D or more things in a given time frame.

                    # Expectation of the number of ligands bound from this source cell
                    prob_array.append(probability*self.N_E*self.para_SASP_fraction)
                    
            # expectation of the number of ligands bound from all source cells within one unit of time

            lamda = np.sum(prob_array)
            # probability of binding more than the number of ligands required for senescence in one unit of time
            #pois_prob = 1 - stats.binom.cdf(self.N_D, self.N_E, lamda/self.N_E)
            pois_prob = 1 - stats.poisson.cdf(self.N_D, lamda)
            self.probs[i] = pois_prob

            ################################################################################
            ### attempting a normal approx to poisson
            ################################################################################
            # if max(prob_array) < 5:
            #     lamda = np.sum(prob_array)
            #     # probability of binding more than the number of ligands required for senescence in one unit of time
            #     pois_prob = 1 - stats.poisson.cdf(self.N_D, lamda)
            #     self.probs[i] = pois_prob
            # else:
            #     mu = np.sum([el for el in prob_array if el > 5])
            #     dist = stats.norm(loc=mu, scale=mu**0.5)
            #     norm_prob = 1 - dist.cdf(self.N_D)
            #     self.probs[i] = norm_prob
            
        return
        
    def prop_calc(self):
        """ Calculate the sum of propensities at a given time
        """
        
        a0 = np.sum(self.probs)
        return a0
    
    
    def time_calc(self):
        """ Calculate the time at which something next happens to the cells
        
        assume that the diffusion time is basically instantaneous - steady state, won't hold if things change
        end up working in terms of the N_E and N_D, the number emitted in a given time, at the moment that's
        a pretty arbitrary time, might as well call it a second
        """
        
        draw = np.random.random()
        rate_sum = np.sum(self.probs)
        if rate_sum == 0:
            time = np.inf
        else:
            time = (1/rate_sum)*np.log(1/draw)
        return time
    
    
    def F_calc(self, at):
        """ Function to calculate F, the cumulative distribution function of time
        
        Parameters
        ----------
        at: float
            a_0(t)*T_1, where a_0 is the sum of propensities at time t and T_1 is the time that the first reaction
             completes, here the time that a cell becomes fully senescent
        
        Returns
        -------
        f: float
            projection of probability of event happening into uniform interval for comparison to random number generated
        """
        f = (1 - np.exp(-at))
        return f

    def cell_calc(self):
        """Determine which cell becomes senescent
        
        Returns
        -------
        sen_cell: tuple
            tuple containing the coordinates of the cell which becomes senescent
        """
        draw = np.random.random()
        # cumulatively add the probabilities and normalise so that they sum to 1. Because some of the values 
        # so small there are some issues with precision here, there are several values that are just "1.0". 
        # However, I don't think this will affect the simulation too much - as the likelihood of choosing any of
        # these values is vanishingly small. 
        cumulative = np.cumsum(self.probs)

        total = np.sum(self.probs)
        normalised = cumulative/total # This should be a list of probabilities between 0 and 1
        cell_no = np.where(normalised > draw)
        try:
            sen_cell = self.sim_list[int(cell_no[0][0])]

        except IndexError:
            print(self.sim_time)
            print("Probabilities too low, no senescence spread at all")
            sen_cell = None
        
        return sen_cell

    def update_juxt_cells(self, x, y, time):
        """ Update the dictionary of Juxtacrine cells. A coordinate of a cell capable of inducing juxtacrine senescence
         in the cells around it is passed to this function and the number of juxtacrine cells updated correctly.
        
        Parameters
        ----------
        x: int
            x coordinate of cell capable of inducing juxtacrine senescent
        y: int
            y coordinate of cell capable of inducing juxtacrine senescent
        time: float
            time at which surrounding cells will become juxtacrine senescent
        """
        # juxtocrine cells surround primary senescent cells, at the moment the nearest sites have a j_prob chance
        # of becoming senescent
                    
        # return the integer above which the coordinate values lies above
        i, j = int(x/self.r_cell), int(y/self.r_cell)
        region = [-3, -2, -1, 0, 1, 2, 3]
        # steps through nearest bins looking for cells
        for di in steps:
            for dj in region:
                c = self.boxes.get((i+di, j+dj))
                # if there is a cell in the bin then check the distance from this cell
                if c:
                    ix, iy = c
                    ix = ix * self.r_cell
                    iy = iy * self.r_cell
                    draw = np.random.random()
                    # can change the if draw < 1.0 parameter to vary the contact Juxtacrine sensence condition
                    if draw < self.j_prob and (ix, iy) not in self.scells and (ix, iy) not in self.jux_scells  and (ix, iy) not in self.paracrine:
                        self.jux_scells[(ix, iy)] = time
                        
        return
                        

    def update_juxt_cells_delay(self, cell_list, update_time):
        """ Update the dictionary of Tstruct with Juxtacrine cells when there is a time delay. A coordinate of a
        cell capable of inducing juxtacrine senescence in the cells around it is passed to this function and the
         number of juxtacrine cells updated correctly.
        
        Parameters
        ----------
        cell_list: list
            a list of cells capable of inducing senescence in the cells around them
        update_time: float
            time at which surrounding cells will become juxtacrine senescent - not simply self.delay as is relative
            to current point in the simulation
        """

        # juxtocrine cells surround primary senescent cells, at the moment the nearest sites have a j_prob chance
        # of becoming senescent
                    
        coords = []            
        # return the integer above which the coordinate values lies above
        for cell in cell_list:
            x = cell[0]
            y = cell[1]
            i, j = int(x/self.r_cell), int(y/self.r_cell)
            region = [-3, -2, -1, 0, 1, 2, 3]
           
            # steps through nearest bins looking for cells
            for di in steps:
                for dj in region:
                    c = self.boxes.get((i+di, j+dj))
                    # if there is a cell in the bin then check the distance from this cell
                    if c:
                        ix, iy = c
                        ix = ix *self.r_cell
                        iy = iy *self.r_cell
                        draw = np.random.random()
                        # can change the j_prob parameter to vary the contact Juxtacrine sensence condition
                        if draw < self.j_prob and (ix, iy) not in self.scells and (ix, iy) not in self.jux_scells \
                                and (ix, iy) not in self.paracrine and (ix, iy) not in self.cells_that_have_waited:
                            coords.append((ix, iy))
                            self.cells_that_have_waited.append((ix, iy))
        if len(coords) > 0:
            self.Tstruct.update({update_time: (coords, 2)})

        return

    def one_run(self, verbose):
        """ One run through the simulation process, the time at which the next cell becomes senescent is calculated and 
        the cell is found. The cells which can become juxtacrine senescent due to this cell are updated.
        
        Parameters
        ----------
        verbose: bool
            If true the simulation prints parameters as it progresses, used for debugging
            
        Returns
        -------
        sim_time: float
            simulation time in hours of the current senescent event occurring
        r: float
            The distance of the newly senescent cell from the origin in units of meters.
            
        """
        
        # update the list of cells which have not yet become senescent
        self.update_senescent_list()
        # calculate the probability of each cell becoming senescent
        self.assign_prob()
        # work out the time at which the next cell becomes senescent.
        # this is complicated as in this time other cells may finish the process
        # of becoming senescent and start to emit SASP
        if self.delay != 0:
            # calculate propensities at this time point
            time_change = 0

            additional_time_to_event = self.time_calc()

            if verbose:
                print("new time step")
                print("TStruct", self.Tstruct)

            if len(self.Tstruct) == 0:
                # there are no cells in the process of becoming senescent so can use usual Gillespie
                if verbose:
                    print("No cells waiting")
            else:
                i = 0
                # T1 = the time at which the first cell in the list would finish becoming senescent
                T1 = list(self.Tstruct)[0]
                if verbose:
                    print("T1", T1)
                    print("additional time to event", additional_time_to_event)

                # start by calculating the time and then update the relevant structures after.
                while additional_time_to_event > T1:
                    # a reaction completed and a cell became senescent before anything else happened

                    # calculate a_t, this involves re-calculating the probabilities including the SASP from
                    # the new cell that has become fully senescent. To do this need position of new cell and the
                    # type of senescence.

                    delayed_reaction_time = list(self.Tstruct)[i]  # time after initial_time cell becomes senescent
                    time_change = delayed_reaction_time

                    if i > 0:
                        self.sim_time = self.sim_time + (list(self.Tstruct)[i] - list(self.Tstruct)[i-1])
                    else:
                        self.sim_time = self.sim_time + delayed_reaction_time

                    coords = self.Tstruct[delayed_reaction_time][0]
                    cell_type = self.Tstruct[delayed_reaction_time][1]

                    # based on the cell type cell is either paracrine or juxtacrine senescent - different consequences
                    # for each of these outcomes
                    if cell_type == 1:
                        # cell is paracrine senescent, add to paracrine cell list
                        for co in coords:
                            self.paracrine.update({co: self.sim_time})
                            if verbose:
                                print("Added paracrine")
                         # deleted inclusion of juxtacrine cells here as want to create juxtacrine cells as soon as cells
                         # are "infected" with senescence this is the main change from "Senescence_grid_delay.py"

                    elif cell_type == 2:
                        # cell is juxtacrine senescent, add to juxtacrine cell list
                        for co in coords:
                            if co not in self.jux_scells and co not in self.paracrine:
                                self.jux_scells.update({(co[0], co[1]): self.sim_time})
                                if verbose:
                                    print("Added juxtacrine")

                    # update the list of cells which can become senescent
                    self.update_senescent_list()
                    # calculate propensities again
                    self.assign_prob()

                    # generate a new waiting time

                    additional_time_to_event = self.time_calc()
                    if verbose:
                        print("Tstruct before 2nd calculation", self.Tstruct)

                    # find the time until the next waiting event, if there is not a next waiting event = np.inf
                    try:
                        T1 = (list(self.Tstruct)[i+1] - list(self.Tstruct)[i])
                    except:
                        T1 = np.inf

                    if verbose:
                        print("additional time to event", additional_time_to_event)
                        print("T1", T1)
                    i = i+1

                # calculate the time for the new reaction
                if verbose:
                    print("time change", time_change)

            # add additional time to sim time
            self.sim_time = self.sim_time + additional_time_to_event

            # alter the dictionary containing the cells waiting to become senescent so that all the times are correct
            # do this at the end of the above loop to save computation time
            self.Tstruct = {key - time_change - additional_time_to_event: val for key, val in self.Tstruct.items() if key > time_change}

            if verbose:
                print("additional time to event", additional_time_to_event)

            new_coord = self.cell_calc()
            
            if new_coord is None:
                self.sim_time = np.nan
            
            # add new reaction to Tstruct
            self.Tstruct.update({self.delay: ([new_coord], 1)})
            # contains all cells that are in the process of becoming senescent, or have become senescent so that they
            # easily be used elsewhere
            self.cells_that_have_waited.append(new_coord)

            ###############################################################################
            # add the juxtacrine cells created due to this to struct
            # first need to calculate which cells will be juxtacrine

            self.update_juxt_cells_delay([new_coord], self.delay + 0.01)

            ##############################################################################

            if verbose:
                print("addded paracrine to Tstruct")
                print(self.Tstruct)

            if new_coord is None:
                r = np.nan
            else:
                r = (new_coord[0]**2 + new_coord[1]**2)**0.5

        else:
            # no delay

            # update the list of cells which have not yet become senescent
            self.update_senescent_list()

            # calculate the probability of each cell becoming senescent
            self.assign_prob()
            # calculate the time at which the next cell becomes senescent
            t = self.time_calc()
            self.sim_time = self.sim_time+t

            # determine which cell becomes senescent
            coord = self.cell_calc()
            # add this cell to the list of paracrine senescent cells
            self.update_paracrine(coord, self.sim_time)

            if coord is None:
                r = np.nan
            else:
                r = (coord[0] ** 2 + coord[1] ** 2) ** 0.5
                self.update_juxt_cells(coord[0], coord[1], self.sim_time)

        return self.sim_time, r

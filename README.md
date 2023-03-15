# Senescence_Spread
Code to simulate the spread of senescence from a lesion of primary senescent cells


## Simulation of delayed senescence spread 

This is a Gillespie simulation with a time delay. We have the probability that a cell at a given distance from a senescent cell which produces a SASP will beome senescent in a given amount of time. The simulation is initialised on a grid of a given size, 100X100 is commonly used as larger grids will increase the computational time. Cells are created so that they populate this plane with the specified density.

The time delay in this simulation is due to the time taken for a cell to become fully senescent after the start of senescence induction. This time can be specified in the "parameters" section of the simulation. Only once a cell has become fully senescent can it produce a SASP.

Primary senescent cells are populated in the centre of the grid at time 0. The simulation proceeds by keeping track of the cells on the grid which are senescent, are in the process of becoming senescent. A cell becomes paracrine senescent when it binds to more than  $𝑁_𝐷$ ligands within one hour.

From the postition of the cell in relation to the cells on the grid which emit a SASP, the probability of a cell becoming senescent in an hour can be calculated.

This probability can then be used to give the "rate" of senescence induction for that cell. If the probability is 0.1 then the rate of senescence induction for a cell in that position is 0.1 per hour. Each cell in the grid which is not yet senescent has a rate of senescence induction. By summing these rates the time at which the next cell becomes paracrine secondary senescent can be calculated. In the usual method used in Gillespie simulations, the event which actually happens is calculated, here this corresponds to calculating which cell becomes senescent.

The situation is complicated by the time delay, it may be that cells have finished becoming fully senescent before the next cell becomes paracrine senescent, in which case the SASP that they produce must also be accounted for. We deal with this by recording the cells which have started to become senescent but are not yet fully senescent and so cannot pass on this senescence to surrounding cells (in a python dictionary).

If there are no cells waiting to become senescent then the simulation progresses as a normal Gillespie simulation would. However, if a cell is waiting to become fully senescent then the simulation must decide if this happens before or after senescent is passed on to a new cell. If the time for the next senescence induction event is after a cell would finish becoming senescent, the simulation progresses to the time at which a cell becomes fully senescent, adds the contribution of that cells SASP to the simulation, and proceeds from there. 

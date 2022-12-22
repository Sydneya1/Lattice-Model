'''
Author: Sydney Ackermann
Code for fragmentation modes on a lattice 
The purpose of this code is to simulate competition between fragmentation modes 
and to see which mode dominates the lattice at the end of the simulation
given the trade-off between propagule size and the potential for vertical transmission of cancer
...

Description of the simulation: This simulation contains: cancer; 
birth and death rates that are size-independant; group sizes of up to 6 cells
'''

# import statements
from random import random, seed, randrange
from numpy import arange, zeros, ones, random, linspace, exp
from matplotlib.pyplot import imshow, colorbar
from pylab import plot, show, legend, subplot, xlabel, ylabel, tight_layout, figure, xlim, ylim, title, grid, pause, draw, legend
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
from math import exp

# -----------------------------------------------------------------------------|
# PARAMETERS
# -----------------------------------------------------------------------------|

# -----------------------------------------------------------------------------|
# PARAMETERS
# -----------------------------------------------------------------------------|
Num_of_gen = 1000  # number of iterations in the simulation
L = 25  # dimension of grid
N0 = 500  # initial population size
N = 6  # maximum group size // group size at fragmentation
dt = 1  # more accurate as this approaches 0
B = 0.3

# a matrix that keeps tracks of the locations of groups that have reached max size (l) and are about to fragment
Locations_of_parents = zeros((L, L), int)

# death rates depending on group size
alpha_d = 1
g2 = zeros(N, float)
deaths = zeros(N, float)

for i in range(N):
    deaths[i] = 0.03

# probability of death depending on group size
dps1 = 1 - exp(-deaths[0] * dt)
dps2 = 1 - exp(-deaths[1] * dt)
dps3 = 1 - exp(-deaths[2] * dt)
dps4 = 1 - exp(-deaths[3] * dt)
dps5 = 1 - exp(-deaths[4] * dt)
dps6 = 1 - exp(-deaths[5] * dt)

# probability of mutating
probm = 0.1  # 0.1

phi = 1

# size dependence A=1, size independence A=0
A = 0

# replicative advantage of cancer cells
R = 2

global = 1 # set this equal to 1 for global dispersal, otherwise the simulation will use local dispersal

# -----------------------------------------------------------------------------|
# define a function for flipping a coin
def coin_toss():
    num = randrange(1, 3)  # will return either 1 or 2
    return num

# -----------------------------------------------------------------------------|
# define a function for rolling a n-1 sided die
def roll_dice(n):
    '''
    This function takes in a number n which dictates how many sides the die has 
    and returns the outcome of rolling an n-1 sided die
    values can range from 1,2,3,...,n-1
    '''

    num = randrange(1, n)  # will return either 1,2,3,4,5, ..., n-1
    return num

# -----------------------------------------------------------------------------|
# DEFINE MATRICES
# -----------------------------------------------------------------------------|
# create N0 randomly scattered individuals on an L by L grid
Grid = zeros((L, L), int)  # 1 represents a single cell, 0 represents no individual
grid_vals = []

random_x = randrange(1, L + 1)
random_y = randrange(1, L + 1)

grid_vals.append([random_x, random_y])

while len(grid_vals) < N0 + 1:

    new = [random_x, random_y]

    if new in grid_vals:
        random_x = randrange(1, L + 1)
        random_y = randrange(1, L + 1)

    else:
        Grid[random_x - 1][random_y - 1] = 1
        grid_vals.append(new)

# plot
fig0 = figure(0)
for i in range(len(Grid)):
    for j in range(len(Grid)):
        if Grid[i][j] == 1:
            plot(i, j, marker='o', color='gray', ms=6)
ax2 = fig0.add_subplot(1, 1, 1)
ax2.set_yticks(linspace(0, L - 1, L))
ax2.set_xticks(linspace(0, L - 1, L))
ax2.grid(True)
title('grid')
show()

# -----------------------------------------------------------------------------|
# create a matrix that records the max group size of each individual
max_size = zeros((L, L), int)  # max group size can be 2,3,4,5 or 6

# -----------------------------------------------------------------------------|
# fragmentation genotype for points on the lattice that are inhabited by cells
Mode = zeros((L, L), int) 

for i in range(L):
    for j in range(L):

        if Grid[i][j] == 1:  # if there is an individual at that spot on the lattice
            num = roll_dice(24)
            Mode[i][j] = num
            if num == 1:
                max_size[i][j] = 2
            elif 2 <= num <= 3:
                max_size[i][j] = 3
            elif 4 <= num <= 7:
                max_size[i][j] = 4
            elif 8 <= num <= 13:
                max_size[i][j] = 5
            elif num >= 14:
                max_size[i][j] = 6

made_colors = ['red', 'orange', 'blue', 'brown', 'pink', 'green', 'purple', 'peru', 'cyan', 'lime', 'teal', 'palegreen', 'magenta', 'rosybrown', 'gold', 'tan', 'navy', 'olive', 'coral', 'darkseagreen', 'bisque', 'lavender', 'darkviolet']
marker_sizes = [0, 6, 8, 10, 12, 14]

# plot
fig1 = figure(1)
ax2 = fig1.add_subplot(1, 1, 1)

for i in range(len(Grid)):
    for j in range(len(Grid)):

        ax2.plot(i, j, marker='o', color=made_colors[Mode[i][j] - 1], ms=marker_sizes[Grid[i][j]])

ax2.set_yticks(linspace(0, L - 1, L))
ax2.set_xticks(linspace(0, L - 1, L))
ax2.grid(True)
a = mlines.Line2D([], [], color='red', marker='o', markersize=6, label='1+1')
b = mlines.Line2D([], [], color='orange', marker='o', markersize=6, label='1+1+1')
c = mlines.Line2D([], [], color='blue', marker='o', markersize=6, label='2+1')
d = mlines.Line2D([], [], color='brown', marker='o', markersize=6, label='1+1+1+1')
e = mlines.Line2D([], [], color='pink', marker='o', markersize=6, label='2+1+1')
f = mlines.Line2D([], [], color='green', marker='o', markersize=6, label='2+2')
g = mlines.Line2D([], [], color='purple', marker='o', markersize=6, label='3+1')

h = mlines.Line2D([], [], color='peru', marker='o', markersize=6, label='1+1+1+1+1')
i = mlines.Line2D([], [], color='cyan', marker='o', markersize=6, label='2+1+1+1')
j = mlines.Line2D([], [], color='lime', marker='o', markersize=6, label='2+2+1')
k = mlines.Line2D([], [], color='teal', marker='o', markersize=6, label='3+1+1')
l = mlines.Line2D([], [], color='palegreen', marker='o', markersize=6, label='3+2')
m = mlines.Line2D([], [], color='magenta', marker='o', markersize=6, label='4+1')

n = mlines.Line2D([], [], color='rosybrown', marker='o', markersize=6, label='1+1+1+1+1+1')
o = mlines.Line2D([], [], color='gold', marker='o', markersize=6, label='2+1+1+1+1')
pp = mlines.Line2D([], [], color='tan', marker='o', markersize=6, label='2+2+1+1')
q = mlines.Line2D([], [], color='navy', marker='o', markersize=6, label='2+2+2')
r = mlines.Line2D([], [], color='olive', marker='o', markersize=6, label='3+1+1+1')
s = mlines.Line2D([], [], color='coral', marker='o', markersize=6, label='3+2+1')
t = mlines.Line2D([], [], color='darkseagreen', marker='o', markersize=6, label='3+3')
u = mlines.Line2D([], [], color='bisque', marker='o', markersize=6, label='4+1+1')
v = mlines.Line2D([], [], color='lavender', marker='o', markersize=6, label='4+2')
w = mlines.Line2D([], [], color='darkviolet', marker='o', markersize=6, label='5+1')

legend(handles=[a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, pp, q, r, s, t, u, v, w], bbox_to_anchor=(1, 1), loc=2, borderaxespad=0.2)
title('fragmentation mode')
show()

# -----------------------------------------------------------------------------|
# FUNCTIONS
# -----------------------------------------------------------------------------|
def vt(num_of_cancer_cells, off):
    '''
    This function takes in the number of cancer cells in the parent 
    (num_of_cancer_cells) and a list of offspring sizes (off)
    and returns a list of the amount of cancer each offspring is to inherit
    '''

    if num_of_cancer_cells == 0:  # if there are no cancer cells, the inheritance distribution will just be [0,0,0]
        spread = []
        for i in range(len(off)):
            spread.append(0)

    else:
        spread = []  # this list is a distribution of the number of cancer cells that is inherited by each offspring
        for i in range(len(off)):  # for however many offspring there are ex: 3
            num = randrange(0, num_of_cancer_cells + 1)  # roll the die 3 times, it can have a value between 0 and the number of cells in a group (l)
            spread.append(num)
        total_cells = sum(spread)
        if total_cells != num_of_cancer_cells:  # check that the sum of the values of the 3 rolls adds up to the amount of cancer in the parent
            total_cells = 0
        for k in range(len(off)):  # check that the number of cancer cells each offspring inherits does not overshoot it's size
            if off[k] < spread[k]:
                total_cells = 0

        while total_cells != num_of_cancer_cells:  # untill we reach the number of cancer cells in parent

            # find three numbers that add to three and don't exceed offspring size
            spread = []
            for i in range(len(off)):
                upper = num_of_cancer_cells + 1
                num = randrange(0, upper)  # 0,1,2,3,4
                spread.append(num)
            total_cells = sum(spread)

            if total_cells != num_of_cancer_cells:  # check that the sum of the values of the 3 rolls adds up to the amount of cancer in the parent
                total_cells = 0

            # check that the number of cancer cells doesn't exceed the number of offspring cells, otherwise, keep trying
            for k in range(len(off)):
                if off[k] < spread[k]:
                    total_cells = 0
                    
    return spread

# -----------------------------------------------------------------------------|
def place_offspring(num_of_cancer_cells, offspring, i, j, inheritance_distrib, offs):
    '''
    This function takes in the number of cancer cells in the parent, the size of
    the offspring, the coordinates of the spot on the lattice, a list that 
    determines the number of cancer cells each offspring inherits and the index 
    of the offspring we are currently considering in the list of offspring.
    The function then randomly assigns the offspring to a spot on the lattice 
    given by (ival, jval) and then calls the winner function to see if it gets 
    to stay there.
    '''

    if global == 1: # global dispersal
        ival = randrange(0, L)  # 0,1,2, ...,L-1
        jval = randrange(0, L)  # 0,1,2, ...,L-1 # never roll a zero or an L so can't fall off the grid

        winner(num_of_cancer_cells, offspring, ival, jval, inheritance_distrib, offs)
    
    else: # local dispersal
        num = randrange(1, 6)  # will return either 1,2,3,4,5, ..., n-1

        if num == 1 and i < L - 1:  # place to the right of parent
            ival = i + 1
            jval = j
            winner(num_of_cancer_cells, offspring, ival, jval, inheritance_distrib, offs)
    
        elif num == 2 and i > 0:  # place to the left of parent
            ival = i - 1
            jval = j
            winner(num_of_cancer_cells, offspring, ival, jval, inheritance_distrib, offs)
    
        elif num == 3 and j < L - 1:  # place above parent
            ival = i
            jval = j + 1
            winner(num_of_cancer_cells, offspring, ival, jval, inheritance_distrib, offs)
    
        elif num == 4 and j > 0:  # place below parent
            ival = i
            jval = j - 1
            winner(num_of_cancer_cells, offspring, ival, jval, inheritance_distrib, offs)
    
        elif num == 5:  # place on the parent spot
            ival = i
            jval = j
            winner(num_of_cancer_cells, offspring, ival, jval, inheritance_distrib, offs)
    
# -----------------------------------------------------------------------------|
def winner(num_of_cancer_cells, offspring, a, b, inheritance_distrib, offs):
    '''
    This function takes in the number of cancer cells in the parent, the size of
    the offspring, the coordinates of the spot on the lattice, a list that 
    determines the number of cancer cells each offspring inherits and the index 
    of the offspring we are currently considering in the list of offspring.
    The function then determines wether or not the offspring wins the spot. If 
    it does win the spot (depending on it's size) then the matrix that records 
    size, the matrix that records fragmentation mode, and the dictionary that 
    records the number of cancer cells at each spot on the lattice all get updated.
    '''
    
    if Grid[a][b] > offspring:
        pass
        #print(Grid[a][b], ' > ', offspring, 'so offspring dies')

    elif Grid[a][b] == offspring:  # if occupied by another the same size then flip a coin to see who stays
        numberr = coin_toss()
        
        if numberr == 1:
            #print(Grid[a][b], ' == ', offspring, 'offspring wins')
            Grid[a][b] = offspring
            Mode[a][b] = gene
            spot_to_num_of_cancer[(a, b)] = inheritance_distrib[offs]
            max_size[a][b] = fragsize

    elif Grid[a][b] < offspring:
        #print(Grid[a][b], '<', offspring, 'so offspring fills an empty spot')
        Grid[a][b] = offspring
        Mode[a][b] = gene
        spot_to_num_of_cancer[(a, b)] = inheritance_distrib[offs]
        max_size[a][b] = fragsize

# -----------------------------------------------------------------------------|
def fragment(off, spot_to_num_of_cancer):
    '''
    This function takes the offspring list; the cancer dictionary
    and obtains the number of cancer cells; calculates the inheritance distribution;
    place's the offspring; and finally, updates the Grid, Mode, and the cancer dictionary
    '''
    
    num_of_cancer_cells = spot_to_num_of_cancer[(i, j)]
    
    # randomly disrtibute cancer cells among offspring
    inheritance_distrib = vt(num_of_cancer_cells, off)  # returns for ex: [1,0,0] meaning first offspring inherits one cancer cell and the second and third offspring don't
    
    for offs in range(len(off)):
        offspring = off[offs]
        place_offspring(num_of_cancer_cells, offspring, i, j, inheritance_distrib, offs)
    
    # update parent's spot
    if Grid[i][j] == 0:
        Mode[i][j] = 0
        spot_to_num_of_cancer[(i, j)] = 0
        max_size[i][j] = 0
        
# -----------------------------------------------------------------------------|
def birth(group_size, numcancer, i, j):  # make sure this works and the number of cancer cells is compatible
    '''
    This function takes in the group size and the number of cancer cells
    
    and calculates the number of cooperators cells; 
    chooses the number of cancer cells that divides and updates the grid;
    chooses the number of cooperator cells that divides and updates the grid;
    chooses the number of cooperator cells that mutate
    calculates the new number of cooperators by subtracting those that mutated and adding those that divided
    calculates the new number of cancer cells by adding the new mutations and those that divided
    then updates the cancer dictionary
    and finally makes sure that the group size (and number of cancer cells) doesn't exceed l
    '''

    g1 = zeros(N, float)
    for k in range(1, N):  # works
        g1[k] = ((k - 1) / (N - 2))

    # calc number of cooperator cells
    numcoop = group_size - numcancer

    # cancer  cell division
    probb = R * ((B) + (A * g1[numcoop - 1]))
    probc = 1 - exp(-numcancer * probb)
    y = random.binomial(1, probc, size=None)
    Grid[i][j] = Grid[i][j] + y

    group_size = group_size + y

    # make sure the number of cells doesn't exceed the max group size
    if max_size[i][j] == 4:
        if Grid[i][j] > 4:
            Grid[i][j] = 4
            group_size = 4
    elif max_size[i][j] == 3:
        if Grid[i][j] > 3:
            Grid[i][j] = 3
            group_size = 3
    elif max_size[i][j] == 2:
        if Grid[i][j] > 2:
            Grid[i][j] = 2
            group_size = 2
    elif max_size[i][j] == 5:
        if Grid[i][j] > 5:
            Grid[i][j] = 5
            group_size = 5
    elif max_size[i][j] == 6:
        if Grid[i][j] > 6:
            Grid[i][j] = 6
            group_size = 6
    


    group_size = Grid[i][j]
    numcancer = numcancer + y
    numcoop = group_size - numcancer

  # cooperator (SIZE-DEPENDENT) cell division
    probb = ((B) + (A * g1[numcoop - 1]))
    probc = 1 - exp(-numcoop * probb)
    x = random.binomial(1, probc, size=None)
    Grid[i][j] += x

    if max_size[i][j] == 4:
        if Grid[i][j] > 4:
            Grid[i][j] = 4
    elif max_size[i][j] == 3:
        if Grid[i][j] > 3:
            Grid[i][j] = 3
    elif max_size[i][j] == 2:
        if Grid[i][j] > 2:
            Grid[i][j] = 2
    elif max_size[i][j] == 5:
        if Grid[i][j] > 5:
            Grid[i][j] = 5
    elif max_size[i][j] == 6:
        if Grid[i][j] > 6:
            Grid[i][j] = 6

    # mutate
    num_of_mutations = random.binomial(x, probm, size=None)

    # new number of cooperator cells
    numcoop = numcoop - num_of_mutations + x

    # update cancer dictionary
    spot_to_num_of_cancer[(i, j)] = numcancer + num_of_mutations

    # but make sure it doesn't exceed l
    if max_size[i][j] == 4:
        if (numcancer + num_of_mutations) > 4:
            spot_to_num_of_cancer[(i, j)] = 4
    elif max_size[i][j] == 3:
        if (numcancer + num_of_mutations) > 3:
            spot_to_num_of_cancer[(i, j)] = 3
    elif max_size[i][j] == 2:
        if (numcancer + num_of_mutations) > 2:
            spot_to_num_of_cancer[(i, j)] = 2
    elif max_size[i][j] == 5:
        if (numcancer + num_of_mutations) > 5:
            spot_to_num_of_cancer[(i, j)] = 5
    elif max_size[i][j] == 6:
        if (numcancer + num_of_mutations) > 6:
            spot_to_num_of_cancer[(i, j)] = 6


# -----------------------------------------------------------------------------|
def death(proportion_of_cancer):
    '''
    This function takes the proportion of cancer in the group
    and kills groups depending on this number
    '''
    
    if Grid[i][j] == 0:  # spot not occupied
        pass
    
    elif proportion_of_cancer == 0:  # no cancer
        pass

    else:
        random_number = random.random_sample()  # pick a random number between 0 and 1

        death_rate = proportion_of_cancer**phi  # nonlinear cancer cost for phi != 0
        death_probability = 1 - exp(-death_rate * dt)  # convert to continuous-time approximation

        if random_number < death_probability:  # if the proportion of cancer is greater than that random number then group dies
            Grid[i][j], Mode[i][j], spot_to_num_of_cancer[(i, j)], max_size[i][j] = 0, 0, 0, 0       

# -----------------------------------------------------------------------------|
def death_size(new_size):
    '''
    This function takes in the size of the group
    and kills groups depending on their size
    '''

    if Grid[i][j] == 0:  # spot not occupied
        pass

    else:
        dpss = [dps1, dps2, dps3, dps4, dps5, dps6]
        if random.random_sample() < dpss[new_size - 1]:
            Grid[i][j], Mode[i][j], spot_to_num_of_cancer[(i, j)], max_size[i][j] = 0, 0, 0, 0

# -----------------------------------------------------------------------------|
# SET UP ITERATION
# -----------------------------------------------------------------------------|

# create dictionary that records the number of cancer cells at each spot on the lattice (key)
spot_to_num_of_cancer = {}
for i in range(L):
    for j in range(L):
        spot_to_num_of_cancer[(i, j)] = 0

# lists to keep track of abuncance of each mode for plotting 

mode1 = []
mode2 = []
mode3 = []
mode4 = []
mode5 = []
mode6 = []
mode7 = []
mode8 = []
mode9 = []
mode10 = []
mode11 = []
mode12 = []
mode13 = []
mode14 = []
mode15 = []
mode16 = []
mode17 = []
mode18 = []
mode19 = []
mode20 = []
mode21 = []
mode22 = []
mode23 = []

size1 = []
size2 = []
size3 = []
size4 = []
size5 = []
size6 = []

cancer5 = []

# -----------------------------------------------------------------------------|
# ITERATE
# -----------------------------------------------------------------------------|
for p in range(Num_of_gen):


    # add cells each iteration ------------------------------------------------|
    for i in range(L):
        for j in range(L):

            if Grid[i][j] != 0:  # if there is an individual (/group) at that spot

                # calculate group size
                group_size = Grid[i][j]

                # calculate number of cancer cells
                numcancer = spot_to_num_of_cancer[(i, j)]

                # call function
                birth(group_size, numcancer, i, j)

                # death based on group size
                new_size = Grid[i][j]
                death_size(new_size)
                
            if Grid[i][j] != 0:
                # death based on proportion of cancer
                proportion_of_cancer = spot_to_num_of_cancer[(i, j)] / Grid[i][j]
                death(proportion_of_cancer)

                # kill adults
                if max_size[i][j] == 4:
                    if Grid[i][j] == 4:  # if the group has reached max size
                        Locations_of_parents[i][j] = 1
                        Grid[i][j] = 0  # kill all fragmenting groups: set size equal to zero but keep fragmentation mode information
                        
                elif max_size[i][j] == 3:
                    if Grid[i][j] == 3:
                        Locations_of_parents[i][j] = 1
                        Grid[i][j] = 0
    
                elif max_size[i][j] == 2:
                    if Grid[i][j] == 2:
                        Locations_of_parents[i][j] = 1
                        Grid[i][j] = 0

                elif max_size[i][j] == 5:
                    if Grid[i][j] == 5:
                        Locations_of_parents[i][j] = 1
                        Grid[i][j] = 0
                        
                elif max_size[i][j] == 6:
                    if Grid[i][j] == 6:
                        Locations_of_parents[i][j] = 1
                        Grid[i][j] = 0

    # iterate through the grid and make groups fragment -----------------------|
    for i in range(L):
        for j in range(L):

            # if there is a parent at that spot
            if Locations_of_parents[i][j] == 1:

                    # obtain genotype
                    gene = Mode[i][j]

                    # obtain size at fragmentation
                    fragsize = max_size[i][j]
                    
                    # fragment
                    offs = [[1,1], [1,1,1], [2, 1], [1,1,1,1], [2, 1, 1], [2, 2], [3, 1], [1, 1, 1, 1, 1], [2, 1, 1, 1], [2, 2, 1], [3,1,1], [3, 2], [4, 1], [1, 1, 1, 1, 1, 1], [2, 1, 1, 1, 1], [2, 2, 1, 1], [2, 2, 2], [3, 1, 1, 1], [3, 2, 1], [3, 3], [4, 1, 1], [4, 2], [5, 1]]
                    fragment(offs[gene - 1], spot_to_num_of_cancer)
                    Locations_of_parents[i][j] = 0

    # at each timestep, update the abundances ---------------------------------|
    counter1 = 0
    counter2 = 0
    counter3 = 0
    counter4 = 0
    counter5 = 0
    counter6 = 0
    counter7 = 0
    counter8 = 0
    counter9 = 0
    counter10 = 0
    counter11 = 0
    counter12 = 0
    counter13 = 0
    counter14 = 0
    counter15 = 0
    counter16 = 0
    counter17 = 0
    counter18 = 0
    counter19 = 0
    counter20 = 0
    counter21 = 0
    counter22 = 0
    counter23 = 0
    
    s1 = 0
    s2 = 0
    s3 = 0
    s4 = 0
    s5 = 0
    s6 = 0

    countercancer = 0

    for i in range(L):
        for j in range(L):

            # modes
            if Mode[i][j] == 0:
                pass
            
            elif Mode[i][j] == 1:
                counter1 += 1
                
            elif Mode[i][j] == 2:
                counter2 += 1
                
            elif Mode[i][j] == 3:
                counter3 += 1
                
            elif Mode[i][j] == 4:
                counter4 += 1

            elif Mode[i][j] == 5:
                counter5 += 1

            elif Mode[i][j] == 6:
                counter6 += 1

            elif Mode[i][j] == 7:
                counter7 += 1
                
            elif Mode[i][j] == 8:
                counter8 += 1
                
            elif Mode[i][j] == 9:
                counter9 += 1
                
            elif Mode[i][j] == 10:
                counter10 += 1
                
            elif Mode[i][j] == 11:
                counter11 += 1

            elif Mode[i][j] == 12:
                counter12 += 1
                
            elif Mode[i][j] == 13:
                counter13 += 1

            elif Mode[i][j] == 14:
                counter14 += 1

            elif Mode[i][j] == 15:
                counter15 += 1

            elif Mode[i][j] == 16:
                counter16 += 1

            elif Mode[i][j] == 17:
                counter17 += 1
                
            elif Mode[i][j] == 18:
                counter18 += 1

            elif Mode[i][j] == 19:
                counter19 += 1

            elif Mode[i][j] == 20:
                counter20 += 1

            elif Mode[i][j] == 21:
                counter21 += 1
                
            elif Mode[i][j] == 22:
                counter22 += 1

            elif Mode[i][j] == 23:
                counter23 += 1            

            # cancer
            if spot_to_num_of_cancer[(i, j)] != 0:  # counting the number of individuals with cancer
                countercancer += 1

            # size
            if Grid[i][j] == 1:
                s1 += 1
            elif Grid[i][j] == 2:
                s2 += 1
            elif Grid[i][j] == 3:
                s3 += 1
            elif Grid[i][j] == 4:
                s4 += 1
            elif Grid[i][j] == 5:
                s5 += 5
            elif Grid[i][j] == 6:
                s6 += 6            

    mode1.append(counter1)
    mode2.append(counter2)
    mode3.append(counter3)
    mode4.append(counter4)
    mode5.append(counter5)
    mode6.append(counter6)
    mode7.append(counter7)
    mode8.append(counter8)
    mode9.append(counter9)
    mode10.append(counter10)
    mode11.append(counter11)
    mode12.append(counter12)
    mode13.append(counter13)
    mode14.append(counter14)
    mode15.append(counter15)
    mode16.append(counter16)
    mode17.append(counter17)
    mode18.append(counter18)
    mode19.append(counter19)
    mode20.append(counter20)
    mode21.append(counter21)
    mode22.append(counter22)
    mode23.append(counter23)
    
    cancer5.append(countercancer)

    size1.append(s1)
    size2.append(s2)
    size3.append(s3)
    size4.append(s4)
    size5.append(s5)
    size6.append(s6)    
    
# -----------------------------------------------------------------------------|

# plot
fig3 = figure(3)
ax4 = fig3.add_subplot(1, 1, 1)
ax4.set_yticks(linspace(0, L - 1, L))
ax4.set_xticks(linspace(0, L - 1, L))
ax4.grid(True)
title('fragmentation mode')

for i in range(len(Grid)):
    for j in range(len(Grid)):

        ax4.plot(i, j, marker='o', color=made_colors[Mode[i][j] - 1], ms=marker_sizes[Grid[i][j]])

        # mutation -----------------------------|
        cancer = spot_to_num_of_cancer[(i, j)]
        if cancer == 0:
            pass
        elif cancer == 1:
            plot(i, j, marker='o', color='black', ms=4)
        elif cancer == 2:
            plot(i, j, marker='o', color='black', ms=6)
        elif cancer == 3:
            plot(i, j, marker='o', color='black', ms=8)
        elif cancer == 4:
            plot(i, j, marker='o', color='black', ms=10)
        elif cancer == 5:
            plot(i, j, marker='o', color='black', ms=12)
        elif cancer == 6:
            plot(i, j, marker='o', color='black', ms=14)       

a = mlines.Line2D([], [], color='red', marker='o', markersize=6, label='1+1')
b = mlines.Line2D([], [], color='orange', marker='o', markersize=6, label='1+1+1')
c = mlines.Line2D([], [], color='blue', marker='o', markersize=6, label='2+1')
d = mlines.Line2D([], [], color='brown', marker='o', markersize=6, label='1+1+1+1')
e = mlines.Line2D([], [], color='pink', marker='o', markersize=6, label='2+1+1')
f = mlines.Line2D([], [], color='green', marker='o', markersize=6, label='2+2')
g = mlines.Line2D([], [], color='purple', marker='o', markersize=6, label='3+1')

h = mlines.Line2D([], [], color='peru', marker='o', markersize=6, label='1+1+1+1+1')
i = mlines.Line2D([], [], color='cyan', marker='o', markersize=6, label='2+1+1+1')
j = mlines.Line2D([], [], color='lime', marker='o', markersize=6, label='2+2+1')
k = mlines.Line2D([], [], color='teal', marker='o', markersize=6, label='3+1+1')
l = mlines.Line2D([], [], color='palegreen', marker='o', markersize=6, label='3+2')
m = mlines.Line2D([], [], color='magenta', marker='o', markersize=6, label='4+1')

n = mlines.Line2D([], [], color='rosybrown', marker='o', markersize=6, label='1+1+1+1+1+1')
o = mlines.Line2D([], [], color='gold', marker='o', markersize=6, label='2+1+1+1+1')
pp = mlines.Line2D([], [], color='tan', marker='o', markersize=6, label='2+2+1+1')
q = mlines.Line2D([], [], color='navy', marker='o', markersize=6, label='2+2+2')
r = mlines.Line2D([], [], color='olive', marker='o', markersize=6, label='3+1+1+1')
s = mlines.Line2D([], [], color='coral', marker='o', markersize=6, label='3+2+1')
t = mlines.Line2D([], [], color='darkseagreen', marker='o', markersize=6, label='3+3')
u = mlines.Line2D([], [], color='bisque', marker='o', markersize=6, label='4+1+1')
v = mlines.Line2D([], [], color='lavender', marker='o', markersize=6, label='4+2')
w = mlines.Line2D([], [], color='darkviolet', marker='o', markersize=6, label='5+1')
z = mlines.Line2D([], [], color='black', marker='o', markersize=6, label='cancer cells')

legend1 = legend(handles=[a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, pp, q, r, s, t, u, v, w, z], bbox_to_anchor=(1, 1), loc=2, borderaxespad=0.2)
ax4 = plt.gca().add_artist(legend1)
aa = mlines.Line2D([], [], color='grey', marker='o', markersize=6, label='1 cell')
bb = mlines.Line2D([], [], color='grey', marker='o', markersize=8, label='2 cell')
cc = mlines.Line2D([], [], color='grey', marker='o', markersize=10, label='3 cell')
dd = mlines.Line2D([], [], color='grey', marker='o', markersize=12, label='4 cell')
legend(handles=[aa, bb, cc, dd], bbox_to_anchor=(1, 0.2), loc=2, borderaxespad=0.2)

show()

# plot time series of abundance of each mode ----------------------------------|
figure(4)
x_axis = zeros(p + 1)
plot(mode1, color='red')
plot(mode2, color='orange')
plot(mode3, color='blue')
plot(mode4, color='brown')
plot(mode5, color='pink')
plot(mode6, color='green')
plot(mode7, color='purple')

plot(cancer5, color='black')

plot(mode8, color='peru')
plot(mode9, color='cyan')
plot(mode10, color='lime')
plot(mode11, color='teal')
plot(mode12, color='palegreen')
plot(mode13, color='magenta')

plot(mode14, color='rosybrown')
plot(mode15, color='gold')
plot(mode16, color='tan')
plot(mode17, color='navy')
plot(mode18, color='olive')
plot(mode19, color='coral')
plot(mode20, color='darkseagreen')
plot(mode21, color='bisque')
plot(mode22, color='lavender')
plot(mode23, color='darkviolet')

a = mlines.Line2D([], [], color='red', marker='o', markersize=6, label='1+1')
b = mlines.Line2D([], [], color='orange', marker='o', markersize=6, label='1+1+1')
c = mlines.Line2D([], [], color='blue', marker='o', markersize=6, label='2+1')
d = mlines.Line2D([], [], color='brown', marker='o', markersize=6, label='1+1+1+1')
e = mlines.Line2D([], [], color='pink', marker='o', markersize=6, label='2+1+1')
f = mlines.Line2D([], [], color='green', marker='o', markersize=6, label='2+2')
g = mlines.Line2D([], [], color='purple', marker='o', markersize=6, label='3+1')

hh = mlines.Line2D([], [], color='peru', marker='o', markersize=6, label='1+1+1+1+1')
i = mlines.Line2D([], [], color='cyan', marker='o', markersize=6, label='2+1+1+1')
j = mlines.Line2D([], [], color='lime', marker='o', markersize=6, label='2+2+1')
k = mlines.Line2D([], [], color='teal', marker='o', markersize=6, label='3+1+1')
l = mlines.Line2D([], [], color='palegreen', marker='o', markersize=6, label='3+2')
m = mlines.Line2D([], [], color='magenta', marker='o', markersize=6, label='4+1')

n = mlines.Line2D([], [], color='rosybrown', marker='o', markersize=6, label='1+1+1+1+1+1')
o = mlines.Line2D([], [], color='gold', marker='o', markersize=6, label='2+1+1+1+1')
pp = mlines.Line2D([], [], color='tan', marker='o', markersize=6, label='2+2+1+1')
q = mlines.Line2D([], [], color='navy', marker='o', markersize=6, label='2+2+2')
r = mlines.Line2D([], [], color='olive', marker='o', markersize=6, label='3+1+1+1')
s = mlines.Line2D([], [], color='coral', marker='o', markersize=6, label='3+2+1')
t = mlines.Line2D([], [], color='darkseagreen', marker='o', markersize=6, label='3+3')
u = mlines.Line2D([], [], color='bisque', marker='o', markersize=6, label='4+1+1')
v = mlines.Line2D([], [], color='lavender', marker='o', markersize=6, label='4+2')
w = mlines.Line2D([], [], color='darkviolet', marker='o', markersize=6, label='5+1')
z = mlines.Line2D([], [], color='black', marker='o', markersize=6, label='cancer cells')

legend(handles=[a, b, c, d, e, f, g, hh, i, j, k, l, m, n, o, pp, q, r, s, t, u, v, w, z], bbox_to_anchor=(1, 1), loc=2, borderaxespad=0.2)

xlabel('generation')
ylabel('abundance of each mode')
title('Time Series of Fragmentation Modes')
show()



#!./venv/bin/python

"""

File: synaptogenesis.py
Ported by: Sean Gallogly
Original Author: Joe Whitley-Casto

Description: This file aims to simulate the process of connection formation 
             between Golgi cell neurons. The neurons are visualized on a grid,
             and the user has a choice of which Golgi cell to look at in order
             to see which of its neighbors it has connected with.

             This file was written initially in Matlab by Joe Whitley-Casto,
             and was ported to Python by Sean Gallogly. Sole credit for the 
             algorithm and visualization goes to the original author.

"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# ============================= UTILITY FUNCTIONS =============================


def random_permute_range(range):
    return np.random.choice(np.arange(range), range, False)


# ================================ MAIN ALGORITHM ==============================

def connect_go_go(go_go_input, go_go_num_input, go_go_output, go_go_num_output,
        num_go_x, num_go_y, num_go, num_con, win_width):

    # total number of possible connections in the window
    num_p_con = (win_width + 1) ** 2
    
    # boolean connectivity array for checking double connections
    con_bool = np.zeros((num_go, num_go), dtype=bool)

    # index arrays representing a window of dimensions (WIN_WIDTH x WIN_WIDTH) 
    win_width_arr_x = np.array([i - win_width // 2 for i in np.arange(win_width + 1)])
    win_width_arr_y = np.array([i - win_width // 2 for i in np.arange(win_width + 1)])

    # repeat index arrays for as many connections as we will make 
    poss_x_coords = np.array([win_width_arr_x[i % (win_width + 1)] for i in np.arange(num_p_con)])
    poss_y_coords = np.array([win_width_arr_y[i // (win_width + 1)] for i in np.arange(num_p_con)])

    # shuffles the indices of every Golgi cell
    go_index = random_permute_range(num_go) 

    # build connections for every Golgi cell
    for i in np.arange(num_go):

        # locate our current cell 
        connector_go_x = go_index[i] % num_go_x
        connector_go_y = go_index[i] // num_go_y

        # randomize the connection order and assign 
        
        all_connected_go_x = poss_x_coords[random_permute_range(num_p_con)]
        all_connected_go_y = poss_y_coords[random_permute_range(num_p_con)]

        # loop through each potential connection within our window
        for j in np.arange(num_p_con):

            # initial offset of this golgi from our 'base' Golgi cell
            this_connected_go_x = connector_go_x + all_connected_go_x[j]
            this_connected_go_y = connector_go_y + all_connected_go_y[j]

            # ensures wrapping within the window and all golgi cells
            this_connected_go_x = (this_connected_go_x % (num_go_x ** 2)) % num_go_x 
            this_connected_go_y = (this_connected_go_y % (num_go_y ** 2)) % num_go_y

            # set the index of jth golgi cell in golgi array
            golgi_j_index = this_connected_go_y * num_go_x + this_connected_go_x 

            # test for self connections, already-established connections,
            # and range checks for number of connections
            if (golgi_j_index != go_index[i]
                and not con_bool[go_index[i], golgi_j_index]
                and num_go_go_input[golgi_j_index] < num_con 
                and num_go_go_output[go_index[i]] < num_con):

                go_go_output_arr[go_index[i], num_go_go_output[go_index[i]]] = golgi_j_index
                num_go_go_output[go_index[i]] += 1 

                go_go_input_arr[golgi_j_index, num_go_go_input[golgi_j_index]] = go_index[i] 
                num_go_go_input[golgi_j_index] += 1 

                con_bool[go_index[i], golgi_j_index] = 1        


# ================================= PLOTTING ===================================


def plot_rel_connections(connector_go_index, num_go_x, num_go_y, num_go, 
        go_go_output, go_go_num_output):

    # encodify 1-D Golgi array into 2D Array
    go_x_coord = np.array([i % num_go_x for i in np.arange(num_go)])
    go_y_coord = np.array([i // num_go_y for i in np.arange(num_go)])

    connector_go_num_con = go_go_num_output[connector_go_index] # extremely likely to be num_con

    # coordinates for our arbitrary point
    connector_go_x = connector_go_index % num_go_x
    connector_go_y = connector_go_index % num_go_y

    # create the connections arrays for visualization
    connecteds_go_x = go_x_coord[go_go_output[connector_go_index, np.arange(connector_go_num_con)]] 
    connecteds_go_y = go_y_coord[go_go_output[connector_go_index, np.arange(connector_go_num_con)]] 

    # plotting time 
    plt.figure(facecolor='k')
    axis = plt.axes(xlim=(-1, num_go_x), ylim=(-1, num_go_y)) 
    axis.set_facecolor('k')
    plt.scatter(go_x_coord, go_y_coord, 20, 'w' )
    plt.scatter(connector_go_x, connector_go_y, 30, 'g')
    plt.scatter(connecteds_go_x, connecteds_go_y, 30, 'r')
    plt.show()


# ================================== MAIN =====================================


if __name__ == '__main__':

    # main variables
    NUM_GO_X      = 32
    NUM_GO_Y      = 32
    NUM_CON       = 16 
    WIN_WIDTH     = 16
   
    NUM_GO = NUM_GO_X * NUM_GO_Y

    # input and output arrays
    go_go_output_arr = -1 * np.ones((NUM_GO, NUM_CON), dtype=np.intc)
    num_go_go_output = np.zeros(NUM_GO, dtype=np.intc)
    go_go_input_arr = -1 * np.ones((NUM_GO, NUM_CON), dtype=np.intc)
    num_go_go_input = np.zeros(NUM_GO, dtype=np.intc)

    # Make the connections
    connect_go_go(go_go_input_arr, num_go_go_input, go_go_output_arr, num_go_go_output,
            NUM_GO_X, NUM_GO_Y, NUM_GO, NUM_CON, WIN_WIDTH)

    # Prepare the golgi x and y coords for plotting
    go_x_coord = np.array([i % NUM_GO_X for i in np.arange(NUM_GO)])
    go_y_coord = np.array([i // NUM_GO_Y for i in np.arange(NUM_GO)])

    # prepare the figure
    fig = plt.figure(facecolor='k')
    ax  = plt.axes(xlim=(-1, NUM_GO_X), ylim=(-1, NUM_GO_Y))

    # declare this function here so that globals in main's scope can be updated.
    # there's a way to delcare functions (as in a statically-typed language) but
    # I've hacked on this enough for the time being.
    def animate_rel_connections(i):
        
        # arrays we'll be updating every fram
        connector_go_num_con = num_go_go_output[i]
        connector_go_x = i % NUM_GO_X 
        connector_go_y = i // NUM_GO_Y
        connecteds_go_x = go_x_coord[go_go_output_arr[i, np.arange(connector_go_num_con)]] 
        connecteds_go_y = go_y_coord[go_go_output_arr[i, np.arange(connector_go_num_con)]] 

        # give us some output while generating frames
        if i % (2 ** 5) == 0:
            print(f'Animating Golgi {i}.')

        ax.clear()
        ax.set_facecolor('k')
        ax.set_xlim([-1, NUM_GO_X])
        ax.set_ylim([-1, NUM_GO_Y])
        ax.scatter(go_x_coord, go_y_coord, 20, 'w' )
        ax.scatter(connector_go_x, connector_go_y, 30, 'g')
        ax.scatter(connecteds_go_x, connecteds_go_y, 30, 'r')

    # Make your own anime
    anime = FuncAnimation(
        fig      = fig,                      # the figure we're plotting upon
        func     = animate_rel_connections,  # function called for every frame
        frames   = NUM_GO,                   # total frames to draw
        interval = 200                       # delay (in ms) between frames
    )
    
    # Save your anime to show to the world.
    anime.save(f'synaptogenesis_{NUM_GO}_{NUM_CON}_{WIN_WIDTH}.gif')


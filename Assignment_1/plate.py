import numpy as np
import numpy.linalg as alg
import matplotlib.pyplot as plt


M = 10      #number of cells in the x direction
X = 0.1     #lenght in m of the x direction
dx =X/M     #lenght of each cell in x

N = 20      #number of cells in the y direction
Y = 0.2     #lenght in m of the y direction
dy =Y/N     #lenght of each cell in y

k = 45      #thermal convectivity coefficient

T_air = 20  #temperature of the surrounding air
T_top = 300 #temperature of the top furnace wall
T_bot = 80  #temperature of the bottom furnace wall

d = 0.1     #depth of support plate

def define_A_b(N, M, T_air, T_top, T_bot, h_left, h_right, k):
    """Creates the matrix A containing how each cell conducts heat, and the corresponding vector b with the constants of the heat transfer"""
    
    #creating a zero matrix and vector for values to be set into
    A = np.zeros((N*M,N*M), float)
    b = np.zeros((N*M,1), float)
    #setting a general offset for each iteration
    offset = 0
    
    #top left corner node
    A[offset][offset] = -3-h_left*dx/k
    A[offset][offset+1] = 1
    A[offset][offset+M] = 1

    b[offset] = -(h_left*dx/k)*T_air-T_top
    offset += 1
    
    
    #top edge nodes (excluding corners)
    for i in range(0, M-2):
        A[offset][offset-1] = 1
        A[offset][offset] = -4
        A[offset][offset+1] = 1
        A[offset][offset+M] = 1
        
        b[offset] = -T_top
        offset+=1


    #top right coner node
    A[offset][offset-1] = 1
    A[offset][offset] = -3-h_right*dx/k
    A[offset][offset+M] = 1
    
    b[offset] = -(h_right*dx/k)*T_air-T_top
    offset+=1
    
    
    #looping through each row
    for x in range(0, N-2):
        #left edge node
        A[offset][offset-M] = 1
        A[offset][offset] = -3 - h_left*dx/k
        A[offset][offset+1] = 1
        A[offset][offset+M] = 1
        
        b[offset] = -(h_left*dx/k)*T_air
        offset+=1
        
        
        for y in range(0, M-2):
            #internal node
            A[offset][offset-M] = 1
            A[offset][offset-1] = 1
            A[offset][offset] = -4
            A[offset][offset+1] = 1
            A[offset][offset+M] = 1
            offset+=1
            
            
        #right edge node
        A[offset][offset-M] = 1
        A[offset][offset-1] = 1
        A[offset][offset] = -3 - h_right*dx/k
        A[offset][offset+M] = 1
        
        b[offset] = -(h_right*dx/k)*T_air
        offset+=1
        
    
    #bottom left corner node
    A[offset][offset-M] = 1
    A[offset][offset] = -3 - h_left*dx/k
    A[offset][offset+1] = 1
    
    b[offset] = -(h_left*dx/k)*T_air-T_bot
    offset+=1
    
    #bottom edge nodes
    for i in range(0, M-2):
        A[offset][offset-M] = 1
        A[offset][offset-1] = 1
        A[offset][offset] = -4
        A[offset][offset+1] = 1
        
        b[offset] = -T_bot
        offset+=1
        
    #bottom right coner node
    A[offset][offset-M] = 1
    A[offset][offset-1] = 1
    A[offset][offset] = -3 - h_right*dx/k
    
    b[offset] = -(h_right*dx/k)*T_air-T_bot
        
    return(A, b)


def create_temp_matrix(A, b, N, M):
    """Creates a temperature matrix containing the temperatures for each cell"""
    #solving for the temperature of each cell
    temp_vector = alg.solve(A, b)
    temp_vector = np.flip(temp_vector)

    #putting each temperature into matrix form
    temp_matrix = np.zeros((N,M), float)

    i = 0
    #looping through each temperature in the temp_vector and putting it in the temp_matrix
    #where each value is put in each row and wraps to next rows.
    for a in range(N):
        for b in range(M):
            temp_matrix[a][b] = temp_vector[i]
            i+=1
    return temp_matrix


def plot_heatmap(temp_matrix, N, M):
    """Plots the heatmap using the generated temperature matrix to show how heat flows through the system"""
    
    #setting up grids to define where each cell sits in the heatmap
    X_grid,Y_grid = np.meshgrid(range(0,M), range(0,N))

    X_grid = dx/2+dx*X_grid
    Y_grid = Y-dy/2+dy*Y_grid

    #plotting the heatmap
    axs = plt.axes()

    heatmap = plt.pcolor(X_grid, Y_grid, temp_matrix, edgecolors='k', linewidths=0.25, shading = 'auto')
    plt.colorbar(heatmap)
    axs.set_aspect("equal")
    plt.show()


def heat_conducted_from_side(d, k, T_row, T_side, M):
    """Calculates the heat transfer to the support plate on the top row (in contact with the hot furnace wall) from the top furnace wall"""

    
    #loops through each node in the top row and adds togeather the Q_dot values
    Q_dot = 0
    for i in range(0, M):
        Q_dot += -k*d*(T_row[i]-T_side)
    return(Q_dot)

def heat_convected_to_air(d, h_left, h_right, T_matrix, T_air, N, M, dy):
    """Calculates the heat transfer from the plate (in contact with the air) to the air"""
    left_temps = []
    for i in range(0, N):
        left_temps.append(T_matrix[i][0])

    right_temps = []
    for i in range(0, N):
        right_temps.append(T_matrix[i][M-1])

    Q_dot_left = 0
    Q_dot_right = 0
    for i in range(0, N):
        Q_dot_left +=   -h_left*d*(T_air - left_temps[i])*dy
        Q_dot_right +=  -h_right*d*(T_air - right_temps[i])*dy
    return(Q_dot_left, Q_dot_right)

def run_program():
    #returns matrix A and vector b
    A, b = define_A_b(N, M, T_air, T_top, T_bot, h_left, h_right, k)

    #returns A and b to get temperatures and forms it into a matrix
    temp_matrix = create_temp_matrix(A,b, N, M)

    #returns the heat transfer from the top furnace wall to the cells in contact with the top wall
    Q_dot_of_top_plate = heat_conducted_from_side(d, k, temp_matrix[N-1],T_top, M)

    #returns the heat transfer from the bottom furnace wall to the cells in contact with the bottom wall
    Q_dot_of_bottom = heat_conducted_from_side(d, k, temp_matrix[0],T_bot, M)
    
    #returns the heat transfer from the cells in contact with air to the air
    Q_dot_of_left, Q_dot_right = heat_convected_to_air(d, h_left,h_right, temp_matrix, T_air, N, M, dy)
    
    return Q_dot_of_top_plate, Q_dot_of_bottom, Q_dot_of_left, Q_dot_right

#Allow anyone to select what they want to see/load (makes the plotting not intrusive when you dont want it :D)
to_see = " "

#allow user to set the h left and right values
h_left = int(input("\nWhat is the left h value: "))     #thermal conductivity coefficient for the left side
h_right = int(input("What is the right h value: "))     #thermal conductivity coefficient for the right side

while to_see != "":
    #ask the user what they want to see
    to_see = input("""\nWhat would you like to see?\n
Top Wall Heat Transfer [T]
Bottom Wall Heat Transfer [B]
Find Lowest value of h having heat flow into system from bottom [H]
Side Air Heat Transfer[S]
Plot the Heatmap with Balanced h values [P]
Reset h left and right values [R]
Leave Empty to Shutdown: """)
    to_see = to_see.lower()

    #runs define_A_b, create_temp_matrix, heat_conducted_from_top, heat_convected_from_bottom and heat_convected_to_air to find the Q vaules of each side of the plate
    Q_dot_of_top_plate, Q_dot_of_bottom, Q_dot_of_left, Q_dot_right = run_program()
    
    #allows the user to reset the h values
    if to_see == "r":
        h_left = int(input("\nWhat is the left h value: "))     #thermal conductivity coefficient for the left side
        h_right = int(input("What is the right h value: "))     #thermal conductivity coefficient for the right side
    
    #prints the value of the heat transfer from the top furnace wall into the plate
    if to_see == "t":
        print(f"Heat Transfer via conduction from top furnace wall into support plate is {Q_dot_of_top_plate:.2F}.")

    #prints the value of the heat transfer from the plate into the bottom furnace wall
    elif to_see == "b":
        print(f"Heat Transfer via conduction from the bottom furnace wall to the support plate is {Q_dot_of_bottom:.2F}.")

    elif to_see == "s":        
        print(f"Heat Transfer via convection from the support plate into the air is {Q_dot_of_left:.2F} on the left side and {Q_dot_right:.2F} on the right. Summing to {Q_dot_of_left+Q_dot_right:.2F}.")
    #high difference between T of the plate column nodes and T_air, so high heat transfer
    
    #find and print the correct steady state value of Q_dot_sides
        print(f"The steady state value of Heat transfer via convection from the support plate into the air is {Q_dot_of_top_plate+Q_dot_of_bottom:.2F} total.")

    #plots the heatmap of the temperatures of the balanced h's
    elif to_see == "p":
        plot_heatmap(temp_matrix, N, M)

    elif to_see == "h":
        #Loop through values of h to find the lowest value of h that gives a positive Q_dot_of_bottom.
        for i in range(0, 1000):
            h_left = i
            h_right = i
            #returns matrix A and vector b
            A, b = define_A_b(N, M, T_air, T_top, T_bot, h_left, h_right, k)
            #returns A and b to get temperatures and forms it into a matrix
            temp_matrix = create_temp_matrix(A,b, N, M)
            Q_dot_of_bottom = heat_conducted_from_side(d, k, temp_matrix[0],T_bot, M)
            if Q_dot_of_bottom > 0:
                break
        print(f"The lowest whole value of h where heat is transfered from the lower temperature furnace wall is {i}.")




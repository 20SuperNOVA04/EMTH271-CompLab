import numpy as np
import numpy.linalg as alg
import matplotlib.pyplot as plt


M = 5      #number of cells in the x direction
X = 0.1     #lenght in m of the x direction
dx =X/M     #lenght of each cell in x

N = 5      #number of cells in the y direction
Y = 0.1     #lenght in m of the y direction
dy =Y/N     #lenght of each cell in y

density = 7800  #denisty of the support plate
Cp = 460        #specific heat of the support plate
k = 45          #thermal conductivity coefficient

T_air = 20      #temperature of the surrounding air

d = 0.1         #depth of support plate

t0 = 0          #initial time
tfinal = 10000  #final time
dt = 25   #time between each simulation

T_past = np.ones((N,M), float)*20

def define_A_b(N, M, T_air, T_top, T_bot, T_past, h_left, h_right, k):
    """Creates the matrix A containing how each cell conducts heat, and the corresponding vector b with the constants of the heat transfer"""
    #creating a zero matrix and vector for values to be set into
    A = np.zeros((N*M,N*M), float)
    b = np.zeros((N*M,1), float)
    #setting a general offset for each iteration
    offset = 0
        
    #top left corner node
    A[offset][offset] = -( 3 + h_left*dx/k + (density*dx**2*Cp)/(k*dt) )
    A[offset][offset+1] = 1
    A[offset][offset+M] = 1

    b[offset] = -(h_left*dx/k)*T_air - T_top - ( (density*dx**2*Cp)/(k*dt) )*T_past[offset][offset]
    offset += 1
    
    #top edge nodes (excluding corners)
    for i in range(0, M-2):
        A[offset][offset-1] = 1
        A[offset][offset] = -(4 + (density*dx**2*Cp)/(k*dt))
        A[offset][offset+1] = 1
        A[offset][offset+M] = 1
        
        b[offset] = -T_top - ((density*dx**2*Cp) / (k*dt) ) * T_past[0][offset]
        offset+=1
        
    #top right coner node
    A[offset][offset-1] = 1
    A[offset][offset] = -(3 + h_right*dx/k + (density*dx**2*Cp)/(k*dt))
    A[offset][offset+M] = 1
    
    b[offset] = -(h_right*dx/k)*T_air - T_top - ((density*dx**2*Cp)/(k*dt))*T_past[0][offset]
    offset+=1

    
    #looping through each row
    for x in range(0, N-2):
        #left edge node
        
        T_past_offset_rows = x+1
        A[offset][offset-M] = 1
        A[offset][offset] = -(3 + h_right*dx/k + (density*dx**2*Cp)/(k*dt))
        A[offset][offset+1] = 1
        A[offset][offset+M] = 1
        
        b[offset] = -(h_left*dx/k)*T_air - ((density*dx**2*Cp)/(k*dt))*T_past[T_past_offset_rows][0]
        offset+=1
        
        for y in range(0, M-2):
            #internal node
            
            T_past_offset_cols = y+1
            A[offset][offset-M] = 1
            A[offset][offset-1] = 1
            A[offset][offset] = -(4 + (density*dx**2*Cp)/(k*dt))
            A[offset][offset+1] = 1
            A[offset][offset+M] = 1
            
            b[offset] = ((density*dx**2*Cp)/(k*dt))*T_past[T_past_offset_rows][T_past_offset_cols]
            offset+=1
            
        #right edge node
        A[offset][offset-M] = 1
        A[offset][offset-1] = 1
        A[offset][offset] = -(3 + h_right*dx/k + (density*dx**2*Cp)/(k*dt))
        A[offset][offset+M] = 1
        
        b[offset] = -(h_right*dx/k)*T_air - ((density*dx**2*Cp)/(k*dt))*T_past[T_past_offset_rows][M-1]
        offset+=1
    
    #bottom left corner node
    A[offset][offset-M] = 1
    A[offset][offset] = -(3 + h_left*dx/k + (density*dx**2*Cp)/(k*dt))
    A[offset][offset+1] = 1
    
    b[offset] = -(h_left*dx/k)*T_air - T_bot - ((density*dx**2*Cp)/(k*dt))*T_past[N-1][0]
    offset+=1
    
    #bottom edge nodes
    for i in range(0, M-2):
        A[offset][offset-M] = 1
        A[offset][offset-1] = 1
        A[offset][offset] = -4 - (density*dx**2*Cp)/(k*dt)
        A[offset][offset+1] = 1
        
        b[offset] = -T_bot - ((density*dx**2*Cp) / (k*dt) ) * T_past[N-1][i+1]
        offset+=1
        
    #bottom right coner node
    A[offset][offset-M] = 1
    A[offset][offset-1] = 1
    A[offset][offset] = -(3 + h_left*dx/k + (density*dx**2*Cp)/(k*dt))
    
    b[offset] = -(h_right*dx/k)*T_air - T_bot - ((density*dx**2*Cp) / (k*dt) ) * T_past[N-1][M-1]
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

def plot_heatmap(temp_matrix, N, M, title):
    """Plots the heatmap using the generated temperature matrix to show how heat flows through the system"""
    
    #setting up grids to define where each cell sits in the heatmap
    X_grid,Y_grid = np.meshgrid(range(0,M), range(0,N))
    X_grid = dx/2+dx*X_grid
    Y_grid = Y-dy/2+dy*Y_grid

    #plotting the heatmap
    axs = plt.axes()
    heatmap = plt.pcolor(X_grid, Y_grid, temp_matrix, edgecolors='k', linewidths=0.25, shading = 'auto', vmin = 20, vmax = 300)
    plt.colorbar(heatmap)
    axs.set_title(title)
    axs.set_xlabel("x position within the plate (mm)")
    axs.set_ylabel("y position within the plate (mm)")
    axs.set_aspect("equal")
    plt.show()
    
def plot_T_distrabuiton(data_list, t0, tfinal, dt):
    all_temps = get_left_temps(data_list)
    # Q_left = []
    # for i in range(0, len(all_temps)):
    #     Q_left.append(data_list[i][3])
        
    t = np.arange(t0, tfinal+1, dt)
    plt.plot(t, all_temps)
    # plt.set_xlabel("time (s)")
    # plt.set_ylabel("")
    plt.show()
    
def get_left_temps(data_list):
    temp_matricies = []
    all_left_temps = []
    for i in range(len(data_list)):
        temp_matricies.append(data_list[i][5])
    for matrix_num in range(len(temp_matricies)):
        temp_matrix = temp_matricies[matrix_num]
        left_temps = []
        for row in range(len(temp_matrix)):
            left_temps.append(temp_matrix[row][0])
        all_left_temps.append(left_temps)
    return(all_left_temps   )


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

def get_time_of_interest():
    stepofinterest = int(input(f"\nWhat timestep in sets of {dt} seconds would you like to see?\nPossible Values are 0 to {int((tfinal-t0)/dt)}\nInput: "))
    timeofinterest = stepofinterest*dt
    return stepofinterest, timeofinterest

def run_program():
    #returns matrix A and vector b
    A, b = define_A_b(N, M, T_air, T_top, T_bot, T_past, h_left, h_right, k)
    #returns A and b to get temperatures and forms it into a matrix
    temp_matrix = create_temp_matrix(A,b, N, M)
    #returns the heat transfer from the top furnace wall to the cells in contact with the top wall
    Q_dot_of_top_plate = heat_conducted_from_side(d, k, temp_matrix[N-1],T_top, M)
    #returns the heat transfer from the bottom furnace wall to the cells in contact with the bottom wall
    Q_dot_of_bottom = heat_conducted_from_side(d, k, temp_matrix[0],T_bot, M)
    #returns the heat transfer from the cells in contact with air to the air
    Q_dot_of_left, Q_dot_right = heat_convected_to_air(d, h_left,h_right, temp_matrix, T_air, N, M, dy)
    return Q_dot_of_top_plate, Q_dot_of_bottom, Q_dot_of_left, Q_dot_right, temp_matrix

#Allow anyone to select what they want to see/load (makes the plotting not intrusive when you dont want it :D)
to_see = " "
#allow user to set the h left and right values
h_left = int(input("\nWhat is the left h value: "))     #thermal conductivity coefficient for the left side
h_right = int(input("What is the right h value: "))     #thermal conductivity coefficient for the right side
print("Solving, Please wait")

data_list = []
for t in range(t0, tfinal+1, dt):
    #setting up top and bottom furnace wall temperatures for given t
    if t < 1000:
        T_top = 20
    else:
        T_top = 20 + 280*(1-np.e**(-(t-1000)/750))  #temperature of the top furnace wall
    T_bot = 20 + 60*(1-np.e**(-t/200))              #temperature of the bottom furnace wall
        
    #runs define_A_b, create_temp_matrix, heat_conducted_from_top, heat_convected_from_bottom and heat_convected_to_air to find the Q vaules of each side of the plate
    Q_dot_of_top_plate, Q_dot_of_bottom, Q_dot_of_left, Q_dot_right, temp_matrix = run_program()
    data_list.append([t, Q_dot_of_top_plate, Q_dot_of_bottom, Q_dot_of_left,Q_dot_right,temp_matrix])

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
      
    
    #allows the user to reset the h values
    if to_see == "r":
        h_left = int(input("\nWhat is the left h value: "))     #thermal conductivity coefficient for the left side
        h_right = int(input("What is the right h value: "))     #thermal conductivity coefficient for the right side
        
    #plots the heatmap of the temperatures of the balanced h's
    elif to_see == "p":
        stepofinterest, timeofinterest = get_time_of_interest()
        temp_matrix = data_list[stepofinterest][5]
        plot_heatmap(temp_matrix, N, M, f"Temperatures at {timeofinterest}s")

    elif to_see == "l":
        plot_T_distrabuiton(data_list, t0, tfinal, dt)

    # elif to_see == "h":
    #     #Loop through values of h to find the lowest value of h that gives a positive Q_dot_of_bottom.
    #     for i in range(0, 1000):
    #         h_left = i
    #         h_right = i
    #         #returns matrix A and vector b
    #         A, b = define_A_b(N, M, T_air, T_top, T_bot, h_left, h_right, k)
    #         #returns A and b to get temperatures and forms it into a matrix
    #         temp_matrix = create_temp_matrix(A,b, N, M)
    #         Q_dot_of_bottom = heat_conducted_from_side(d, k, temp_matrix[0],T_bot, M)
    #         if Q_dot_of_bottom > 0:
    #             break
    #     print(f"The lowest whole value of h where heat is transfered from the lower temperature furnace wall is {i}.")
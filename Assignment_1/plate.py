import numpy as np
import numpy.linalg as alg
import matplotlib.pyplot as plt


M = 10      #number of cells in the x direction
X = 0.1     #lenght in m of the x direction
dx =X/M     #lenght of each cell in x

N = 20      #number of cells in the y direction
Y = 0.2     #lenght in m of the y direction
dy =Y/N     #lenght of each cell in y

h = 25      #thermal conductivity coefficient
k = 45      #thermal convectivity coefficient

T_air = 20  #temperature of the surrounding air
T_top = 300 #temperature of the top furnace wall
T_bot = 80  #temperature of the bottom furnace wall



def create_temp_matrix(N, M, T_air, T_top, T_bot):
    A = np.zeros((N*M,N*M), int)
    b = np.zeros((N*M,1), float)
    #setting a general offset for each iteration
    offset = 0
    print(offset)
    
    #top left corner node
    A[offset][offset] = -3-h*dx/k
    A[offset][offset+1] = 1
    A[offset][offset+M] = 1

    b[offset] = -(h*dx/k)*T_air-T_top
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
    A[offset][offset] = -3-h*dx/k
    A[offset][offset+M] = 1
    
    b[offset] = -(h*dx/k)*T_air-T_top
    offset+=1
    
    
    #looping through each row
    for x in range(0, N-2):
        #left edge node
        A[offset][offset-M] = 1
        A[offset][offset] = -3 - h*dx/k
        A[offset][offset+1] = 1
        A[offset][offset+M] = 1
        
        b[offset] = -(h*dx/k)*T_air
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
        A[offset-1][offset] = 1
        A[offset][offset-1] = 1
        A[offset][offset] = -3 - h*dx/k
        A[offset][offset+M] = 1
        
        b[offset] = -(h*dx/k)*T_air
        offset+=1
        
    
    #bottom left corner node
    A[offset][offset-M] = 1
    A[offset][offset] = -3 - h*dx/k
    A[offset][offset+1] = 1
    
    b[offset] = -(h*dx/k)*T_air-T_bot
    offset+=1
    
    #bottom edge nodes
    for i in range(0, M-2):
        if offset+M < M*N:
            A[offset][offset+M] = 1
        A[offset][offset-1] = 1
        A[offset][offset] = -4
        A[offset+1][offset+1] = 1
        
        b[offset] = -T_bot
        offset+=1
        
    #bottom right coner node
    A[offset][offset-M] = 1
    A[offset][offset-1] = 1
    A[offset][offset] = -3 - h*dx/k
    
    b[offset] = -(h*dx/k)*T_air-T_bot
        
    
    
    return(A, b)

#creating matrix A and vector b
A, b = create_temp_matrix(N, M, T_air, T_top, T_bot)

print(A, "\n")
print(b, "\n")

#solving for the temperature of each cell
temp_vector = alg.solve(A, b)

#putting each temperature into matrix form
temp_matrix = np.zeros((N,M), float)

i = 0
for a in range(N):
    for b in range(M):
        temp_matrix[a][b] = temp_vector[i]
        i+=1

X_grid,Y_grid = np.meshgrid(range(0,M), range(0,N))

X_grid = dx/2+dx*X_grid
Y_grid = Y-dy/2+dy*Y_grid

axs = plt.axes()

heatmap = plt.pcolor(X_grid, Y_grid, temp_matrix, edgecolors='k', linewidths=0.25, shading = 'auto')
plt.colorbar(heatmap)
plt.show()



import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate
from scipy.optimize import fsolve
import scipy

# list1 = [50.19, 27.04, 28.70, 28.286, 33.79]

# avgv = sum(list1)/len(list1)
# numv = []

# for i in range(len(list1)):
#     avglist = [avgv,avgv,avgv,avgv,avgv]
#     numv.append((list1[i]-avglist[i])**2)
# SDv = np.sqrt(sum(numv)/(len(list1)-1))
# print(avgv, SDv)

# mlistPb = [list1[0]*11.34,list1[1]*11.34,list1[2]*11.34,list1[3]*11.34,list1[4]*11.34]
# mlistSn = [(100-list1[0])*7.29,(100-list1[1])*7.29,(100-list1[2])*7.29,(100-list1[3])*7.29,(100-list1[4])*7.29]


# mfraction = []
# numm=[]

# for i in range(len(list1)):
#     mfraction.append(mlistPb[i]/(mlistPb[i]+mlistSn[i]))

# avgmfraction = sum(mfraction)/len(mfraction)

# for i in range(len(list1)):
#     avglistm = [avgmfraction,avgmfraction,avgmfraction,avgmfraction,avgmfraction]
#     numm.append((mfraction[i]-avglistm[i])**2)
# SDm = np.sqrt(sum(numm)/(len(list1)-1))
# print(avgmfraction, SDm)

# xs = np.linspace(-2,2, 1000)

# ys1 = xs**3+1
# ys2 = 1-np.sin(4*xs+np.pi/3)

# solns = []
# for x0 in range(-2*100,2*100):
#     x0 = x0/100
#     soln = float(fsolve(lambda x: x**3+1- 1-np.sin(4*x+np.pi/3), x0))
    
#     if round(soln,4) not in solns:
#         solns.append(round(soln,4))
# print(solns)

# fig, ax = plt.subplots()
# ax.plot(xs, ys1);
# ax.plot(xs, ys2, color = "r");
# plt.show()

# A = 800
# b = 0.003
# t0 = 3
# K = 0.7
# Ts = 50
# y = 170
# h = 0.5

# t = t0

# def f(y, t):
#     return(A/((t-t0)**2+b)-K*(y-Ts))

# while t <= 10:
#     tnew = t+h
#     m = f(y,t)
#     ynew = y + (tnew-t)*m
#     print(ynew)
#     t = tnew
#     y = ynew


# fn = lambda y,t:A/((t-t0)**2+b)-K*(y-Ts)

# print(scipy.integrate.solve_ivp(fn, [3, 10], [170], method="RK23"))

T_past = np.array([[1,2,3],[4,5,6],[7,8,9]])
print(T_past)

N=3
M=5

A = np.zeros((N,M), float)

print(A)
"""
@author: Gloria Bazargan
"""

import math
import numpy as np

#Lattice parameters
print('Enter lattice parameters (in units of angstroms and degrees) separated by commas: a,b,c,alpha,beta,gamma')
cell = input("-->").split(',')
a,b,c,alpha,beta,gamma = float(cell[0]),float(cell[1]),float(cell[2]),math.radians(float(cell[3])),math.radians(float(cell[4])),math.radians(float(cell[5]))
#File containing cartesian coordianates 
f = input('Enter name of XYZ file:  ')
xyz = open(f)
N = int(xyz.readline())
header = xyz.readline()
atom_symbol, cart_coords = ([] for i in range (2))
for line in xyz:
    atom,x,y,z = line.split()
    atom_symbol.append(atom)
    cart_coords.append([float(x),float(y),float(z)])
xyz.close()
C = np.asarray(cart_coords)
#Transformation matrix
T = np.zeros([3,3])
v = (a*b*c) * np.sqrt(1 - np.cos(alpha)**2 - np.cos(beta)**2 - np.cos(gamma)**2 + 2*np.cos(alpha)*np.cos(beta)*np.cos(gamma))
T[0] = 1/a, -np.cos(gamma)/(a*np.sin(gamma)), b*c*((np.cos(beta)*np.cos(gamma)-np.cos(beta))/(v*np.sin(gamma)))
T[1] = 0, 1/(b*np.sin(gamma)), a*c*(np.cos(beta)*np.cos(gamma)-np.cos(alpha)/(v*np.sin(gamma)))
T[2] = 0,0, (a*b*np.sin(gamma))/v
#Coordinate transformation cartesian to fractional 
F = np.zeros([N,3])
for i in range (N):
    F[i,0] = T[0,0]*C[i,0] + T[0,1]*C[i,1] + T[0,2]*C[i,2]
    F[i,1] = T[1,0]*C[i,0] + T[1,1]*C[i,1] + T[1,2]*C[i,2]
    F[i,2] = T[2,0]*C[i,0] + T[2,1]*C[i,1] + T[2,2]*C[i,2]
#Write output file containing fracitonal coordinates
a = atom_symbol
b = F.tolist()
f = open('fractional_coordinates.txt','w+')
f.write(str(N)+'\n')
f.write('a,b,c,alpha,beta,gamma: ' + str(cell) +'\n')
out = [i + '    ' + str(j) for i,j in zip (a,b)]
out = ([i.replace('[', '     ') for i in out])
out = ([i.replace(']', '     ') for i in out])
out = ([i.replace(',', '     ') for i in out])
for i in out:
    f.write('\n' + i)
f.close()



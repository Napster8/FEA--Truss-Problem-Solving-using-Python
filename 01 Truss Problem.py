# Dependencies
import copy
import math
import numpy as np
import matplotlib.pyplot as plt

# Constants
# Young's Modulus for Steel in N/m^2
E=200 * 10 ** 9

# Area in (m^2)
A=0.005

# Element 1 (node 1 to 2)
theta=0  # in radians
L=3  # in metres

# 01 Top left quadrant
e11=math.cos(theta) ** 2
e12=math.cos(theta) * math.sin(theta)
e21=math.cos(theta) * math.sin(theta)
e22=math.sin(theta) ** 2

k11_12=(E * A / L) * np.array([[e11, e12], [e21, e22]])

# 02 Top right quadrant
e11=-math.cos(theta) ** 2
e12=-math.cos(theta) * math.sin(theta)
e21=-math.cos(theta) * math.sin(theta)
e22=-math.sin(theta) ** 2

k12_12=(E * A / L) * np.array([[e11, e12], [e21, e22]])

# 03 Bottom left quadrant
el1=-math.cos(theta) ** 2
el2=-math.cos(theta) * math.sin(theta)
e21=-math.cos(theta) * math.sin(theta)
e22=-math.sin(theta) ** 2

k21_12=(E * A / L) * np.array([[e11, e12], [e21, e22]])

# 04 Bottom right quadrant
e11=math.cos(theta) ** 2
e12=math.cos(theta) * math.sin(theta)
e21=math.cos(theta) * math.sin(theta)
e22=math.sin(theta) ** 2

k22_12=(E * A / L) * np.array([[e11, e12], [e21, e22]])

# Joining/Concatenating Top quadrants and Bottom quadrants
top=np.concatenate((k11_12, k12_12), axis=1)

# Joining/Concatenating bottom quadrants and Bottom quadrants
bottom=np.concatenate((k21_12, k22_12), axis=1)

# Joining/Concatenating 'top' and 'bottom' along horizontal axis
# Global Stiffness Matrix for Element 1
k1g=np.concatenate((top, bottom), axis=0)

print("\nGlobal Stiffness matrix for 1: \n")
print(np.round(k1g, 3))

# Element 2 (node 1 to 2)
theta=2.2143  # in radians
L=5  # in metres

# 01 Top left quadrant
e11=math.cos(theta) ** 2
e12=math.cos(theta) * math.sin(theta)
e21=math.cos(theta) * math.sin(theta)
e22=math.sin(theta) ** 2

k11_23=(E * A / L) * np.array([[e11, e12], [e21, e22]])

# 02 Top right quadrant
e11=-math.cos(theta) ** 2
e12=-math.cos(theta) * math.sin(theta)
e21=-math.cos(theta) * math.sin(theta)
e22=-math.sin(theta) ** 2

k12_23=(E * A / L) * np.array([[e11, e12], [e21, e22]])

# 03 Bottom left quadrant
el1=-math.cos(theta) ** 2
el2=-math.cos(theta) * math.sin(theta)
e21=-math.cos(theta) * math.sin(theta)
e22=-math.sin(theta) ** 2

k21_23=(E * A / L) * np.array([[e11, e12], [e21, e22]])

# 04 Bottom right quadrant
e11=math.cos(theta) ** 2
e12=math.cos(theta) * math.sin(theta)
e21=math.cos(theta) * math.sin(theta)
e22=math.sin(theta) ** 2

k22_23=(E * A / L) * np.array([[e11, e12], [e21, e22]])

# Joining/Concatenating Top quadrants and Bottom quadrants
top=np.concatenate((k11_23, k12_23), axis=1)

# Joining/Concatenating bottom quadrants and Bottom quadrants
bottom=np.concatenate((k21_23, k22_23), axis=1)

# Joining/Concatenating 'top' and 'bottom' along horizontal axis
# Global Stiffness Matrix for Element 1
k2g=np.concatenate((top, bottom), axis=0)

print("\nGlobal Stiffness matrix for 2: \n")
print(np.round(k2g, 3))

# Primary Stiffness matrix for the structure

k11=k11_12
k12=k12_12
k13=np.zeros([2, 2])

k21=k21_12
k22=k22_12 + k11_23
k23=k12_23

k31=np.zeros([2, 2])
k32=k21_23
k33=k22_23

r1=np.concatenate((k11, k12, k13), axis=1)
r2=np.concatenate((k11, k22, k23), axis=1)
r3=np.concatenate((k31, k32, k33), axis=1)
kp=np.concatenate((r1, r2, r3), axis=0)
print(kp)

# Reduce to Structure stiffness matrix and kp is copied using copy library

kp_red=copy.copy(kp)

# Impose the fact that U_xl = 0

kp_red[:, 0]=0
kp_red[0, :]=0
kp_red[0, 0]=1
print(kp_red)

# Impose the fact that U_yl = 0

kp_red[:, 1]=0
kp_red[1, :]=0
kp_red[1, 1]=1
print(kp_red)

# Impose the fact that U_x3 = 0

kp_red[:, 4]=0
kp_red[4, :]=0
kp_red[4, 4]=1
print(kp_red)

# Extracting the relevant portion from the structure stiffness matrix
ks=kp[2:4, 2:4]
ks=np.matrix(ks)
u2=ks.I * np.array([[0], [-150000]])
u_x2=u2[0, 0]
u_y2=u2[1, 0]

one=u_x2

print(f"\nThe horizontal displacement at node 2 is {one} m (to the left)")
print(f"The vertical displacement at node 2 is {one} m (downwards)")

## Determining the reaction forces

UG=np.array([[0, 0, u_x2, u_y2, 0, 0]])
FG=np.matmul(kp, UG.T)
print(np.round(FG, 3))

# Forces: Element A(nodes 1 to 2)
theta=0
L=3

# Transformation matrix
c=math.cos(theta)
s=math.sin(theta)
T=np.array([[c, s, 0, 0], [0, 0, c, s]])
disp=np.array([[0, 0, u_x2, u_y2]]).T
disp_local=np.matmul(T, disp)
F12=(E * A / L) * (disp_local[1] - disp_local[0]).item()
print("The force in element A is {one} kN".format(one=round(F12 / 1000, 1)))

## Forces: Element B(nodes 2 to 3)


theta=2.2143
L=5

# Transformation matrix
c = math.cos(theta)
s = math.sin(theta)
T = np.array([[c, s, 0, 0], [0, 0, c, s]])
disp = np.array([[u_x2, u_y2, 0, 0]]).T
disp_local = np.matmul(T, disp)
F12 = (E * A / L) * (disp_local[1] - disp_local[0]).item()
print("The force in element B is {one} kN".format(one=round(F12 / 1000, 1)))

# Visualizations

xfactor = 100

label = "Deflection scale factor: " + str(xfactor) + "\nux2 = " + str(round(u_x2, 6)) + " m \nuy2 = " + str(round(u_y2, 6)) + " m"

fig = plt.figure()
axes = fig.add_axes([0,0,2,2])
# For the figure to be in the correct 1:1 aspect ratio
fig.gca().set_aspect('equal', adjustable='box')

axes.plot([0,3],[0,0],'b')
axes.plot([0,3],[4,0],'b')

axes.plot([0,3+u_x2*xfactor],[0,u_y2],'--r')
axes.plot([0,3+u_x2*xfactor],[4,u_y2],'--r')

axes.plot([0],[0],'bo')
axes.plot([3],[0],'bo')
axes.plot([0],[4],'bo')

plt.text(1,4, label, fontsize=14, verticalalignment="top")
axes.grid()
axes.set_xlabel("Distance (m)")
axes.set_ylabel("Distance (m)")
axes.set_title("Deflected shape")

plt.show()
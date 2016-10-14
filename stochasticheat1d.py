import numpy as np
import matplotlib.pyplot as plt
from pdetools import *
import math

tBegin = 0
tEnd = 5
xBegin=0
xEnd=10
f=lambda x: math.sin(x)
l=lambda x: 0.0
r=lambda x: 0.0
Nt=1000
Nx=50

dt = float(tEnd-tBegin)/Nt
dx=float(xEnd-xBegin)/Nx

theta = 1
mu = 1.2
sigma = 0.3

sqrtdt = np.sqrt(dt)

IC=[]
for i in range(0,Nx+1):
	x=xBegin+i*dx
	IC.append(f(x))


Solution=[]
Solution.append(IC)

xr=[]
t=[]

for i in range(0,Nt+1):
	tg=tBegin+i*dt
	t.append(tg)

for i in range(0,Nx+1):
	xg=xBegin+i*dx
	xr.append(xg)



c=1
for i in range(1,Nt+1):
	temp=[]
	pastsol=Solution[-1]
	for k in range(0,Nx+1):
		if k==0:
			temp.append(l(t[i]))
		elif k==Nx:
			temp.append(r(t[i]))
		else:
			currentRan=np.random.gamma(1.0,3.0)
			temp.append(pastsol[k] + (c + currentRan )*(dt/dx**2)*(pastsol[k-1] - 2*pastsol[k] + pastsol[k+1]))

	Solution.append(temp)

makeMovie(xr,Solution,t,'Still trying to find f hat','x','u')


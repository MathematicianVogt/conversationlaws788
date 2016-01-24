import math

#flux function for project one
def projectICRarefaction(self,x,pl):
	if(x < 0):
		return pl
	else:
		return 1
def projectICShock(self,x,pr):
	if(x < 0):
		return 1
	else:
		return pr




def projflux(u):
	return u*(1-u)

#Numerical flux for lax fredrich scheme
def laxf(flux,a,b,deltax,deltat):
	return .5*(flux(a) +flux(b)) - (deltax/(2.0*deltat))*(b-a)
#Numerical flux for godonov scheme
def god(flux,a,b,deltax,deltat):
	pass
#Numerical flux for eng scheme
def eng(flux,a,b,deltax,deltat):
	return .5*(flux(a) +flux(b)) -.5*((b-a)/6.0)*(math.fabs(derv(flux,a,deltax)) +4.0*derv(flux,(a+b)/2.0,deltax) + derv(flux,b,deltax) )
#Numerical flux for lax wendroff scheme
def laxw(flux,a,b,deltax,deltat):
	return .5*(flux(a) +flux(b)) - (deltax/(2.0*deltat))*derv(flux,(a+b)/2.0,deltax)*(flux(b)-flux(a))

#centeral difference approximation for first derivative 
def derv(f,x, deltax):
	return (f(x+deltax) -f(x-deltax))/(2*deltax)


#code variables
Nx=100
Nt=100
a=-1
b=1
t0=0
t1=1
deltax=(b-a)/float(Nx)
deltat=(t1-t0)/float(Nt)

'''
Python Code for Solving 1-D conservation law with arbitary flux function using
classical numerical methods.

Usage:
-Define Flux function
-Define x range [a,b], t range [t0,t1] and respective steps taken Nx, Ny
-Define initial condition

'''


tmat=[]
xmat=[]

for i in range(0,Nx+1):
	xmat.append(a+i*deltax)

for i in range(0,Nt+1):
	tmat.append(t0+i*deltat)













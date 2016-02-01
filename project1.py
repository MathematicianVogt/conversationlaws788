# This example uses a MovieWriter directly to grab individual frames and
# write them to a file. This avoids any event loop integration, but has
# the advantage of working with even the Agg backend. This is not recommended
# for use in an interactive setting.
# -*- noplot -*-
import math
import matplotlib.pyplot as plt
import time
import matplotlib.animation as manimation
import numpy as np
import matplotlib


#flux function for project one
def projectICShock(x,pl):
	if(x < 0):
		return pl
	else:
		return 1
def projectICRarefaction(x,pr):
	if(x < 0):
		return 1
	else:
		return pr

def projflux(u):
	return u*(1-u)

def trueSolShock(x,t,pl):
	if(t< (-x/pl)):
		return pl
	if(t>(-x/pl)):
		return 1
	if(t==(-x/pl)):
		return pl
def trueSolRare(x,t,ul,ur):
	if(t<=x/(1-2*ul)):
		return ul
	elif(t>x/(1-2*ul) and t>x/(1-2*ur)):
		return (-x/(2*t)) +(1.0/2.0)

	elif(t<=x/(1-2*ur)):
		return ur

	else:
		print 'erros'



#Numerical flux for lax fredrich scheme
def laxf(flux,a,b,deltax,deltat):
	return .5*(flux(a) +flux(b)) - (deltax/(2.0*deltat))*(b-a)
#Numerical flux for godonov scheme
def god(flux,a,b,deltax,deltat):
	if(a<=b and a<= (1.0/2.0) and b<=(1.0/2.0)):
		return flux(a)
	elif(a<=b and a>= (1.0/2.0) and b>=(1.0/2.0)):
		return flux(b)
	elif(b<=a and a<= (1.0/2.0) and b<=(1.0/2.0)):
		return flux(a)
	elif (b<=a and a>= (1.0/2.0) and b>=(1.0/2.0)):
		return flux(b)
	elif (b<=a and b<=(1.0/2.0) and a>=(1.0/2.0)):
		
		return flux(1.0/2.0)
	elif(a<=b and a<=(1.0/2.0) and b>=(1.0/2.0)):
		val=min(a,1-b)
		return flux(val)


#Numerical flux for eng scheme
def eng(flux,a,b,deltax,deltat):
	return .5*(flux(a) +flux(b)) -.5*((b-a)/6.0)*(math.fabs(derv(flux,a,deltax)) +4.0*math.fabs(derv(flux,(a+b)/2.0,deltax)) + math.fabs(derv(flux,b,deltax)) )
#Numerical flux for lax wendroff scheme
def laxw(flux,a,b,deltax,deltat):
	return .5*(flux(a) +flux(b)) - (deltat/(2.0*deltax))*(1-2*((a+b)/2.0))*(flux(b)-flux(a))

#centeral difference approximation for first derivative 
def derv(f,x, deltax):
	return (f(x+deltax) -f(x-deltax))/(2*deltax)

class trafficconservationlaw:
	'''
	Python Code for Solving 1-D conservation law with arbitary flux function using
	classical numerical methods.

	Usage:
	-Define Flux function
	-Define x range [a,b], t range [t0,t1] and respective steps taken Nx, Ny
	-Define initial condition

	'''

	def __init__(self,a,b,t0,t1,Nx,Nt,flux):
		self.Nx=Nx
		self.Nt=Nt
		self.a=a
		self.b=b
		self.t0=t0
		self.t1=t1
		self.deltax=(b-a)/float(Nx)
		self.deltat=(t1-t0)/float(Nt)
		self.trafficflux=flux
		self.tmat=[]
		self.xmat=[]
		self.pl=.25
		self.pr=.25
		self.trueshock=[]
		self.truerare=[]
		self.shocklaxf=[]
		self.shockgod=[]
		self.shockeng=[]
		self.shocklaxw=[]
		self.rarelaxf=[]
		self.raregod=[]
		self.rareeng=[]
		self.rarelaxw=[]


		for i in range(0,Nx+1):
			self.xmat.append(a+i*self.deltax)

		for i in range(0,Nt+1):
			self.tmat.append(t0+i*self.deltat)

		#set up IC.
		uint1=[]
		uint2=[]
		for x in self.xmat:
			uint1.append(projectICShock(x,self.pl))
			uint2.append(projectICRarefaction(x,self.pr))
		self.trueshock.append(uint1)
		self.truerare.append(uint2)
		self.shocklaxf.append(uint1)
		self.shockgod.append(uint1)
		self.shockeng.append(uint1)
		self.shocklaxw.append(uint1)
		self.rarelaxf.append(uint2)
		self.raregod.append(uint2)
		self.rareeng.append(uint2)
		self.rarelaxw.append(uint2)

		#generate true solution plot. 
		for i in range(1,len(self.tmat)):
			tempsolshock=[]
			tempsolrare=[]
			
			for j in range(0,len(self.xmat)):
				tempsolshock.append(trueSolShock(self.xmat[j],self.tmat[i],self.pl))
				tempsolrare.append(trueSolRare(self.xmat[j],self.tmat[i],1.0,self.pr))
			#print tempsolrare
			self.trueshock.append(tempsolshock)
			self.truerare.append(tempsolrare)
			
		#Compute Solution for Shock Problem using different schemes
		self.computeSolution(self.shocklaxf,self.pl,1.0,"Solution - Lax Fred - Shock", laxf)
		self.computeSolution(self.shockgod,self.pl,1.0,"Solution - Godonov - Shock", god)
		self.computeSolution(self.shockeng,self.pl,1.0,"Solution - Eng Osh - Shock", eng)
		self.computeSolution(self.shocklaxw,self.pl,1.0,"Solution - Lax Wend - Shock", laxw)

		#Compute Solution for Rarefaction Problem using different schemes
		self.computeSolution(self.rarelaxf,1.0,self.pr,"Solution - Lax Fred - Rarefaction", laxf)
		self.computeSolution(self.raregod,1.0,self.pr,"Solution - Godonov- Rarefaction", god)
		self.computeSolution(self.rareeng,1.0,self.pr,"Solution - Eng Osh- Rarefaction", eng)
		self.computeSolution(self.rarelaxw,1.0,self.pr,"Solution - Lax Wend- Rarefaction", laxw)

	def computeSolution(self, solutionList, ul, ur,title,numericalFlux):
		u=solutionList[0]
		unext=[]
		for i in range(1,len(self.tmat)):
			unext=[]
			sol=0
			for j in range(0,len(self.xmat)):
				if(j==0):
					sol = u[j] - (self.deltat/self.deltax)*(numericalFlux(projflux,u[j],u[j+1],self.deltax,self.deltat) - numericalFlux(projflux,ul,u[j],self.deltax,self.deltat) )
				elif(j==len(self.xmat)-1):
					sol = u[j] - (self.deltat/self.deltax)*(numericalFlux(projflux,u[j],ur,self.deltax,self.deltat) - numericalFlux(projflux,u[j-1],u[j],self.deltax,self.deltat) )
				else:
					sol = u[j] - (self.deltat/self.deltax)*(numericalFlux(projflux,u[j],u[j+1],self.deltax,self.deltat) - numericalFlux(projflux,u[j-1],u[j],self.deltax,self.deltat) )
				unext.append(sol)
			u=unext
			solutionList.append(unext)
		self.makeMovie(self.xmat,solutionList,self.tmat,title)
		print "done"

	

	def plotshock(self):
		self.makeMovie(self.xmat,self.trueshock,self.tmat,"trueshockvogt")
	def plotrare(self):
		self.makeMovie(self.xmat,self.truerare,self.tmat,"truerarefacetionvogt")

	def solvePDE(self):
		pass
	def makeMovie(self,xList, solList,tList,title):
		FFMpegWriter = manimation.writers['ffmpeg']
		metadata = dict(title=title, artist='Matplotlib',
		                comment='Movie support!')
		writer = FFMpegWriter(fps=1, metadata=metadata)

		fig = plt.figure()
		l, = plt.plot([], [], 'k-o')

		plt.ylim(0, 2)
		plt.xlim(self.a-1,self.b+1)

		with writer.saving(fig, title+ ".mp4", 100):
			for i in range(0,len(tList)):

				x0 = self.xmat
				y0 =solList[i]

				plt.xlim([-10,10])
				plt.ylim([0,2])
				plt.plot(x0,y0)
				plt.title(title + " Time = " + str(tList[i]))
				plt.xlabel('x')
				plt.ylabel(r'$\rho$')
				writer.grab_frame()
				plt.clf()


sim=trafficconservationlaw(-5,5,0,10,100,200,projflux)
sim.plotshock()
sim.plotrare()















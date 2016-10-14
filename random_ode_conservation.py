#solve dx/dt = x , x(0)~uni(0,2) - Find the PDF

import math
import matplotlib.pyplot as plt
import time
import matplotlib.animation as manimation
import numpy as np
import matplotlib

class solz:
	def __init__(self,a,b,t0,t1,Nx,Nt):
			self.Nx=Nx
			self.Nt=Nt
			self.a=a
			self.b=b
			self.t0=t0
			self.t1=t1
			self.deltax=(b-a)/float(Nx)
			self.deltat=(t1-t0)/float(Nt)
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
			self.unia=0
			self.unib=2


			for i in range(0,Nx+1):
				self.xmat.append(a+i*self.deltax)

			for i in range(0,Nt+1):
				self.tmat.append(t0+i*self.deltat)

			self.intial_conditon=[]
			for i in range(0,len(self.xmat)):
				if(self.unia<=self.xmat[i]<=self.unib):
					self.intial_conditon.append(1.0/float(self.unib-self.unia))
				else:
					self.intial_conditon.append(0.0)

			total_sol=[]
			total_sol.append(self.intial_conditon)

			for i in range(1,len(self.tmat)):
				current_sol=[]
				prev_sol=total_sol[-1]
				print len(self.xmat)
				print len(prev_sol)
				for j in range(0,len(self.xmat)):
					print j
					if j==0:
						current_sol.append(0.0)
					elif j==len(self.xmat)-1:
						current_sol.append(0.0)
					else:
						current_sol.append(prev_sol[j] - (self.deltat/self.deltax)*(laxf(prev_sol[j],prev_sol[j+1],self.deltax,self.deltat,self.xmat[j]) - laxf(prev_sol[j-1],prev_sol[j],self.deltax,self.deltat,self.xmat[j])  ))
				total_sol.append(current_sol)

			makeMovie(self.xmat,total_sol,self.tmat,"random_sol")


def laxf(a,b,deltax,deltat,xpos):
	return .5*(xpos* a+ xpos*b) - (deltax/(2.0*deltat))*(b-a)


def makeMovie(xList, solList,tList,title):
		FFMpegWriter = manimation.writers['ffmpeg']
		metadata = dict(title=title, artist='Matplotlib',
		                comment='Movie support!')
		writer = FFMpegWriter(fps=1, metadata=metadata)

		fig = plt.figure()
		l, = plt.plot([], [], 'k-o')

		

		with writer.saving(fig, title+ ".mp4", 100):
			for i in range(0,len(tList)):

				x0 = xList
				y0 =solList[i]

				plt.xlim([-10,10])
				plt.ylim([0,2])
				plt.plot(x0,y0)
				plt.title(title + " Time = " + str(tList[i]))
				plt.xlabel('x')
				plt.ylabel(r'$f_x$')
				writer.grab_frame()
				plt.clf()






x=solz(0,10,0,5,100,300)
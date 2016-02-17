from  pdetools import *
import math
import numpy as np
import time

class wave:
	def __init__(self,a,b,t0,t1,L,R,IC,ICC,Nx,Nt):
		self.xint=onedinterval(a,b,float(b-a)/Nx,Nx)
		self.tint=onedinterval(t0,t1,float(t1-t0)/Nt,Nt)
		self.w1Left=[]
		self.w1Right=[]
		self.w1IC=[]
		self.w2Left=[]
		self.w2Right=[]
		self.w2IC=[]
		self.deltah=float(b-a)/Nx
		self.deltat=float(t1-t0)/Nt
		self.t0=t0
		self.a=a
		self.b=b
		self.L=L
		self.R=R
		self.IC=IC
		self.ICC=ICC
		print str(self.tint)

		for i in range(0,Nx+1):
			x=self.a+i*self.deltah
			#print x
			self.w1IC.append( (self.ICC(x) - derv(self.IC,x,self.deltah))/2.0 )
			self.w2IC.append( (self.ICC(x) + derv(self.IC,x,self.deltah))/2.0 )



		#Find solution for w1
		w1Sol=[]
		w1Sol.append(self.w1IC)

		for i in range(1,len(self.tint)):
			tempsol=[]
			for k in range(0,len(self.xint)):
				if k==0:
					tempsol.append(derv(self.L,self.tint[i],self.deltat)/2.0)
				elif k==len(self.xint)-1:
					tempsol.append(derv(self.R,self.tint[i],self.deltat)/2.0)
				elif(self.tint[i]>=self.xint[k]-self.a):
					tempsol.append(derv(self.L,self.tint[i]-self.xint[k]+self.a,self.deltat)/2.0)
				else:
					tempsol.append((self.ICC(self.xint[k]-self.tint[i]) - derv(self.IC,self.xint[k]-self.tint[i],self.deltah))/2.0)
			w1Sol.append(tempsol)

		#Find solution for w2
		w2Sol=[]
		w2Sol.append(self.w2IC)
		for i in range(1,len(self.tint)):
			tempsol=[]
			for k in range(0,len(self.xint)):
				if k==0:
					tempsol.append(derv(self.L,self.tint[i],self.deltat)/2.0)
				elif k==len(self.xint)-1:
					tempsol.append(derv(self.R,self.tint[i],self.deltat)/2.0)
				elif(self.tint[i]>=self.b-self.xint[k]):
					tempsol.append(derv(self.R,self.tint[i]+self.xint[k]-self.b,self.deltat)/2.0)
				else:
					tempsol.append((self.ICC(self.xint[k]+self.tint[i]) + derv(self.IC,self.xint[k]+self.tint[i],self.deltah))/2.0)
			w2Sol.append(tempsol)

		makeMovie(self.xint,w1Sol,self.tint,'w1 sol','x','w1')
		makeMovie(self.xint,w2Sol,self.tint,'w2 sol','x','w2')
		u1Sol=[]
		u2Sol=[]
		for i in range(0,len(self.tint)):
			tempsol1=[]
			tempsol2=[]
			w1part=w1Sol[i]
			w2part=w2Sol[i]

			for k in range(0,len(w1Sol[0])):
			
				tempsol1.append(w2part[k]-w1part[k])
				tempsol2.append(w1part[k]+w2part[k])
			u1Sol.append(tempsol1)
			u2Sol.append(tempsol2)



		
		
		makeMovie(self.xint,u1Sol,self.tint,'u1 sol','x','u1')
		makeMovie(self.xint,u2Sol,self.tint,'u2 sol','x','u2')
		
		mat=np.zeros((len(u1Sol[0]),len(u1Sol[0])))

		mat[0,0]=1.0
		mat[0,1]=0
		for i in range(1,len(u1Sol[0])-1):
			mat[i,i]=-1.0/(self.deltah)
			mat[i,i+1]=1.0/(self.deltah)
		mat[len(u1Sol[0])-1,len(u1Sol[0])-1]=1.0
		mat[len(u1Sol[0])-1,len(u1Sol[0])-2]=0

		truesol=[]
		for i in range(0,len(self.tint)):
			
			x=np.linalg.solve(mat,u1Sol[i])
			truesol.append(x)

		makeMovie(self.xint,truesol,self.tint,'u sol -using del x','x','u')
		
		mat=np.zeros((len(u1Sol[0]),len(u1Sol[0])))

		mat[0,0]=1.0
		mat[0,1]=0
		for i in range(1,len(u1Sol[0])-1):
			mat[i,i]=-1.0/(self.deltat)
			mat[i,i+1]=1.0/(self.deltat)
		mat[len(u1Sol[0])-1,len(u1Sol[0])-1]=1.0
		mat[len(u1Sol[0])-1,len(u1Sol[0])-2]=0
		truesol=[]
		for i in range(0,len(self.tint)):
			
			x=np.linalg.solve(mat,u2Sol[i])
			truesol.append(x)

		makeMovie(self.xint,truesol,self.tint,'u sol- using del t','x','u')




L = lambda x: 0
R= lambda x: 0
IC = lambda x: math.sin(x)
ICC = lambda x: 0

sol=wave(-4.0*math.pi,4.0*math.pi,0,4*math.pi,L,R,IC,ICC,50,100)





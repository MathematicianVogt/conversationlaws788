import matplotlib.pyplot as plt
import time
import matplotlib.animation as manimation
import numpy as np
import matplotlib
import math

def onedinterval(a,b,h,n):
	interval = []
	for i in range(0,n+1):
		interval.append(a+i*h)
	return interval
def derv(f,x,h):
	return (f(x+h)-f(x-h))/(2*h)



def makeMovie(xList, solList,tList,title,xlab,ylab):
		FFMpegWriter = manimation.writers['ffmpeg']
		metadata = dict(title=title, artist='Matplotlib',
		                comment='Movie support!')
		writer = FFMpegWriter(fps=1, metadata=metadata)

		fig = plt.figure()
		l, = plt.plot([], [], 'k-o')

		plt.ylim(-2, 2)
		plt.xlim(xList[0]-1,xList[len(xList)-1]+1)

		with writer.saving(fig, title+ ".mp4", 100):
			for i in range(0,len(tList)):

				x0 = xList
				y0 =solList[i]

				plt.xlim([0,10])
				plt.ylim([-2,2])
				plt.plot(x0,y0)
				plt.title(title + " Time = " + str(tList[i]))
				plt.xlabel(xlab)
				plt.ylabel(ylab)
				writer.grab_frame()
				plt.clf()

import os
import clawpack
from matplotlib import animation
import matplotlib.pyplot as plt
from clawpack.visclaw.JSAnimation import IPython_display
from clawpack import pyclaw
from clawpack import riemann
import numpy as np

path = os.path.abspath(clawpack.__file__)
print path

def fplot(frame_number):
    frame = claw.frames[frame_number]
    density = frame.q[0,:,:]
    (vx,vy) = np.gradient(density)
    vs = np.sqrt(vx**2 + vy**2)
    im.set_data(vs.T)
    #print im
    return im,



claw = pyclaw.Controller()
claw.tfinal = 0.6
claw.num_output_times = 40

riemann_solver = riemann.euler_4wave_2D
claw.solver = pyclaw.ClawSolver2D(riemann_solver)
claw.solver.all_bcs = pyclaw.BC.extrap

grid_size = (300,300)
domain = pyclaw.Domain( (0.,0.), (1.,1.), grid_size)

claw.solution = pyclaw.Solution(claw.solver.num_eqn,domain)
gam = 1.4
claw.solution.problem_data['gamma']  = gam

# Set initial data
q = claw.solution.q
xx,yy = domain.grid.p_centers
l = xx<0.5; r = xx>=0.5; b = yy<0.5; t = yy>=0.5
q[0,...] = 2.*l*t + 1.*l*b + 1.*r*t + 3.*r*b
q[1,...] = 0.75*t - 0.75*b
q[2,...] = 0.5*l  - 0.5*r
q[3,...] = 0.5*q[0,...]*(q[1,...]**2+q[2,...]**2) + 1./(gam-1.)

claw.keep_copy = True       # Keep solution data in memory for plotting
claw.output_format = None   # Don't write solution data to file

status = claw.run()


fig = plt.figure(figsize=[4,4])

frame = claw.frames[0]
density = frame.q[0,:,:]
(vx,vy) = np.gradient(density)
vs = np.sqrt(vx**2 + vy**2)

x, y = frame.state.grid.c_centers    

# This essentially does a pcolor plot, but it returns the appropriate object
# for use in animation.  See http://matplotlib.org/examples/pylab_examples/pcolor_demo.html.
# Note that it's necessary to transpose the data array because of the way imshow works.
im = plt.imshow(vs.T, cmap='Greys', vmin=vs.min(), vmax=vs.max()/20,
           extent=[x.min(), x.max(), y.min(), y.max()],
           interpolation='nearest', origin='lower')
b=animation.FuncAnimation(fig, fplot, frames=len(claw.frames), interval=20)
FFwriter = animation.FFMpegWriter()
b.save('basic_animation.mp4', writer = FFwriter, fps=50)

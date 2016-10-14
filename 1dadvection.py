#!/usr/bin/env python
# encoding: utf-8
import time
import numpy as np
from clawpack import riemann
from advectionsolver import adsolver


def setup(use_petsc=0,kernel_language='Python',outdir='./_output',solver_type='classic'):

    if use_petsc:
        import clawpack.petclaw as pyclaw
    else:
        from clawpack import pyclaw

    if kernel_language == 'Python': 
        riemann_solver = adsolver
    elif kernel_language == 'Fortran':
        riemann_solver = riemann.burgers_1D

    if solver_type=='sharpclaw':
        solver = pyclaw.SharpClawSolver1D(riemann_solver)
    else:
        solver = pyclaw.ClawSolver1D(riemann_solver)
        solver.limiters = pyclaw.limiters.tvd.vanleer

    solver.kernel_language = kernel_language
        
    solver.bc_lower[0] = pyclaw.BC.periodic
    #solver.user_bc_lower=custom_bc_one
    solver.bc_upper[0] = pyclaw.BC.periodic
    #solver.user_bc_upper=custom_bc_two
    solver.num_waves=1
    solver.num_eqn=1
    x = pyclaw.Dimension(0.0,1.0,500,name='x')
    domain = pyclaw.Domain(x)
    num_eqn = 1
    state = pyclaw.State(domain,num_eqn)

    xc = state.grid.x.centers
    state.q[0,:] = np.sin(np.pi*2*xc) + 0.50
    state.problem_data['efix']=True
    state.problem_data['a']=1.0

    claw = pyclaw.Controller()
    claw.tfinal = 2.0
    claw.solution = pyclaw.Solution(state,domain)
    claw.solver = solver
    claw.outdir = outdir
    claw.setplot = setplot
    claw.keep_copy = True

    return claw


def setplot(plotdata):
    """ 
    Plot solution using VisClaw.
    """ 
    plotdata.clearfigures()  # clear any old figures,axes,items data

    # Figure for q[0]
    plotfigure = plotdata.new_plotfigure(name='q[0]', figno=0)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = [-1., 20.]
    plotaxes.title = 'q[0]'

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='1d')
    plotitem.plot_var = 0
    plotitem.plotstyle = '-o'
    plotitem.color = 'b'
    
    return plotdata


if __name__=="__main__":
    from clawpack.pyclaw.util import run_app_from_main
    output = run_app_from_main(setup,setplot)
#!/usr/bin/env python
# encoding: utf-8
r"""
2D shallow water: flow over a sill
==================================

Solve the 2D shallow water equations with diagonal zero width barrier and
variable bathymetry:

.. :math:
    h_t + (hu)_x + (hv)_y & = 0 \\
    (hu)_t + (hu^2 + \frac{1}{2}gh^2)_x + (huv)_y & = -g h b_x \\
    (hv)_t + (huv)_x + (hv^2 + \frac{1}{2}gh^2)_y & = -g h b_y.

The bathymetry is either a sloping beach or flat.
The BCs are outflow on the side of wave's direction and wall BC on the side of the initial dam break gradient.
"""
from __future__ import absolute_import
from clawpack import riemann
from clawpack import pyclaw
from clawpack.riemann.shallow_roe_with_efix_2D_constants import depth, x_momentum, y_momentum, num_eqn
import numpy as np
from clawpack.pyclaw.plot import plot
import matplotlib.pyplot as plt
from clawpack.visclaw.colormaps import make_colormap
from matplotlib import cm

## barrier / grid information:
bar_height = 1.5
my =  300 # num x cells
mx = my  # num y cells
mbc = 2 # num ghost cells
front = int(my/10)

def bathymetry(x,y):
    # sloping beach
    r = np.zeros((mx,my))
    for i in range(mx):
        for j in range(my):
            r[i,j] = 0.01*j -0.01*i
            r[i,j] = 0.01*j - 0.01*i
    r = r-2

    # making flat near barrier
    r[range(mx-1),range(1,my)] = r[mx-1,my-1]
    r[range(1,mx),range(my-1)] = r[mx-1,my-1]
    r[range(mx-2),range(2,my)] = r[mx-1,my-1]
    r[range(mx-2),range(2,my)] = r[mx-1,my-1]
    r[range(2,mx),range(my-2)] = r[mx-1,my-1]
    r[range(2,mx),range(my-2)] = r[mx-1,my-1]

    # Flat bathymetry (uncomment the following to change to this bathymetric setting)
    r = -2

    return r

# the setup of example
def setup(kernel_language='Fortran', solver_type='classic', use_petsc=False,
          outdir='./_output'):

    solver = pyclaw.ClawSolver2D(riemann.sw_aug_2D)#hallow_bathymetry_fwave_2D)##
    solver.dimensional_split = True # No transverse solver available
    solver.order = 1
    # solver.limiters =4# pyclaw.limiters.tvd.minmod_limiter(r,cfl)
    solver.bc_lower[0] = pyclaw.BC.wall
    solver.bc_upper[0] = pyclaw.BC.wall
    solver.bc_lower[1] = pyclaw.BC.wall
    solver.bc_upper[1] = pyclaw.BC.extrap

    solver.aux_bc_lower[0] = pyclaw.BC.wall
    solver.aux_bc_upper[0] = pyclaw.BC.wall
    solver.aux_bc_lower[1] = pyclaw.BC.wall
    solver.aux_bc_upper[1] = pyclaw.BC.extrap

    x = pyclaw.Dimension(0.,1.,mx,name='x')
    y = pyclaw.Dimension(0.,1.,my,name='y')
    domain = pyclaw.Domain([x,y])
    state = pyclaw.State(domain,num_eqn,num_aux=1) # use aux array to do the small cell business. the second and third aux var are small cell values
    X, Y = state.p_centers
    # solver.limiters = pyclaw.limiters.tvd.minmod_limiter(state,cfl)

    state.aux[0,:,:] =bathymetry(X,Y)
    state.auxu[0,:,:] =bathymetry(X,Y) # bathymetry value for upper side of barrier (currently obsolete as they are same as the lower side bathymetry)

    # lower small cell state values in the cut cells
    state.q[depth,:,:] =1.2 # -0.55-state.aux[0,:,:]
    m1 = int(mx/2)
    m2 = int(mx/4)
    m3 = mx- m2
    state.q[depth,:,:front] += 0.8
    state.q[x_momentum,:,:] = 0.
    state.q[y_momentum,:,:] = 0.

    # upper small-cell-replaced values in cut cells
    state.q2[depth,:,:] =1.2 # -0.55 -state.auxu[0,:,:]
    state.q2[depth,:,:front] += 0.8
    state.q2[x_momentum,:,:] = 0.
    state.q2[y_momentum,:,:] = 0.

    state.problem_data['grav'] = 1.0
    state.problem_data['dry_tolerance'] = 1.e-8
    state.problem_data['sea_level'] = 0.
    state.problem_data['bar_ht'] = bar_height
    state.problem_data['method'] = 'hbox'
    state.problem_data['BC2'] = "we"  # Boundary conditions for array q2, the upper half grid. This sets BC at left and top side.

    claw = pyclaw.Controller()
    claw.tfinal = 1.4
    claw.solution = pyclaw.Solution(state,domain)
    claw.solver = solver
    claw.num_output_times = 14
    claw.output_style =1

    # Gauge for outputting 1D slice:
    state.grid.add_gauges([(0.5,0.8),(0.5,0.39)])#(0.125,0.875),(0.24,0.76),(0.375,0.625),(0.5,0.5),(0.625,0.375),(0.75,0.25),(0.875,0.125)])
    solver.compute_gauge_values = gauge_height
    state.keep_gauges = True

    claw.setplot = setplot
    claw.keep_copy = True

    claw.write_aux_always = True

    claw.run()

    plot(setplot=setplot,outdir='./_output',plotdir='./_plots_slanted15',iplot=False,htmlplot=True)
    #return claw

# functions for plotting:
def barrier_draw(current_data):
    x_1 = 1.0 #0.995#
    x_0 = 0.0
    y_0 =  0.3#0.393#0.005#0.345
    y_1 = 0.6533#639#0.776 #.755# 1.0#
    axis = plt.gca()
    axis.plot([x_0,x_1],[y_0,y_1],'chartreuse',linewidth=1.5)
    return

def barrier_draw_1d(current_data):
    x_1 = 1.0
    x_0 = 0.0
    y_1 = 1.0
    y_0 = 0.0
    bar_loc = 0.50
    b = bathymetry(current_data.x,current_data.y)
    aux_wall = -2
    axis = plt.gca()
    axis.plot([0.45,0.45],[aux_wall,aux_wall+bar_height],'g',linewidth=1.5)
    # axis.plot(np.linspace(x_0,x_1,mx),b[range(mx),range(my-1,-1,-1)],'k-') # comment this out for flat bathy

    return

def surface_height(current_data):
    h = current_data.q[0,:,:]
    b = bathymetry(current_data.x,current_data.y)
    return h-1.2

def gauge_height(q,aux):
    h = q[0]
    return h

def gauge_spots(current_data):
    barrier_draw(current_data)
    gauge_points = [(0.5,0.8),(0.5,0.39)]#[(0.125,0.875),(0.24,0.76),(0.375,0.625),(0.5,0.5),(0.625,0.375),(0.75,0.25),(0.875,0.125)]
    axis = plt.gca()
    # x_0 = 0.0 ; x_1 = 1.0
    # axis.plot([x_0,x_1],[x_0,x_1],'g',linewidth=1.5)
    for i in range(len(gauge_points)):
        axis.plot(gauge_points[i][0],gauge_points[i][1],'k*')
        axis.annotate(str(i+1),(gauge_points[i][0],gauge_points[i][1]))
    return

def momentum_x(current_data):
    hu = current_data.q[1,:,1]
    return hu
def height_x(current_data):
    h = current_data.q[0,10,:]
    b = bathymetry(current_data.x,current_data.y)
    # b2 = b[range(mx),range(my-1,-1,-1)] # comment this out for flat bathymetry setting
    return h+b#2 # h+b for flat bathymetry setting


def setplot(plotdata):
    from clawpack.visclaw import colormaps

    plotdata.clearfigures()  # clear any old figures,axes,items data

    # Figure for q[0]
    plotfigure = plotdata.new_plotfigure(name='Water height', figno=0)
    plotfigure.kwargs = {'facecolor': '#FFFFFF'}

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'Surface height with barrier height ' + str(bar_height)
    plotaxes.scaled = False

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_contourf')
    plotitem.plot_var =  surface_height
    plotitem.add_colorbar = True
    plotitem.colorbar_ticks = [-1.0,-0.8,-0.6,-0.4,-0.2,0.0,0.2,0.4,0.6,0.8,1.0]
    cmap2 = cm.get_cmap('bwr')

    plotitem.fill_cmap = cmap2
    plotitem.contour_min = -1.0
    plotitem.contour_max = 1.0
    # plotitem.contour_colors = 'b'
    plotaxes.afteraxes = barrier_draw
    plotaxes.xlimits = [0,1]
    plotaxes.ylimits = [0,1]
    plotaxes.afteraxes = gauge_spots # for gauge points

    plotfigure = plotdata.new_plotfigure(name="Momentum",figno=1)
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = "Momentum in x"
    plotitem = plotaxes.new_plotitem(plot_type='1d')
    plotitem.plot_var = momentum_x


    plotfigure = plotdata.new_plotfigure(name="Height",figno=2)
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = "Diagonal cross section"
    plotaxes.xlimits = [0.,1.]
    plotaxes.ylimits = [-1.1,0.1]
    plotitem = plotaxes.new_plotitem(plot_type='1d')
    plotitem.plot_var = height_x
    # plotaxes.afteraxes = barrier_draw_1d

    return plotdata

if __name__=="__main__":
    setup()

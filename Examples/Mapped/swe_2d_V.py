#!/usr/bin/env python
# encoding: utf-8
r"""
Two-dimensional SWE on a mapped grid with V shaped barrier
===============================================================

Solve the shallow water equations in 2D with a V zero width barrier:

"""
from __future__ import absolute_import
import numpy as np
from six.moves import range
import SWE_Vmap
from clawpack.pyclaw.plot import plot
import matplotlib.pyplot as plt
from matplotlib import cm


m1 = 0.616  # slope of one edge of the V barrier
def m(y):
    # slope function with respect to position y under the V barrier
    y_low = y<=0.72
    y_abo = y>0.72
    slope = ((m1/0.72 * y)*y_low + (-m1/0.28 * y + m1 + (m1/0.28)*0.72)*y_abo)*(y<0.9) + 0*(y>0.9)
    return slope
def Vmap(xc,yc):
    # the V mapping function
    xp = xc.copy()
    yp = yc.copy()
    yp = yc - np.absolute(np.multiply(m(yc),np.minimum(1-xc,xc)))
    return xp,yp

def compute_geometry(grid):
    r"""Computes
        a_x  normal vector x comp
        a_y   "      "     y comp
        length_ratio_left
        b_x  normal vector x comp
        b_y  normal vector y comp
        length_ratio_bottom
        cell_area
    """

    dx, dy = grid.delta
    area_min = 1.e6
    area_max = 0.0

    x_corners, y_corners = grid.p_nodes

    lower_left_y, lower_left_x = y_corners[:-1,:-1], x_corners[:-1,:-1]
    upper_left_y, upper_left_x = y_corners[:-1,1: ], x_corners[:-1,1: ]
    lower_right_y, lower_right_x = y_corners[1:,:-1], x_corners[1:,:-1]
    upper_right_y, upper_right_x = y_corners[1:,1: ], x_corners[1:,1: ]

    a_x =   upper_left_y - lower_left_y  #upper left and lower left
    a_y = -(upper_left_x - lower_left_x)
    anorm = np.sqrt(a_x**2 + a_y**2)
    a_x, a_y = a_x/anorm, a_y/anorm
    length_ratio_left = anorm/dy

    b_x = -(lower_right_y - lower_left_y)  #lower right and lower left
    b_y =   lower_right_x - lower_left_x
    bnorm = np.sqrt(b_x**2 + b_y**2)
    b_x, b_y = b_x/bnorm, b_y/bnorm
    length_ratio_bottom = bnorm/dx

    area = 0*grid.c_centers[0]
    area += 0.5 * (lower_left_y+upper_left_y)*(upper_left_x-lower_left_x)
    area += 0.5 * (upper_left_y+upper_right_y)*(upper_right_x-upper_left_x)
    area += 0.5 * (upper_right_y+lower_right_y)*(lower_right_x-upper_right_x)
    area += 0.5 * (lower_right_y+lower_left_y)*(lower_left_x-lower_right_x)
    area = area/(dx*dy)
    area_min = min(area_min, np.min(area))
    area_max = max(area_max, np.max(area))

    return a_x, a_y, length_ratio_left, b_x, b_y, length_ratio_bottom, area

def incoming_square_wave(state,dim,t,qbc,auxbc,num_ghost):
    """
    Incoming square wave at left boundary.
    """
    if t<0.05:
        s = 1.0
    else:
        s = 0.0

    for i in range(num_ghost):
        reflect_ind = 2*num_ghost - i - 1
        alpha = 0#auxbc[0,:,i]
        beta  = -1#auxbc[1,:,i]
        u_normal = alpha*qbc[1,:,reflect_ind] + beta*qbc[2,:,reflect_ind]
        u_tangential = -beta*qbc[1,:,reflect_ind] + alpha*qbc[2,:,reflect_ind]
        u_normal = 2.0*s - u_normal
        qbc[0,:,i] = qbc[0,:,reflect_ind]
        qbc[1,:,i] = alpha*u_normal - beta*u_tangential
        qbc[2,:,i] = beta*u_normal + alpha*u_tangential


def setup(kernel_language='Fortran', use_petsc=False, outdir='./_output',
          solver_type='classic', time_integrator='SSP104',
          num_output_times=14, disable_output=False, num_cells=200):
    from clawpack import riemann

    if use_petsc:
        import clawpack.petclaw as pyclaw
    else:
        from clawpack import pyclaw
    # riemann_solver = riemann.shallow_bathymetry_fwave_2D

    if solver_type=='classic':
        solver=pyclaw.ClawSolver2D(SWE_Vmap)
        solver.dimensional_split=True
        solver.order = 1
        # solver.limiters = pyclaw.limiters.tvd.MC
    elif solver_type=='sharpclaw':
        solver=pyclaw.SharpClawSolver2D(riemann_solver)
        solver.time_integrator=time_integrator

    solver.bc_lower[0]=pyclaw.BC.wall
    solver.bc_upper[0]=pyclaw.BC.wall
    solver.bc_lower[1]=pyclaw.BC.wall
    solver.bc_upper[1]=pyclaw.BC.wall
    # solver.user_bc_upper = incoming_square_wave


    solver.aux_bc_lower[0]=pyclaw.BC.wall
    solver.aux_bc_upper[0]=pyclaw.BC.wall
    solver.aux_bc_lower[1]=pyclaw.BC.wall
    solver.aux_bc_upper[1]=pyclaw.BC.wall

    x = pyclaw.Dimension(0.,1.0,num_cells,name='x')
    y = pyclaw.Dimension(0.,1.0,num_cells,name='y')
    domain = pyclaw.Domain([x,y])

    num_eqn = 3
    num_waves = 3
    num_aux = 8 # geometry (7), impedance, sound speed
    state = pyclaw.State(domain,num_eqn,num_aux)
    state.problem_data['grav'] = 1.0
    state.grid.mapc2p = Vmap

    a_x, a_y, length_left, b_x, b_y, length_bottom, area = compute_geometry(state.grid)
    state.aux[0,:,:] = -2.0 # bathymetry
    state.aux[1,:,:] = a_x
    state.aux[2,:,:] = a_y
    state.aux[3,:,:] = length_left
    state.aux[5,:,:] = b_x
    state.aux[4,:,:] = b_y
    state.aux[6,:,:] = length_bottom
    state.aux[7,:,:] = area
    state.index_capa = 7 # aux[7,:,:] holds the capacity function

    grid = state.grid
    xp, yp = grid.p_centers

    # for i, circle in enumerate(circles):
    #     # Set impedance and sound speed in each inclusion
    #     radius = circle[0][0]
    #     x0, y0 = circle[1]
    #     distance = np.sqrt( (xp-x0)**2 + (yp-y0)**2 )
    #     in_circle = np.where(distance <= radius)
    #     state.aux[7][in_circle] = impedance[i]
    #     state.aux[8][in_circle] = sound_speed[i]

    # Set initial condition
    state.q[0,:,:] = 1.2 + (yp>0.9)*0.8*np.ones(state.q[0,:,:].shape)# - 2.0 *np.ones(state.q[0,:,:].shape)*(yp>1.0)

    state.q[1,:,:] = 0.
    state.q[2,:,:] = 0.

    claw = pyclaw.Controller()
    claw.keep_copy = True
    if disable_output:
        claw.output_format = None
    claw.solution = pyclaw.Solution(state,domain)
    claw.solver = solver
    claw.outdir = outdir
    claw.tfinal = 0.9
    claw.num_output_times = num_output_times
    claw.write_aux_init = True
    claw.setplot = setplot
    if use_petsc:
        claw.output_options = {'format':'binary'}
    claw.run()
    plot(setplot=setplot,outdir='./_output',plotdir='./_plots_mapped_V',iplot=False,htmlplot=True)


    # return claw
def surface_height(current_data):
    h = current_data.q[0,:,:]

    return h-1.2

def setplot(plotdata):
    """
    Plot solution using VisClaw.

    This example shows how to mark an internal boundary on a 2D plot.
    """

    from clawpack.visclaw import colormaps

    plotdata.clearfigures()  # clear any old figures,axes,items data
    plotdata.mapc2p = Vmap

    # Figure for pressure
    plotfigure = plotdata.new_plotfigure(name='Height', figno=0)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'Height'
    plotaxes.scaled = True      # so aspect ratio is 1
    plotaxes.afteraxes = barrier_draw



    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_contourf')
    plotitem.plot_var = surface_height
    plotitem.add_colorbar = True
    plotitem.colorbar_ticks = [-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5]
    cmap2 = cm.get_cmap('bwr')
    plotitem.fill_cmap = cmap2
    plotitem.contour_min = -0.5
    plotitem.contour_max = 0.5
    plotaxes.xlimits = [0,1]
    plotaxes.ylimits = [0,1]

    # Figure for x-velocity plot
    plotfigure = plotdata.new_plotfigure(name='x-momentum', figno=1)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'hu'
    plotaxes.afteraxes = barrier_draw

    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = 1
    plotitem.pcolor_cmap = 'RdBu'
    plotitem.add_colorbar = True
    plotitem.pcolor_cmin = -1.0
    plotitem.pcolor_cmax=   1.0

    return plotdata

def barrier_draw(current_data):
    # x_1 = 1.0 #0.995#
    # x_0 = 0.0
    # y_1 =  0.393#0.005#0.345
    # y_0 = 0.776 #.755# 1.0#
    x_0 = 0.0
    y_0 = .72 #0.719
    x_e = 0.5
    y_e = 0.412
    x_1 = x_e
    x_2 = 1
    y_1 = y_e
    y_2 = y_0
    axis = plt.gca()
    axis.plot([x_0,x_e],[y_0,y_e],'chartreuse',linewidth=1.5)
    axis.plot([x_1,x_2],[y_1,y_2],'chartreuse',linewidth=1.5)
    return


if __name__=="__main__":
    # import sys
    # from clawpack.pyclaw.util import run_app_from_main
    # output = run_app_from_main(setup,setplot)
    setup()

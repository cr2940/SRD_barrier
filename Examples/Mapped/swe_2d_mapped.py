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


# m1 = 0.3639  # slope of one edge of the V barrier
bar_loc = 0.65333# 0.65333
m1 = (bar_loc-0.3)
straight = 0.0#0.1   LH: (0,0) -> (0.5,0.4)  y=(0.4)/(0.5)x  UH: (0.5,0.4) -> (1,0)  y = 0.4-(0.4/0.5)(x-0.5)
mx= 450
ycen = (bar_loc+0.3)/2
print(int(mx*bar_loc))
def m(y):
    # slope function with respect to position y under the V barrier
    y_low = y<=bar_loc
    y_abo = y>bar_loc
    y_01 = y>=straight
    y_01l = y<straight
    # slope = ((-m1/(straight-bar_loc) * 0.72 + m1 + (m1/(straight-bar_loc))*bar_loc)*y_low + (-m1/(straight-bar_loc) * y + m1 + (m1/(straight-bar_loc))*bar_loc)*y_abo)#*y_09 + 0*y_09p
    # slope = ((m1/bar_loc * y)*y_low + (m1)*y_abo)#*y_09 + 0*y_09p

    slope = ((m1/(bar_loc-straight) * y - straight*m1/(bar_loc-straight))*y_low + (-m1/(1-bar_loc) * y + (m1/(1-bar_loc)))*y_abo)*y_01 + 0*y_01l
    return slope
def Vmap(xc,yc):
    # the V mapping function
    # [0.3,0.7] --> [0.7,1.0]
    # [0,0.3] <-- [0,0.7]
    xp = xc.copy()
    y_star = m1*xc + 0.3 # the point where things shift from being crunched to even
    dist_to_top = 1-y_star
    abo = yc.copy()# > y_star
    bel = yc.copy()# <= y_star
    bel[:,:] = 0
    bel[:,:int(mx*ycen)] = 1# el[1,:]
    abo[:,:] = 1
    abo[:,:int(mx*ycen)] = 0
    # abo[0,:] = abo[1,:]
    yp = (np.multiply(y_star/ycen, yc) )* bel + (np.multiply(dist_to_top/ycen, yc) + y_star - (dist_to_top))*abo #(0.7-0.4*xc)/0.5 * yc /(0.7/0.5) * abo #

    # print("xc: ", xc)
    # print("here 1:" , np.multiply(y_star/0.5, yc) )
    # print("here 2:", (np.multiply(dist_to_top/0.5, yc) + y_star - (dist_to_top)))
    # print("yc: ", yc)

    # print("Y_STAR: ", y_star)
    # print("dist to top ", dist_to_top)
    # print("above:: ", abo )
    # print("below:: ",bel)
    # print("Yp: ", yp)

    # yp = yc - np.absolute(np.multiply(m(yc),1-xc))#*yc -np.absolute(np.multiply(m(yc),xc))*yc>=0.5.all()
    return xp,yp
def Vmap_inv(xp,yp):
    xc = xp.copy()
    yc = yp -np.absolute(np.multiply(m(yp),np.minimum(1-xp,xp)))
    return xc,yc

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
    # print("Xcorn",x_corners,"Ycorn",y_corners)

    lower_left_y, lower_left_x = y_corners[:-1,:-1], x_corners[:-1,:-1]
    upper_left_y, upper_left_x = y_corners[:-1,1: ], x_corners[:-1,1: ]
    lower_right_y, lower_right_x = y_corners[1:,:-1], x_corners[1:,:-1]
    upper_right_y, upper_right_x = y_corners[1:,1: ], x_corners[1:,1: ]

    x_center = 0.5*lower_right_x + 0.5*lower_left_x
    y_center = 0.25*(lower_right_y + lower_left_y + upper_right_y + upper_left_y)

    y_edge_center = 0.5*(lower_left_y + upper_left_y)
    x_edge_center = 0.5*(lower_left_x + lower_right_x)

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


    print(np.amin(upper_right_y-lower_right_y)/dy)
    return a_x, a_y, length_ratio_left, b_x, b_y, length_ratio_bottom, area

def wave_maker_bc(state,dim,t,qbc,auxbc,num_ghost):
    "Generate waves at left boundary as if there were a moving wall there."
    amplitude = 0.1  # Height of incoming wave
    t_bdy = 0.1       # Stop sending in waves at this time
    # if dim.on_upper_boundary:

    t=state.t
    if t <= t_bdy:
        qbc[0,:,-2]=2.0
        qbc[0,:,-1] = 2.0
        #     vwall = -amplitude*(np.sin(t*np.pi/1.5))
    else:
        qbc[0,:,-2] = qbc[0,:,-3]
        qbc[0,:,-1] = qbc[0,:,-2]

        #     vwall=0.
        # for ibc in range(num_ghost-1):
        #     qbc[2,:,-(num_ghost-ibc-1)] = 2*vwall + qbc[2,:,-(num_ghost+ibc)]


def setup(kernel_language='Fortran', use_petsc=False, outdir='./_output_Lmapped',
          solver_type='classic', time_integrator='SSP104',
          num_output_times=14, disable_output=False, num_cells=mx):
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
    solver.bc_upper[1]=pyclaw.BC.extrap
    # solver.user_bc_upper = wave_maker_bc


    solver.aux_bc_lower[0]=pyclaw.BC.wall
    solver.aux_bc_upper[0]=pyclaw.BC.wall
    solver.aux_bc_lower[1]=pyclaw.BC.wall
    solver.aux_bc_upper[1]=pyclaw.BC.wall

    x = pyclaw.Dimension(0.,1.0,num_cells,name='x')
    y = pyclaw.Dimension(0.0,1.0,num_cells,name='y')
    domain = pyclaw.Domain([x,y])

    num_eqn = 3
    num_aux = 8 # geometry (7), impedance, sound speed
    bar_ht = 1.5
    state = pyclaw.State(domain,num_eqn,num_aux)
    state.problem_data['grav'] = 1.0
    state.problem_data['wall_height']= bar_ht
    state.grid.mapc2p = Vmap


    grid = state.grid
    xp, yp = grid.p_centers
    a_x, a_y, length_left, b_x, b_y, length_bottom, area = compute_geometry(state.grid)
    state.aux[0,:,:] = -2.0 #+ (yp>1.0) *5*np.ones(state.q[0,:,:].shape)# bathymetry
    state.aux[1,:,:] = a_x
    state.aux[2,:,:] = a_y
    state.aux[3,:,:] = length_left
    state.aux[4,:,:] = b_y
    state.aux[5,:,:] = b_x
    state.aux[6,:,:] = length_bottom
    state.aux[7,:,:] = area
    state.index_capa = 7 # aux[7,:,:] holds the capacity function


    # for i, circle in enumerate(circles):
    #     # Set impedance and sound speed in each inclusion
    #     radius = circle[0][0]
    #     x0, y0 = circle[1]
    #     distance = np.sqrt( (xp-x0)**2 + (yp-y0)**2 )
    #     in_circle = np.where(distance <= radius)
    #     state.aux[7][in_circle] = impedance[i]
    #     state.aux[8][in_circle] = sound_speed[i]


    # Set initial condition
    state.q[0,:,:] = 1.2 + (yp<0.1)*0.8*np.ones(state.q[0,:,:].shape)# - 2.0 *np.ones(state.q[0,:,:].shape)*(yp>1.0)

    # state.q[0,:,180:] += 1.0
    state.q[1,:,:] = 0.
    state.q[2,:,:] = 0.

    state.grid.add_gauges([(0.5,0.785),(0.5,0.789),(0.5,0.788),(0.5,0.783),(0.5,0.75),(0.5,0.74),(0.5,0.73),(0.5,0.72),(0.5,0.41),(0.5,0.42),(0.5,0.43),(0.5,0.45),(0.5,0.35),(0.5,0.25),(0.5,0.15)])#,(0.5,0.82),(0.5,0.83),(0.5,0.87),(0.5,0.869),(0.5,0.861),(0.5,0.545),(0.5,0.863),(0.5,0.864),(0.5,0.865),(0.5,0.866),(0.5,0.867),(0.5,0.868),(0.5,0.94),(0.5,0.95),(0.5,0.96),(0.5,0.97),(0.5,0.98)])#,(0.25,0.375),(0.25,0.73),(0.75,0.73),(0.75,0.375)])#,(0.24,0.76),(0.375,0.625),(0.5,0.5),(0.625,0.375),(0.75,0.25),(0.875,0.125)])
    solver.compute_gauge_values = gauge_height
    state.keep_gauges = True

    claw = pyclaw.Controller()
    claw.keep_copy = True
    if disable_output:
        claw.output_format = None
    claw.solution = pyclaw.Solution(state,domain)
    claw.solver = solver
    claw.outdir = outdir
    claw.tfinal = 1.4
    claw.num_output_times = num_output_times
    claw.write_aux_init = True
    claw.setplot = setplot
    if use_petsc:
        claw.output_options = {'format':'binary'}
    claw.run()
    plot(setplot=setplot,outdir='./_output_Lmapped',plotdir='./_plots_mapped',iplot=False,htmlplot=True)


    # return claw
def surface_height(current_data):
    h = current_data.q[0,:,:]

    return h-1.2
def gauge_height(q,aux):
    h = q[0]
    return h

def gauge_spots(current_data):
    gauge_points = [(0.5,0.8),(0.5,0.39)]#[(0.5,0.4),(0.5,0.41),(0.5,0.42),(0.5,0.43),(0.5,0.44),(0.5,0.45),(0.5,0.46),(0.5,0.47),(0.5,0.48),(0.5,0.49),(0.5,0.5)]#(0.5,0.8),(0.5,0.4)] #[(0.5,0.5)] #
    axis = plt.gca()
    barrier_draw(current_data)
    # x_0 = 0.0 ; x_1 = 1.0
    # axis.plot([x_0,x_1],[x_0,x_1],'g',linewidth=1.5)
    for i in range(len(gauge_points)):
        axis.plot(gauge_points[i][0],gauge_points[i][1],'k*')
        axis.annotate(str(i+1),(gauge_points[i][0],gauge_points[i][1]))
    return

def setplot(plotdata):
    """
    Plot solution using VisClaw.

    This example shows how to mark an internal boundary on a 2D plot.
    """

    from clawpack.visclaw import colormaps


    # x = pyclaw.Dimension(0.,1.0,num_cells,name='x')
    # y = pyclaw.Dimension(0.,1.0,num_cells,name='y')
    # domain = pyclaw.Domain([x,y])

    # num_eqn = 3
    # num_waves = 2
    # num_aux = 8 # geometry (7), impedance, sound speed
    # state = pyclaw.State(domain,num_eqn,num_aux)
    # state.problem_data['grav'] = 1.0
    # state.grid.mapc2p = Vmap

    # a_x, a_y, length_left, b_x, b_y, length_bottom, area= compute_geometry(state.grid)
    # grid = state.grid
    # xp, yp = grid.p_centers

    plotdata.clearfigures()  # clear any old figures,axes,items data
    plotdata.mapc2p = Vmap

    # Figure for pressure
    plotfigure = plotdata.new_plotfigure(name='Height', figno=0)
    plotfigure.kwargs = { 'facecolor': '#FFFFFF'}


    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'Height'
    plotaxes.scaled = False      # so aspect ratio is 1



    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_contourf')
    plotitem.plot_var = surface_height
    plotitem.add_colorbar = True
    plotitem.colorbar_ticks = [-1.0,-0.8,-0.6,-0.4,-0.2,0.0,0.2,0.4,0.6,0.8,1.0]#[-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5]
    cmap2 = cm.get_cmap('bwr')
    plotitem.fill_cmap = cmap2
    plotitem.contour_min = -1
    plotitem.contour_max = 1
    plotaxes.xlimits = [0,1]
    plotaxes.ylimits = [0,1]
    plotaxes.afteraxes = gauge_spots # for gauge points


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

    # Figure for y-velocity plot
    plotfigure = plotdata.new_plotfigure(name='y-momentum', figno=2)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'hv'
    plotaxes.afteraxes = barrier_draw

    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = 2
    plotitem.pcolor_cmap = 'RdBu'
    plotitem.add_colorbar = True
    plotitem.pcolor_cmin = -1.0
    plotitem.pcolor_cmax=   1.0

    # mapped grid
    plotfigure = plotdata.new_plotfigure(name='mapped grid', figno=3)
    plotfigure.kwargs = { 'facecolor': '#FFFFFF'}

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'xp,yp'
    plotaxes.afteraxes = barrier_draw
    plotaxes.title_with_t = False


    plotitem = plotaxes.new_plotitem(plot_type='2d_patch')
    # plotitem.plot_var = 1
    # plotitem.pcolor_cmap = 'RdBu'
    # plotitem.add_colorbar = True
    # plotitem.pcolor_cmin = -1.0
    # plotitem.pcolor_cmax=   1.0

    return plotdata

def barrier_draw(current_data):
    # x_1 = 1.0 #0.995#
    # x_0 = 0.0
    # y_1 =  0.393#0.005#0.345
    # y_0 = 0.776 #.755# 1.0#
    x_0 = 0.0
    y_0 = .3 #0.719
    x_e = 1.0
    y_e = bar_loc#0.7#
    # x_1 = x_e
    # x_2 = 1
    # y_1 = y_e
    # y_2 = y_0
    axis = plt.gca()
    axis.plot([x_0,x_e],[y_0,y_e],'chartreuse',linewidth=1.5)
    # axis.plot([x_1,x_2],[y_1,y_2],'chartreuse',linewidth=1.5)
    return


if __name__=="__main__":
    # import sys
    # from clawpack.pyclaw.util import run_app_from_main
    # output = run_app_from_main(setup,setplot)
    setup()

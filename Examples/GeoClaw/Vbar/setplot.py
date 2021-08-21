
"""
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.

"""


from __future__ import absolute_import
from __future__ import print_function

from clawpack.visclaw import gaugetools

from clawpack.visclaw import particle_tools
from clawpack.visclaw import legend_tools
from clawpack.visclaw import colormaps
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy
green = [0.0,1.0,0.0];
dark_green = [0.1,0.4,0.0];
light_green = [0.8,1.0,0.5];
tan = [0.9,0.8,0.2];
white = [1.0,1.0,1.0]
red = [1.0,0.0,0.0]
blue = [0.0,0.0,1.0];
land2_colormap = colormaps.make_colormap({-2:dark_green,
                                         -1:green,
                                          0:light_green,
                                          5:tan})
#--------------------------
def setplot(plotdata=None):
#--------------------------

    """
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of pyclaw.plotters.data.ClawPlotData.
    Output: a modified version of plotdata.

    """

    if plotdata is None:
        from clawpack.visclaw.data import ClawPlotData
        plotdata = ClawPlotData()


    from clawpack.visclaw import colormaps, geoplot

    plotdata.clearfigures()  # clear any old figures,axes,items data

    plotdata.format = 'ascii'                # Format of output

    print('Reading all gauges...')
    gauge_solutions = particle_tools.read_gauges(gaugenos='all',
                                                 outdir=plotdata.outdir)
    #
    # gaugenos_lagrangian = [k for k in gauge_solutions.keys() \
    #             if gauge_solutions[k].gtype=='lagrangian']
    gaugenos_stationary = [k for k in gauge_solutions.keys() \
                if gauge_solutions[k].gtype=='stationary']

    #print('+++ gaugenos_lagrangian: ',gaugenos_lagrangian)

    def add_particles(current_data):
        t = current_data.t

        # # plot recent path:
        # t_path_length = 0.5   # length of path trailing particle
        # kwargs_plot_path = {'linewidth':1, 'color':'k'}
        # particle_tools.plot_paths(gauge_solutions,
        #                           t1=t-t_path_length, t2=t,
        #                           gaugenos=gaugenos_lagrangian,
        #                           kwargs_plot=kwargs_plot_path)

        # # plot current location:
        # kwargs_plot_point = {'marker':'o','markersize':3,'color':'k'}
        # particle_tools.plot_particles(gauge_solutions, t,
        #                               gaugenos=gaugenos_lagrangian,
        #                               kwargs_plot=kwargs_plot_point)
        #
        # # plot any stationary gauges:
        gaugetools.plot_gauge_locations(current_data.plotdata, \
             gaugenos=gaugenos_stationary, format_string='kx', add_labels=False)
        kwargs={'loc':'upper left'}
        legend_tools.add_legend(['Lagrangian particle','Stationary gauge'],
                linestyles=['',''], markers=['o','x'],
                loc='lower right', framealpha=0.5, fontsize=10)

    def land(current_data):
       """
       Return a masked array containing the surface elevation only in dry cells.
       """
       drytol = 1.e-3
       q = current_data.q
       aux = current_data.aux
       h = q[0,:,:]
       beta = aux[0,:,:]
       land = numpy.ma.masked_where(h>drytol, beta)
       return land
    def water(current_data):
        from pylab import sqrt, where, zeros
        from numpy.ma import masked_where, allequal
        q = current_data.q
        h = q[0,:,:]
        b= q[3,:,:]-h
        hs = sqrt(q[1,:,:]**2 + q[2,:,:]**2)
        where_hpos = (h > 1e-3)
        s = zeros(h.shape)
        s[where_hpos] = hs[where_hpos]/h[where_hpos]
        s = masked_where(h<1e-3, h) # if you want 0's masked out
        #s = s * 1.94384  # convert to knots
        return h-1.2

    speed_cmap = colormaps.make_colormap({0:[0,1,1], 0.5:[1,1,0], 1:[1,0,0]})

    def gauge_spots(current_data):
        gauge_points = [(0.25,0.6),(0.75,0.6),(0.25,0.3),(0.75,0.3)] #[(0.5,0.5)]]
        axis = plt.gca()
        barrier_draw(current_data)
        # x_0 = 0.0 ; x_1 = 1.0
        # axis.plot([x_0,x_1],[x_0,x_1],'g',linewidth=1.5)
        for i in range(len(gauge_points)):
            axis.plot(gauge_points[i][0],gauge_points[i][1],'k*')
            axis.annotate(str(i+1),(gauge_points[i][0],gauge_points[i][1]))
        return
    def barrier_draw(current_data):
        # x_1 = 1.0 #0.995#
        # x_0 = 0.0
        # y_1 =  0.393#0.005#0.345
        # y_0 = 0.776 #.755# 1.0#
        x_0 = 0.0
        y_0 = .72
        x_e = 0.5
        y_e = 0.412
        x_1 = x_e
        x_2 = 1
        y_1 = y_e
        y_2 = y_0
        axis = plt.gca()
        axis.plot([x_0,x_e],[y_0,y_e],'green',linewidth=1.5)
        axis.plot([x_1,x_2],[y_1,y_2],'green',linewidth=1.5)
        return
    #-----------------------------------------
    # Figure for pcolor plot
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='pcolor', figno=0)
    plotfigure.kwargs = {'facecolor': '#FFFFFF'}

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('pcolor')
    plotaxes.title = 'Surface height'
    plotaxes.scaled = False
    plotaxes.xlimits = [0,1]
    plotaxes.ylimits = [0,1]

    # Water
    plotitem = plotaxes.new_plotitem(plot_type='2d_contourf')
    plotitem.plot_var = water
    cmap2 = cm.get_cmap('bwr')
    plotitem.fill_cmap = cmap2
    plotitem.add_colorbar = True
    # plotitem.pcolor_cmap = "Blues"
    plotitem.colorbar_ticks = [-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5]
    plotitem.contour_min = -0.5
    plotitem.contour_max = 0.5
    plotaxes.afteraxes = gauge_spots


    # plotitem.pcolor_cmin = 0.
    # plotitem.pcolor_cmax = 10
    # plotitem.colorbar_label = 'm'
    # plotitem.amr_celledges_show = [0,0,0]
    # plotitem.amr_patchedges_show = [1]
    # plotitem.amr_patchedges_color = ['m','g','w']

    # Land
    # plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    # #plotitem.show = False
    # plotitem.plot_var = land
    # plotitem.pcolor_cmap = land2_colormap
    # plotitem.pcolor_cmin = -1
    # plotitem.pcolor_cmax = 6
    # plotitem.add_colorbar = False
    # plotitem.amr_celledges_show = [0,0,0]


    # Add contour lines of topography:
    plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
    plotitem.show = False
    plotitem.plot_var = geoplot.topo
    from numpy import arange, linspace
    plotitem.contour_levels = arange(-75,75,10)
    #plotitem.contour_nlevels = 10
    plotitem.amr_contour_colors = ['g']  # color on each level
    plotitem.kwargs = {'linestyles':'solid'}
    # plotitem.amr_contour_show = [1,1,1]  # show contours only on finest level
    # plotitem.celledges_show = 0
    plotitem.patchedges_show = 0


    #-----------------------------------------
    # Figures for gauges
    #-----------------------------------------
    # plotfigure = plotdata.new_plotfigure(name='Surface', figno=300, \
    #                 type='each_gauge')
    #
    # plotfigure.clf_each_gauge = True
    #
    # # Set up for axes in this figure:
    # plotaxes = plotfigure.new_plotaxes()
    # plotaxes.xlimits = 'auto'
    # plotaxes.ylimits = [-25,50]
    # plotaxes.title = 'Surface'  # reset in fix_gauge
    #
    # # Plot surface as blue curve:
    # plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    # plotitem.plot_var = 3
    # plotitem.plotstyle = 'b-'
    #
    # def fix_gauge(current_data):
    #     from pylab import plot, title
    #     t = current_data.t
    #     plot(t, 0*t, 'k')
    #     gaugeno = current_data.gaugeno
    #     if gaugeno in gaugenos_stationary:
    #         title('Surface elevation at stationary gauge %s' % gaugeno)
    #     else:
    #         title('Surface elevation at lagrangian gauge %s' % gaugeno)
    #
    # plotaxes.afteraxes = fix_gauge

    #-----------------------------------------

    # Parameters used only when creating html and/or latex hardcopy
    # e.g., via pyclaw.plotters.frametools.printframes:

    plotdata.printfigs = True                # print figures
    plotdata.print_format = 'png'            # file format
    plotdata.print_framenos = range(40)
    # plotdata.print_gaugenos = [15,25]        # list of gauges to print
    plotdata.print_fignos = 'all'            # list of figures to print
    plotdata.html = True                     # create html files of plots?
    plotdata.html_homelink = '../README.html'   # pointer for top of index
    plotdata.latex = True                    # create latex file of plots?
    plotdata.latex_figsperline = 2           # layout of plots
    plotdata.latex_framesperline = 1         # layout of plots
    plotdata.latex_makepdf = False           # also run pdflatex?
    plotdata.parallel = True                 # make multiple frame png's at once
    plotdata.html_movie_width = 700         # width used in JSAnimation

    return plotdata

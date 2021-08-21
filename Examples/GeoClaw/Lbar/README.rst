
.. _geoclaw_examples_Linear_barrier:

Flow past a linear single cell barrier
==========================================


To run the code:

   insert this example folder into your repo clone of GeoClaw, under `examples` directory

   and do : make all


In this code, :math:`x` and :math:`y` are in meters (coordinate_system=1
in `setrun.py`).

In `setrun.py` 2 gauges are set up.

In `setaux.f90` the bathymetry of the single cell barrier is set up.

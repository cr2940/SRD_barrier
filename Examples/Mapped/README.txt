To run the mapped examples, FIRST run:

f2py3 -c redist_module.f90 rpn2_shallow_Vmapped.f90 -m SWE_Vmap

in your example directory (where the example file is).

NEXT go to static.py in riemann/riemann/ and insert into list 

'SWE_Vmap' : 3

under both num_eqns and num_waves.

FINALLY run 

python swe_*.py

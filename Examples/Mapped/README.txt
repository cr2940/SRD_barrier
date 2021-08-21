To run the mapped examples, FIRST run:

f2py3 -c redist_module.f90 rpn2_shallow_Vmapped.f90 -m SWE_Vmap

in your example directory (where the example file is).

THEN run python swe_*.py

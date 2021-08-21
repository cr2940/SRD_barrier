To run these: 

  go to PyClaw_SRD_Codes
  put the superclass files into your PyClaw src directory
  put the classic files (Linear or V, corresponding to the example you want to run) into your PyClaw classic directory 
  run `python3 setup.py install` 
  update the .so compiled file with the newly compiled one (usually in the build/lib*/python..2d*.so)
  Also bring the appropriate cut-cell data file from `cut_cell_data` directory (organized by cell number) to the directory (where the examples are)
  Then run `python3 barrier_*.py`

# JWFL_lbm
LBM code for porous media applications

This implementation is fully parallel and use boost serialisation and CGNS export format.
The code cannot work in serial.
The code was tested with:
- gcc and intel compiler
- openmpi: 1.8.3 and 1.8.4
- mpich: 3.1.3 and 3.1.4
- cgns: 3.2.1 with hdf5 and parallel
- boost: 1.59.0

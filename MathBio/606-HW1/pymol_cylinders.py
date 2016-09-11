
# (Based on a program found at: http://www.rubor.de/bioinf/tips_python.html)
from pymol.cgo import *
from pymol import cmd
#
# CYLINDER, 0., 0., 0., 10., 0., 0., 0.2, 1.0, 1.0, 1.0, 1.0, 0.0, 0.,
#           --------    ----------  ----  -------------  -----------       
#            xyz 1        xyz 2      rad   rgb 1          rgb2
obj = [
  CYLINDER,   61.871,  -68.268,  -12.848,   62.489,  -68.061,  -10.957,  0.10, 1.000, 0.500, 0.250,  1.000, 0.500, 0.250 ,
  CYLINDER,   59.273,  -67.761,  -10.649,   58.310,  -69.508,  -10.786,  0.10, 1.000, 0.500, 0.250,  1.000, 0.500, 0.250 ,
  CYLINDER,   56.650,  -68.283,   -8.887,   56.421,  -67.495,   -7.064,  0.10, 1.000, 0.500, 0.250,  1.000, 0.500, 0.250 
]

# load it into PyMOL
cmd.load_cgo(obj,'axes')


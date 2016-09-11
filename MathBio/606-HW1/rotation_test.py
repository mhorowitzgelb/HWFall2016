#!/usr/bin/env python
#
howtousemessage ='''

To run this program you must provide one piece of information on the command line:

   python rotation_test.py  "degrees_to_rotate_by"

For example, to use a test angle of 45 degrees,

   python rotation_test.py  45

'''
import sys
import math
from numpy import *


def set_versor( angle, axis):
   '''
   NOTE: The axis vector must be normalized before passing to this function!
   Define a rotation matrix to rotate about vector "axis" by "angle" radians.
   Very cool, read all about "versors" in the American classic:
   "Vector Analysis", J. Willard Gibbs / E. B. Wilson
   also at http://www.archive.org/stream/117714283#page/338/mode/2up/search/versor

   NOTE: THIS IS EQUIVALENT TO THE USE OF Axis-Angle Representation in Section 2.3.2
   from class text. (page 36 or so)
   '''
   a1 = axis[0]
   a2 = axis[1]
   a3 = axis[2]
   cosq = math.cos(angle)
   cosqm1 = 1.0 - cosq
   sinq = math.sin(angle)
   versor=[]
   versor.append([])
   versor.append([])
   versor.append([])
   versor[0].append( a1*a1*cosqm1 + cosq )    
   versor[0].append( a1*a2*cosqm1 - a3*sinq ) 
   versor[0].append( a1*a3*cosqm1 + a2*sinq ) 
   versor[1].append( a1*a2*cosqm1 + a3*sinq ) 
   versor[1].append( a2*a2*cosqm1 + cosq )    
   versor[1].append( a2*a3*cosqm1 - a1*sinq ) 
   versor[2].append( a1*a3*cosqm1 - a2*sinq ) 
   versor[2].append( a2*a3*cosqm1 + a1*sinq ) 
   versor[2].append( a3*a3*cosqm1 + cosq )    
   return versor

def normalize_3vector( v ):
   vv = v[0]*v[0] + v[1]*v[1] + v[2]*v[2]
   if vv > 0. :
      vnorm = math.sqrt(vv)
      scale = 1. / vnorm
   else:
      return False

   return [ scale*v[0], scale*v[1], scale*v[2] ]


if len( sys.argv ) < 2:
   print "Insufficient number of parameters used on command line."
   print howtousemessage
   sys.exit(1)   
else:
   try:
      print sys.argv[1]
      Degrees = float(sys.argv[1])
      print "Will use:", "degrees for testing rotation matrices."
   except:
      print "Command line argument 1, number of degrees to rotate by, must be a number.  Please try again."
      sys.exit(1)

rad2deg = 180.0 / pi
deg2rad = pi / 180.0

############# HOEMWORK: TRY OUT SOME DIFFERENT VECTORS #####################################
testvecset = [[2.,2.,2.],[3.,4.,5.], [ 1., 0., 0.], [ -1., 0., 0.], [ 0., 0., 1.], [ 1., 0., 1.], [ 1., 1., 1.], [ 1., 2., 3.] ]

for rotaxis in testvecset:
  print
  print "Rotation axis (unnormalized):", rotaxis,  "Angle of rotation:", Degrees
  axis = normalize_3vector(rotaxis)
  rmat = set_versor( deg2rad*Degrees, axis)
  # cast to the Python list type object into a Numpy array type...
  r = array(rmat)
  print r
  print "Rotation matrix determinant:", linalg.det(r)


  # Check to see if matrix return the correct value for the angle:
  trace = r.trace()
  # check angle using Eq (2.27), page 30
  theta =  rad2deg * math.acos(max(-1,min(1.,(trace - 1.) / 2.)))
  print "... checking the rotation angle of matrix using Eq (2.27):", theta

  # Check to see if matrix returns the correct value for the axis:
  evals, evecs = linalg.eig(r)
  print "... eigenvalues and eigenvectors of the matrix (look for all real case):"
  print '  ', evals[0], evecs[:,0]
  print '  ', evals[1], evecs[:,1]
  print '  ', evals[2], evecs[:,2]


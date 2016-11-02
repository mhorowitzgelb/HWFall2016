import sys
import math
import numpy

#
# Overall size of figure, adjust for your screen.
# Uncomment your choice and comment out your non-choice:
#

bigscreen = True   # For Imacs and comparable sized displays
bigscreen = False  # For most Laptops

################################################
# Default Parameters for Van der Waals Potential
#        See Equation 3.11:
################################################
sigmaVdW = 3  # note: this uses the same parameters for all atom types
welldepthVdW = 0.1 # note: this uses the same parameters for all atom types
ConVdW = welldepthVdW
topval_VderWaals = 100000. # note: this sets max VdW so it's not huge
##########################################


def help_page():
   help = '''
   You must use 9 command line arguments:
       1   built-in model number for Fixed molecule
       2   electrostatic energy function number 1 for Coulomb, 2 for attenuated Coulomb 
       3   alpha value
       4   epsilon value (add to 1/r to give 1 / (r+epsilon) ) use small nonnegative value
       5   grid lowest value
       6   grid highest value
       7   grid step size
       8   plot option: 0 for No Plot, else Plot
       9   built-in model number for Moved molecule
       10  x position of minimum found by gradient descent
       11  y position of minimum found by gradient descent
     [12]  sigmaVdW, where 6-12 pot is 0.0 (default = 3.0)
     [13]  max well depth for VdW (default = 0.1)
     [14]  set all VdW greater than this value to this value (default = 100000.)

   The last three are optional Van der Waals potential pars, need to give all three or none.

   Examples:

   python energyAB2.py 12 1 .75 .1 -5 5 .25 1 11

   python energyAB2.py 12 1 .75 .1 -5 5 .25 1 11 2.5 .5 100.

   python energyAB2.py 19 1 .75 .1 -5 5 .5 1 18 1.2 .5 100.

   '''
   print help

if len(sys.argv) < 12:
   print
   print "Error: not enough command line arguments used."
   help_page()
   sys.exit(2)

modelA  = int(sys.argv[1])
efunc  = int(sys.argv[2])
alpha = float(sys.argv[3])
eps   = float(sys.argv[4])
lb    = float(sys.argv[5])
ub    = float(sys.argv[6])
step  = float(sys.argv[7])
plotoption = int(sys.argv[8])
modelB  = int(sys.argv[9])
grad_x = float(sys.argv[10])
grad_y = float(sys.argv[11])

if len(sys.argv) == 15:
   sigmaVdW = float(sys.argv[12])
   welldepthVdW = float(sys.argv[13])
   topval_VderWaals = float(sys.argv[14])
   ConVdW = welldepthVdW

print 'using fixed molecule model:', modelA
print 'using moved molecule model:', modelB
print 'using energy function:', efunc
print 'plotoption:', plotoption,
if plotoption == 0:
   print '(Skip all plots.)'
else:
   print '(Will draw plots.)'

print

model_descr = { 1: " 1CRN-THR-1/2 ",
                2: " two CA++ verti",
                3: " CA++/O-- verti",
                4: " two O-- verti",     
                5: " two CA++ horiz",
                6: " CA++/O-- horiz",
                7: " two O-- horiz",
                8: " 1CRN-GLU-23",
                18: " 1CRN-side-GLU-23",
                9: " 1CRN-ARG 10",
                19: " 1CRN-side-ARG-10",
                10: " one CA++",
                11: " one O--",
                12: " L-shape-3xCA++"}

maintitle = "A: " + sys.argv[1] + model_descr[modelA]+ "; "
maintitle += "B: " + sys.argv[9] + model_descr[modelB]+ "; "
maintitle += "lb: " + sys.argv[5] + "; "
maintitle += "ub: " + sys.argv[6] + "; "
maintitle += "step: " + sys.argv[7] +  "; "
maintitle += "Eqq type: " + sys.argv[2] + "; "
maintitle += "alpha: " + sys.argv[3] + "; "
maintitle += "eps: " + sys.argv[4] + "; "
maintitle += "sigVdW: " + str(sigmaVdW) + "; "
maintitle += "wellVdW: " + str(welldepthVdW) + "; "
maintitle += "topVdW: " + str(topval_VderWaals)

def get_molecule_model(model):
   if model == 1:
      # THREONINE 1 and 2 from PDB 1crn
      # ( x, y, z, charge, radius )
      atoms = [
      (17.047, 14.099, 3.625, 0.1812 , 1.824),
      (16.967, 12.784, 4.338, 0.0034 , 1.908),
      (15.685, 12.755, 5.133, 0.6163 , 1.908),
      (15.268, 13.825, 5.594, -0.5722, 1.6612),
      (18.17 , 12.703, 5.337, 0.4514 , 1.908),
      (19.334, 12.829, 4.463, -0.6764, 1.721),
      (18.15 , 11.546, 6.304, -0.2554, 1.908),
      (16.944, 12.070, 3.646, 0.1087 , 1.1),
      (18.142, 13.494, 5.933, -0.0323, 1.387),
      (19.089, 12.928, 3.520, 0.407  , 0),
      (17.308, 11.590, 6.859, 0.0627 , 1.487),
      (18.163, 10.677, 5.806, 0.0627 , 1.487),
      (18.957, 11.601, 6.911, 0.0627 , 1.487),
      (16.122, 14.421, 3.452, 0.1934, 0.6),
      (17.532, 13.966, 2.767, 0.1934, 0.6),
      (17.536, 14.744, 4.203, 0.1934, 0.6),
      (15.115, 11.555, 5.265, -0.4157, 1.824),
      (13.856, 11.469, 6.066, -0.0389, 1.908),
      (14.164, 10.785, 7.379, 0.5973, 1.908),
      (14.993,  9.862, 7.443, -0.5679, 1.6612),
      (12.732, 10.711, 5.261, 0.3654, 1.908),
      (13.308,  9.439, 4.926, -0.6761, 1.721),
      (12.484, 11.442, 3.895, -0.2438, 1.908),
      (15.511, 10.776, 4.852, 0.2719, 0.6),
      (13.548, 12.399, 6.243, 0.1007, 1.387),
      (11.957, 10.566, 5.874, 0.0043, 1.387),
      (13.421,  9.330, 3.945, 0.4102, 0),
      (12.166, 12.364, 4.050, 0.0642, 1.487),
      (13.345, 11.456, 3.383, 0.0642, 1.487),
      (11.803, 10.942, 3.372, 0.0642, 1.487),
      ]  
      imodx = 1
      imody = 2

   elif model == 2:
      # two CA++
      atoms = [
      (0., 2., 0., 2., 2.),
      (0., -2., 0., 2., 2.)
      ]  
      imodx = 1
      imody = 2

   elif model == 3:
      # CA++  and O--
      atoms = [
      (0., 2., 0., 2., 2.),
      (0., -2., 0., -2., 2.)
      ]  
      imodx = 1
      imody = 2

   elif model == 4:
      # two O--      
      atoms = [
      (0., 2., 0., -2., 2.),
      (0., -2., 0., -2., 2.)
      ]  
      imodx = 1
      imody = 2

   elif model == 5:
      # two CA++
      atoms = [
      (0., 2., 0., 2., 2.),
      (0., -2., 0., 2., 2.)
      ]  
      imodx = 2
      imody = 1

   elif model == 6:
      # CA++  and O--
      atoms = [
      (0., 2., 0., 2., 2.),
      (0., -2., 0., -2., 2.)
      ]  
      imodx = 2
      imody = 1

   elif model == 7:
      # two O--      
      atoms = [
      (0., 2., 0., -2., 2.),
      (0., -2., 0., -2., 2.)
      ]  
      imodx = 2
      imody = 1

   elif model == 8:
      # GLU 23 from PDB 1crn
      atoms = [
      (   -2.397,   -3.805,    0.719, -0.5163 , 1.8240 ),   # N   
      (   -1.905,   -3.393,   -0.621,  0.0397 , 1.9080 ),   # CA  
      (   -1.462,   -4.562,   -1.444,  0.5366 , 1.9080 ),   # C   
      (   -1.705,   -4.570,   -2.679, -0.5819 , 1.6612 ),   # O   
      (   -0.704,   -2.421,   -0.466,  0.0560 , 1.9080 ),   # CB  
      (   -1.149,   -1.022,    0.033,  0.0136 , 1.9080 ),   # CG  
      (    0.000,    0.000,    0.000,  0.8054 , 1.9080 ),   # CD  
      (    1.125,   -0.395,   -0.000, -0.8188 , 1.6612 ),   # OE1 
      (   -0.384,    1.196,    0.000, -0.8188 , 1.6612 ),   # OE2 
      (   -2.032,   -3.462,    1.539,  0.2936 , 0.6000 ),   # H   
      (   -2.646,   -2.889,   -1.092,  0.1105 , 1.3870 ),   # HA  
      (   -0.108,   -2.791,    0.204, -0.0173 , 1.4870 ),   # HB2 
      (   -0.304,   -2.303,   -1.342, -0.0173 , 1.4870 ),   # HB3 
      (   -1.848,   -0.699,   -0.565, -0.0425 , 1.4870 ),   # HG2 
      (   -1.438,   -1.110,    0.960, -0.0425 , 1.4870 )    # HG3 
      ]  
      imodx = 0
      imody = 1

   elif model == 18:
      # GLU 23 from PDB 1crn
      atoms = [
      (   -0.704,   -2.421,   -0.466,  0.0560 , 1.9080 ),   # CB  
      (   -1.149,   -1.022,    0.033,  0.0136 , 1.9080 ),   # CG  
      (    0.000,    0.000,    0.000,  0.8054 , 1.9080 ),   # CD  
      (    1.125,   -0.395,   -0.000, -0.8188 , 1.6612 ),   # OE1 
      (   -0.384,    1.196,    0.000, -0.8188 , 1.6612 ),   # OE2 
      (   -0.108,   -2.791,    0.204, -0.0173 , 1.4870 ),   # HB152 
      (   -0.304,   -2.303,   -1.342, -0.0173 , 1.4870 ),   # HB3 
      (   -1.848,   -0.699,   -0.565, -0.0425 , 1.4870 ),   # HG2 
      (   -1.438,   -1.110,    0.960, -0.0425 , 1.4870 )    # HG3 
      ]  
      imodx = 0
      imody = 1

   elif model == 9:
      # ARG 10 from PDB 1crn
      atoms = [
      (    1.638,    4.502,   -3.554, -0.3479 , 1.8240 ),   # N  
      (    0.753,    4.211,   -2.410, -0.2637 , 1.9080 ),   # CA 
      (   -0.150,    5.428,   -2.074,  0.7341 , 1.9080 ),   # C  
      (   -1.338,    5.246,   -1.773, -0.5894 , 1.6612 ),   # O  
      (    1.599,    3.836,   -1.171, -0.0007 , 1.9080 ),   # CB 
      (    0.758,    3.599,    0.078,  0.0390 , 1.9080 ),   # CG 
      (   -0.196,    2.448,    0.037,  0.0486 , 1.9080 ),   # CD 
      (    0.613,    1.196,   -0.010, -0.5295 , 1.8240 ),   # NE 
      (    0.000,   -0.000,    0.000,  0.8076 , 1.9080 ),   # CZ 
      (   -1.338,   -0.118,    0.000, -0.8627 , 1.8240 ),   # NH1
      (    0.790,   -1.082,   -0.000, -0.8627 , 1.8240 ),   # NH2
      (    2.621,    4.514,   -3.460,  0.2747 , 0.6000 ),   # H  
      (    0.160,    3.437,   -2.648,  0.1560 , 1.3870 ),   # HA 
      (    2.091,    3.002,   -1.378,  0.0327 , 1.4870 ),   # HB2
      (    2.222,    4.582,   -0.990,  0.0327 , 1.4870 ),   # HB3
      (    1.387,    3.467,    0.871,  0.0285 , 1.4870 ),   # HG2
      (    0.227,    4.449,    0.276,  0.0285 , 1.4870 ),   # HG3
      (   -0.733,    2.438,    0.854,  0.0687 , 1.3870 ),   # HD2
      (   -0.735,    2.499,   -0.777,  0.0687 , 1.3870 ),   # HD3
      (    1.580,    1.256,   -0.049,  0.3456 , 0.6000 ),   # HE 
      (   -1.894,    0.712,   -0.028,  0.4478 , 0.6000 ),   #HH11
      (   -1.758,   -1.015,    0.026,  0.4478 , 0.6000 ),   #HH12
      (    1.777,   -1.007,   -0.027,  0.4478 , 0.6000 ),   #HH21
      (    0.363,   -1.996,    0.027,  0.4478 , 0.6000 )    #HH22
      ]  
      imodx = 0
      imody = 1

   elif model == 19:
      # ARG 10 from PDB 1crn
      atoms = [
      (    1.599,    3.836,   -1.171, -0.0007 , 1.9080 ),   # CB 
      (    0.758,    3.599,    0.078,  0.0390 , 1.9080 ),   # CG 
      (   -0.196,    2.448,    0.037,  0.0486 , 1.9080 ),   # CD 
      (    0.613,    1.196,   -0.010, -0.5295 , 1.8240 ),   # NE 
      (    0.000,   -0.000,    0.000,  0.8076 , 1.9080 ),   # CZ 
      (   -1.338,   -0.118,    0.000, -0.8627 , 1.8240 ),   # NH1
      (    0.790,   -1.082,   -0.000, -0.8627 , 1.8240 ),   # NH2
      (    2.091,    3.002,   -1.378,  0.0327 , 1.4870 ),   # HB2
      (    2.222,    4.582,   -0.990,  0.0327 , 1.4870 ),   # HB3
      (    1.387,    3.467,    0.871,  0.0285 , 1.4870 ),   # HG2
      (    0.227,    4.449,    0.276,  0.0285 , 1.4870 ),   # HG3
      (   -0.733,    2.438,    0.854,  0.0687 , 1.3870 ),   # HD2
      (   -0.735,    2.499,   -0.777,  0.0687 , 1.3870 ),   # HD3
      (    1.580,    1.256,   -0.049,  0.3456 , 0.6000 ),   # HE 
      (   -1.894,    0.712,   -0.028,  0.4478 , 0.6000 ),   #HH11
      (   -1.758,   -1.015,    0.026,  0.4478 , 0.6000 ),   #HH12
      (    1.777,   -1.007,   -0.027,  0.4478 , 0.6000 ),   #HH21
      (    0.363,   -1.996,    0.027,  0.4478 , 0.6000 )    #HH22
      ]  
      imodx = 0
      imody = 1

   elif model == 10:
      # one CA++      
      atoms = [
      (0., 0., 0., 2., 2.),
      ]  
      imodx = 2
      imody = 1

   elif model == 11:
      # one O--      
      atoms = [
      (0., 0., 0., -2., 2.),
      ]  
      imodx = 2
      imody = 1

   elif model == 12:
      # L-shape, three CA++      
      atoms = [
      (3., 0., 0., 2., 2.),
      (0., 0., 0., 2., 2.),
      (0., 2., 0., 2., 2.),
      ]  
      imodx = 0
      imody = 1

   else:
      print "Error, no model defined for number:", model
      sys.exit(1)
   
   return (atoms, imodx, imody)


def gridMol( mol, lb, step ):
   '''
   Create new set of coordintes so that all atom positions are exactly on
   grid points.
   '''
   (atoms, imodx, imody) = mol

   gatoms = []
   for i in range(0,len(atoms)):
      # round to nearby grid point
      xf = int(round( ( atoms[i][imodx] - lb ) / step ))
      yf = int(round( ( atoms[i][imody] - lb ) / step ))
      xfg = lb + xf*step
      yfg = lb + yf*step
      (x,y,z,q,r) = atoms[i]
      gatoms.append( (xfg,yfg,z,q,r) )
      # make sure to use imodx = 0 and imody = 1
   return gatoms


def EnergyPair( alpha, eps, efunc,  q1, q2,  x1, x2, rAB ):
   '''
   Evaluate pair energy of two atoms located at x1 and x2.
   x1 and x2 must be Numpy arrays.
   '''
   if efunc < 1 or efunc > 5:
      print "Error in EnergyPair: no energy function for efunc =", efunc
      sys.exit(3)

   x12 = abs( x2 - x1 + rAB )
   d2 = x12**2
   r12 = math.sqrt( numpy.sum(d2) )
   E12 = 0.0

   if efunc == 1 or efunc == 4:
      E12 = 332.0636 * alpha * q1*q2 / ( r12 + eps )

   elif efunc == 2 or efunc == 5:
      E12 = 332.0636 * q1*q2 * numpy.exp( -alpha*r12 ) / ( r12 + eps )

   if efunc == 3 or efunc == 4 or efunc == 5:  # these functions include both electrostatic and vdW terms
      # Equation 3.11  (with eps added in)
      r6term = ( sigmaVdW / ( r12 + eps ) )**6
      Evdw = ConVdW * r6term * ( r6term - 1. )
      Evdw = min( Evdw, topval_VderWaals )
      E12 += Evdw

   return E12


def EnergyInterMol(molA, molB, efunc, alpha, eps,  rAB ):
   '''
   Total interaction energy of two molecules, A and B.
   rAB is a Numpy array giving displacement of B initial coordinates from A.
   '''
   (atomsA, imodxA, imodyA) = molA
   (atomsB, imodxB, imodyB) = molB

   E_AB = 0.0
   for i in range(0,len(atomsA)):
      xi = numpy.array( [ atomsA[i][imodxA], atomsA[i][imodyA]] )
      qi = atomsA[i][3]
      for j in range(0,len(atomsB)):
         xj = numpy.array( [ atomsB[j][imodxB], atomsB[j][imodyB]] )
         qj = atomsB[j][3]
         E_AB += EnergyPair( alpha, eps, efunc,  qi, qj,  xi, xj, rAB )
   return E_AB


def fourierTwoPots( mol, efunc, alpha, eps, ub, lb, step):
   '''
   Given molecule specs in mol, return 2D FFT for two potenial
   functions, 1: Coulombic-like and  2: Van der Waals - like
   Still use efunc to select type of Coulombic-like to use.
   Return the two FFTs needed to use these potentials for given molecule.
   '''
   #
   # Mess with these numbers to get nice region
   #
   print "fourierTwoPots options set to:"
   print '  efunc:', efunc
   print '  alpha:', alpha
   print '  eps:  ', eps
   print '  lb:   ', lb
   print '  ub:   ', ub
   print '  step: ', step
   print

   if abs( lb + ub ) > 0.00001:
       print "WARNING from fourierTwoPots: lb and ub should have same magnitude for this program."

   n = ( (ub-lb)/step ) / 2
   if n <= 0:
      print "Error in fourierTwoPots: ub and/or lb values yielded n<=0."
      sys.exit()

   #############################
   # Check bounds of Moleucle...
   #############################

   (atoms, imodx, imody) = mol

   #
   # Check all atoms, will they all fit within the specified grid?
   #
   maxx = -1000000000000.
   minx =  1000000000000.
   maxy = -1000000000000.
   miny =  1000000000000.
   for i in range(0,len(atoms)):
       # round to nearby grid point
       xat = atoms[i][imodx]
       yat = atoms[i][imody]
       maxx = max(maxx,xat)
       minx = min(minx,xat)
       maxy = max(maxy,yat)
       miny = min(miny,yat)

   print 'Atom locations in are in region ', minx,'<= x <=',maxx, ' : ',  miny,'<= y <=',maxy
   if maxx > ub or minx < lb or maxy > ub or miny < lb:
      print "Error in fourierTwoPots: specified grid is too small for molecule, please use larger grid."
      sys.exit(3)

   ############################################
   # SET-UP GRID...
   ############################################
   xgrid = numpy.arange( lb, ub, step )
   ygrid = numpy.arange( lb, ub, step )
   x,y = numpy.meshgrid( xgrid, ygrid )

   r = numpy.sqrt( x**2 + y**2 )

   #############################
   # Define Van der Waals Field
   #############################
   # Equation 3.11  (with eps added in)
   r6term = ( sigmaVdW / (r + eps) )**6
   g = ConVdW * r6term * ( r6term - 1. )

   igg = 0
   for gg in g:
      iv = 0
      for v in gg:
         #print igg,iv, g[igg,iv]
         if v > topval_VderWaals:
            g[igg,iv] = topval_VderWaals
         iv += 1
      igg += 1

   gmin = g.min()
   gmax = g.max()
   print "gmin, gmax:",gmin, gmax
   if math.isnan(gmin) or  math.isnan(gmax):
      print "ERROR in fourierTwoPots: g VderWaals function is flat, try using different parameters."
      print
      sys.exit(6)
   elif math.isinf(gmin) or  math.isinf(gmax):
      print "ERROR in fourierTwoPots: g VderWaals function is infinite, try using different parameters."
      print
      sys.exit(6)

   if gmin > 0.0:
      print "WARNING from fourierTwoPots: VderWaals potential is strictly positive, no obvious well."

   debug = False
   if debug:
     for gg in g:
         print "gg", gg

   ################################################
   # define Van der Waals source terms for Molecule
   ################################################

   f = numpy.zeros(numpy.shape(x))

   for i in range(0,len(atoms)):
       # round to nearby grid point
       xf = int(round( ( atoms[i][imodx] - lb ) / step ))
       yf = int(round( ( atoms[i][imody] - lb ) / step ))
       #f[xf,yf] = f[xf,yf] + 1
       f[yf,xf] = f[yf,xf] + 1
       if f[yf,xf] > 1:
          #print "Warning from fourierTwoPots: course grid, more than one atom at same grid point:",xf,yf,f[xf,yf] 
          print "Warning from fourierTwoPots: course grid, more than one atom at same grid point:",xf,yf,f[yf,xf] 

   # Calculate 2d FFT's for Van der Waals-like function...
   fft2g = numpy.fft.fft2(g)

   # Calculate 2d FFT's for atom locations in Molecule...
   fft2f = numpy.fft.fft2(f)
   #fft2fgVdW = numpy.conj(fft2f)*fft2g
   fft2fgVdW = fft2f*numpy.conj(fft2g)

   ############################################
   # Define Coulombic Field
   ############################################

   # (add "eps" to r in 1/r to avoid singularity and keep from getting too large for plot)
   # doing electrostatics?

   if efunc == 1:
       g = 332.0636 * alpha / ( r + eps )
       elec = True

   elif efunc == 2:
       g = 332.0636 * numpy.exp( -alpha*r ) / ( r + eps )
       elec = True
   else:
      print "Error fourierTwoPots: no electrostatic energy function for efunc =", efunc
      sys.exit(3)

   gmin = g.min()
   gmax = g.max()
   print "gmin, gmax:",gmin, gmax
   if math.isnan(gmin) or  math.isnan(gmax):
      print "ERROR fourierTwoPots: g Coul function is flat, try using different parameters."
      print
      sys.exit(6)
   elif math.isinf(gmin) or  math.isinf(gmax):
      print "ERROR fourierTwoPots: g Coul function is infinite, try using different parameters."
      print
      sys.exit(6)

   debug = False
   if debug:
     for gg in g:
         print "gg", gg

   ############################################
   # define Coulombic source terms for Molecule
   ############################################
   f = numpy.zeros(numpy.shape(x))

   # Note the difference in how individual elements of multi-dimensional
   #   Numpy arrays (f) are addressed compared to standard Python lists (atoms)...

   for i in range(0,len(atoms)):
      # round to nearby grid point
      xf = int(round( ( atoms[i][imodx] - lb ) / step ))
      yf = int(round( ( atoms[i][imody] - lb ) / step ))
      f[yf,xf] = f[yf,xf] + atoms[i][3] # for electrostatics only

   # Calculate 2d FFT's for Coulombic function...
   fft2g = numpy.fft.fft2(g)

   # Calculate 2d FFT's for Stationary Molecule Charges
   fft2f = numpy.fft.fft2(f)

   fft2fgCoul = fft2f*numpy.conj(fft2g)

   return x, y, f, fft2fgCoul, fft2fgVdW


def fourierSources( mol, efunc, alpha, eps, ub, lb, step):
   '''
   Given molecule specs in mol, return 2D FFT for two molecule
   source terms to use with the two potenial functions,
   1: Coulombic-like and  2: Van der Waals - like
   Return the two source FFTs needed to use these potentials
   in combination with the FFTs for a given fixed molecule. (defined elsewhere) 
   The FFT's from this function will be used to generate the correlation functions
   for the two energy terms using the corresponding FFTs for the other molecule.
   '''
   #
   # Mess with these numbers to get nice region
   #
   print "fourierSources options set to:"
   print '  efunc:', efunc
   print '  alpha:', alpha
   print '  eps:  ', eps
   print '  lb:   ', lb
   print '  ub:   ', ub
   print '  step: ', step
   print

   if abs( lb + ub ) > 0.00001:
       print "WARNING: fourierSources, lb and ub should have same magnitude for this program."

   n = ( (ub-lb)/step ) / 2
   if n <= 0:
      print "Error fourierSources, ub and/or lb values yielded n<=0."
      sys.exit()

   ############################################
   # Check bounds of Moleucle...
   ############################################

   (atoms, imodx, imody) = mol

   #
   # Check all atoms, will they all fit within the specified grid?
   #
   maxx = -1000000000000.
   minx =  1000000000000.
   maxy = -1000000000000.
   miny =  1000000000000.
   for i in range(0,len(atoms)):
       # round to nearby grid point
       xat = atoms[i][imodx]
       yat = atoms[i][imody]
       maxx = max(maxx,xat)
       minx = min(minx,xat)
       maxy = max(maxy,yat)
       miny = min(miny,yat)

   print 'fourierSources: Atom locations in are in region ', minx,'<= x <=',maxx, ' : ',  miny,'<= y <=',maxy
   if maxx > ub or minx < lb or maxy > ub or miny < lb:
      print "Error fourierSources: specified grid is too small for molecule, please use larger grid."
      sys.exit(3)

   ############################################
   # SET-UP GRID...
   ############################################
   xgrid = numpy.arange( lb, ub, step )
   ygrid = numpy.arange( lb, ub, step )
   x,y = numpy.meshgrid( xgrid, ygrid )

   ################################################
   # define Van der Waals source terms for Molecule
   ################################################

   f = numpy.zeros(numpy.shape(x))

   for i in range(0,len(atoms)):
       # round to nearby grid point
       xf = int(round( ( atoms[i][imodx] - lb ) / step ) )
       yf = int(round( ( atoms[i][imody] - lb ) / step ) )
       f[yf,xf] = f[yf,xf] + 1
       if f[yf,xf] > 1:
          print "Warning fourierSources: course grid, more than one atom at same grid point:",xf,yf,f[yf,xf] 

   # Calculate 2d FFT's for atom locations in Molecule...
   fft2fVdW = numpy.fft.fft2(f)
   ############################################
   # define Coulombic source terms for Molecule
   ############################################
   f = numpy.zeros(numpy.shape(x))

   # Note the difference in how individual elements of multi-dimensional
   #   Numpy arrays (f) are addressed compared to standard Python lists (atoms)...

   for i in range(0,len(atoms)):
      # round to nearby grid point
      xf = int(round( ( atoms[i][imodx] - lb ) / step ) )
      yf = int(round( ( atoms[i][imody] - lb ) / step ) )
      f[yf,xf] = f[yf,xf] + atoms[i][3] # for electrostatics only

   # Calculate 2d FFT's for Moving Molecule Charges
   fft2fCoul = numpy.fft.fft2(f)
   return f, fft2fCoul, fft2fVdW


def net_charge_mol( mol ):
   (atoms, imodx, imody) = mol
   sigmaQ = 0.
   for i in range(0,len(atoms)):
      sigmaQ += atoms[i][3]
   return sigmaQ

mA = get_molecule_model(modelA)
mB = get_molecule_model(modelB)
QtotA = net_charge_mol( mA )
QtotB = net_charge_mol( mB )
print "Total charge molecule A:", QtotA
print "Total charge molecule B:", QtotB

x, y, fA, fft2fgCoul_A, fft2fgVdW_A = fourierTwoPots( mA, efunc, alpha, eps, ub, lb, step)
fB, fft2fCoulSrc_B, fft2fVdWSrc_B = fourierSources( mB, efunc, alpha, eps, ub, lb, step)

fft2ECoul = numpy.conj(fft2fCoulSrc_B)*fft2fgCoul_A
ECoul = numpy.fft.ifft2(fft2ECoul)

fft2EVdW = numpy.conj(fft2fVdWSrc_B)*fft2fgVdW_A
EVdW = numpy.fft.ifft2(fft2EVdW)


Etot = ECoul + EVdW

CheckEnergy = False
CheckEnergy = True

if CheckEnergy:

   # To get exact agreement between discrete sum Energy and FFT Energy
   #  need to reposition the atoms to sit exactly at grid points.

   PrintAllPoints = True
   PrintAllPoints = False
   UseGridAtoms = False
   UseGridAtoms = True

   if UseGridAtoms:
      gatoms = gridMol( mA, lb, step )
      mA = (gatoms, 0, 1)
      gatoms = gridMol( mB, lb, step )
      mB = (gatoms, 0, 1)

   ixAB = 0
   minEAB = 1000000000000.
   maxEAB = -1000000000000.
   minEfft = 1000000000000.
   maxEfft = -1000000000000.
   for xAB in numpy.arange( lb, ub, step ):
      iyAB = 0
      for yAB in numpy.arange( lb, ub, step ):
         rAB = numpy.array( [ xAB, yAB ] )
         EsumAB = EnergyInterMol( mA, mB, efunc+3, alpha, eps,  rAB )
         if minEAB > EsumAB:
            minEAB = EsumAB
            min_xAB = xAB
            min_yAB = yAB
         if maxEAB < EsumAB:
            maxEAB = EsumAB
            max_xAB = xAB
            max_yAB = yAB

         Efft = numpy.real(Etot[iyAB][ixAB])
         if minEfft > Efft:
            minEfft = Efft
            min_xfft = xAB
            min_yfft = yAB
         if maxEfft < Efft:
            maxEfft = Efft
            max_xfft = xAB
            max_yfft = yAB

         if PrintAllPoints:
            print 'ix (x) | iy (y):', ixAB, '(', xAB, ') |', iyAB, '(', yAB, ') || Efft, Esum:', Efft, ',', EsumAB
         iyAB += 1
      ixAB += 1

   print "Min AB:", minEAB, 'at', min_xAB, min_yAB
   print "Max AB:", maxEAB, 'at', max_xAB, max_yAB
   print "Min fft:", minEfft, 'at', min_xfft, min_yfft
   print "Max fft:", maxEfft, 'at', max_xfft, max_yfft


print "Min and Max of FFT energies:"
lowCoul = numpy.real(ECoul).min()
hiCoul = numpy.real(ECoul).max()
print "  Electrostatic Energy: %12.3f %12.3f" %(lowCoul, hiCoul)
lowVdW = numpy.real(EVdW).min()
hiVdW = numpy.real(EVdW).max()
print "  Van der Waals Energy: %12.3f %12.3f" %(lowVdW, hiVdW)
lowE = numpy.real(Etot).min()
hiE = numpy.real(Etot).max()
print "   Total  FFT  Energy : %12.3f %12.3f" %(lowE, hiE)

if plotoption != 0:
   from matplotlib import cm
   import matplotlib.pyplot as plt
   import mpl_toolkits.mplot3d.axes3d as p3

   #
   # Overall size of figure, adjust for your screen:
   #
   if bigscreen:
      fig = plt.figure(figsize=(18.,6.6))
   else:
      fig = plt.figure(figsize=(15.,7.8))

   fig.suptitle(maintitle, fontsize=12)

   ax = fig.add_subplot(2, 3, 1)
   ax.set_title("Molecule A: Qtot = " + str(QtotA) )
   surf1 = ax.contourf(x,y,fA, cmap=cm.RdBu, linewidth=0, antialiased=True)
   fig.colorbar(surf1, shrink=0.5, aspect=10)

   ax = fig.add_subplot(2, 3, 2)
   ax.set_title("Molecule B: Qtot = " + str(QtotB) )
   surf2 = ax.contourf(x,y,fB, cmap=cm.RdBu, linewidth=0, antialiased=True)
   fig.colorbar(surf2, shrink=0.5, aspect=10)

   ax = fig.add_subplot(2, 3, 4)
   ax.set_title("Corr_E_Coul")
   surf4 = ax.contourf(x,y,numpy.real(ECoul), cmap=cm.RdBu, linewidth=0, antialiased=True)
   fig.colorbar(surf4, shrink=0.5, aspect=10)

   ax = fig.add_subplot(2, 3, 5)
   ax.set_title("Corr_E_VdW")
   surf5 = ax.contourf(x,y,numpy.real(EVdW), cmap=cm.RdBu, linewidth=0, antialiased=True)
   fig.colorbar(surf5, shrink=0.5, aspect=10)

   ax = fig.add_subplot(2, 3, 6)
   ax.set_title("Corr_E_tot")
   surf6 = ax.contourf(x,y,numpy.real(Etot), cmap=cm.RdBu, linewidth=0, antialiased=True)
   fig.colorbar(surf6, shrink=0.5, aspect=10)


   # Top threshhold value for plot...
   iee = 0

   top_Evalplot = 10.
   for ee in Etot:
      iv = 0
      for v in ee:
         #print iee,iv, Etot[iee,iv]
         if v > top_Evalplot:
            Etot[iee,iv] = top_Evalplot
         iv += 1
      iee += 1

   ax = fig.add_subplot(2, 3, 3)
   surf3 = ax.contourf(x,y,numpy.real(Etot), cmap=cm.RdBu, linewidth=0, antialiased=True)

   if CheckEnergy:
      ax.set_title("Direct_E_tot")
      # Plot FFT min point
      ax.plot(min_xfft,min_yfft,'go',markersize=20)
      print min_xfft, min_yfft, min_xAB, min_yAB

      # Plot discrete E min point
      ax.plot(min_xAB,min_yAB,'yo',markersize=10)

      ax.plot(grad_x, grad_y, 'bo', markersize=10)


   else:
      ax.set_title("real(E_tot)")

   fig.colorbar(surf3, shrink=0.5, aspect=10)

# print out a PDB file to view min E config if using two AA models...

pymolcase=(8,9,18,19)
if CheckEnergy and modelA in pymolcase and modelB in pymolcase:
   pdblineparts9 = '''
ATOM    121  N   ARG    10       1.638   4.502  -3.554  -0.3479 1.8240/
ATOM    122  CA  ARG    10       0.753   4.211  -2.410  -0.2637 1.9080/
ATOM    123  C   ARG    10      -0.150   5.428  -2.074   0.7341 1.9080/
ATOM    124  O   ARG    10      -1.338   5.246  -1.773  -0.5894 1.6612/
ATOM    125  CB  ARG    10       1.599   3.836  -1.171  -0.0007 1.9080/
ATOM    126  CG  ARG    10       0.758   3.599   0.078   0.0390 1.9080/
ATOM    127  CD  ARG    10      -0.196   2.448   0.037   0.0486 1.9080/
ATOM    128  NE  ARG    10       0.613   1.196  -0.010  -0.5295 1.8240/
ATOM    129  CZ  ARG    10       0.000  -0.000   0.000   0.8076 1.9080/
ATOM    130  NH1 ARG    10      -1.338  -0.118   0.000  -0.8627 1.8240/
ATOM    131  NH2 ARG    10       0.790  -1.082  -0.000  -0.8627 1.8240/
ATOM    132  H   ARG    10       2.621   4.514  -3.460   0.2747 0.6000/
ATOM    133  HA  ARG    10       0.160   3.437  -2.648   0.1560 1.3870/
ATOM    134  HB2 ARG    10       2.091   3.002  -1.378   0.0327 1.4870/
ATOM    135  HB3 ARG    10       2.222   4.582  -0.990   0.0327 1.4870/
ATOM    136  HG2 ARG    10       1.387   3.467   0.871   0.0285 1.4870/
ATOM    137  HG3 ARG    10       0.227   4.449   0.276   0.0285 1.4870/
ATOM    138  HD2 ARG    10      -0.733   2.438   0.854   0.0687 1.3870/
ATOM    139  HD3 ARG    10      -0.735   2.499  -0.777   0.0687 1.3870/
ATOM    140  HE  ARG    10       1.580   1.256  -0.049   0.3456 0.6000/
ATOM    141 HH11 ARG    10      -1.894   0.712  -0.028   0.4478 0.6000/
ATOM    142 HH12 ARG    10      -1.758  -1.015   0.026   0.4478 0.6000/
ATOM    143 HH21 ARG    10       1.777  -1.007  -0.027   0.4478 0.6000/
ATOM    144 HH22 ARG    10       0.363  -1.996   0.027   0.4478 0.6000/'''
   pdblineparts19 = '''
ATOM    125  CB  ARG    10       1.599   3.836  -1.171  -0.0007 1.9080/
ATOM    126  CG  ARG    10       0.758   3.599   0.078   0.0390 1.9080/
ATOM    127  CD  ARG    10      -0.196   2.448   0.037   0.0486 1.9080/
ATOM    128  NE  ARG    10       0.613   1.196  -0.010  -0.5295 1.8240/
ATOM    129  CZ  ARG    10       0.000  -0.000   0.000   0.8076 1.9080/
ATOM    130  NH1 ARG    10      -1.338  -0.118   0.000  -0.8627 1.8240/
ATOM    131  NH2 ARG    10       0.790  -1.082  -0.000  -0.8627 1.8240/
ATOM    134  HB2 ARG    10       2.091   3.002  -1.378   0.0327 1.4870/
ATOM    135  HB3 ARG    10       2.222   4.582  -0.990   0.0327 1.4870/
ATOM    136  HG2 ARG    10       1.387   3.467   0.871   0.0285 1.4870/
ATOM    137  HG3 ARG    10       0.227   4.449   0.276   0.0285 1.4870/
ATOM    138  HD2 ARG    10      -0.733   2.438   0.854   0.0687 1.3870/
ATOM    139  HD3 ARG    10      -0.735   2.499  -0.777   0.0687 1.3870/
ATOM    140  HE  ARG    10       1.580   1.256  -0.049   0.3456 0.6000/
ATOM    141 HH11 ARG    10      -1.894   0.712  -0.028   0.4478 0.6000/
ATOM    142 HH12 ARG    10      -1.758  -1.015   0.026   0.4478 0.6000/
ATOM    143 HH21 ARG    10       1.777  -1.007  -0.027   0.4478 0.6000/
ATOM    144 HH22 ARG    10       0.363  -1.996   0.027   0.4478 0.6000/'''
   pdblineparts8 = '''
ATOM    322  N   GLU    23      -2.397  -3.805   0.719  -0.5163 1.8240/
ATOM    323  CA  GLU    23      -1.905  -3.393  -0.621   0.0397 1.9080/
ATOM    324  C   GLU    23      -1.462  -4.562  -1.444   0.5366 1.9080/
ATOM    325  O   GLU    23      -1.705  -4.570  -2.679  -0.5819 1.6612/
ATOM    326  CB  GLU    23      -0.704  -2.421  -0.466   0.0560 1.9080/
ATOM    327  CG  GLU    23      -1.149  -1.022   0.033   0.0136 1.9080/
ATOM    328  CD  GLU    23       0.000   0.000   0.000   0.8054 1.9080/
ATOM    329  OE1 GLU    23       1.125  -0.395  -0.000  -0.8188 1.6612/
ATOM    330  OE2 GLU    23      -0.384   1.196   0.000  -0.8188 1.6612/
ATOM    331  H   GLU    23      -2.032  -3.462   1.539   0.2936 0.6000/
ATOM    332  HA  GLU    23      -2.646  -2.889  -1.092   0.1105 1.3870/
ATOM    333  HB2 GLU    23      -0.108  -2.791   0.204  -0.0173 1.4870/
ATOM    334  HB3 GLU    23      -0.304  -2.303  -1.342  -0.0173 1.4870/
ATOM    335  HG2 GLU    23      -1.848  -0.699  -0.565  -0.0425 1.4870/
ATOM    336  HG3 GLU    23      -1.438  -1.110   0.960  -0.0425 1.4870/'''
   pdblineparts18 = '''
ATOM    326  CB  GLU    23      -0.704  -2.421  -0.466   0.0560 1.9080/
ATOM    327  CG  GLU    23      -1.149  -1.022   0.033   0.0136 1.9080/
ATOM    328  CD  GLU    23       0.000   0.000   0.000   0.8054 1.9080/
ATOM    329  OE1 GLU    23       1.125  -0.395  -0.000  -0.8188 1.6612/
ATOM    330  OE2 GLU    23      -0.384   1.196   0.000  -0.8188 1.6612/
ATOM    333  HB2 GLU    23      -0.108  -2.791   0.204  -0.0173 1.4870/
ATOM    334  HB3 GLU    23      -0.304  -2.303  -1.342  -0.0173 1.4870/
ATOM    335  HG2 GLU    23      -1.848  -0.699  -0.565  -0.0425 1.4870/
ATOM    336  HG3 GLU    23      -1.438  -1.110   0.960  -0.0425 1.4870/'''

   pdblookup = { 8:pdblineparts8,  18:pdblineparts18,  9:pdblineparts9,  19:pdblineparts19 }

#ATOM    325  O   GLU    23      -1.705  -4.570  -2.679  -0.5819 1.6612
#0123456789012345678901234567890123456789012345678901234567890123456789

   fpdbh = open("AB.pdb","w")
   fpdbffth = open("ABfft.pdb","w")
   for l in pdblookup[modelA].split('/'):
      l = l.strip()
      if l == '': continue
      print >> fpdbh, l
      print >> fpdbffth, l

   (atomsB, imodxB, imodyB) = mB
   ithatm = 0
   for l in pdblookup[modelB].split('/'):
      l = l.strip()
      if l == '': continue
      pdbcoords = []
      pdbcoords.append(l[30:38])
      pdbcoords.append(l[38:46])
      pdbcoords.append(l[46:54])
      xold= atomsB[ithatm][imodxB]
      yold= atomsB[ithatm][imodyB]

      # Use discrete energy min
      xnew = "%8.3f" %(xold + min_xAB)
      ynew = "%8.3f" %(yold + min_yAB)
      pdbcoords[imodxB] = xnew
      pdbcoords[imodyB] = ynew
      newcoords = ''.join( pdbcoords )
      print >> fpdbh, l[0:30] + newcoords + l[53:]

      # Use FFT energy min
      xnew = "%8.3f" %(xold + min_xfft)
      ynew = "%8.3f" %(yold + min_yfft)
      pdbcoords[imodxB] = xnew
      pdbcoords[imodyB] = ynew
      newcoords = ''.join( pdbcoords )
      print >> fpdbffth, l[0:30] + newcoords + l[53:]

      ithatm += 1

   print "Wrote PDB file AB.pdb to view min discrete E config with Pymol."
   print "Wrote PDB file ABfft.pdb to view min FFT E config with Pymol."
   fpdbh.close()
   fpdbffth.close()

if plotoption != 0:
   plt.show()

sys.exit()


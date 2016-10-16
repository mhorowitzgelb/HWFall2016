import sys
import math
import numpy

#
# Overall size of figure, adjust for your screen.
# Uncomment your choice and comment out your non-choice:
#
#bigscreen = False  # For most Laptops
bigscreen = True   # For Imacs and comparable sized displays

def help_page():
   help = '''
   You must use 8 command line arguments:
       1   built-in molecular model number
       2   energy function number 1 for Coulomb, 2 for attenuated Coulomb, 3 for LJ6-12
       3   alpha value
       4   epsilon value (add to 1/r to give 1 / (r+epsilon) ) use small nonnegative value
       5   grid lowest value
       6   grid highest value
       7   grid step size
       8   plot option: 1 for full 3D, 2 for contour lines, 3 for filled contour

   For example:

   python fouriercorr.py 8 2 1. .01 -10 10 .5 3

   '''
   print help


if len(sys.argv) < 9:
   print
   print "Error: not enough command line arguments used."
   help_page()
   sys.exit(2)

model  = int(sys.argv[1])
efunc  = int(sys.argv[2])
alpha0 = float(sys.argv[3])
eps0   = float(sys.argv[4])
lb0    = int(sys.argv[5])
ub0    = int(sys.argv[6])
step0  = float(sys.argv[7])
plotoption = int(sys.argv[8])

maintitle =  "molecule model: " + sys.argv[1] + ", "
maintitle += "energy function: " + sys.argv[2] + ", "
maintitle += "alpha: " + sys.argv[3] + ", "
maintitle += "eps: " + sys.argv[4] + ", "
maintitle += "lb: " + sys.argv[5] + ", "
maintitle += "ub: " + sys.argv[6] + ", "
maintitle += "step: " + sys.argv[7]

print 'using molecule model:', model
print 'using energy function:', efunc
print 'plotoption: ', plotoption
print

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
      (   -2.397,   -3.805,    0.719, -0.5163 , 1.8240 ),
      (   -1.905,   -3.393,   -0.621,  0.0397 , 1.9080 ),
      (   -1.462,   -4.562,   -1.444,  0.5366 , 1.9080 ),
      (   -1.705,   -4.570,   -2.679, -0.5819 , 1.6612 ),
      (   -0.704,   -2.421,   -0.466,  0.0560 , 1.9080 ),
      (   -1.149,   -1.022,    0.033,  0.0136 , 1.9080 ),
      (    0.000,    0.000,    0.000,  0.8054 , 1.9080 ),
      (    1.125,   -0.395,   -0.000, -0.8188 , 1.6612 ),
      (   -0.384,    1.196,    0.000, -0.8188 , 1.6612 ),
      (   -2.032,   -3.462,    1.539,  0.2936 , 0.6000 ),
      (   -2.646,   -2.889,   -1.092,  0.1105 , 1.3870 ),
      (   -0.108,   -2.791,    0.204, -0.0173 , 1.4870 ),
      (   -0.304,   -2.303,   -1.342, -0.0173 , 1.4870 ),
      (   -1.848,   -0.699,   -0.565, -0.0425 , 1.4870 ),
      (   -1.438,   -1.110,    0.960, -0.0425 , 1.4870 )
      ]  
      imodx = 0
      imody = 1

   elif model == 9:
      # ARG 10 from PDB 1crn
      atoms = [
      (    1.638,    4.502,   -3.554, -0.3479 , 1.8240 ),
      (    0.753,    4.211,   -2.410, -0.2637 , 1.9080 ),
      (   -0.150,    5.428,   -2.074,  0.7341 , 1.9080 ),
      (   -1.338,    5.246,   -1.773, -0.5894 , 1.6612 ),
      (    1.599,    3.836,   -1.171, -0.0007 , 1.9080 ),
      (    0.758,    3.599,    0.078,  0.0390 , 1.9080 ),
      (   -0.196,    2.448,    0.037,  0.0486 , 1.9080 ),
      (    0.613,    1.196,   -0.010, -0.5295 , 1.8240 ),
      (    0.000,   -0.000,    0.000,  0.8076 , 1.9080 ),
      (   -1.338,   -0.118,    0.000, -0.8627 , 1.8240 ),
      (    0.790,   -1.082,   -0.000, -0.8627 , 1.8240 ),
      (    2.621,    4.514,   -3.460,  0.2747 , 0.6000 ),
      (    0.160,    3.437,   -2.648,  0.1560 , 1.3870 ),
      (    2.091,    3.002,   -1.378,  0.0327 , 1.4870 ),
      (    2.222,    4.582,   -0.990,  0.0327 , 1.4870 ),
      (    1.387,    3.467,    0.871,  0.0285 , 1.4870 ),
      (    0.227,    4.449,    0.276,  0.0285 , 1.4870 ),
      (   -0.733,    2.438,    0.854,  0.0687 , 1.3870 ),
      (   -0.735,    2.499,   -0.777,  0.0687 , 1.3870 ),
      (    1.580,    1.256,   -0.049,  0.3456 , 0.6000 ),
      (   -1.894,    0.712,   -0.028,  0.4478 , 0.6000 ),
      (   -1.758,   -1.015,    0.026,  0.4478 , 0.6000 ),
      (    1.777,   -1.007,   -0.027,  0.4478 , 0.6000 ),
      (    0.363,   -1.996,    0.027,  0.4478 , 0.6000 )
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
      # three CA++      
      atoms = [
      (4., 0., 0., 2., 2.),
      (0., 0., 0., 2., 2.),
      (0., 2., 0., 2., 2.),
      ]  
      imodx = 0
      imody = 1


   else:
      print "Error, no model defined for number:", model
      sys.exit(1)
   
   return atoms, imodx, imody



def fouriercorr(efunc, alpha, eps, ub, lb, step):

   #
   # Mess with these numbers to get nice region
   #
   print "fouriercorr options set to:"
   print '  efunc:', efunc
   print '  alpha:', alpha
   print '  eps:  ', eps
   print '  lb:   ', lb
   print '  ub:   ', ub
   print '  step: ', step
   print

   #
   # Check bounds...
   #
   n = ( (ub-lb)/step ) / 2
   if n <= 0:
      print "Error, ub and/or lb values yielded n<=0."
      sys.exit()

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

   print 'Atom locations in model are in region ', minx,'<= x <=',maxx, ' : ',  miny,'<= y <=',maxy
   if maxx > ub or minx < lb or maxy > ub or miny < lb:
      print "Error: specified grid is too small for molecule, please use larger grid."
      sys.exit(3)

   xgrid = numpy.arange( lb, ub, step )
   ygrid = numpy.arange( lb, ub, step )
   x,y = numpy.meshgrid( xgrid, ygrid )
   r = numpy.sqrt( x**2 + y**2 )

   # define fundamental function
   # (add "eps" to r in 1/r to avoid singularity and keep from getting too large for plot)
   # doing electrostatics?

   if efunc ==1:
       g = alpha / (r+eps)
       elec = True

   elif efunc ==2:
       g = numpy.exp(-alpha*r) / (r + eps)
       elec = True

   elif efunc ==3:
       r6term = (alpha/r+eps)**6
       g = 2.* r6term * ( r6term - 2)

       igg = 0
       for gg in g:
          iv = 0
          for v in gg:
             #print igg,iv, g[igg,iv]
             if v > 10.0:
                g[igg,iv] = 10.0
             iv += 1
          igg += 1
       elec = False

   else:
      print "Error: no energy function for efunc =", efunc
      sys.exit(3)

   gmin = g.min()
   gmax = g.max()
   print "gmin, gmax:",gmin, gmax
   if math.isnan(gmin) or  math.isnan(gmax):
      print "ERROR: g function is flat, try using different parameters."
      print
      sys.exit(6)
   elif math.isinf(gmin) or  math.isinf(gmax):
      print "ERROR: g function is infinite, try using different parameters."
      print
      sys.exit(6)

   if efunc ==3 and gmin > 0.0:
      print "WARNING: LJ 6-12 potential is strictly positive, no obvious well."

   debug = False
   if debug:
     for gg in g:
         print "gg", gg

   # define source terms
   # Note, next line: original Matlab code read: f = zeros(size(x));
   #  so Numpy requires "shape" in place of "size" ..

   f = numpy.zeros(numpy.shape(x))

   # Note the difference in how individual elements of multi-dimensional
   #   Numpy arrays (f) are addressed compared to standard Python lists (atoms)...

   for i in range(0,len(atoms)):
       # round to nearby grid point
       xf = int(round( ( atoms[i][imodx] - lb ) / step ))
       yf = int(round( ( atoms[i][imody] - lb ) / step ))
       if elec:
           f[xf,yf] = f[xf,yf] + atoms[i][3] # for electrostatics only
       else:
           f[xf,yf] = f[xf,yf] + 1

       print 'atomi/xf/yf/f/x/y',i,xf,yf,f[xf,yf],( atoms[i][imodx] - lb ) / step,( atoms[i][imody] - lb ) / step

   print "fmin,fmax",f.min(), f.max()

   # calculate correlation
   fft2f = numpy.fft.fft2(f)
   fft2g = numpy.fft.fft2(g)
   fft2fg = fft2f*numpy.conj(fft2g)
   ifftfg = numpy.fft.ifft2(fft2fg)
   ifftfg = numpy.fft.fftshift(ifftfg)

   # fix offset issue for discrete FFT
   # note that numpy has a built-in function for this (fftshift)

#   fft2f = numpy.fft.fftshift(fft2f)
#   fft2g = numpy.fft.fftshift(fft2g)
#   fft2fg = numpy.fft.fftshift(fft2fg)
#   ifftfg = numpy.fft.fftshift(ifftfg)

   return x, y, f, g, fft2f, fft2g, fft2fg, ifftfg

atoms, imodx, imody = get_molecule_model(model)

x,y,f,g,fft2f,fft2g,fft2fg,ifftfg = fouriercorr( efunc, alpha0, eps0, ub0, lb0, step0 )

from matplotlib import cm
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
#
# Overall size of figure, adjust for your screen:
#
if bigscreen:
   fig = plt.figure(figsize=(18.,9.6))
else:
   fig = plt.figure(figsize=(15.,7.8))

fig.suptitle(maintitle, fontsize=13)

if plotoption == 1:
   #---- First subplot
   ax = fig.add_subplot(2, 3, 1, projection='3d')
   ax.set_title("f")
   #ax.view_init(90,0)
   # rstride=1, cstride=1

   surf1 = ax.plot_surface(x,y,f, cmap=cm.RdBu, linewidth=0, antialiased=True)
   #maxf = f.max()
   #ax.set_zlim3d(-maxf, maxf )
   ## no! ax.set_xticks([-80,-40,0,40,80])
   fig.colorbar(surf1, shrink=0.5, aspect=10)

   #---- Second subplot
   ax = fig.add_subplot(2, 3, 2, projection='3d')
   #ax.view_init(90,0)
   ax.set_title("g")
   surf2 = ax.plot_surface(x,y,g, cmap=cm.RdBu, linewidth=0, antialiased=True)
   #maxg = g.max()
   #ax.set_zlim3d(-maxg, maxg)
   fig.colorbar(surf2, shrink=0.5, aspect=10)

   #---- Third subplot
   ax = fig.add_subplot(2, 3, 3, projection='3d')
   #ax.view_init(90,0)
   ax.set_title("real(inverse_fft2fg)")
   surf3 = ax.plot_surface(x,y,numpy.real(ifftfg), cmap=cm.RdBu, linewidth=0, antialiased=True)
   #maxifftfg = numpy.real(ifftfg).max()
   #ax.set_zlim3d(-maxifftfg, maxifftfg)
   fig.colorbar(surf3, shrink=0.5, aspect=10)

   # Row 2

   #---- Fourth subplot

   ax = fig.add_subplot(2, 3, 4, projection='3d')
   #ax.view_init(90,0)
   ax.set_title("abs(fft2f)")
   surf4 = ax.plot_surface(x,y,numpy.abs(fft2f), cmap=cm.spectral, linewidth=0, antialiased=True)
   #maxfft2f = numpy.abs(fft2f).max()
   #ax.set_zlim3d(-maxfft2f, maxfft2f)
   fig.colorbar(surf4, shrink=0.5, aspect=10)

   #---- Fifth subplot
   ax = fig.add_subplot(2, 3, 5, projection='3d')
   #ax.view_init(90,0)
   ax.set_title("abs(fft2g)")
   surf5 = ax.plot_surface(x,y,numpy.abs(fft2g), cmap=cm.spectral, linewidth=0, antialiased=True)
   #maxfft2g = numpy.abs(fft2g).max()
   #ax.set_zlim3d(-maxfft2g, maxfft2g)
   fig.colorbar(surf5, shrink=0.5, aspect=10)

   #---- Sixth subplot
   ax = fig.add_subplot(2, 3, 6, projection='3d')
   #ax.view_init(90,0)
   ax.set_title("abs(fft2f*conj(fft2g))" )
   surf6 = ax.plot_surface(x,y,numpy.abs(fft2fg), cmap=cm.spectral, linewidth=0, antialiased=True)
   #maxfft2fg = numpy.abs(fft2fg).max()
   #ax.set_zlim3d(-maxfft2fg, maxfft2fg)
   fig.colorbar(surf6, shrink=0.5, aspect=10)

elif plotoption == 2 or plotoption == 3:
   #---- First subplot
   ax = fig.add_subplot(2, 3, 1)
   ax.set_title("f")
   if plotoption == 2:
      surf1 = ax.contour(x,y,f, cmap=cm.RdBu, linewidth=0, antialiased=True)
   else:
      surf1 = ax.contourf(x,y,f, cmap=cm.RdBu, linewidth=0, antialiased=True)
   #maxf = f.max()
   #ax.set_zlim3d(-maxf, maxf )
   ## no! ax.set_xticks([-80,-40,0,40,80])
   fig.colorbar(surf1, shrink=0.5, aspect=10)

   #---- Second subplot
   ax = fig.add_subplot(2, 3, 2)
   ax.set_title("g")
   if plotoption == 2:
      surf2 = ax.contour(x,y,g, cmap=cm.RdBu, linewidth=0, antialiased=True)
   else:
      surf2 = ax.contourf(x,y,g, cmap=cm.RdBu, linewidth=0, antialiased=True)
   #maxg = g.max()
   #ax.set_zlim3d(-maxg, maxg)
   fig.colorbar(surf2, shrink=0.5, aspect=10)

   #---- Third subplot
   ax = fig.add_subplot(2, 3, 3)
   ax.set_title("real(inverse_fft2fg)")
   if plotoption == 2:
      surf3 = ax.contour(x,y,numpy.real(ifftfg), cmap=cm.RdBu, linewidth=0, antialiased=True)
   else:
      surf3 = ax.contourf(x,y,numpy.real(ifftfg), cmap=cm.RdBu, linewidth=0, antialiased=True)
   #maxifftfg = numpy.real(ifftfg).max()
   #ax.set_zlim3d(-maxifftfg, maxifftfg)
   fig.colorbar(surf3, shrink=0.5, aspect=10)

   # Row 2

   #---- Fourth subplot
   ax = fig.add_subplot(2, 3, 4)
   ax.set_title("abs(fft2f)")
   if plotoption == 2:
      surf4 = ax.contour(x,y,numpy.abs(fft2f), cmap=cm.spectral, linewidth=0, antialiased=True)
   else:
      surf4 = ax.contourf(x,y,numpy.abs(fft2f), cmap=cm.spectral, linewidth=0, antialiased=True)
   #maxfft2f = numpy.abs(fft2f).max()
   #ax.set_zlim3d(-maxfft2f, maxfft2f)
   fig.colorbar(surf4, shrink=0.5, aspect=10)

   #---- Fifth subplot
   ax = fig.add_subplot(2, 3, 5)
   ax.set_title("abs(fft2g)")
   if plotoption == 2:
      surf5 = ax.contour(x,y,numpy.abs(fft2g), cmap=cm.spectral, linewidth=0, antialiased=True)
   else:
      surf5 = ax.contourf(x,y,numpy.abs(fft2g), cmap=cm.spectral, linewidth=0, antialiased=True)
   #maxfft2g = numpy.abs(fft2g).max()
   #ax.set_zlim3d(-maxfft2g, maxfft2g)
   fig.colorbar(surf5, shrink=0.5, aspect=10)

   #---- Sixth subplot
   ax = fig.add_subplot(2, 3, 6)
   ax.set_title("abs(fft2f*conj(fft2g))" )
   if plotoption == 2:
      surf6 = ax.contour(x,y,numpy.abs(fft2fg), cmap=cm.spectral, linewidth=0, antialiased=True)
   else:
      surf6 = ax.contourf(x,y,numpy.abs(fft2fg), cmap=cm.spectral, linewidth=0, antialiased=True)
   #maxfft2fg = numpy.abs(fft2fg).max()
   #ax.set_zlim3d(-maxfft2fg, maxfft2fg)
   fig.colorbar(surf6, shrink=0.5, aspect=10)

'''
# Below code attempts to usr slider widget like the Matlab program did to vary alpha value..,
# PARK THIS CODE FOR FUTURE RELEASE OF matplotlib WHICH SUPPORTS set_zdata FUNCTIONALITY IN 3D PLOTS.
# ( will need to add comma after variable name in all above 'surfN =' to read  'surfN, =' to use this)

from matplotlib.widgets import Slider, Button, RadioButtons
from pylab import *

axcolor = 'lightgoldenrodyellow'
axfreq = axes([0.25, 0.1, 0.65, 0.03], axisbg=axcolor)
axamp  = axes([0.25, 0.15, 0.65, 0.03], axisbg=axcolor)

newalp = Slider(axfreq, 'Alpha', 0.01, 200.0, valinit=alpha0)
neweps = Slider(axamp, 'Eps', 0.0, 1.0, valinit=eps0)

def update(val):
   alpha = newalp.val
   eps = neweps.val
   x, y, f, g, fft2f, fft2g, fft2fg, ifftfg = fouriercorr(alpha, eps)

   surf1.set_zdata(f)
   surf2.set_zdata(g)
   surf3.set_zdata(numpy.real(ifftfg))
   surf4.set_zdata(numpy.abs(fft2f))
   surf5.set_zdata(numpy.abs(fft2g))
   surf6.set_zdata(numpy.abs(fft2fg))
   draw()

newalp.on_changed(update)
neweps.on_changed(update)
'''

plt.show()

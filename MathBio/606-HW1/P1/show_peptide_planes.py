
def cross3vecs( v1, v2 ):
    '''
    Given two vectors, v1, v2, stored in lists,
    return the normalized cross-product or die if 0.
    '''
    from math import sqrt
###############################################
# HOMEWORK: Finish the code for C1 and C2
###############################################
    c0 =  v1[1]*v2[2] - v2[1]*v1[2]
    c1 =  v1[2]*v2[0] - v1[0]*v2[2]
    c2 =  v1[0]*v2[1] - v1[1]*v2[0]

    csq = c0**2 + c1**2 + c2**2
    if not csq > 0:
       print "ERROR: cross product is Null vector."
       exit(1)

    zz = 1.0 / sqrt( csq )
    c = []
    c.append( zz*c0 )
    c.append( zz*c1 )
    c.append( zz*c2 )
    return c

def bondvec( a1, a2 ):
    b = []
    b.append( a2[0] - a1[0] )
    b.append( a2[1] - a1[1] )
    b.append( a2[2] - a1[2] )
    return b

# not used in this program...
def distance3( x, y ):
    from math import sqrt
    d2 = ( y[0] - x[0] )**2 + ( y[1] - x[1] )**2 + ( y[2] - x[2] )**2
    if d2 > 0:
       d = sqrt(d2)
    else:
       d = False
    return d

#ATOM     32  CB  VAL A   6      54.977 -70.796  -8.926  1.00 17.11           C  
#ATOM     33  CG1 VAL A   6      55.702 -71.127 -10.195  1.00 16.13           C  
#01234567890123456789012345678901234567890123456789012345678901234567890123456789
#         1         2         3         4         5         6         7
#          atnam13:16           x30:38  y38:46  z46:54

def get_coords_from_PDB_line( line ):
   x = float( line[30:38] )
   y = float( line[38:46] )
   z = float( line[46:54] )
   return (x, y, z)

#
# read in file and determine normal vector to peptide plane
#
fh = open("tetrapeptide.pdb","r")
lines = fh.readlines()
fh.close()

N = []; CA = []; C = []; O = []

for l in lines:
   if l[0:6] == 'ATOM  ':
      print l[13:16], "**", l
      if l[13:16].strip() == "N" :
         N.append( get_coords_from_PDB_line( l ) )
      elif l[13:16].strip() == "CA" :
         CA.append( get_coords_from_PDB_line( l ) )
      elif l[13:16].strip() == "C" :
         C.append( get_coords_from_PDB_line( l ) )
      elif l[13:16].strip() == "O" :
         O.append( get_coords_from_PDB_line( l ) )
      atmnolast = int( l[6:11] )

debug = True
if debug:
   print "N ", N
   print "CA", CA
   print "C ", C
   print "O ", O
   print "atmnolast:", atmnolast

# We will now calculate the normals and write out a new file
# including all of the original file with the extra fake atoms defining
# the peptide normal vectors appended...

debug = False
if debug:
   for i in range(0,3):
      print i,  C[i], O[i], N[i+1], CA[i+1]

atmno = atmnolast + 1
resnum = 9000
newpdb = open("tetrapepnew.pdb","w")

# re-write all of the previous lines into a new PDB file...
for l in lines:
   print >> newpdb, l.strip()
#

# Add fake HETATMS to new pdb representing midpoint of CN bond and
#   endpoint of normal vector filere-write all of the previous lines
#    into a new PDB file...
#
# Also, save the atom pairs for drawing normal vectors for second Pymol
# rendering demo...
#
normal_pts = []

for i in range(0,3):
   # Calculate two approximate normals for each peptide bond...
   #  n1_i = C_i--O_i cross C_i--N_i+i  (will be up )
   #  n2_i = C_i--N_i+i cross N_i+i--CA_i+1  (will be down)
   n1 = cross3vecs( bondvec( C[i], O[i] ), bondvec( C[i], N[i+1] ) )
   n2 = cross3vecs( bondvec( N[i+1], CA[i+1]),  bondvec( C[i], N[i+1] ) )
   #print "n1:", n1
   #print "n2:", n2
   # for drawing normal, chose midpoint of peptide bond as basepoint...
   origx = ( C[i][0] +  N[i+1][0] ) / 2.
   origy = ( C[i][1] +  N[i+1][1] ) / 2.
   origz = ( C[i][2] +  N[i+1][2] ) / 2.
   s = "HETATM%5d   X1  UNK %4d    %8.3f%8.3f%8.3f  1.00 50.00           X" %(atmno,resnum,origx,origy,origz)
   print >> newpdb, s
   atmno += 1
   scale = 2.0
   x = scale*(n1[0] + n2[0]) / 2. +  origx
   y = scale*(n1[1] + n2[1]) / 2. +  origy
   z = scale*(n1[2] + n2[2]) / 2. +  origz
   s = "HETATM%5d   X2  UNK %4d    %8.3f%8.3f%8.3f  1.00 50.00           X" %(atmno,resnum,x,y,z)
   print >> newpdb, s
   atmno += 1
   resnum += 1
   normal_pts.append( ( origx, origy, origz, x,y,z ) )
newpdb.close()


# 
# Write out a small Python program to be run within a PyMol .pml file to
#  render the normal vectors as cylinders.
#

fhpy = open("pymol_cylinders.py", "w" )

part1 = '''
# (Based on a program found at: http://www.rubor.de/bioinf/tips_python.html)
from pymol.cgo import *
from pymol import cmd
#
# CYLINDER, 0., 0., 0., 10., 0., 0., 0.2, 1.0, 1.0, 1.0, 1.0, 0.0, 0.,
#           --------    ----------  ----  -------------  -----------       
#            xyz 1        xyz 2      rad   rgb 1          rgb2
obj = [
'''
print >>fhpy, part1,

color1 = ( 1.0, 0.5, 0.25 )
color2 = ( 1.0, 0.5, 0.25 )

cylrad = 0.1 # Angstroms
first = True
for pairs in normal_pts:

   if not first:
      print >>fhpy, ','
   else:
      first = False

   coords_string = "  CYLINDER, %8.3f, %8.3f, %8.3f, %8.3f, %8.3f, %8.3f," %pairs
   print coords_string
   cylrad_string = " %5.2f," %cylrad
   colors_string = " %5.3f, %5.3f, %5.3f, " %color1 + " %5.3f, %5.3f, %5.3f" %color2
   print >>fhpy,  coords_string + cylrad_string + colors_string,

part2 = '''
]

# load it into PyMOL
cmd.load_cgo(obj,'axes')
'''

print >>fhpy, part2
fhpy.close()


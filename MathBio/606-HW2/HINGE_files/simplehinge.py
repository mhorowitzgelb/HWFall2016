#!/usr/bin/env python
#
howtousemessage ='''

To run this program you must provide four pieces of information on the command line:

   python simplehinge.py  "pdb_file_name"  "residue_number_to_rotate_phi"  "AngleDegreesStepSize" "IncludHetAtms"

For example, using molecule 4cln.pdb, rotate around residue 67 Phi dihedral angle by increments of
 20 degrees and include heteroatoms in output.

   python simplehinge.py 4cln.pdb 67 20 T

All of the residues after specified residue number will be rotated and all the resdiues before
that residue will be left unchanged.   You can also select which HETATM', to leave fixed or to rotate
but this feature is hardwired into the code.

'''
import sys
from PDButils import *
from numpy import *
from math import sin, cos, acos, pi, sqrt

# set shorthand ptrs

_rnums = Res_["resnums"]
_ratms = Res_["atoms"]
_rnams = Res_["names"]
_rptrs = Res_["resptrs"]

_aname = Atm_["name"]
_arnum = Atm_["resnum"]
_xyz   = Atm_["coords"]
_B     = Atm_["biso"]

_N  = 0
_CA = 1
_C  = 2
_O  = 3
_CB = 4

# Energy constants...

qq_epscon = 331.878
vdw_epsw = 20.
vdw_rijmin_all = 2.0

def Eqq_pair( rij_inv, Q1, Q2):
   '''
   Pair term in Eq. (3.8) page 43.
   compute Coulombic energy between two charges, given positive value for distance
   '''
   return qq_epscon * Q1 *Q2 * rij_inv


def Vdw_pair( rij_inv, vdwpar ):
   '''
   Pair term in Eq. (3.13) page 45.
   compute approximate Van der Waals energy for atom pair
   '''
   rij_inv6 =  (vdwpar * rij_inv)**6
   return vdw_epsw* (  rij_inv6 * ( rij_inv6 - 2.0 ) )


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
   cosq = cos(angle)
   cosqm1 = 1.0 - cosq
   sinq = sin(angle)
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

if len( sys.argv ) < 5:
   print "Insufficient number of parameters used on command line."
   print howtousemessage
   sys.exit(1)   
else:
   pdbfile = sys.argv[1]

# next line should read all the data from the PDB and store it in memory...
prot1 = PDBin( pdbfile )

if not prot1.stat:
   print "Failed to open file:", pdbfile
   sys.exit(1)
else:
   print "Opened and read file:", pdbfile

(model_id, chain_id) = prot1.chainlist[0]
if model_id == '':
   print "The name of the first MODEL is: (MODEL records not used in this PDB file.)"
else:
   print "The name of the first MODEL is:", model_id

print "The name of the first CHAIN in MODEL is:", chain_id

# check/set other command line arguments...

Good = True
try:
   RotatePhiResNum = int(sys.argv[2])
   HingeStart = RotatePhiResNum - 2
   HingeEnd = RotatePhiResNum + 2
   print "Will rotate Phi for residue number:", RotatePhiResNum
except:
   print "Command line argument 2, residue number, must br an integer.  Please try again."
   Good = False

try:
   StepDegrees = float(sys.argv[3])
   print "Will rotate Phi by in steps of:", StepDegrees, "degrees"
except:
   print "Command line argument 3, number of degrees step by, must be a number.  Please try again."
   Good = False

if sys.argv[4].upper() == 'T':
   do_hetatms = True
elif sys.argv[4].upper() == 'F':
   do_hetatms = False
else:
   print 'For HETATM inclusion you must specify either "T" or "F" as commnad line argument 3.'
   Good = False

if not Good:
   print ''
   sys.exit(1)

########################################################################################
# NOTE, CODE IS NOT  GENERAL HERE::::

# Also select the which chain to use... (this is not general but specific to
#  this particular project - you really need to know what's in the PDB file beforehand
#  before you can make all these selections.  This could be specified as an additional
#  command line argument but has not been implemented for this demo.  If you run this program
#  for a different PDB file you will have to edit this section to define which chains
#  to use.

# Select the polymer (protein or nucleic acid) chain to use...
chain1 = prot1.chains[ (model_id, chain_id) ]
chainlab = prot1.chainlist[0][1]

# Also select the HETATMs chain to use... (this also is not general but somewhat custom to
#  this particular project - you really need to know what's in the PDB file beforehand
#  before you can make all these selections.

hetchain1 = prot1.hetchains[ (model_id, chain_id) ][0]
hetchainlab = prot1.hetchainlist[0][1]

# END NOTE. (there are more cases of application specific stuff below...)
########################################################################################


rad2deg = 180.0 / pi
deg2rad = pi / 180.0
iresnumprevious = -100000
resnumsinchain = chain1[_rptrs].keys()

if RotatePhiResNum in resnumsinchain:
   iresnumprev = RotatePhiResNum - 1
   if iresnumprev not in resnumsinchain:
      print "Unable to rotate Phi for residue number", RotatePhiResNum,". Previous residue is not in PDB file."
      sys.exit(2)
else:           
      print "Unable to rotate Phi for residue number", RotatePhiResNum,". The residue is not in PDB file."
      sys.exit(3)

j1 = chain1[_rptrs][iresnumprev]
j2 = chain1[_rptrs][RotatePhiResNum]

C1  = array( chain1[_ratms][j1][_C][_xyz] )
N2  = array( chain1[_ratms][j2][_N][_xyz] )
CA2 = array( chain1[_ratms][j2][_CA][_xyz] )
C2  = array( chain1[_ratms][j2][_C][_xyz] )

bnd_C1_N2  =  N2 - C1
bnd_N2_CA2 =  CA2 - N2
bnd_CA2_C2 =  C2 - CA2

#-----------------------------------------------------------------------
# Phi: dihedral angle of C(i-1)-N(i)---Ca(i)-C(i)
#-----------------------------------------------------------------------

n_C1_N2_CA2 = cross( bnd_C1_N2, bnd_N2_CA2 )
n_C1_N2_CA2 = n_C1_N2_CA2 / linalg.norm( n_C1_N2_CA2 )      

n_N2_CA2_C2 = cross( bnd_N2_CA2, bnd_CA2_C2 )
n_N2_CA2_C2 = n_N2_CA2_C2 / linalg.norm( n_N2_CA2_C2 )      

Phi = rad2deg * acos( dot( n_C1_N2_CA2, n_N2_CA2_C2 ) )

n1Xn2 = cross( n_C1_N2_CA2, n_N2_CA2_C2 )
if dot(bnd_N2_CA2,n1Xn2) < 0: Phi = -Phi

print "Current Phi:", Phi

'''

COMPARISON OF TRANSFORMATION CODE USING Numpy TO CODE NOT USING Numpy

# Matrix multiply with numpy

def atrans( armat, avec, atrans ):
   # using numpy arrays...
   print "armat:", armat
   print "avec:", avec
   print "atrans:", atrans
   # all in one line!!!, function not needed.
   newvec = dot( armat, avec-atrans ) + atrans
   return newvec


# Matrix multiply NOT using numpy, using standard Python lists...

rmat = set_versor( deg2rad*Degrees, rotaxis)
print rmat

def trans( rmat, vec, trans ):
   # using standard python lists...
   print "rmat:", rmat
   print "vec:", vec
   print "trans:", trans
   newvec = []
   newvec.append(rmat[0][0]*(vec[0]-trans[0])+ rmat[0][1]*(vec[1]-trans[1])+rmat[0][2]*(vec[2]-trans[2]))
   newvec.append(rmat[1][0]*(vec[0]-trans[0])+ rmat[1][1]*(vec[1]-trans[1])+rmat[1][2]*(vec[2]-trans[2]))
   newvec.append(rmat[2][0]*(vec[0]-trans[0])+ rmat[2][1]*(vec[1]-trans[1])+rmat[2][2]*(vec[2]-trans[2]))
   newvec[0] += trans[0]
   newvec[1] += trans[1] 
   newvec[2] += trans[2]
   return newvec
'''


# Rotate about: bnd_N2_CA2

#create a clone of PDB: (to leave original unchanged for possible future use )
hinge = PDBin( pdbfile )
chain1rot = hinge.chains[ (model_id, chain_id) ]
hetchain1rot = hinge.hetchains[ (model_id, chain_id) ][0]

# First rotate all atoms in the residue with the Phi angle dihedral...
# (Ca will stay fixed but should be the origin of the rotation...)

origin = CA2
rotaxis = bnd_N2_CA2 / linalg.norm( bnd_N2_CA2 )      

nsteps = int( 360. / StepDegrees + 0.5 )

import time
now = time.asctime(time.localtime())
key = (model_id, chain_id)

hinge.add_to_TITLE( "homework 2" )
hinge.add_to_TITLE( "rotate about Phi for residue number " + str(RotatePhiResNum) + "." )
hinge.add_to_TITLE( "use angle step size of " + str(StepDegrees) + "." )
hinge.add_to_TITLE( "this file created at: " + now )
hinge.rewrite_pdb_header_records( [ "REMARK" ] )
hinge.open_pdboutfile( "hinge_models_res" + str(RotatePhiResNum) + ".pdb" )
hinge.write_PDB_header()

# will use this to save inforamtion about the intger-domain interaction energy for plotting later...
savedata = []

for iang in range( 0, nsteps + 1):

   # Create Rotation Matrix to Rotate about bnd_N2_CA2 

   Degrees = iang*StepDegrees
   rotmat = array( set_versor( deg2rad*Degrees, rotaxis) )

   jatm = 0
   residue = chain1[_ratms][j2]
   for atom in residue:
      # print "atom id:", jatm
      if jatm not in ( _N, _CA ):
         oldcoords = array( atom[_xyz] )
         newcoords = dot( rotmat, oldcoords - origin ) + origin
         chain1rot[_ratms][j2][jatm][_xyz][0] = newcoords[0]
         chain1rot[_ratms][j2][jatm][_xyz][1] = newcoords[1]
         chain1rot[_ratms][j2][jatm][_xyz][2] = newcoords[2]
      jatm += 1

   transpart = origin - dot(rotmat, origin) 

   jres = j2 + 1
   for residue in chain1[_ratms][j2+1:]:
      jatm = 0
      for atom in residue:
         oldcoords = array( atom[_xyz] )
         newcoords = dot( rotmat, oldcoords) + transpart
         #xyz = newcoords.tolist()
         chain1rot[_ratms][jres][jatm][_xyz][0] = newcoords[0]
         chain1rot[_ratms][jres][jatm][_xyz][1] = newcoords[1]
         chain1rot[_ratms][jres][jatm][_xyz][2] = newcoords[2]
         jatm += 1
      jres += 1


   if do_hetatms:

      # for better runtime efficiency, some of these variables should really be defined above,
      # outside of rotation loop, but we leave them here to improve readability of code...

      # HERE AGAIN WE HARDWRIRE SOME PROBLEM SPECIFIC INFORMATION ABOUT PARTICULAR PDB FILE
      # WE ARE USING, NAMELY, WE WANT TO EXAMINE CALCIUM IONS...

      het_names_to_skip_in_print = ["HOH"]
      het_resnums_to_skip_rotation = [149,150]  # fixed Calciums
      het_resnums_moving_Calciums = [151, 152]  # moving Calciums
      resnumber = hetchain1[_rnums]
      resname = hetchain1[_rnams]
      hetjres = 0

      for residue in hetchain1[_ratms]:
         if resname[hetjres] not in het_names_to_skip_in_print:
            if resnumber[hetjres] not in het_resnums_to_skip_rotation:
               jatm = 0
               for atom in residue:
                  oldcoords = array( atom[_xyz] )
                  newcoords = dot( rotmat, oldcoords) + transpart
                  #xyz = newcoords.tolist()
                  hetchain1rot[_ratms][hetjres][jatm][_xyz][0] = newcoords[0]
                  hetchain1rot[_ratms][hetjres][jatm][_xyz][1] = newcoords[1]
                  hetchain1rot[_ratms][hetjres][jatm][_xyz][2] = newcoords[2]
                  jatm += 1

         hetjres += 1

   # COMPUTE INTER-DOMAIN ENERGY... (very approximately)

   nltcut = 0
   ngtcut = 0
   jres = j2 + 1
   # ------------------------------------
   # Loop over atoms of moving domain...
   # ------------------------------------
   for residue in chain1[_ratms][j2+1:]:
      if chain1[_rnums][jres] > HingeEnd:
         jatm = 0
         for atom in residue:
            newcoords_domain2 = array( chain1rot[_ratms][jres][jatm][_xyz] )
            jres_stat = 0
            # ------------------------------------------------
            # Loop over residues/atoms in stationary domain...
            # ------------------------------------------------
            for residue_stat in chain1[_ratms][0:j2]:
               if chain1[_rnums][jres_stat] < HingeStart:
                  for atomstat in residue_stat:
                     oldcoords_domain1 = array( atomstat[_xyz] )
                     d2 = ( ( newcoords_domain2 - oldcoords_domain1 )**2 ).sum()
                     if d2 > 0.0:
                        d = sqrt(d2)
                        if d > 3.:
                           ngtcut += 1
                        else:
                           nltcut += 1
               jres_stat += 1
         jatm += 1
      jres += 1

   # Get Calcium ion inter-domain energy...

   index_CA149 = hetchain1[_rptrs][149]
   index_CA150 = hetchain1[_rptrs][150]
   index_CA151 = hetchain1[_rptrs][151]
   index_CA152 = hetchain1[_rptrs][152]
   r149 = array( hetchain1[_ratms][index_CA149][0][_xyz] )
   r150 = array( hetchain1[_ratms][index_CA150][0][_xyz] )
   r151 = array( hetchain1rot[_ratms][index_CA151][0][_xyz] )
   r152 = array( hetchain1rot[_ratms][index_CA152][0][_xyz] )

   print 'r149:', r149
   print 'r150:', r150
   print 'r151:', r151
   print 'r152:', r152


   CalciumE = 0.0
   d2 = ( ( r151 - r149 )**2 ).sum()
   if d2 > 0.0:
      dinv = 1./sqrt(d2)
      CalciumE += Eqq_pair( dinv, 2., 2.)
   else:
      CalciumE += 1000000000000000.

   d2 = ( ( r151 - r150 )**2 ).sum()
   if d2 > 0.0:
      dinv = 1./sqrt(d2)
      CalciumE += Eqq_pair( dinv, 2., 2.)
   else:
      CalciumE += 1000000000000000.

   d2 = ( ( r152 - r149 )**2 ).sum()
   if d2 > 0.0:
      dinv = 1./sqrt(d2)
      CalciumE += Eqq_pair( dinv, 2., 2.)
   else:
      CalciumE += 1000000000000000.

   d2 = ( ( r152 - r150 )**2 ).sum()
   if d2 > 0.0:
      dinv = 1./sqrt(d2)
      CalciumE += Eqq_pair( dinv, 2., 2.)
   else:
      CalciumE += 1000000000000000.

   print 'Degrees/CalciumE/n>dcut/n<dcut:', Degrees, CalciumE, ngtcut, nltcut
   savedata.append( (Degrees, CalciumE, ngtcut, nltcut) )
   hinge.write_pdb_modelstart()
   hinge.write_PDB_atoms( key )
   hinge.write_PDB_hets( key, het_names_to_skip_in_print )
   hinge.write_pdb_modelend()

hinge.write_pdb_end()


# write out data to file for post processing...

fhdat = open("interdomain_interaction_hinge_models_res" + str(RotatePhiResNum) + ".txt", 'w' )
print >> fhdat, '# Degrees, CalciumE, ngtcut, nltcut'
for data in savedata:
   (Degrees, CalciumE, ngtcut, nltcut) = data
   print >> fhdat, Degrees, CalciumE, ngtcut, nltcut

fhdat.close()

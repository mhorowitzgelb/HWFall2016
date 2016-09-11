#
# A set of demos for learning how to use the PDButils module, and a basic
#  introduction to Numpy and Matplotlib.
#
# To run this, you will need a PDB file which exists on your
#  computer and specify its name on the command line:
#
#   python ramachandran.py "pdb_filename"
#
# Directions for calculating the dihedral angles for the peptide backbone
# may be viewed at  the bottom of this file.

# Next two lines import modules.  Module "sys" is a standard module which
# comes with Python and is needed here in order to use the command line
# argument.  "PDButils" is a module for reading in PDB files and making the
# data contained in the files accessible to your Python programs.

import sys
from PDButils import *

# set shorthand for easier access to Residue lists:
# ( You can use any names you a want for the LHS, but pick something
#   short and unique looking to avoid name collsions, that's why I begin
#   the name with an underscore (_).  Note also that RHS terms have been
#   defined in the DATA section in the PDButils module. )

_rnums = Res_["resnums"]
_ratms = Res_["atoms"]
_rnams = Res_["names"]
_rptrs = Res_["resptrs"]

# set shorthand pointers for easier access to ATOM parameters:

_aname = Atm_["name"]
_arnum = Atm_["resnum"]
_xyz   = Atm_["coords"]
_B     = Atm_["biso"]

# some more, this time the first 5 atoms in amino acids, except
# glycine (GLY) will not have a Cbeta (CB).

_N  = 0
_CA = 1
_C  = 2
_O  = 3
_CB = 4

#
# The examples below are intended to demonstrate how to use the PDB data read
# using the PDButils module.  The locally defined shorthand items defined above
# will be used.
#

#-------------------------------------------------------------------------------------
# I. CREATE A SINGLE OBJECT CONTAINING ALL THE COORDINATES FROM A PDB FILE...
#-------------------------------------------------------------------------------------
print; print "( *** below output is from Demo Part I. *** )"; print

# pdbfile = "4cln.pdb"   # this file must exist in you current directory
if len( sys.argv ) < 2:
   print "Please give the name of a PDB file on the command line."
else:
   pdbfile = sys.argv[1]

# next line should read all the data from the PDB and store it in memory...
prot1 = PDBin( pdbfile )

if not prot1.stat:
   print "Failed to open file:", pdbfile
   sys.exit(1)
else:
   print "Opened and read file:", pdbfile



#-------------------------------------------------------------------------------------
# II. SELECT FIRST CHAIN OF FIRST MODEL...
#-------------------------------------------------------------------------------------
print; print "( *** below output is from Demo Part II. *** )"; print

(model_id, chain_id) = prot1.chainlist[0]
if model_id == '':
   print "The name of the first MODEL is: (MODEL records not used in this PDB file.)"
else:
   print "The name of the first MODEL is:", model_id

print "The name of the first CHAIN in MODEL is:", chain_id



#-------------------------------------------------------------------------------------
# III. DEFINE "chain1" TO REFER ALL THE DATA FOR THE FIRST CHAIN IN THE FIRST MODEL...
#-------------------------------------------------------------------------------------
print; print "( *** below output is from Demo Part III. *** )"; print

chain1 = prot1.chains[ (model_id, chain_id) ]

print len(chain1), "residues are in this chain."

print "coords of atom N in first residue of first chain:", chain1[_ratms][0][_N][_xyz]

#                                                                    ^    ^   ^    ^
#                                                                    |    |   |    |
#          select the residue list describing all atoms in  this chain    |   |    |
#                                                                         |   |    |
#                                        select the first residue from list   |    |
#                                                                             |    |
#               select block of info about atom with label "N" for this residue    |
#                                                                                  |
#                                             retrieve the coordinates for this atom



# 
#-------------------------------------------------------------------------------------
# IV. GET COORDINATES for ALL BACKBONE ATOMS FOR ONE RESIDUE, RESIDUE 0...
#-------------------------------------------------------------------------------------
print; print "( *** below output is from Demo Part IV. *** )"; print

print "backbone atoms for residue:", chain1[_rnams][0],  chain1[_rnums][0]
for atomid in ( _N, _CA, _C, _O ):
   print chain1[_ratms][0][atomid][_aname], chain1[_ratms][0][atomid][_xyz]



# 
#-------------------------------------------------------------------------------------
# V. LOOK FOR CHAIN GAPS (MISSING RESIDUES IN PDB MODEL)
#-------------------------------------------------------------------------------------
print; print "( *** below output is from Demo Part V. *** )"; print

resnumsinchain = chain1[_rptrs].keys()
iresnumprevious = resnumsinchain[0] - 1
Ngaps = 0
for iresnum in resnumsinchain:
   if iresnum-1 != iresnumprevious:
      iprev = chain1[_rptrs][iresnumprevious]
      i = chain1[_rptrs][iresnum]
      print "Sequence gap in model between residues",  chain1[_rnams][iprev], iresnumprevious, "and", \
            chain1[_rnams][i], iresnum
      Ngaps += 1
   iresnumprevious = iresnum

if Ngaps == 0:
   print "No gaps found in this chain."



# 
#-------------------------------------------------------------------------------------
# VI. COMPUTE PHI-PSI ANGLES FOR BACKBONE...
#-------------------------------------------------------------------------------------
print; print "( *** below output is from Demo Part V. *** )"; print


chainlab = prot1.chainlist[0][1]

from numpy import *
from math import acos, pi
         
Phiprev = -999.0
rad2deg = 180.0 / pi
iresnumprevious = -100000

# Save angles for plotting...
phiang = []
psiang = []
for iresnum in resnumsinchain:

   if iresnum-1 != iresnumprevious:
      # To calculate Phi(j) requires AA(j) and AA(j-1) ==> not defined for N-term
      # To calculate Psi(j) requires AA(j) and AA(j+1) ==> not defined for C-term
      # To calculate Omega(j) requires AA(j) and AA(j+1) ==> not defined for C-term

      j1 = chain1[_rptrs][iresnum]
      if Phiprev != -999.0:
         Psi = -999.0; Omega = -999.0
         print "Chain/Res/Phi/Psi/Omega: %s %s%4d %7.1f %7.1f  %8.2f" \
               %( chainlab, chain1[_rnams][j1], iresnumprevious, Phiprev, Psi, Omega )

      Phiprev = -999.0

   else:

      j2 = chain1[_rptrs][iresnum]

      N1  = array( chain1[_ratms][j1][_N][_xyz] )
      CA1 = array( chain1[_ratms][j1][_CA][_xyz] )
      C1  = array( chain1[_ratms][j1][_C][_xyz] )
      N2  = array( chain1[_ratms][j2][_N][_xyz] )
      CA2 = array( chain1[_ratms][j2][_CA][_xyz] )
      C2  = array( chain1[_ratms][j2][_C][_xyz] )
      
      bnd_N1_CA1 =  CA1 - N1
      bnd_CA1_C1 =  C1 - CA1
      bnd_C1_N2  =  N2 - C1
      bnd_N2_CA2 =  CA2 - N2
      bnd_CA2_C2 =  C2 - CA2

      #-----------------------------------------------------------------------
      # Phi: dihedral angle of C(i-1)-N(i)---Ca(i)-C(i)
      #-----------------------------------------------------------------------
      # ( this phi is actually for the NEXT residue )
      n_C1_N2_CA2 = cross( bnd_C1_N2, bnd_N2_CA2 )
      n_C1_N2_CA2 = n_C1_N2_CA2 / linalg.norm( n_C1_N2_CA2 )      

      n_N2_CA2_C2 = cross( bnd_N2_CA2, bnd_CA2_C2 )
      n_N2_CA2_C2 = n_N2_CA2_C2 / linalg.norm( n_N2_CA2_C2 )      

      Phi = rad2deg * acos( dot( n_C1_N2_CA2, n_N2_CA2_C2 ) )

      n1Xn2 = cross( n_C1_N2_CA2, n_N2_CA2_C2 )
      if dot(bnd_N2_CA2,n1Xn2) < 0: Phi = -Phi

      #-----------------------------------------------------------------------
      # Psi: dihedral angle of N(i)-Ca(i)---C(i)-N(i+1)
      #-----------------------------------------------------------------------

      n_CA1_C1_N2 = cross( bnd_CA1_C1, bnd_C1_N2 )
      n_CA1_C1_N2 = n_CA1_C1_N2 / linalg.norm( n_CA1_C1_N2 )     

      n_N1_CA1_C1 = cross(bnd_N1_CA1, bnd_CA1_C1)
      n_N1_CA1_C1 = n_N1_CA1_C1 / linalg.norm(n_N1_CA1_C1)
      #
      #  HOME WORK PROBLEM:  WRITE CODE TO CALCULATE Psi
      #
      ################################################################################

      Psi = -777.0  # set a temporary default value
      Psi = rad2deg * acos(dot(n_N1_CA1_C1,n_CA1_C1_N2))

      #-----------------------------------------------------------------------
      # Omega: dihedral angle of Ca(i)-C(i)===N(i+1)-Ca(i+1)
      #-----------------------------------------------------------------------
      # n_C1_N2_CA2 already defined for Phi
      # n_CA1_C1_N2 already defined for Psi

      Omega = rad2deg * acos( dot( n_CA1_C1_N2, n_C1_N2_CA2 ) )

      n1Xn2 = cross( n_CA1_C1_N2, n_C1_N2_CA2 )
      if dot(bnd_CA1_C1,n1Xn2) < 0: Omega = -Omega

      if abs(Omega) < 90.0 :
         print "Cis-peptide between:",  chain1[_rnams][j1], iresnumprevious, "and", \
               chain1[_rnams][j2], iresnum

      print "Chain/Res/Phi/Psi/Omega: %s %s%4d %7.1f %7.1f  %8.2f" \
            %( chainlab, chain1[_rnams][j1], iresnumprevious, Phiprev, Psi, Omega )


      # Save angles for plottin...
      phiang.append( Phiprev )
      psiang.append( Psi )

      j1 = j2
      Phiprev = Phi

   iresnumprevious = iresnum

print; print "( *** ALL Demos completed. *** )"; print


#
# Draw a Ramachandran plot of the data...
#
# NOTE: you will need Matplotlib installed to make the plot
#
x = []
y = []
i = 0
for phi in phiang:
   if phi > -181. and  psiang[i] > -181.:
      x.append(phi)
      y.append(psiang[i])
   i+= 1

from pylab import *
phipsiplot = subplot(111)
phipsiplot.scatter( array(x), array(y), s=40, c='b', marker='s', edgecolors='none' )
phipsiplot.set_xlim( -190.0, 190.0 )
phipsiplot.set_ylim( -190.0, 190.0 )
phipsiplot.set_xticks( [-180,-120,-60, 0,  60, 120, 180] )
phipsiplot.set_yticks( [-180,-120,-60, 0,  60, 120, 180] )
if model_id == '':
   phipsiplot.set_title("Ramachandran Plot for: " + pdbfile + "\n Chain " + chain_id )
else:
   phipsiplot.set_title("Ramachandran Plot for: " + pdbfile +
                         "\n Model " + str(model_id) + ", Chain " + chain_id )
grid(True)
xlabel(r"$\Phi$", fontsize=14)
ylabel(r"$\Psi$", fontsize=14)

show()

raw_input("Press Return key to exit.")


helppage = '''

NOTES FOR PHI, PSI, OMEGA:

 N-term                                                         
 
      i    No Phi(i)
     N    /             O                                   
      \  /   Psi(i)   //                         
       \      |      //                          
        \     V     //                           
         Ca---------C                              
        /            \               H            
       /   peptide -> \             /             
      H     bond       \ i+1       / 
           Omega(i)     N---------Ca   Psi(i+1)
                       /    ^      \  /            H
                      /     |       \             /
                     H    Phi(i+1)   \           /   
                                      C---------N i+2             O                          
                                     //    ^     \               //                           
                                    //     |      \ <-Phi(i+2)  //                            
                                   //   peptide    \           //                             
                                  O      bond      Ca---------C                               
                                      Omega(i+1)   /           \               H            
                                                  /  peptide -> \             /             
                                                 H    bond       \ i+3       /              
                                                    Omega(i+2)    N---------Ca              
                                                                  /    ^     \              
                                                                 /     |      \  <-- No Psi(i+3)             
                                                                H   Phi(i+3)   \            
                                                                                C---------O
                                                                               //           
                                                                              //            
                                                                             //      C-term       
                                                                             O               
ANGLE DEFINITIONS:

Phi(i): dihedral angle of C(i-1)-N(i)---Ca(i)-C(i)

        looking down N(i)---Ca(i) bond, Phi is the angle between the 
        projection of C(i-1):N(i) bond and the projection of Ca(i):C(i) bond
        (right hand rotation is positive, i.e., clockwise)

        To calculate Phi(i) requires C(i-1) ==> not defined for N-term.


Psi(i): dihedral angle of N(i)-Ca(i)---C(i)-N(i+1)

        looking down Ca(i):C(i) bond, Psi is the angle between the 
        projection of C(i-1):N(i) bond and the projection of Ca(i):C(i) bond
        (right hand rotation is positive, i.e., clockwise)

        To calculate Psi(i) requires N(i+1) ==> not defined for C-term.


Omega(i): dihedral angle of Ca(i)-C(i)===N(i+1)-Ca(i+1)
          Omega ~ 0 almost always,  cis-peptides do occur and will have ~180.
          Omega not defined for C-term, no peptide bond.


Use convention:  ( -180 < Angle <= 180 )

Use the four bonded atoms to calcualte dihedral about the center bond.
Define two normal vectors by taking the cross-products of middle bond against
each end bond.

Use the dot product of these two normal vectors to get |Angle| of the dihedral.

Determine sign of Angle using dot product of rotation bond against the vector
formed by thecross-product of the two normal vectors

'''

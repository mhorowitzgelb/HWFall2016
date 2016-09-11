#-------------------------------------------------------------------------------------------
# what: PDButils Python module
# for: Biochem 606 - Math Methods - Structural Biology
# when: Fall Term, 2010
# by: Gary Wesenberg, UW-Madison, BACTER Institute. (gary@biochem.wisc.edu)
#-------------------------------------------------------------------------------------------
# This file is contains Python code and is intended to be used as a "module" by other
# Python programs. By using the "import" command in your Python program, you may utitlize
# the data and tools contained in this module.  For example, by including the following
# line in your program:
#
#    from PDButils import *
#
# you may directly use the objects defined in this module.
#
#-------------------------------------------------------------------------------------------


#-------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------
# DATA: This section of this module defines some data objects.
#-------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------

# Below, aacodes1, is an example of a Python "dictionary".  The braces ( { } ) enclose the
# entries of the dictionary.  Notice the entries are separated by commas.  Each entry
# consists of a pair of objects separated by a colon ( : ). This is bascially a lookup table,
# with the "keys" being specified by left member of the pair, and the item which is retrieved
# being the right side member.  To use the dictionary to lookup the single letter amino acid
# code, say for example, tryptophan, you would write:
#
#     aacodes1['TRP']
#
# Also, note that a long line of code may be broken into multiple lines to improve
# readability using the backslash character. (no white space allowed after \)

aacodes1 = { 'ALA':'A', \
             'HIS':'H', \
             'TRP':'W', \
             'SER':'S', \
             'ASP':'D', \
             'VAL':'V', \
             'PRO':'P', \
             'MET':'M', \
             'THR':'T', \
             'GLU':'E', \
             'LEU':'L', \
             'PHE':'F', \
             'CYS':'C', \
             'ASN':'N', \
             'LYS':'K', \
             'ILE':'I', \
             'TYR':'Y', \
             'GLY':'G', \
             'GLN':'Q', \
             'ARG':'R', \
             'CSH':'c', \
             'MSE':'m' }

aacodes3 = { 'A':'ALA', \
             'H':'HIS', \
             'W':'TRP', \
             'S':'SER', \
             'D':'ASP', \
             'V':'VAL', \
             'P':'PRO', \
             'M':'MET', \
             'T':'THR', \
             'E':'GLU', \
             'L':'LEU', \
             'F':'PHE', \
             'C':'CYS', \
             'N':'ASN', \
             'K':'LYS', \
             'I':'ILE', \
             'Y':'TYR', \
             'G':'GLY', \
             'Q':'GLN', \
             'R':'ARG', \
             'c':'CSH', \
             'm':'MSE' }

#
# AtomPars is a Python "tuple", a simple ordered container of items.
# The contents of Python tuples can not be changed once defined. To retrieve
# the value of the third item you would write: AtomPars[2]
# 
AtomPars = ( "serialid", "atomnam", "altloc", "insertcode", "x", "y", "z", \
              "occ", "biso", "segid", "element", "charge" )
#
# These are the items stored in the atom list. Note that present implementaion
# does not allow you to change anything except 'coords'.  
#
#    ( SerialId, AtomNam, AltLoc, seqPOS, insertCode, coords,
#      occ, Biso, SegId, Element, Charge ) = atom

Atm_ = { "id":0, "name":1, "alt":2, "resnum":3, "insertcode":4, "coords":5, \
         "occ":6, "biso":7, "segid":8, "element":9, "charge":10 }

Res_ = { "resnums":0, "atoms":1, "names":2, "resptrs":3 }

HelixType = ( "Unknown", \
              "Right-handed alpha", \
              "Right-handed omega", \
              "Right-handed pi", \
              "Right-handed gamma", \
              "Right-handed 3-10", \
              "Left-handed alpha", \
              "Left-handed omega", \
              "Left-handed gamma", \
              "2-7 ribbon/helix", \
              "Polyproline" )


#-------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------
# FUNCTIONS: Here are some "functions" included in this "module".  These will be accessible
# from your program when you import this module.
#-------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------

def open_and_read_a_file( filename, whoshallisayiscalling, verbose = True ):
   '''
   Generic file open and reader with some fault checking.
   [ Here is my attempt to demonstrate how good code should check for error conditions
     and report back to the human in front the terminal what when wrong should there
     be a problem.  This is a bit tedious, but it is often well worth the investment if
     the code is likely to be used often.  Try to trap bad things in your code and let
     it be known when something unexpected or unwanted occurs. ]
   '''
   from os.path import isfile

   if type(filename).__name__ != 'str':
      if verbose:
         print '!ERR* filename passed to "open_and_read_a_file" was not a string.' 
         print "!ERR* [called from:", whoshallisayiscalling + ".]"
      return False
   else:
      f = filename.strip()
      if len(f) == 0:
         if verbose:
            print '!ERR* filename string passed to "open_and_read_a_file" was empty.' 
            print "!ERR* [called from:", whoshallisayiscalling + ".]"
         return False
      else:
         if not isfile(f):
            if verbose:
               print '!ERR* File not found. filename passed to "open_and_read_a_file" was:', f 
               print "!ERR* [called from:", whoshallisayiscalling + ".]"
            return False
         else:
            try:
               fh = open(f,'r')
               try:
                  lines = fh.readlines()
                  fh.close()
                  return lines
               except:
                  if verbose:
                     print "!ERR* Problem reading lines from file:", f
                     print '!ERR* ["open_and_read_a_file" called from:', whoshallisayiscalling + ".]"
                  return False
            except:
               if verbose:
                  print "!ERR* Unable to open file:", f
                  print '!ERR* ["open_and_read_a_file" called from:', whoshallisayiscalling + ".]"
                  print "!ERR* (check access permissions)"
               return False


def open_file_for_write( filename, whoshallisayiscalling, verbose = True ):
   '''
   Generic file open for writing with some fault checking.
   Note that the file will be over-written if it already exists.
   '''
   if type(filename).__name__ != 'str':
      if verbose:
         print '!ERR* filename passed to "open_file_for_write" was not a string.' 
         print "!ERR* [called from:", whoshallisayiscalling + ".]"
      return False
   else:
      f = filename.strip()
      if len(f) == 0:
         if verbose:
            print '!ERR* filename string passed to "open_file_for_write" was empty.' 
            print "!ERR* [called from:", whoshallisayiscalling + ".]"
         return False
      else:
         try:
            fh = open(f,'w')
            return fh
         except:
            if verbose:
               print "!ERR* Unable to open file for writing:", f
               print '!ERR* ["open_file_for_write" called from:', whoshallisayiscalling + ".]"
               print "!ERR* (check access permissions)"
            return False


def write_PDB( fileout_name, headlines, chainlist, chains, hetchainlist, hetchains, option = "ALL" ):
   '''
   In order to use this function the PDB data must be contained in the data structures
   defined in the PDButils class PDBin. Note that some of the elements of these structures
   may be "immutable", though the coordinates are not.
   '''

   fhnew = open_file_for_write( fileout_name, "write_PDB" )
   if not fhnew:
      print 'ERROR: function "write_PDB": Unable to open file for writing:', fileout_name
      return False

   for l in headlines:
      #if l[0:4] == "ATOM" or  l[0:5] == "MODEL": break
      #print >> fhnew, l.rstrip()
      print >> fhnew, l

   if  len( chainlist ) > 0:
      hitModel = False
      LastModel = "None"

      for chainid in chainlist:
     
         ( model_Id, chainLabel ) = chainid
         if model_Id != '':
            if hitModel:
               if LastModel != model_Id:
                  print >> fhnew, "ENDMDL"
                  print >> fhnew, "MODEL     %4d" %model_Id
            else:
               print >> fhnew, "MODEL     %4d" %model_Id
               hitModel = True
            LastModel = model_Id         

         ( SEQpositions, resAtoms, resNames, SEQpositions_ptr ) = chains[chainid]
         ith = 0
         for resnum in SEQpositions:
            resname = resNames[ith]
            for atom in resAtoms[ith]:
     
               ( SerialId, AtomNam, AltLoc, seqPOS, insertCode, coords, \
                 occ, Biso, SegId, Element, Charge ) = atom
     
               ATOMfrmtPDB = "ATOM  %5d %-4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s%2s%2s" 
     
               ATOMdata = ( SerialId, AtomNam, AltLoc, resname, chainLabel, \
                            seqPOS, insertCode, coords[0],coords[1],coords[2], \
                            occ, Biso, SegId, Element, Charge )
     
               print >> fhnew, ATOMfrmtPDB %ATOMdata
     
            ith += 1
         TERdata = ( SerialId+1, resname, chainLabel,seqPOS, insertCode )
         TERfrmtPDB = "TER   %5d      %3s %1s%4d%1s" 
         print >> fhnew, TERfrmtPDB %TERdata
     
      if hitModel:
         print >> fhnew, "ENDMDL"
  
   if len( hetchainlist ) > 0:
      hitModel = False
      LastModel = "None"
      for chainid in hetchainlist:
         (model_Id, chainLabel) = chainid
     
         if model_Id != '':
            if hitModel:
               if LastModel != model_Id:
                  print >> fhnew, "ENDMDL"
                  print >> fhnew, "MODEL     %4d" %model_Id
            else:
               print >> fhnew, "MODEL     %4d" %model_Id
               hitModel = True
            LastModel = model_Id         
     
         for ( SEQpositions, resAtoms, resNames, SEQpositions_ptr ) in hetchains[chainid]:
            #print "CHAIN", chainid
            #print ( SEQpositions, resAtoms, resNames, SEQpositions_ptr )
            ith = 0
            for resnum in SEQpositions:
               resname = resNames[ith]
               #print 'resname:', resname, "ith:", ith
               #print 'resAtoms[ith]:', resAtoms[ith]
               for atom in resAtoms[ith]:
                  print >> fhnew,  'resAtoms[ith]:', resAtoms[ith]

                  ( SerialId, AtomNam, AltLoc, seqPOS, insertCode, coords, \
                    occ, Biso, SegId, Element, Charge ) = atom
     
                  HETATMfrmtPDB = "HETATM%5d %-4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s%2s%2s" 
                  HETATMdata = ( SerialId, AtomNam, AltLoc, resname, chainLabel, \
                                 seqPOS, insertCode, coords[0],coords[1],coords[2], \
                                 occ, Biso, SegId, Element, Charge )
     
                  print >> fhnew, HETATMfrmtPDB %HETATMdata
     
               ith += 1

      if hitModel:
         print >> fhnew, "ENDMDL"
   #
   # WARNING NOTE! There are likely more lines toward the bottom of the original PDB
   # file which are currently being left out... these might be important in some
   # cases.
   #
   print >> fhnew, "END"

   fhnew.close()
   return True


def define_versor( angle, axis):
   '''
   Define a rotation matrix to rotate about vector "axis" by "angle" radians.
   Very cool, read all about "versors" in :
   "Vector Analysis", J. Willard Gibbs / E. B. Wilson, Dover
   also at http://www.archive.org/stream/117714283#page/338/mode/2up/search/versor
   '''
   import math
   cosq = math.cos(ang)
   cosqm1 = 1.0 - cosq
   sinq = math.sin(ang)
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

def define_versor( angle, axis):
   '''
   Define a rotation matrix to rotate about vector "axis" by "angle" radians.
   Very cool, read all about "versors" in :
   "Vector Analysis", J. Willard Gibbs / E. B. Wilson, Dover
   also at http://www.archive.org/stream/117714283#page/338/mode/2up/search/versor
   '''
   import math
   cosq = math.cos(ang)
   cosqm1 = 1.0 - cosq
   sinq = math.sin(ang)
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


def getAtomPars( mylist, atomsdata ):
   '''
   ( SerialId, AtomNam, AltLoc, seqPOS, insertCode, coords,
     occ, Biso, SegId, Element, Charge ) = atom
   '''
   orderAtomPars = { "serialid":0, "atomnam":1, "altloc":2, "insertcode":3, "seqpos":4, \
                     "coords":5, "occ":6, "biso":7, "segid":8, "element":9, "charge":10 }

   allowed = orderAtomPars.keys()

   if len(mylist) == 0: 
         print 'Error in call to "getAtomPars", EMPTY "mylist" set.'
         return False

   coords = atomsdata[orderAtomPars["coords"]]
   print "len atomdata:", len(atomsdata)
   wantlow = []
   for want in mylist:
      wantlowercase = want.lower()
      wantlow.append(wantlowercase)
      if wantlowercase in allowed:
         command = wantlowercase + ' = atomsdata[orderAtomPars["' + wantlowercase + '"]]'
         #print command
         exec( command )
      elif "x" in wantlowercase: x = coords[0]
      elif "y" in wantlowercase: y = coords[1]
      elif "z" in wantlowercase: z = coords[2]

      else:
         print 'Error in call to "getAtomPars", unknown atom item set to UNDEFINED:', want
         command = wantlowercase + ' = "UNDEFINED"'
         exec( command )

   command = "goods = ["
   first = True
   for w in wantlow:
      if first:
         command += w
         first = False
      else:
         command += ', ' + w
   command += ']'

   #print command
   
   exec(command)
   return goods

#-------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------
# CLASSES: Class definitions begin here.  A class my be "instantiated" from another piece
# of code which "imports" this module.  
#-------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------

class ModelBasic(object):
   def __init__( self, modelnam ):
      self.name = modelnam

class ChainBasic(object):
   def __init__( self, chainnam ):
      self.name = chainnam

class AtomBasic(object):
   def __init__( self, atnam, coords, Bisotropic ):
      self.name = atnam
      self.x = coords[0]
      self.y = coords[1]
      self.z = coords[2]
      self.B = Bisotropic

class ResBasic(object):
   def __init__( self, resnum, resnam ):
      self.snum = resnum
      self.name = resnam
      self.atom = []

class PDBin(object):
   '''
   PURPOSE: For reading in standard PDB files.  This class attempts to follow the PDB 
            format specified at: http://www.wwpdb.org/documentation/format32/v3.2.html
            Some other programs will write out files in PDB-like format, but they may
            not strictly adhere to the PDB format.  In such cases, the methods in this 
            class may produce undesired results.  Please beware.
   '''

   def __init__( self, filename=None, Options="All" ):
      '''
      Check it out! Get all the goods from a PDB file into data structures.
      '''
      # save name of file for possible future use...
      self.filename = filename

      # this is just a fancy way to get the name of 'this' location for debugging...
      _classname = self.__class__.__name__

      # if all goes well, lines will return as an array containig ALL the lines in the file... 
      self.lines = open_and_read_a_file( filename, _classname )

      if self.lines:
         self.stat = True
         if Options == "All":
            self._get_pdb_header_records( )
            self._secondary_structure( )
            self._read_polymer_atoms( )
            self._read_hetatoms( )
            #self._sequences( )
      else:
         self.stat = False

   def _read_polymer_atoms( self ):
      '''
      Note the strict formatting of PDB file!  

      Keywrd  Atm Atm  Res Chn Res       x       y      z     occ    B            Atm
              Num Name Nam Nam Number                                             Type


      ATOM    803 HH21 ARG A  42       1.502  19.627  10.115  1.00 20.00           H  
      ATOM    804 HH22 ARG A  42       1.720  20.400   8.585  1.00 20.00           H  
      012345678901234567890123456789012345678901234567890123456789012345678901234567890
                                    123456781234567812345678123456123456
      ATOM    805  N   GLY A  43       6.741  14.655  13.200  1.00 16.96           N  
      ATOM    806  CA  GLY A  43       7.998  14.953  13.849  1.00 16.03           C  
      ATOM    807  C   GLY A  43       9.077  14.061  13.293  1.00 15.55           C  
      ATOM    808  O   GLY A  43       8.765  13.064  12.655  1.00 12.92           O  
      ATOM    809  H   GLY A  43       6.685  13.836  12.649  1.00 20.00           H  
      ATOM    810  N   ALA A  44      10.332  14.398  13.575  1.00 14.66           N  
      ATOM    811  CA  ALA A  44      11.467  13.621  13.100  1.00 15.80           C  
      ATOM    812  C   ALA A  44      12.677  13.951  13.924  1.00 16.35           C  
      ATOM    813  O   ALA A  44      12.825  15.065  14.438  1.00 17.58           O  
      ATOM    814  CB  ALA A  44      11.770  13.925  11.614  1.00 14.67           C  
      ATOM    815  H   ALA A  44      10.512  15.200  14.111  1.00 20.00           H  
      ATOM    816  N   VAL A  45      13.558  12.975  14.032  1.00 14.61           N  
      ATOM    817  CA  VAL A  45      14.792  13.145  14.766  1.00 17.80           C  
      '''   
      _atomsinres = ''
      _SEQpositions = ''
      _SEQpositions_ptr = ''
      _chainKEYlast = 'INIT'      
      _model_id = ''
      _model_id_last = ''

      self.chains={}
      self.chainlist=[]

      # read in all ATOMS

      for _l in self.lines:

         if _l[0:5] == "MODEL":
            _model_id = int(_l[5:].strip())

         elif _l[0:4] == "ATOM":

            # Just hit an ATOM record.  Must now:
            # 1) parse line with current atom info
            # 2) check chain label and see if it's a new one for this model
            #     2a) if new chain label
            #         3i) if have previous saved atom,  must add it to previous resAtoms
            #         3ii) must add last residue to previous chain
            #         3iii) initialize variables for next chain

            _SerialId =  int(_l[6:11])
            _AtomNam = _l[12:16].strip()
            _AltLoc  = _l[16].strip()
            _resnam = _l[17:20].strip()
            _chain = _l[21]
            _chainKEY = ( _model_id, _chain )
            _seqPOS = int(_l[22:26])
            _insertCode = _l[26]
            _coords = [ float(_l[30:38]), float(_l[38:46]), float(_l[46:54]) ]
            _occ = float(_l[54:60])
            _Biso = float(_l[60:66])
            _SegId = _l[72:76].strip()
            _Element = _l[76:78].strip()
            _Charge = _l[78:80].strip()

            # Check for change in Chain label.  If new chain started, save last one, if needed.
            if _chainKEY not in self.chainlist:
               
               # save new chainkey
               self.chainlist.append(_chainKEY)
            
               if _chainKEYlast != 'INIT' and _chainKEYlast not in self.chains:
                  if _atomsinres != '':
                     _resAtoms.append( _atomsinres )
                  self.chains[_chainKEYlast] =  ( _SEQpositions, _resAtoms, _resNames, _SEQpositions_ptr )

               _NewChain = True
               _chainKEYlast = _chainKEY
               _atomsinres = ''

               # Three arrays indexed to same residue within a chain
               # (Keep in mind that there may be chain breaks.)
               _resNames = []
               _resAtoms = []
               _SEQpositions = []

               # Could always derive this dictionary if needed, but included here now, for simplicity.
               # (Use this to directly access a particular residue on the fly.)
               _SEQpositions_ptr = {}
               _iSEQptr = 0
      
            if _seqPOS not in _SEQpositions:
               if _NewChain: 
                   _resAtoms = []
                   _NewChain = False
               else:
                   _resAtoms.append( _atomsinres )

               _SEQpositions.append( _seqPOS )
               _SEQpositions_ptr[_seqPOS] = _iSEQptr
               _iSEQptr += 1
               _resNames.append( _resnam )
               _atomsinres = [ ( _SerialId, _AtomNam, _AltLoc, _seqPOS, \
                             _insertCode, _coords, _occ, _Biso, _SegId, _Element, _Charge ) ] 
            else:
               _atomsinres.append(( _SerialId, _AtomNam, _AltLoc, _seqPOS, \
                             _insertCode, _coords, _occ, _Biso, _SegId, _Element, _Charge ))

         elif _l[0:6] == "ENDMDL":
            _resAtoms.append( _atomsinres )
            self.chains[_chainKEY] =  ( _SEQpositions, _resAtoms, _resNames, _SEQpositions_ptr )

      if _model_id == '':
         _resAtoms.append( _atomsinres )
         self.chains[_chainKEY] =  ( _SEQpositions, _resAtoms, _resNames, _SEQpositions_ptr )


   def _read_hetatoms( self ):
      '''
      This method is modeled after "_read_polymer_atoms".

      Note the strict formatting of PDB file!  

      Keywrd  Atm Atm  Res Chn Res       x       y      z     occ    B            Atm
              Num Name Nam Nam Number                                             Type

      HETATM 7173  C1  NAG A 253      43.384 -12.964  33.458  0.50  9.44           C
      HETATM 7174  C2  NAG A 253      43.585 -14.438  33.196  0.50 15.04           C
      HETATM 7175  C3  NAG A 253      44.987 -14.832  33.638  0.50 14.32           C
      HETATM 7176  C4  NAG A 253      46.043 -13.945  33.025  0.50 12.96           C
      HETATM 7177  C5  NAG A 253      45.684 -12.459  33.246  0.50 15.18           C
      HETATM 7178  C6  NAG A 253      46.593 -11.482  32.536  0.50 19.74           C
      ...
      HETATM 7186  O7  NAG A 253      41.763 -16.326  32.220  0.50 14.18           O
      HETATM 7187 MN    MN A 254      14.452 -16.518  50.661  1.00 36.68          MN
      HETATM 7188 CA    CA A 255      17.001 -19.887  52.040  1.00 27.83          CA
      012345678901234567890123456789012345678901234567890123456789012345678901234567890
                                    123456781234567812345678123456123456
      '''   
      _atomsinres = ''
      _SEQpositions = ''
      _SEQpositions_ptr = ''
      _chainKEYlast = 'INIT'      
      _model_id = ''

      self.hetchains={}
      self.hetchainlist=[]

      # read in all ATOMS
      for _l in self.lines:

         if _l[0:5] == "MODEL":
            _model_id = int(_l[5:].strip())

         elif _l[0:6] == "HETATM":
            _SerialId = int(_l[6:11])
            _AtomNam = _l[12:16].strip()
            _AltLoc = _l[16].strip()
            _resnam = _l[17:20].strip()
            _chain = _l[21]
            _seqPOS = int(_l[22:26])
            _insertCode = _l[26]
            _coords = [ float(_l[30:38]), float(_l[38:46]), float(_l[46:54]) ]
            _occ = float(_l[54:60])
            _Biso = float(_l[60:66])
            _SegId = _l[72:76].strip()
            _Element = _l[76:78].strip()
            _Charge = _l[78:80].strip()

            _chainKEY = ( _model_id, _chain )

            # Check for change in Chain label.  If new chain started, save last one, if needed.

            if _chainKEY not in self.hetchainlist:

               self.hetchainlist.append(_chainKEY)

               if _atomsinres != '':
                  _resAtoms.append( _atomsinres )

               if _chainKEYlast != 'INIT':
                  if _chainKEYlast not in self.hetchains:
                     # In case was the first residue in previous chain...
                     self.hetchains[_chainKEYlast] = [( _SEQpositions, _resAtoms, _resNames, _SEQpositions_ptr )]
                  else:
                     self.hetchains[_chainKEYlast].append(( _SEQpositions, _resAtoms, _resNames, _SEQpositions_ptr ))

               _NewChain = True
               _chainKEYlast = _chainKEY
               _atomsinres = ''

               # Three arrays indexed to same residue within a chain
               # (Keep in mind that there may be chain breaks.)
               _resNames = []
               _resAtoms = []
               _SEQpositions = []

               # Could always derive this dictionary if needed, but included here now, for simplicity.
               # (Use this to directly access a particular residue on the fly.)
               _SEQpositions_ptr = {}
               _iSEQptr = 0

            else:

               if _chainKEY != _chainKEYlast:
                  # we've seen this chain before, so save previous one with last
                  if _atomsinres != '':
                     _resAtoms.append( _atomsinres )

                  if _chainKEYlast not in self.hetchains:
                     # In case was the first residue in previous chain...
                     self.hetchains[_chainKEYlast] = [( _SEQpositions, _resAtoms, _resNames, _SEQpositions_ptr )]
                  else:
                     self.hetchains[_chainKEYlast].append(( _SEQpositions, _resAtoms, _resNames, _SEQpositions_ptr ))

                  _chainKEYlast = _chainKEY
                  _atomsinres = ''
                  # Three arrays indexed to same residue within a chain
                  # (Keep in mind that there may be chain breaks.)
                  _resNames = []
                  _resAtoms = []
                  _SEQpositions = []

                  # Could always derive this dictionary if needed, but included here now, for simplicity.
                  # (Use this to directly access a particular residue on the fly.)
                  _SEQpositions_ptr = {}
                  _iSEQptr = 0
                  _NewChain = True

            if _seqPOS not in _SEQpositions:
               if _NewChain: 
                   _resAtoms = []
                   _NewChain = False
               else:
                   _resAtoms.append( _atomsinres )

               _SEQpositions.append( _seqPOS )
               _SEQpositions_ptr[_seqPOS] = _iSEQptr
               _iSEQptr += 1
               _resNames.append( _resnam )
               _atomsinres = [ ( _SerialId, _AtomNam, _AltLoc, _seqPOS, \
                             _insertCode, _coords, _occ, _Biso, _SegId, _Element, _Charge ) ] 
            else:
               _atomsinres.append(( _SerialId, _AtomNam, _AltLoc, _seqPOS, \
                             _insertCode, _coords, _occ, _Biso, _SegId, _Element, _Charge ))

         elif _l[0:6] == "ENDMDL" and _atomsinres != '':
            _resAtoms.append( _atomsinres )
            self.hetchains[_chainKEY].append(( _SEQpositions, _resAtoms, _resNames, _SEQpositions_ptr ))

      if _model_id == '' and _atomsinres != '':
         _resAtoms.append( _atomsinres )
         if _chainKEY not in self.hetchains:
            self.hetchains[_chainKEY] = [( _SEQpositions, _resAtoms, _resNames, _SEQpositions_ptr )]
         else:
            self.hetchains[_chainKEY].append(( _SEQpositions, _resAtoms, _resNames, _SEQpositions_ptr ))


   def printall( self ):
      '''
      DEBUG: dump all lines read from file
      '''
      for _l in self.lines:
         print _l.rstrip()

   def _get_pdb_header_records( self ):
      self.header_recs = {}
      self.header_keywrd_order = []
      for _l in self.lines:
         _keywrd = _l[0:6]
         if _keywrd == 'ATOM  ': break
         if _keywrd in self.header_recs:
             self.header_recs[_keywrd].append( _l[6:].rstrip() )
         else:
             self.header_recs[_keywrd] = [ _l[6:].rstrip() ]
             self.header_keywrd_order.append(_keywrd)

   def rewrite_pdb_header_records(self, skip_key_list ):
      self.headerlines = []
      for _key in self.header_keywrd_order:
         if _key in skip_key_list: continue
         for _l in self.header_recs[_key]:
            self.headerlines.append( _key + _l )

   def add_to_TITLE( self, titleplus ):
      _Nchars = len( titleplus )
      # Sorry, only room for 80 columns in PDB format
      _Nchar = min(_Nchars, 69 )
      if "TITLE " not in self.header_recs:
         # currently there is no TITLE recores in PDB header, will add one...
         if len(  self.header_keywrd_order ) == 0:
            self.header_keywrd_order.append("TITLE ")
         elif "TITLE " not in self.header_keywrd_order:
            self.header_keywrd_order.insert(1,"TITLE ")
         self.header_recs["TITLE "].append('    ' + titleplus[:_Nchars])
      else:
         _next = len(self.header_recs["TITLE "]) + 1
         if _next > 99: return 1
         _continuation_str = "  %2d" %_next
         self.header_recs["TITLE "].append( _continuation_str + titleplus[:_Nchars] )

      
   def showChains( self ):
      print "CHAINS from file", self.filename
      for _chainid in self.chainlist:
         (_model_Id, _chainLabel) = _chainid
         if _model_Id != '':
            print "Model:", _model_Id,",  Polymer Chain:", _chainLabel
         else:
            print "Model: (Null),  Polymer Chain:", _chainLabel


   def _helix_parse( self, _l ):
      '''
      Extract the Helix SSE from PDB header.  Here's what it looks like:

      HELIX    1   1 THR A   17  LYS A   19  5                                   3    
      HELIX    2   2 GLU A   47  ASP A   49  5                                   3    
      01234567890123456789012345678901234567890
      HELIX    3   3 LYS A   52  ARG A   63  1                                  12    
      SHEET    1   A 2 LYS A   3  TYR A   8  0                                        
      SHEET    2   A 2 GLU A  11  ASP A  16 -1  N  VAL A  15   O  VAL A   4           
      01234567890123456789012345678901234567890123456789012345678901234567890
      SHEET    1   B 3 VAL A  23  VAL A  26  0                                        
      SHEET    2   B 3 MET A  29  ASP A  36 -1  N  SER A  31   O  TRP A  24           
      SHEET    3   B 3 LYS A  39  SER A  46 -1  N  VAL A  45   O  VAL A  30           

      NOTE, some PDB files will not include H-bond info for sheets...
      '''   
      _hnum      = int(_l[7:10])
      _hname     = _l[11:13].strip()

      _aa1resnam = _l[15:18]
      _aa1chain  = _l[19].strip()
      _aa1resn   = int(_l[21:25])

      _aa2resnam = _l[27:30].strip()
      _aa2chain  = _l[31].strip()
      _aa2resn   = int(_l[33:37])

      _type      = int( _l[38:40] )
      _comment   = _l[40:70].strip()

      return  _hnum, (_hname, _aa1resnam, _aa1chain, _aa1resn, _aa2resnam, \
              _aa2chain, _aa2resn, _type, _comment)

   def _strand_parse( self, _l ):
      _strandid        = _l[7:10].strip()
      _sheetid         = _l[11:14].strip()    
      _nstrandsinsheet = int(_l[14:16])
      _aa1name         = _l[28:31].strip()
      _aa1chain        = _l[21].strip()
      _aa1resn         = int(_l[22:26])
      _aa2name         = _l[28:31].strip()
      _aa2chain        = _l[32].strip()
      _aa2resn         = int(_l[33:37])
      _strandsense     = int(_l[38:40])
      return _sheetid, (_strandid, _nstrandsinsheet, _aa1name, _aa1chain, _aa1resn, \
             _aa2name, _aa2chain, _aa2resn, _strandsense)

   def _secondary_structure( self ):

      self.helices = {}
      self.helixids = []
      self.sheets = {}
      self.sheetids = []

      for _l in self.lines:
         if _l[0:5] == "HELIX":
            _hnum, _helix = self._helix_parse( _l )
            if _hnum not in self.helices:
               self.helices[_hnum] = _helix
               self.helixids.append( _hnum )
         
         elif _l[0:5] == "SHEET":
            _sheetid, _strand = self._strand_parse( _l )
            if _sheetid not in self.sheets:
               self.sheets[_sheetid] = []
               self.sheetids.append(_sheetid)

            self.sheets[_sheetid].append( _strand )

         elif _l[0:4] == "ATOM":
            return

   def showsecondary( self ):
      self.showSecondary( )

   def showSecondary( self ):

      _Nsheets = len( self.sheetids )
      if _Nsheets == 0:
         print "No SHEET records found in PDB file."
      else:
         print "SHEETS:",_Nsheets, "defined in PDB SHEET header records."

         for _sheetid in self.sheetids:
            print "  sheet", _sheetid
            for _strand in self.sheets[_sheetid]:
               _strandid, _nstrandsinsheet, _aa1name, _aa1chain, _aa1resn, \
                 _aa2name, _aa2chain, _aa2resn, _strandsense = _strand
               print "      strand", _strandid,":",
               if _aa1name != '':
                  print _aa1name, _aa1chain, _aa1resn,  _aa2name, _aa2chain, _aa2resn, _strandsense
               else:
                  print _aa1name, _aa1chain, _aa1resn, _strandsense

      _Nhelices = len( self.helices )
      if _Nhelices == 0:
         print "No HELIX records found in PDB file header records."
      else:
         print "HELICES:",_Nhelices, "defined in PDB HELIX header records."
         for _helixid in self.helixids:
            (_hname, _aa1name, _aa1chain, _aa1resn, _aa2name, _aa2chain, _aa2resn, _type, _comment) = self.helices[_helixid]
            print "  helix", _helixid, _aa1name, _aa1chain, _aa1resn,  _aa2name, \
                   _aa2chain, _aa2resn, HelixType[_type], _comment
            # save list of residues in helices...
            #if aa1chain not in helix_all:
            #    helix_all[aa1chain] = [(aa1resn,aa2resn)]
            #else:
            #    helix_all[aa1chain].append([(aa1resn,aa2resn)])

   def set_atom_coords( self, key, nthres, nthatmm, x, y, z ):
      pass


#-------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------
# MAIN: The code below demonstratesis a niffty Python way to test all the code contained
# in this module. Typically you would never directly run a module, but you can!  The "if"
# statement below detects the case when you directly pass this file to the Python interpreter
# and any subsequent code which follows will be executed.  On the otherhand, if this
# file is included in another program using the "import" mechanism, then the code
# is ignored.  This is an ideal set-up for developing, debugging, and testing new modules.
#-------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------

if __name__ == "__main__":
    
   import sys
   A = PDBin(sys.argv[1])
   if not A.stat:
       print "!!! ERROR: Failed to open/read PDB file."
       sys.exit(1)
  
   print "opened and read file:", sys.argv[1]
   print "test of writing out new PDB file from data just read..."
   print " will write to file: test_PDButil.pdb"

   A.add_to_TITLE("this is a revised test write of pdb from class PDButil")
   A.rewrite_pdb_header_records( ["REMARK"] )
   write_PDB( "test_PDButil.pdb", A.headerlines, A.chainlist, A.chains, \
               A.hetchainlist, A.hetchains )
  

#-------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------
# MISCELLANEOUS DOCUMENTAION:
#-------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------


'''
Details about ATOM records in PDB
-----------------------------------------------------------------------------

see: http://www.wwpdb.org/documentation/format32/v3.2.html
also see: http://deposit.rcsb.org/format-faq-v1.html

ATOM record fields form PDB

COLUMNS        DATA TYPE       CONTENTS                            
--------------------------------------------------------------------------------
 1 -  6        Record name     "ATOM  "                                            
 7 - 11        Integer         Atom serial number.                   
13 - 16        Atom            Atom name.                            
17             Character       Alternate location indicator.         
18 - 20        Residue name    Residue name.                         
22             Character       Chain identifier.                     
23 - 26        Integer         Residue sequence number.              
27             AChar           Code for insertion of residues.       
31 - 38        Real(8.3)       Orthogonal coordinates for X in Angstroms.                       
39 - 46        Real(8.3)       Orthogonal coordinates for Y in Angstroms.                            
47 - 54        Real(8.3)       Orthogonal coordinates for Z in Angstroms.                            
55 - 60        Real(6.2)       Occupancy.                            
61 - 66        Real(6.2)       Temperature factor (Default = 0.0).                   
73 - 76        LString(4)      Segment identifier, left-justified.   
77 - 78        LString(2)      Element symbol, right-justified.      
79 - 80        LString(2)      Charge on the atom.     



see: http://deposit.rcsb.org/adit/docs/pdb_atom_format.html )
also: http://deposit.rcsb.org/format-faq-v1.html

COLUMNS        DATA TYPE       CONTENTS                            
--------------------------------------------------------------------------------
 1 -  6        Record name     "ATOM  "                                            
 7 - 11        Integer         Atom serial number.                   
13 - 16        Atom            Atom name.                            
17             Character       Alternate location indicator.         
18 - 20        Residue name    Residue name.                         
22             Character       Chain identifier.                     
23 - 26        Integer         Residue sequence number.              
27             AChar           Code for insertion of residues.       
31 - 38        Real(8.3)       Orthogonal coordinates for X in Angstroms.                       
39 - 46        Real(8.3)       Orthogonal coordinates for Y in Angstroms.                            
47 - 54        Real(8.3)       Orthogonal coordinates for Z in Angstroms.                            
55 - 60        Real(6.2)       Occupancy.                            
61 - 66        Real(6.2)       Temperature factor (Default = 0.0).                   
73 - 76        LString(4)      Segment identifier, left-justified.   
77 - 78        LString(2)      Element symbol, right-justified.      
79 - 80        LString(2)      Charge on the atom.     





THIS IS THE BEST:  http://www.wwpdb.org/documentation/format32/v3.2.html



                                     CLASS NUMBER
TYPE OF  HELIX                     (COLUMNS 39 - 40)
--------------------------------------------------------------
Right-handed alpha (default)                1
Right-handed omega                          2
Right-handed pi                             3
Right-handed gamma                          4
Right-handed 3 - 10                         5
Left-handed alpha                           6
Left-handed omega                           7
Left-handed gamma                           8
2 - 7 ribbon/helix                          9
Polyproline                                10


Overview

HELIX records are used to identify the position of helices in the molecule. Helices are named, numbered, and classified by type. The residues where the helix begins and ends are noted, as well as the total length.

Record Format

COLUMNS        DATA  TYPE     FIELD         DEFINITION
-----------------------------------------------------------------------------------
 1 -  6        Record name    "HELIX "
 8 - 10        Integer        serNum        Serial number of the helix. This starts
                                            at 1  and increases incrementally.
12 - 14        LString(3)     helixID       Helix  identifier. In addition to a serial
                                            number, each helix is given an 
                                            alphanumeric character helix identifier.
16 - 18        Residue name   initResName   Name of the initial residue.
20             Character      initChainID   Chain identifier for the chain containing
                                            this  helix.
22 - 25        Integer        initSeqNum    Sequence number of the initial residue.
26             AChar          initICode     Insertion code of the initial residue.
28 - 30        Residue  name  endResName    Name of the terminal residue of the helix.
32             Character      endChainID    Chain identifier for the chain containing
                                            this  helix.
34 - 37        Integer        endSeqNum     Sequence number of the terminal residue.
38             AChar          endICode      Insertion code of the terminal residue.
39 - 40        Integer        helixClass    Helix class (see below).
41 - 70        String         comment       Comment about this helix.
72 - 76        Integer        length        Length of this helix.

Details

    * Additional HELIX records with different serial numbers and identifiers occur if more than one helix is present.
    * The initial residue of the helix is the N-terminal residue.
    * Helices are classified as follows:

Overview

SHEET records are used to identify the position of sheets in the molecule. Sheets are both named and numbered. The residues where the sheet begins and ends are noted.

Record Format

COLUMNS       DATA  TYPE     FIELD          DEFINITION
-------------------------------------------------------------------------------------
 1 -  6        Record name   "SHEET "
 8 - 10        Integer       strand         Strand  number which starts at 1 for each
                                            strand within a sheet and increases by one.
12 - 14        LString(3)    sheetID        Sheet  identifier.
15 - 16        Integer       numStrands     Number  of strands in sheet.
18 - 20        Residue name  initResName    Residue  name of initial residue.
22             Character     initChainID    Chain identifier of initial residue 
                                            in strand. 
23 - 26        Integer       initSeqNum     Sequence number of initial residue
                                            in strand.
27             AChar         initICode      Insertion code of initial residue
                                            in  strand.
29 - 31        Residue name  endResName     Residue name of terminal residue.
33             Character     endChainID     Chain identifier of terminal residue.
34 - 37        Integer       endSeqNum      Sequence number of terminal residue.
38             AChar         endICode       Insertion code of terminal residue.
39 - 40        Integer       sense          Sense of strand with respect to previous
                                            strand in the sheet. 0 if first strand,
                                            1 if  parallel,and -1 if anti-parallel.
42 - 45        Atom          curAtom        Registration.  Atom name in current strand.
46 - 48        Residue name  curResName     Registration.  Residue name in current strand
50             Character     curChainId     Registration. Chain identifier in
                                            current strand.
51 - 54        Integer       curResSeq      Registration.  Residue sequence number
                                            in current strand.
55             AChar         curICode       Registration. Insertion code in
                                            current strand.
57 - 60        Atom          prevAtom       Registration.  Atom name in previous strand.
61 - 63        Residue name  prevResName    Registration.  Residue name in
                                            previous strand.
65             Character     prevChainId    Registration.  Chain identifier in
                                            previous  strand.
66 - 69        Integer       prevResSeq     Registration. Residue sequence number
                                            in previous strand.
70             AChar         prevICode      Registration.  Insertion code in
                                            previous strand.

Details

    * The initial residue for a strand is its N-terminus. Strand registration information is provided in columns 39 - 70. Strands are listed starting with one edge of the sheet and continuing to the spatially adjacent strand.
    * The sense in columns 39 - 40 indicates whether strand n is parallel (sense = 1) or anti-parallel (sense = -1) to strand n-1. Sense is equal to zero (0) for the first strand of a sheet.
    * The registration (columns 42 - 70) of strand n to strand n-1 may be specified by one hydrogen bond between each such pair of strands. This is done by providing the hydrogen bonding between the current and previous strands. No register information should be provided for the first strand.
    * Split strands, or strands with two or more runs of residues from discontinuous parts of the amino acid sequence, are explicitly listed. Detail description can be included in the REMARK 700 .

Relationships to Other Record Types

If the entry contains bifurcated sheets or beta-barrels, the relevant REMARK 700 records must be provided.
See the REMARK section for details.

Examples

         1         2         3         4         5         6         7         8
12345678901234567890123456789012345678901234567890123456789012345678901234567890
SHEET    1   A 5 THR A 107  ARG A 110  0
SHEET    2   A 5 ILE A  96  THR A  99 -1  N  LYS A  98   O  THR A 107
SHEET    3   A 5 ARG A  87  SER A  91 -1  N  LEU A  89   O  TYR A  97
SHEET    4   A 5 TRP A  71  ASP A  75 -1  N  ALA A  74   O  ILE A  88
SHEET    5   A 5 GLY A  52  PHE A  56 -1  N  PHE A  56   O  TRP A  71
SHEET    1   B 5 THR B 107  ARG B 110  0
SHEET    2   B 5 ILE B  96  THR B  99 -1  N  LYS B  98   O  THR B 107
SHEET    3   B 5 ARG B  87  SER B  91 -1  N  LEU B  89   O  TYR B  97
SHEET    4   B 5 TRP B  71  ASP B  75 -1  N  ALA B  74   O  ILE B  88
SHEET    5   B 5 GLY B  52  ILE B  55 -1  N  ASP B  54   O  GLU B  73

The sheet presented as BS1 below is an eight-stranded beta-barrel. This is represented by a nine-stranded sheet in which the first and last strands are identical.

SHEET    1 BS1 9  VAL   13  ILE    17  0                               
SHEET    2 BS1 9  ALA   70  ILE    73  1  O  TRP    72   N  ILE    17  
SHEET    3 BS1 9  LYS  127  PHE   132  1  O  ILE   129   N  ILE    73  
SHEET    4 BS1 9  GLY  221  ASP   225  1  O  GLY   221   N  ILE   130  
SHEET    5 BS1 9  VAL  248  GLU   253  1  O  PHE   249   N  ILE   222  
SHEET    6 BS1 9  LEU  276  ASP   278  1  N  LEU   277   O  GLY   252  
SHEET    7 BS1 9  TYR  310  THR   318  1  O  VAL   317   N  ASP   278  
SHEET    8 BS1 9  VAL  351  TYR   356  1  O  VAL   351   N  THR   318  
SHEET    9 BS1 9  VAL   13  ILE    17  1  N  VAL    14   O  PRO   352  

The sheet structure of this example is bifurcated. In order to represent this feature, two sheets are defined. Strands 2 and 3 of BS7 and BS8 are identical.

SHEET    1 BS7 3  HIS  662  THR   665  0                               
SHEET    2 BS7 3  LYS  639  LYS   648 -1  N  PHE   643   O  HIS   662  
SHEET    3 BS7 3  ASN  596  VAL   600 -1  N  TYR   598   O  ILE   646  
SHEET    1 BS8 3  ASN  653  TRP   656  0                               
SHEET    2 BS8 3  LYS  639  LYS   648 -1  N  LYS   647   O  THR   655  
SHEET    3 BS8 3  ASN  596  VAL   600 -1  N  TYR   598   O  ILE   646  
'''

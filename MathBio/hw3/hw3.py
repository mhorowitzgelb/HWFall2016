import sys
import math
import numpy

#
# Overall size of figure, adjust for your screen.
# Uncomment your choice and comment out your non-choice:
#

bigscreen = True  # For Imacs and comparable sized displays
bigscreen = False  # For most Laptops

################################################
# Default Parameters for Van der Waals Potential
#        See Equation 3.11:
################################################
sigmaVdW = 3  # note: this uses the same parameters for all atom types
welldepthVdW = 0.1  # note: this uses the same parameters for all atom types
ConVdW = welldepthVdW
topval_VderWaals = 100000.  # note: this sets max VdW so it's not huge


##########################################


def help_page():
	help = '''
   You must use 9 command line arguments:
       1   built-in model number for Fixed molecule
       2   alpha value
       3   epsilon value (add to 1/r to give 1 / (r+epsilon) ) use small nonnegative value
       4   built-in model number for Moved molecule
       5  x start position for gradient descent optimization
       6  y start position for gradient descent optimization
       7  step size for gradient ascent
     [7]  sigmaVdW, where 6-12 pot is 0.0 (default = 3.0)
     [8]  max well depth for VdW (default = 0.1)
     [9]  set all VdW greater than this value to this value (default = 100000.)

   The last three are optional Van der Waals potential pars, need to give all three or none.

   '''
	print help

if __name__=='__main__':
	if len(sys.argv) < 8:
		print
		print "Error: not enough command line arguments used."
		help_page()
		sys.exit(2)

	modelA = int(sys.argv[1])
	efunc = 1
	alpha = float(sys.argv[2])
	eps = float(sys.argv[3])
	modelB = int(sys.argv[4])
	x_start = float(sys.argv[5])
	y_start = float(sys.argv[6])
	grad_step = float(sys.argv[7])

	if len(sys.argv) == 11:
		sigmaVdW = float(sys.argv[8])
		welldepthVdW = float(sys.argv[9])
		topval_VderWaals = float(sys.argv[10])
		ConVdW = welldepthVdW

	print 'using fixed molecule model:', modelA
	print 'using moved molecule model:', modelB
	print 'using energy function:', efunc
	print

	model_descr = {1: " 1CRN-THR-1/2 ",
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

	maintitle = "A: " + sys.argv[1] + model_descr[modelA] + "; "
	maintitle += "B: " + sys.argv[4] + model_descr[modelB] + "; "
	maintitle += "alpha: " + sys.argv[2] + "; "
	maintitle += "eps: " + sys.argv[3] + "; "
	maintitle += "sigVdW: " + str(sigmaVdW) + "; "
	maintitle += "wellVdW: " + str(welldepthVdW) + "; "
	maintitle += "topVdW: " + str(topval_VderWaals)


def get_molecule_model(model):
	if model == 1:
		# THREONINE 1 and 2 from PDB 1crn
		# ( x, y, z, charge, radius )
		atoms = [
			(17.047, 14.099, 3.625, 0.1812, 1.824),
			(16.967, 12.784, 4.338, 0.0034, 1.908),
			(15.685, 12.755, 5.133, 0.6163, 1.908),
			(15.268, 13.825, 5.594, -0.5722, 1.6612),
			(18.17, 12.703, 5.337, 0.4514, 1.908),
			(19.334, 12.829, 4.463, -0.6764, 1.721),
			(18.15, 11.546, 6.304, -0.2554, 1.908),
			(16.944, 12.070, 3.646, 0.1087, 1.1),
			(18.142, 13.494, 5.933, -0.0323, 1.387),
			(19.089, 12.928, 3.520, 0.407, 0),
			(17.308, 11.590, 6.859, 0.0627, 1.487),
			(18.163, 10.677, 5.806, 0.0627, 1.487),
			(18.957, 11.601, 6.911, 0.0627, 1.487),
			(16.122, 14.421, 3.452, 0.1934, 0.6),
			(17.532, 13.966, 2.767, 0.1934, 0.6),
			(17.536, 14.744, 4.203, 0.1934, 0.6),
			(15.115, 11.555, 5.265, -0.4157, 1.824),
			(13.856, 11.469, 6.066, -0.0389, 1.908),
			(14.164, 10.785, 7.379, 0.5973, 1.908),
			(14.993, 9.862, 7.443, -0.5679, 1.6612),
			(12.732, 10.711, 5.261, 0.3654, 1.908),
			(13.308, 9.439, 4.926, -0.6761, 1.721),
			(12.484, 11.442, 3.895, -0.2438, 1.908),
			(15.511, 10.776, 4.852, 0.2719, 0.6),
			(13.548, 12.399, 6.243, 0.1007, 1.387),
			(11.957, 10.566, 5.874, 0.0043, 1.387),
			(13.421, 9.330, 3.945, 0.4102, 0),
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
			(-2.397, -3.805, 0.719, -0.5163, 1.8240),  # N
			(-1.905, -3.393, -0.621, 0.0397, 1.9080),  # CA
			(-1.462, -4.562, -1.444, 0.5366, 1.9080),  # C
			(-1.705, -4.570, -2.679, -0.5819, 1.6612),  # O
			(-0.704, -2.421, -0.466, 0.0560, 1.9080),  # CB
			(-1.149, -1.022, 0.033, 0.0136, 1.9080),  # CG
			(0.000, 0.000, 0.000, 0.8054, 1.9080),  # CD
			(1.125, -0.395, -0.000, -0.8188, 1.6612),  # OE1
			(-0.384, 1.196, 0.000, -0.8188, 1.6612),  # OE2
			(-2.032, -3.462, 1.539, 0.2936, 0.6000),  # H
			(-2.646, -2.889, -1.092, 0.1105, 1.3870),  # HA
			(-0.108, -2.791, 0.204, -0.0173, 1.4870),  # HB2
			(-0.304, -2.303, -1.342, -0.0173, 1.4870),  # HB3
			(-1.848, -0.699, -0.565, -0.0425, 1.4870),  # HG2
			(-1.438, -1.110, 0.960, -0.0425, 1.4870)  # HG3
		]
		imodx = 0
		imody = 1

	elif model == 18:
		# GLU 23 from PDB 1crn
		atoms = [
			(-0.704, -2.421, -0.466, 0.0560, 1.9080),  # CB
			(-1.149, -1.022, 0.033, 0.0136, 1.9080),  # CG
			(0.000, 0.000, 0.000, 0.8054, 1.9080),  # CD
			(1.125, -0.395, -0.000, -0.8188, 1.6612),  # OE1
			(-0.384, 1.196, 0.000, -0.8188, 1.6612),  # OE2
			(-0.108, -2.791, 0.204, -0.0173, 1.4870),  # HB152
			(-0.304, -2.303, -1.342, -0.0173, 1.4870),  # HB3
			(-1.848, -0.699, -0.565, -0.0425, 1.4870),  # HG2
			(-1.438, -1.110, 0.960, -0.0425, 1.4870)  # HG3
		]
		imodx = 0
		imody = 1

	elif model == 9:
		# ARG 10 from PDB 1crn
		atoms = [
			(1.638, 4.502, -3.554, -0.3479, 1.8240),  # N
			(0.753, 4.211, -2.410, -0.2637, 1.9080),  # CA
			(-0.150, 5.428, -2.074, 0.7341, 1.9080),  # C
			(-1.338, 5.246, -1.773, -0.5894, 1.6612),  # O
			(1.599, 3.836, -1.171, -0.0007, 1.9080),  # CB
			(0.758, 3.599, 0.078, 0.0390, 1.9080),  # CG
			(-0.196, 2.448, 0.037, 0.0486, 1.9080),  # CD
			(0.613, 1.196, -0.010, -0.5295, 1.8240),  # NE
			(0.000, -0.000, 0.000, 0.8076, 1.9080),  # CZ
			(-1.338, -0.118, 0.000, -0.8627, 1.8240),  # NH1
			(0.790, -1.082, -0.000, -0.8627, 1.8240),  # NH2
			(2.621, 4.514, -3.460, 0.2747, 0.6000),  # H
			(0.160, 3.437, -2.648, 0.1560, 1.3870),  # HA
			(2.091, 3.002, -1.378, 0.0327, 1.4870),  # HB2
			(2.222, 4.582, -0.990, 0.0327, 1.4870),  # HB3
			(1.387, 3.467, 0.871, 0.0285, 1.4870),  # HG2
			(0.227, 4.449, 0.276, 0.0285, 1.4870),  # HG3
			(-0.733, 2.438, 0.854, 0.0687, 1.3870),  # HD2
			(-0.735, 2.499, -0.777, 0.0687, 1.3870),  # HD3
			(1.580, 1.256, -0.049, 0.3456, 0.6000),  # HE
			(-1.894, 0.712, -0.028, 0.4478, 0.6000),  # HH11
			(-1.758, -1.015, 0.026, 0.4478, 0.6000),  # HH12
			(1.777, -1.007, -0.027, 0.4478, 0.6000),  # HH21
			(0.363, -1.996, 0.027, 0.4478, 0.6000)  # HH22
		]
		imodx = 0
		imody = 1

	elif model == 19:
		# ARG 10 from PDB 1crn
		atoms = [
			(1.599, 3.836, -1.171, -0.0007, 1.9080),  # CB
			(0.758, 3.599, 0.078, 0.0390, 1.9080),  # CG
			(-0.196, 2.448, 0.037, 0.0486, 1.9080),  # CD
			(0.613, 1.196, -0.010, -0.5295, 1.8240),  # NE
			(0.000, -0.000, 0.000, 0.8076, 1.9080),  # CZ
			(-1.338, -0.118, 0.000, -0.8627, 1.8240),  # NH1
			(0.790, -1.082, -0.000, -0.8627, 1.8240),  # NH2
			(2.091, 3.002, -1.378, 0.0327, 1.4870),  # HB2
			(2.222, 4.582, -0.990, 0.0327, 1.4870),  # HB3
			(1.387, 3.467, 0.871, 0.0285, 1.4870),  # HG2
			(0.227, 4.449, 0.276, 0.0285, 1.4870),  # HG3
			(-0.733, 2.438, 0.854, 0.0687, 1.3870),  # HD2
			(-0.735, 2.499, -0.777, 0.0687, 1.3870),  # HD3
			(1.580, 1.256, -0.049, 0.3456, 0.6000),  # HE
			(-1.894, 0.712, -0.028, 0.4478, 0.6000),  # HH11
			(-1.758, -1.015, 0.026, 0.4478, 0.6000),  # HH12
			(1.777, -1.007, -0.027, 0.4478, 0.6000),  # HH21
			(0.363, -1.996, 0.027, 0.4478, 0.6000)  # HH22
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


def gridMol(mol, lb, step):
	'''
   Create new set of coordintes so that all atom positions are exactly on
   grid points.
   '''
	(atoms, imodx, imody) = mol

	gatoms = []
	for i in range(0, len(atoms)):
		# round to nearby grid point
		xf = int(round((atoms[i][imodx] - lb) / step))
		yf = int(round((atoms[i][imody] - lb) / step))
		xfg = lb + xf * step
		yfg = lb + yf * step
		(x, y, z, q, r) = atoms[i]
		gatoms.append((xfg, yfg, z, q, r))
		# make sure to use imodx = 0 and imody = 1
	return gatoms


def EnergyPair(alpha, eps, efunc, q1, q2, x1, x2, rAB):
	'''
   Evaluate pair energy of two atoms located at x1 and x2.
   x1 and x2 must be Numpy arrays.
   '''
	if efunc < 1 or efunc > 5:
		print "Error in EnergyPair: no energy function for efunc =", efunc
		sys.exit(3)

	x12 = abs(x2 - x1 + rAB)
	d2 = x12 ** 2
	r12 = numpy.sqrt(numpy.sum(d2))
	E12 = 0.0

	if efunc == 1 or efunc == 4:
		E12 = 332.0636 * alpha * q1 * q2 / (r12 + eps)

	elif efunc == 2 or efunc == 5:
		E12 = 332.0636 * q1 * q2 * numpy.exp(-alpha * r12) / (r12 + eps)

	if efunc == 3 or efunc == 4 or efunc == 5:  # these functions include both electrostatic and vdW terms
		# Equation 3.11  (with eps added in)
		r6term = (sigmaVdW / (r12 + eps)) ** 6
		Evdw = ConVdW * r6term * (r6term - 1.)
		if Evdw > topval_VderWaals:
			Evdw = topval_VderWaals
		E12 += Evdw

	return E12


def GradEnergyPair(alpha, eps, efunc, q1, q2, x1, x2, rAB):
	if efunc < 1 or efunc > 5:
		print "Error in EnergyPairGradient: no energy function for efunc = ", efunc
		sys.exit(3)


	x12 = abs(x2 - x1 + rAB)
	d2 = x12 ** 2
	r12 = math.sqrt(numpy.sum(d2))

	x11 = x1[0]
	x12 = x1[1]

	x21 = x2[0]
	x22 = x2[1]

	r1 = rAB[0]
	r2 = rAB[1]
	eleGrad = numpy.array([
		-332.0636 * alpha * q1 * q2 * (x21 - x11 + r1) / (r12 * (r12 + eps)**2),
		-332.0636 * alpha * q1 * q2 * (x22 - x12 + r2) / (r12 * (r12 + eps)**2)
	])

	r6term = (sigmaVdW / (r12 + eps)) ** 6
	Evdw = ConVdW * r6term * (r6term - 1.)
	if Evdw > topval_VderWaals:
		vdwGrad = numpy.array([0,0])
	else:
		xterm = (x21-x11 + r1)
		yterm = (x22 - x12 + r2)
		vdwGrad = numpy.array([
			ConVdW * (6*sigmaVdW**6*(xterm/(r12*(r12+eps)**7)) - 12*sigmaVdW**12*(xterm/(r12*(r12+eps)**13))),
			ConVdW * (6*sigmaVdW**6*(yterm/(r12*(r12+eps)**7)) - 12*sigmaVdW**12*(yterm/(r12*(r12+eps)**13)))
		])

	return eleGrad + vdwGrad



def EnergyInterMol(molA, molB, efunc, alpha, eps, rAB):
	'''
   Total interaction energy of two molecules, A and B.
   rAB is a Numpy array giving displacement of B initial coordinates from A.
   '''
	(atomsA, imodxA, imodyA) = molA
	(atomsB, imodxB, imodyB) = molB

	E_AB = 0.0
	for i in range(0, len(atomsA)):
		xi = numpy.array([atomsA[i][imodxA], atomsA[i][imodyA]])
		qi = atomsA[i][3]
		for j in range(0, len(atomsB)):
			xj = numpy.array([atomsB[j][imodxB], atomsB[j][imodyB]])
			qj = atomsB[j][3]
			E_AB += EnergyPair(alpha, eps, efunc, qi, qj, xi, xj, rAB)
	return E_AB

def GradEnergyInterMol(molA, molB, efunc, alpha, eps, rAB):
	(atomsA, imodxA, imodyA) = molA
	(atomsB, imodxB, imodyB) = molB

	GradE_AB = 0.0
	for i in range(0, len(atomsA)):
		xi = numpy.array([atomsA[i][imodxA], atomsA[i][imodyA]])
		qi = atomsA[i][3]
		for j in range(0, len(atomsB)):
			xj = numpy.array([atomsB[j][imodxB], atomsB[j][imodyB]])
			qj = atomsB[j][3]
			GradE_AB += GradEnergyPair(alpha, eps, efunc+3, qi, qj, xi, xj, rAB)
	return GradE_AB

def net_charge_mol(mol):
	(atoms, imodx, imody) = mol
	sigmaQ = 0.
	for i in range(0, len(atoms)):
		sigmaQ += atoms[i][3]
	return sigmaQ

if __name__ == '__main__':
	mA = get_molecule_model(modelA)
	mB = get_molecule_model(modelB)
	QtotA = net_charge_mol(mA)
	QtotB = net_charge_mol(mB)
	print "Total charge molecule A:", QtotA
	print "Total charge molecule B:", QtotB

	'''
	Do Grad Descent to find minimum
	'''

	step_size = grad_step
	precision = grad_step / 10.
	steps = 0

	gradX = numpy.array([x_start, y_start])
	trajX = [gradX[0]]
	trajY = [gradX[1]]
	energyOld = 1
	energyNew = 0
	while abs(energyNew - energyOld) > precision:
		energyOld = energyNew
		steps += 1
		energyNew = EnergyInterMol(mA,mB,efunc+3,alpha,eps,gradX)
		print "Running gradient descent step ", steps, " pos" , gradX ,"energy", energyNew
		gradX -= grad_step * GradEnergyInterMol(mA,mB,efunc,alpha,eps,gradX)



	print 'Min GradE: ', energyNew, 'at ', gradX


	sys.exit()

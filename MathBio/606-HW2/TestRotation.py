__author__ = 'mhorowitzgelb'

import numpy as np
import math as meth

vector = np.array([1.,0.])


cos_45 = meth.cos(meth.pi / 4.)
sin_45 = meth.cos(meth.pi / 4.)

rot_mat = np.array([[cos_45, -sin_45],[sin_45, cos_45]])

for i in range(0,8):
	print vector
	vector = np.transpose(np.dot(rot_mat, np.transpose(vector)))

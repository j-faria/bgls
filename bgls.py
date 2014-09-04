# -*- coding: utf-8 -*-
#================================================================================
# Copyright (c) 2014 João Faria, Annelies Mortier
# Distributed under the MIT License.
# (See accompanying file LICENSE or copy at http://opensource.org/licenses/MIT)
#================================================================================

import numpy as np
import sys

try:
	import mpmath  # https://code.google.com/p/mpmath/
except ImportError, e1:
	try:
		from sympy import mpmath  # http://sympy.org/en/index.html
	except ImportError, e2:
		raise e2
	finally:
		raise e1

pi = np.pi

def bgls(t, y, err, plow=0.5, phigh=100, ofac=1):

	f = np.linspace(1./phigh, 1./plow, int(100*ofac))

	omegas = 2. * pi * f

	err2 = err * err
	w = 1./err2
	W = sum(w)

	bigY = sum(w*y)  # Eq. (10)

	p = []
	constants = []
	exponents = []

	for i, omega in enumerate(omegas):
		theta = 0.5 * np.arctan2(sum(w*np.sin(2.*omega*t)), sum(w*np.cos(2.*omega*t)))
		x = omega*t - theta
		cosx = np.cos(x)
		sinx = np.sin(x)
		wcosx = w*cosx
		wsinx = w*sinx

		C = sum(wcosx)
		S = sum(wsinx)

		YCh = sum(y*wcosx)
		YSh = sum(y*wsinx)
		CCh = sum(wcosx*cosx)
		SSh = sum(wsinx*sinx)

		if (CCh != 0 and SSh != 0):
			K = (C*C*SSh + S*S*CCh - W*CCh*SSh)/(2.*CCh*SSh)

			L = (bigY*CCh*SSh - C*YCh*SSh - S*YSh*CCh)/(CCh*SSh)

			M = (YCh*YCh*SSh + YSh*YSh*CCh)/(2.*CCh*SSh)

			constants.append(1./np.sqrt(CCh*SSh*abs(K)))

		elif (CCh == 0):
			K = (S*S - W*SSh)/(2.*SSh)

			L = (bigY*SSh - S*YSh)/(SSh)

			M = (YSh*YSh)/(2.*SSh)

			constants.append(1./np.sqrt(SSh*abs(K)))

		elif (SSh == 0):
			K = (C*C - W*CCh)/(2.*CCh)

			L = (bigY*CCh - C*YCh)/(CCh)

			M = (YCh*YCh)/(2.*CCh)

			constants.append(1./np.sqrt(CCh*abs(K)))

		if K > 0:
			raise RuntimeError('K is positive. This should not happen.')

		exponents.append(M - L*L/(4.*K))

	constants = np.array(constants)
	exponents = np.array(exponents)

	logp = np.log10(constants) + (exponents * np.log10(np.exp(1.)))

	p = [10**mpmath.mpf(x) for x in logp]

	p = np.array(p) / max(p)  # normalize

	p[p < (sys.float_info.min * 10)] = 0
	p = np.array([float(pp) for pp in p])

	return 1./f, p

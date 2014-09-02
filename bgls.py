# from periodogram import *
# from utils import *
from scipy.special import dawsn

import numpy as np
from pylab import *
import random

import mpmath

pi = np.pi

t, y, err = np.loadtxt('/home/joao/Work/OPEN/14her.rv', skiprows=2, unpack=True)
# t, y, err = np.loadtxt('/home/joao/phd/data/extra/24sex-keck.rv', skiprows=0, unpack=True)
# t, y, err = np.loadtxt('/home/joao/phd/data/HD41248_harps_mean.rdb', skiprows=2, unpack=True, usecols=(0,1,2))




def bgls(t, y, err, plow=0.5, phigh=100, ofac=1):


	# f = np.array([0.0005555555555555556])
	# f = np.linspace(1./4000, 1./100., 1000)
	f = np.linspace(1./phigh, 1./plow, int(100*ofac))
	# f = per1.freq

	g1 = min(y)
	g2 = max(y)


	omegas = 2.* pi * f

	err2 = err * err
	W = sum(1./err2)
	w = 1./err2

	yh = sum(w*y)        # Eq. (7)
	# yh = self.y - self._Y          # Subtract weighted mean
	_YY = sum(w * y**2)

	# Unnormalized power
	_upow = zeros(len(omegas))
	_a = zeros(len(omegas))
	_b = zeros(len(omegas))
	_c = zeros(len(omegas))

	KK = []
	p = []

	list_constants = []
	list_exponents = []

	# def bgls_periodogram(omega):
	# 	print omega
	# 	theta = 0.5 * np.arctan2(sum(w*np.sin(2.*omega*t)), sum(w*np.cos(2.*omega*t)))
	# 	x = omega*t - theta
	# 	cosx = cos(x)
	# 	sinx = sin(x)
	# 	wcosx = w*cosx         # attach weights
	# 	wsinx = w*sinx         # attach weights

	# 	C = sum(wcosx)         # Eq. (8)
	# 	S = sum(wsinx)         # Eq. (9)

	# 	YC = sum(y*wcosx)     # Eq. (11)
	# 	YS = sum(y*wsinx)     # Eq. (12)
	# 	CCh = sum(wcosx*cosx)  # Eq. (13)
	# 	# CSh = sum(wcosx*sinx)  # Eq. (15)
	# 	SSh = sum(wsinx*sinx)
	# 	# SSh = 1.-CCh

	# 	# CC = CCh-C*C           # Eq. (13)
	# 	# SS = SSh-S*S           # Eq. (14)
	# 	# CS = CSh-C*S           # Eq. (15)
	# 	# D = CC*SS-CS*CS        # Eq. (6)

	# 	K = (C*C*SSh + S*S*CCh - W*CCh*SSh)/(2.*CCh*SSh)
	# #	KK.append(K)
	# #	if K>0:
	# #	print 'error!!!!!!'
	# 	# else:
	# 	# 	print K
	# 	# K = abs(K)

	# 	L = (yh*CCh*SSh - C*YC*SSh - S*YS*CCh)/(CCh*SSh)

	# 	M = (YC*YC*SSh + YS*YS*CCh)/(2.*CCh*SSh)
	# 	# M = 0.

	# #	list_constants.append(1./np.sqrt(CCh*SSh*abs(K)))
	# #	list_exponents.append(M - L*L/(4.*K))
	# 	prob = 1./np.sqrt(CCh*SSh*abs(K)) * mpmath.exp(M - L*L/(4.*K)) 
	# 	       # (np.exp(g2*g2*K + g2*L) * dawsn((2.*g2*K + L)/(2.*np.sqrt(K))) - np.exp(g1*g1*K + g1*L) * dawsn((2.*g1*K + L)/(2.*np.sqrt(K))))
	# 	# p.append(prob)
	# 	return prob


	for i, omega in enumerate(omegas):
		theta = 0.5 * np.arctan2(sum(w*np.sin(2.*omega*t)), sum(w*np.cos(2.*omega*t)))
		x = omega*t - theta
		cosx = cos(x)
		sinx = sin(x)
		wcosx = w*cosx         # attach weights
		wsinx = w*sinx         # attach weights

		C = sum(wcosx)         # Eq. (8)
		S = sum(wsinx)         # Eq. (9)

		YC = sum(y*wcosx)     # Eq. (11)
		YS = sum(y*wsinx)     # Eq. (12)
		CCh = sum(wcosx*cosx)  # Eq. (13)
		# CSh = sum(wcosx*sinx)  # Eq. (15)
		SSh = sum(wsinx*sinx)
		# SSh = 1.-CCh

		# CC = CCh-C*C           # Eq. (13)
		# SS = SSh-S*S           # Eq. (14)
		# CS = CSh-C*S           # Eq. (15)
		# D = CC*SS-CS*CS        # Eq. (6)

		K = (C*C*SSh + S*S*CCh - W*CCh*SSh)/(2.*CCh*SSh)
		KK.append(K)
		if K>0:
			print 'error!!!!!!'
		# else:
		# 	print K
		# K = abs(K)

		L = (yh*CCh*SSh - C*YC*SSh - S*YS*CCh)/(CCh*SSh)

		M = (YC*YC*SSh + YS*YS*CCh)/(2.*CCh*SSh)
		# M = 0.

		list_constants.append(1./np.sqrt(CCh*SSh*abs(K)))
		list_exponents.append(M - L*L/(4.*K))
		# prob = 1./np.sqrt(CCh*SSh*abs(K)) * mpmath.exp(M - L*L/(4.*K)) 
		       # (np.exp(g2*g2*K + g2*L) * dawsn((2.*g2*K + L)/(2.*np.sqrt(K))) - np.exp(g1*g1*K + g1*L) * dawsn((2.*g1*K + L)/(2.*np.sqrt(K))))
		# p.append(prob)



	# self._a[i] = (YC*SS-YS*CS) / D
	# self._b[i] = (YS*CC-YC*CS) / D
	# self._off[i] = -self._a[i]*C - self._b[i]*S

	# self._upow[i] = (SS*YC*YC + CC*YS*YS - 2.*CS*YC*YS) / (self._YY*D) # Eq. (5) in ZK09

	list_constants = np.array(list_constants)
	list_exponents = np.array(list_exponents)

	mean_exp = np.mean(list_exponents)

	logp = np.log10(list_constants) + (list_exponents * np.log10(np.exp(1.)))

	p = [10**mpmath.mpf(x) for x in logp]
	# normalize
	p = np.array(p) / max(p)

	# min_number = sys.float_info.min * 1
	# print 'minimum number: ', min_number

	p[p < (sys.float_info.min * 10)] = 0
	p = np.array([float(pp) for pp in p])


	return 1./f, p

	# figure()
	# semilogx(1./f, p, 'k-'); show()

	# # KK = np.array(KK)

	# # # p = np.array(p)
	# # p = np.array(p) 
	# # p /= mpmath.mean(p)
	# # pp = np.array([float(x) for x in p])
	# figure()
	# semilogx(1./f, logp, 'k-')
	# # figure()
	# # semilogx(1./per1.freq, per1.power, 'b-')
	# show()

#### bgls
# period_bgls, prob_bgls = bgls(t, y, err, plow=10, phigh=100)
# figure()
# semilogx(period_bgls, prob_bgls, 'b-')
# show()
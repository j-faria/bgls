from periodogram import *
from utils import *
from bgls import bgls

import numpy as np
from pylab import *
import matplotlib as mpl
import random
from copy import copy
import sys


pi = np.pi

# set tick width
mpl.rcParams['xtick.major.size'] = 6
mpl.rcParams['ytick.major.size'] = 6
# mpl.rcParams['xtick.major.width'] = 4
mpl.rcParams['xtick.minor.size'] = 4
mpl.rcParams['ytick.minor.size'] = 4
# mpl.rcParams['xtick.minor.width'] = 2


which_example = int(raw_input('Which example: '))
which_plots = raw_input('Which plots: ')

if which_example == 1:
	#############
	## example 1
	#############


	## Make fake data

	realP = 105.
	f = [1./realP]
	K = 1.0

	t = sorted([random.random()*180. + 10. for i in range(50)])
	y = [K * sin(2.*pi/realP*i) + random.random()*1.4 - 0.7 for i in t]
	err = [0.4 + random.random()*0.1 - 0.05 for i in t]

	fakeP = 5.*realP/11.
	ft = sorted([random.random()*180. + 10. for i in range(50)])
	fy = [4.*K/3. * sin(2.*pi/fakeP*i - 3.0) + random.random()*1.8 - 0.9 for i in ft]
	ferr = [1.1 + random.random()*0.3 - 0.15 for i in ft]


	tlong = [i/2+10. for i in range(360)]
	ylong = [K * sin(2.*pi/realP*i) for i in tlong]

	t2 = [i[0] for i in sorted(zip(t+ft,y+fy,err+ferr))]
	y2 = [i[1] for i in sorted(zip(t+ft,y+fy,err+ferr))]
	err2 = [i[2] for i in sorted(zip(t+ft,y+fy,err+ferr))]


	# freq,power,freqmax,powermax,a_cos,b_sin,c_cte,phi,fNy,timespan = GLS.periodogram(np.array(t2),np.array(y2),np.array(err2),4,2.0)


	t2 = np.asarray(t2, order='F')
	y2 = np.asarray(y2, order='F')
	err2 = np.asarray(err2, order='F')


	tt = np.linspace(0, 200, 500)
	yy = K * sin(2.*pi/realP*tt)

	if '1' in which_plots:
		figure('obs')
		subplot(111)
		errorbar(t2, y2, yerr=err2, fmt='ko', ms=7)
		plot(tt, yy, 'k-', lw=1.5, alpha=0.9)
		xlim([0, 200])
		ylim([-3, 3])
		xlabel('Time (days)')
		ylabel('RV')
		# subplot(212)
		show()

	#sys.exit(0)

	# c = raw_input('Continue (y/n) ')
	# if c == 'n': sys.exit(0)

	TimeSeries = BasicTimeSeries()
	TimeSeries.time = t2
	TimeSeries.val = y2
	TimeSeries.error = err2
	per = LombScargle(TimeSeries, ofac=10)
	per1 = gls(TimeSeries, hifac=500, plow=1.)
	perb = bls(TimeSeries, hifac=500, plow=1.)

	#### bgls
	period_bgls, prob_bgls = bgls(t2, y2, err2, plow=1, ofac=100, phigh=1000)

	# for normalization
	a1 = perb.power.max() / 100.
	# a2 = per.power.max()

	if '2' in which_plots:
		figure('per')
		ax1 = subplot(2, 1, 1)
		ax1.semilogx(1./perb.freq, perb.power/a1, 'r--', label=r'Bayesian LS')
		ax1.semilogx(period_bgls, prob_bgls * 100., 'k', label=r'Bayesian GLS')

		for true in f:
			ax1.axvline(x=1./true, color='g', lw=4, alpha=0.4, label='true')

		ax1.set_ylim([0, 110])
		ax1.set_xlim([5, 1000])
		ax1.legend(frameon=False, loc='upper right')
		ax1.set_ylabel('Probability (%)')

		ax2 = subplot(2, 1, 2, sharex=ax1)
		ax2.semilogx(1./per.freq, per.power, 'r--', label='LS')
		ax2.semilogx(1./per1.freq, per1.power, 'k', label='GLS')
		# ax2.semilogx(1./freq, power, 'k', label='GLS')

		for true in f:
			ax2.axvline(x=1./true, color='g', lw=4, alpha=0.4, label='true')
		ax2.legend(frameon=False, loc='upper right')
		ax2.set_xlim([5, 1000])
		# text(0.03, 0.85, '%s)' % ascii_lowercase[j], transform=ax.transAxes)
		ax2.set_ylabel('Power')
		ax2.set_xlabel('Period (days)')
		show()

	##################################################################################
	# figure with 3 periodogram panels
	##################################################################################

	if '3' in which_plots:
		figure('per_log', figsize=(8, 11))
		ax1 = subplot(3, 1, 1)  # <====
		ax1.semilogx(1./per.freq, per.power, 'r--', lw=1.5, label='LS')
		ax1.semilogx(1./per1.freq, per1.power, 'k', lw=1.5, label='GLS')
		# ax1.semilogx(1./freq, power, 'k', label='GLS')

		for true in f:
			ax1.axvline(x=1./true, color='g', lw=4, alpha=0.4) #, label='true')
		ax1.legend(frameon=False, loc='upper right')
		ax1.set_xlim([5, 1000])
		ax1.set_ylabel('Power')


		ax2 = subplot(3, 1, 2, sharex=ax1)  # <====
		ax2.semilogx(1./perb.freq, perb.power/a1, 'r--', lw=1.5, label=r'Bayesian LS')
		ax2.semilogx(period_bgls, prob_bgls*100, 'k', lw=1.5, label=r'Bayesian GLS')

		# ax2.grid()

		for true in f:
			ax2.axvline(x=1./true, color='g', lw=4, alpha=0.4) #, label='true')

		ax2.set_ylim([0, 110])
		ax2.set_xlim([5, 1000])
		ax2.legend(frameon=False, loc='upper right')
		ax2.set_ylabel('Probability (%)')


		ax3 = subplot(3, 1, 3, sharex=ax1)  # <====
		ax3.loglog(1./perb.freq, perb.power/a1, 'r--', lw=1.5, label=r'Bayesian LS')
		ax3.loglog(period_bgls, prob_bgls*100, 'k', lw=1.5, label=r'Bayesian GLS')

		# ax3.grid()

		for true in f:
			ax3.axvline(x=1./true, color='g', lw=4, alpha=0.4) #, label='true')

		ax3.set_ylim([0, 10000])
		ax3.set_xlim([5, 1000])
		loc = [1e0, 1e-10, 1e-20]
		if (min(prob_bgls*100) < 1e-29): loc.append(1e-30)
		# loc_labels = ['10$^{%i}$' % (int(math.log10(i))) for i in loc]
		ax3.set_yticks(loc)
		ax3.minorticks_on()
		ax3.tick_params(axis='y',which='major')
		# ax3.set_yticklabels(loc_labels)
		# ax3.yaxis.major.formatter.set_locs([1e-3])

		leg = ax3.legend(frameon=False, loc='center right', shadow=True)
		frame = leg.get_frame()
		frame.set_color('white')
		ax3.set_ylabel('Probability (%)')
		#ax2.set_xlabel('Period (days)')



		# text(0.03, 0.85, '%s)' % ascii_lowercase[j], transform=ax.transAxes)
		ax3.set_xlabel('Period (days)')
		show()

elif which_example == 2:

	#############
	## example 2
	#############

	realP = 50.
	f = [1./realP]
	K = 1.0

	t = sorted([random.random()*180. + 10. for i in range(100)])
	# y = [K * sin(2.*pi/realP*i) + random.random()*1.4 - 0.7 for i in t]
	y = [K * sin(2.*pi/realP*i) for i in t]
	err = [0.2 for i in t]


	t2 = np.asarray(t, order='F')
	y2 = np.asarray(y, order='F')
	err2 = np.asarray(err, order='F')

	cond = y2 > -0.3
	t2 = t2[cond]
	y2 = y2[cond]
	err2 = err2[cond]

	tt = np.linspace(0, 200, 500)
	yy = K * sin(2.*pi/realP*tt)

	if '1' in which_plots:
		figure('obs')
		subplot(111)
		errorbar(t2, y2, yerr=err2, fmt='ko', ms=7)
		plot(tt, yy, 'k-', lw=1.5, alpha=0.9)
		xlim([0, 200])
		ylim([-1.5, 1.5])
		xlabel('Time (days)')
		ylabel('RV')
		# subplot(212)
		show()

	print np.mean(y2)
	# sys.exit(0)

	# c = raw_input('Continue (y/n) ')
	# if c == 'n': sys.exit(0)

	TimeSeries = BasicTimeSeries()
	TimeSeries.time = t2
	TimeSeries.val = y2
	TimeSeries.error = err2
	per = LombScargle(TimeSeries, ofac=10)
	per1 = gls(TimeSeries, hifac=500, plow=1.)
	perb = bls(TimeSeries, hifac=500, plow=1.)

	#### bgls
	period_bgls, prob_bgls = bgls(t2, y2, err2, plow=1, ofac=100, phigh=1000)

	# for normalization
	a1 = perb.power.max() / 100.
	# a2 = per.power.max()

	if '2' in which_plots:
		figure('per')
		ax1 = subplot(2, 1, 1)
		ax1.semilogx(1./perb.freq, perb.power/a1, 'r-.', lw=1.5, label=r'Bayesian LS')
		ax1.semilogx(period_bgls, prob_bgls*100, 'k', lw=1.5, label=r'Bayesian GLS')

		for true in f:
			ax1.axvline(x=1./true, color='g', lw=4, alpha=0.4, label='true')

		ax1.set_ylim([0, 110])
		ax1.set_xlim([10, 150])
		ax1.legend(frameon=False, loc='upper right')
		ax1.set_ylabel('Probability (%)')

		ax2 = subplot(2, 1, 2, sharex=ax1)
		ax2.semilogx(1./per.freq, per.power, 'r--', lw=1.5, label='LS')
		ax2.semilogx(1./per1.freq, per1.power, 'k', lw=1.5, label='GLS')
		# ax2.semilogx(1./freq, power, 'k', label='GLS')

		for true in f:
			ax2.axvline(x=1./true, color='g', lw=4, alpha=0.4, label='true')
		ax2.legend(frameon=False, loc='upper right')
		ax2.set_xlim([10, 150])
		# text(0.03, 0.85, '%s)' % ascii_lowercase[j], transform=ax.transAxes)
		ax2.set_ylabel('Power')
		ax2.set_xlabel('Period (days)')
		show()


	##################################################################################
	# figure with 3 periodogram panels
	##################################################################################

	if '3' in which_plots:
		figure('per_log', figsize=(8, 11))
		ax1 = subplot(3, 1, 1)  # <====
		ax1.semilogx(1./per.freq, per.power, 'r--', lw=1.5, label='LS')
		ax1.semilogx(1./per1.freq, per1.power, 'k', lw=1.5, label='GLS')
		# ax1.semilogx(1./freq, power, 'k', label='GLS')

		for true in f:
			ax1.axvline(x=1./true, color='g', lw=4, alpha=0.4) #, label='true')
		ax1.legend(frameon=False, loc='upper right')
		ax1.set_xlim([10, 150])
		ax1.set_ylabel('Power')


		ax2 = subplot(3, 1, 2, sharex=ax1)  # <====
		ax2.semilogx(1./perb.freq, perb.power/a1, 'r--', lw=1.5, label=r'Bayesian LS')
		ax2.semilogx(period_bgls, prob_bgls*100, 'k', lw=1.5, label=r'Bayesian GLS')

		# ax2.grid()

		for true in f:
			ax2.axvline(x=1./true, color='g', lw=4, alpha=0.4) #, label='true')

		ax2.set_ylim([0, 110])
		ax2.set_xlim([10, 150])
		ax2.legend(frameon=False, loc='upper right')
		ax2.set_ylabel('Probability (%)')


		ax3 = subplot(3, 1, 3, sharex=ax1)  # <====
		ax3.loglog(1./perb.freq, perb.power/a1, 'r--', lw=1.5, label=r'Bayesian LS')
		ax3.loglog(period_bgls, prob_bgls*100, 'k', lw=1.5, label=r'Bayesian GLS')

		# ax3.grid()

		for true in f:
			ax3.axvline(x=1./true, color='g', lw=4, alpha=0.4) #, label='true')

		ax3.set_ylim([0, 10000])
		ax3.set_xlim([10, 150])
		loc = [1e0, 1e-10, 1e-20, 1e-30, 1e-40]
		if (min(prob_bgls*100) < 1e-50): loc.append(1e-50)
		# loc_labels = ['10$^{%i}$' % (int(math.log10(i))) for i in loc]
		ax3.set_yticks(loc)
		# ax3.minorticks_off()
		ax3.tick_params(axis='y',which='major')
		# ax3.set_yticklabels(loc_labels)
		# ax3.yaxis.major.formatter.set_locs([1e-3])

		ax3.legend(frameon=False, loc='upper right')
		ax3.set_ylabel('Probability (%)')
		#ax2.set_xlabel('Period (days)')



		# text(0.03, 0.85, '%s)' % ascii_lowercase[j], transform=ax.transAxes)
		ax3.set_xlabel('Period (days)')
		show()
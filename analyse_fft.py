import matplotlib.pyplot as plt
#import plotly.plotly as py
import numpy as np
# Learn about API authentication here: https://plot.ly/python/getting-started
# Find your api_key here: https://plot.ly/settings/api

# comparaision 1
#res = np.loadtxt('chainage.dat')
#res2 = np.loadtxt('couplage.dat')

# comparaision 2
#res = np.loadtxt('coupling_mu_0001_interne.dat')
#res2 = np.loadtxt('coupling_mu_0001_interne.dat')
#res2 = np.loadtxt('coupling_mu_0001_Chen.dat')
#res = np.loadtxt('coupling_mu_0001_Chen_coase.dat')
#res = np.loadtxt('coupling_mu_0001_geo_new.dat')
#res2 = np.loadtxt('coupling_mu_001_Chen.dat')
res = np.loadtxt('Fy.dat')
#res2 = np.loadtxt('coupling2.dat')

# comparaision 3
#res = np.loadtxt('chainage_u5_E6e4.dat')
#res2 = np.loadtxt('couplage_u5_E6e4.dat')



time = res[10000:, 0]
dep = res[10000:, 1]

# comparaision 1 & 3
#dep2 = res2[:, 1]


#for i in range(1, 1000-1):
#	if (res2[i, 1] > res2[i-1, 1]) & ((res2[i, 1] > res2[i+1, 1]) ):
#		print('res :', res2[i, 0], res2[i, 1])

# comparaision 2
#dep2 = res2[:, 5]

# pas de temps
dt = 1.e-4

#print(res)
#print(res[:, 0])

n = len(dep) # length of the signal
k = np.arange(n)
Ts = dt; # sampling interval 
Fs = 1.0/Ts; # sampling rate
T = n/Fs
frq = k/T # two sides frequency range
frq = frq[range(n/2)] # one side frequency range

Y = np.fft.fft(dep)/n # fft computing and normalization
Y = Y[range(n/2)]

#Y2 = np.fft.fft(dep2)/n
#Y2 = Y2[range(n/2)]

fig, ax = plt.subplots(2, 1)
ax[0].plot(time,dep, 'b')
ax[0].set_xlabel('Time')
ax[0].set_ylabel('Amplitude')

#ax[0].plot(time,dep2, 'g')
#ax[0].legend(['', 'coase'])

#ax[0].set_ylim(-0.008, 0.008)
#ax[0].set_xlim(0, 1)

ax[1].plot(frq,abs(Y),'b') # plotting the spectrum
ax[1].set_xlabel('Freq (Hz)')
ax[1].set_ylabel('|Y(freq)|')

f_max = frq[0]
for i in range(0, n/2):
	if (abs(Y[i]) > abs(Y[f_max])):
		f_max = i
		print(f_max)
	
print('f_max', frq[f_max])

#print('f_max', frq(index))
#ax[1].plot(frq,abs(Y2),'g')

plt.show()
#plot_url = py.plot_mpl(fig, filename='mpl-basic-fft')

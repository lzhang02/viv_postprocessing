import numpy as np
from scipy.optimize import leastsq
import pylab as plt

res = np.loadtxt('Fy.dat')

# parameters
rho = 1.2
U = 5.47
d_alu = 0.003902
# cercle inscrit
D = 0.013804*2
mu = 1.813e-5
# dimensionless force
res[:, 1] = res[:, 1]/(0.5*rho*U**2*D)/0.01

t = res[10000:, 0]
data = res[10000:, 1] # create artificial data with noise
f = 41.2 # Hz

# nombres adimensionnels


#N = 1000 # number of data points
#t = np.linspace(0, 4*np.pi, N)
#data = 3.0*np.sin(t+0.005) + 0.5 #+ np.random.randn(N)

#guess_mean = np.mean(data)
#guess_std = 3*np.std(data)/(2**0.5)
#guess_phase = 0

guess_mean = np.mean(data)
guess_std = 0.003#*np.std(data)/(2**0.004)
guess_phase = 0

print(guess_mean)
print(guess_std)
print(guess_phase)
# we'll use this to plot our first estimate. This might already be good enough for you
data_first_guess = guess_std*np.sin(2*np.pi*f*t+guess_phase) + guess_mean

# Define the function to optimize, in this case, we want to minimize the difference
# between the actual data and our "guessed" parameters
optimize_func = lambda x: x[0]*np.sin(2*np.pi*f*t+x[1]) + x[2] - data
est_std, est_phase, est_mean = leastsq(optimize_func, [guess_std, guess_phase, guess_mean])[0]

# recreate the fitted curve using the optimized parameters
data_fit = est_std*np.sin(2*np.pi*f*t+est_phase) + est_mean

print('C', est_mean)
print('C0', est_std)
print('est_phase', est_phase)
print('St', f*D/U)
print('Re', rho*D*U/mu)

fig, ax = plt.subplots()
ax.set_xlabel('Time (s)')
ax.set_ylabel('Dimensionless lift force (-)')

ax.plot(res[:, 0], res[:, 1], '-', label='numerical result')
ax.plot(t, data_fit, label='sinus approximation')
#plt.plot(data_first_guess, label='first guess')

ax.legend(loc = 3, bbox_to_anchor=(0.,0.75))
#ax.set_ylim(0, 1)
#plt.legend()
plt.savefig('Fl.eps', format='eps')#, dpi=1000)
plt.show()

# calcul de RMS
mydata = data**2
print(data)
print(mydata)
print(sum(mydata))
print(np.sqrt(sum(mydata)/t.shape[0]))


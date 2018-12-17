import matplotlib.pyplot as plt
import numpy as np

################################### !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! #####################################################
nb_proc = 30 # nb de processus dans le calcul
period_output = 10 # periode de sauvegarde, definie dans cs_user_extra_operations-global_efforts.f90
nb_digits = 2 # on suppose que moins de 100 processeurs sont utilises dans le calcul
list_ntcabs = np.arange(10000,15000,1) # start stop step
location = '/home/lzhang/Postdoc_Rte/FSI/Saturne/Cas4publication/10_VIV_cylindre_fixe/Turb_k_omega_SST/RESU/Fusion_20181214-1918/'#20180918-1940/'

################################### !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! #####################################################


fig, ax = plt.subplots()
fig2, ax2 = plt.subplots()

nb_cells_large = np.zeros([nb_proc, 1])
for ii in range(0, nb_proc):
	n_digit = len(str(ii))
	str_blank = ''
	for jj in range(0, nb_digits-n_digit):
		str_blank = str_blank+' '

	#print(n_digit)
	
	file_name = "nb_cells"+str_blank+str(ii)+".dat"
	#print(file_name)
	res = int(np.loadtxt(location+file_name))
	#print(res)
	nb_cells_large[ii] = res

Nthetamax = int(np.sum(nb_cells_large)) # nb des cellues dans la direction circonferentielle
# initialisation
res_large = np.zeros([Nthetamax, 4])
res_local = np.zeros([Nthetamax, 4])

nb_points_begin = 0
nb_points_end = 0
for ii in range(0, nb_proc):
	if (nb_cells_large[ii] > 0):	
		# number of digits
		n_digit = len(str(ii))
		#print(n_digit)
		# l'espace dans le nom des fichiers
		str_blank = ''
		for jj in range(0, nb_digits-n_digit):
			str_blank = str_blank+' '

		file_name = "yplus"+str_blank+str(ii)+".dat"
		print(file_name)
	
		res = np.loadtxt(location+file_name)
		#print(res)
		nb_points_end = nb_points_end + nb_cells_large[ii]

		res_large[nb_points_begin:nb_points_end, :] = 0.
		#print(ii, nb_points_begin, nb_points_end)
		for ntcabs in list_ntcabs:
			print(ntcabs)
			res_local[:, :] = 0.
			res_local = res[(ntcabs/period_output-1)*nb_cells_large[ii]:(ntcabs/period_output)*nb_cells_large[ii], :]
			#print(res_large[nb_points_begin:nb_points_end, :])
			#print(res_local[0:nb_cells_large[ii], :])
			res_large[nb_points_begin:nb_points_end, :] = res_large[nb_points_begin:nb_points_end, :] + res_local[0:nb_cells_large[ii], :]
#			print("res_local", res_local)
		

		res_large[nb_points_begin:nb_points_end, :] = res_large[nb_points_begin:nb_points_end, :]/list_ntcabs.shape[0]#/len(list_ntcabs)#.shape[0]
		#print(len(list_ntcabs), list_ntcabs.shape[0])

		nb_points_begin = nb_points_end
		
#print("res_large", nb_points_end, res_large[0:nb_points_end, :])
theta = np.zeros([Nthetamax, 1])
yplus = np.zeros([Nthetamax, 1])

# plot pressure in function of theta
for ii in range(0, Nthetamax):
	x = res_large[ii, 1]
	y = res_large[ii, 2]
	yplus[ii] = res_large[ii, 3]
	r = np.sqrt(x**2 + y**2)

	#dx = 
	
	if (x>0):
		theta[ii] = (np.pi - np.arcsin(y/r))/np.pi*180.
	elif (x<0 and y>0):
		theta[ii] = (np.arcsin(y/r))/np.pi*180.
	else:
		theta[ii] = (2*np.pi + np.arcsin(y/r))/np.pi*180.


cof = 200.
xlim_max = 0
xlim_min = 0
ylim_max = 0
ylim_min = 0
# plot p*vec(n) around the cylinder
for ii in range(0, Nthetamax):
	x = res_large[ii, 1]
	y = res_large[ii, 2]
	yplus_local = yplus[ii, 0]
	dx = yplus_local/np.sqrt((y/x)**2+1.)
	dy = y*dx/x
	#print(dx, dy, yplus_local)
	if (x>0):
		ax2.arrow(cof*x, cof*y, dx, dy, head_width=0.05, head_length=0.1, fc='k', ec='k')
		ax2.plot(cof*x+1/np.sqrt((y/x)**2+1.), cof*y+y*1/np.sqrt((y/x)**2+1.)/x, '.')
		if (cof*x + dx > xlim_max):
			xlim_max = cof*x + dx
		if (cof*y + dy > ylim_max):
			ylim_max = cof*y + dy
		if (cof*x + dx < xlim_min):
			xlim_min = cof*x + dx
		if (cof*y + dy < ylim_min):
			ylim_min = cof*y + dy
	else:
		ax2.arrow(cof*x, cof*y, -dx, -dy, head_width=0.05, head_length=0.1, fc='k', ec='k')
		ax2.plot(cof*x-1/np.sqrt((y/x)**2+1.), cof*y-y*1/np.sqrt((y/x)**2+1.)/x, '.')
		if (cof*x - dx > xlim_max):
			xlim_max = cof*x - dx		
		if (cof*y - dy > ylim_max):
			ylim_max = cof*y - dy
		if (cof*x - dx < xlim_min):
			xlim_min = cof*x - dx
		if (cof*y - dy < ylim_min):
			ylim_min = cof*y - dy
	
lim_max = np.max([xlim_max, ylim_max])
lim_min = np.min([xlim_min, ylim_min])
lim = np.max([lim_max, -lim_min])
lim = lim + 1.
ax2.set_xlim(-lim, lim)
ax2.set_ylim(-lim, lim)
ax2.set_aspect('equal')


ax.plot(theta, yplus, '+')
#ax[0].set_ylim(-0.008, 0.008)
ax.set_xlim(0., 360.)

print('averaged yplus:' + str(np.sum(yplus)/Nthetamax))

plt.show()

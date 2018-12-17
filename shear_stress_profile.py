import matplotlib.pyplot as plt
import numpy as np

################################### !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! #####################################################
nb_proc = 30 # nb de processus dans le calcul
period_output = 10 # periode de sauvegarde, definie dans cs_user_extra_operations-global_efforts.f90
nb_digits = 2 # on suppose que moins de 100 processeurs sont utilises dans le calcul
list_ntcabs = np.arange(10000,20000,1) # start stop step
location = '/home/lzhang/Postdoc_Rte/FSI/Saturne/Cas4publication/10_VIV_cylindre_fixe/Turb_k_omega_SST/RESU/Fusion_20181214-1918/'#20180918-1940/'

################################### !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! #####################################################


fig, ax = plt.subplots()
#fig2, ax2 = plt.subplots()

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
res_large = np.zeros([Nthetamax, 5])
res_local = np.zeros([Nthetamax, 5])

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

		file_name = "force"+str_blank+str(ii)+".dat"
		print(file_name)
	
		res = np.loadtxt(location+file_name)
		#print(res)
		nb_points_end = nb_points_end + nb_cells_large[ii]

		res_large[nb_points_begin:nb_points_end, :] = 0.
		#print(ii, nb_points_begin, nb_points_end)
		for ntcabs in list_ntcabs:
			print(ntcabs)
			res_local[:, :] = 0.
			res_local[0:nb_cells_large[ii], :] = res[(ntcabs/period_output-1)*nb_cells_large[ii]:(ntcabs/period_output)*nb_cells_large[ii], :]
			#print(res_large[nb_points_begin:nb_points_end, :])
			#print(res_local[0:nb_cells_large[ii], :])
			res_large[nb_points_begin:nb_points_end, :] = res_large[nb_points_begin:nb_points_end, :] + res_local[0:nb_cells_large[ii], :]
#			print("res_local", res_local)
		

		res_large[nb_points_begin:nb_points_end, :] = res_large[nb_points_begin:nb_points_end, :]/list_ntcabs.shape[0]#/len(list_ntcabs)#.shape[0]
		#print(len(list_ntcabs), list_ntcabs.shape[0])

		nb_points_begin = nb_points_end
		
#print("res_large", nb_points_end, res_large[0:nb_points_end, :])
theta = np.zeros([Nthetamax, 1])
Fx = np.zeros([Nthetamax, 1])
Fy = np.zeros([Nthetamax, 1])
Ftan = np.zeros([Nthetamax, 1])

# plot pressure in function of theta
for ii in range(0, Nthetamax):
	x = res_large[ii, 1]
	y = res_large[ii, 2]
	Fx[ii] = res_large[ii, 3]
	Fy[ii] = res_large[ii, 4]
	r = np.sqrt(x**2 + y**2)

	#dx = 
	
	if (x>0 and y>0):
		theta[ii] = (np.pi - np.arcsin(y/r))/np.pi*180.
		theta_local = np.arcsin(y/r)
		Ftan[ii] = Fy[ii]*np.cos(theta_local) - Fx[ii]*np.sin(theta_local)
	elif (x>0 and y<0):
		theta[ii] = (np.pi - np.arcsin(y/r))/np.pi*180.
		theta_local = -np.arcsin(y/r)
		Ftan[ii] = -Fy[ii]*np.cos(theta_local) - Fx[ii]*np.sin(theta_local)
	elif (x<0 and y>0):
		theta[ii] = (np.arcsin(y/r))/np.pi*180.
		theta_local = np.arcsin(y/r)
		Ftan[ii] = -Fy[ii]*np.cos(theta_local) - Fx[ii]*np.sin(theta_local)
	else:
		theta[ii] = (2*np.pi + np.arcsin(y/r))/np.pi*180.
		theta_local = -np.arcsin(y/r)
		Ftan[ii] = Fy[ii]*np.cos(theta_local) - Fx[ii]*np.sin(theta_local)




############################# attention le signe ici #######################, ce qu'on veut est la force de la paroi aux fluides, saturne donne la force des fluides a la paroi
ax.plot(theta, -Ftan, '+') 
#ax[0].set_ylim(-0.008, 0.008)
ax.set_xlim(0., 360.)

plt.show()


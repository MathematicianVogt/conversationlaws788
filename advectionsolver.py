num_eqn = 1
num_waves = 1
import time
import numpy as np

def adsolver(q_l,q_r,aux_l,aux_r,problem_data):
	num_rp = q_l.shape[1]
	# Output arrays
	wave = np.empty( (num_eqn, num_waves, num_rp) )
	s = np.empty( (num_waves, num_rp) )
	amdq = np.empty( (num_eqn, num_rp) )
	apdq = np.empty( (num_eqn, num_rp) )
	print q_l
	print q_r
	print "dongs"
	# Basic solve
	wave[0,:,:] = q_r - q_l
	s[0,:] = problem_data['a']
	print problem_data['a']
	s_index = np.zeros((2,num_rp))
	s_index[0,:] = s[0,:]
	
	if(problem_data['a']>=0):
		amdq[0,:] = 0.0* wave[0,0,:]
		apdq[0,:] = problem_data['a']* wave[0,0,:]
	else:
		amdq[0,:] = problem_data['a'] * wave[0,0,:]
		apdq[0,:] = 0.0 * wave[0,0,:]


	return wave, s, amdq, apdq
num_eqn = 1
num_waves = 1

import numpy as np

def solver(q_l,q_r,aux_l,aux_r,problem_data):

	num_rp = q_l.shape[1]
	# Output arrays
	wave = np.empty( (num_eqn, num_waves, num_rp) )
	s = np.empty( (num_waves, num_rp) )
	amdq = np.empty( (num_eqn, num_rp) )
	apdq = np.empty( (num_eqn, num_rp) )

	# Basic solve
	wave[0,:,:] = q_r - q_l
	s[0,:] = problem_data['a']

	s_index = np.zeros((2,num_rp))
	s_index[0,:] = s[0,:]
	amdq[0,:] = np.min(s_index,axis=0) * wave[0,0,:]
	apdq[0,:] = np.max(s_index,axis=0) * wave[0,0,:]
	


	return wave, s, amdq, apdq
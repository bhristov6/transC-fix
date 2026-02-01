from base import *

# --------------------------------------------------------------------------- #
# Do a Random Walk
# --------------------------------------------------------------------------- #	

# Do a random walk analytically
def do_random_walk(a, gamma, start_locs, fname):
	start_time = time.time()
	nbins, pi = len(a), [0]
	
	# create the transitional matrix P
	p = np.zeros((nbins,nbins))
	
	# for every node i
	for i in range(nbins):
		outgoing = sum(a[i])
		
		# find the probability to move to node j
		for j in range(nbins):			
			# you can get via i -> j or you restart at j
			if outgoing > 0:
				p[i][j] = (1-gamma)*a[i][j]/outgoing + gamma*1.0*(j in start_locs)/len(start_locs)
			# special case: all zeros in that row i
			else:
				p[i][j] = gamma*1.0*(j in start_locs)/len(start_locs)

	# sanity test
	for i in range(nbins):
		break
		if not math.isclose(sum(p[i]),1):
			print(f"i = {i}; {sum(p[i])}")
			#sys.exit()
			pass

	
	# find the left eigenvector corresponding to eigenvalue 1
	w, vl = eig(p,left=True, right=False)
	
	# extract it (the real part)
	prec = False
	for i in range(len(w)):
		if math.isclose(w[i].real,1, rel_tol=1e-6):
			pi = vl[:,i]
			prec = True
			break

	# if there was numeric instability during eigenvector decomposition
	if not prec:
		eig_vl, eig_dist = [], []
		for i in range(len(w)):
			if math.isclose(w[i].real,1, rel_tol=1e-2):
				eig_vl.append(i)
				eig_dist.append(abs(w[i].real - 1.0))
		
		if len(eig_dist) > 0:
			idx = np.argsort(eig_dist)[0]
			pi = vl[:,eig_vl[idx]]
			prec = True
		else:
			print(f"couldn't carry the eigenvector decomposition\n")
			sys.exit(1)
			
	
	# normalize the vector
	pi = [float(i.real)/sum(pi).real for i in pi]
	
	write_pi(fname, pi)
	return pi
	
	
def write_pi(fname, pi):
	fo = open(f"{fname}/bin_scores.txt", "w")
	for i in reversed(np.argsort(pi)):
		fo.write(f"{i}\t{pi[i]}\n")
	fo.close()


if __name__ == "__main__":
	pass
	
	

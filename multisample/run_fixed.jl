

using HDF5


include("../pseudogp_student_fixed.jl"); # load in inference & plotting methods;


function run_fixed(i; h5infile = "../data/multisample.h5", 
	h5outfile = "/net/isi-project/CW010_CAMPBELL_SCNGEN/data/GP/multisample/multitrace.h5",
	niter = 2e5, tau_init_file = "../data/tauchain.h5")

	sample_i = string("sample_", i)
	X = h5read(h5infile, string(sample_i, "/X"))

	srand(123)

	tau_chain = h5read(tau_init_file, "tauchain")
	bp = int(size(tau_chain)[1] / 2)

	include("../pseudogp_student_fixed.jl");
	srand(123)

	N, P = size(X)


	n_iter = int(5e5)
	burn = n_iter / 2

	thin = int(n_iter / 1000)

	# pseudotime parameters
	eps = 1e-6
	t =  rand(Uniform(.5 - eps, .5 + eps), N) # t_gt #  
	tvar = 6.5e-3

	# GP function
	s = [30, 30]
	tau = fill(0.0, size(X))

	for j in 1:P
	    for i in 1:N
	        tau[i,j] = mean(vec(tau_chain[bp:end,j,j]))
	    end
	end

	# kernel parameters
	lambda = [50, 10]

	r = 1e-3 # repulsion parameter

	return_burn = true # should the burn period be returned in the traces?

	gamma = 50000.0

	mhf = pseudogp_student_fixed(X, n_iter, burn, thin, tau, 
		t, tvar, lambda, s, r = r, gamma = gamma)

	tchain = mhf["tchain"]
	h5write(h5outfile, string(sample_i, "/tchain"), tchain)
end

i = ARGS[1]

run_fixed(i)

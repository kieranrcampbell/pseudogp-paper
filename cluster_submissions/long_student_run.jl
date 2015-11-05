
# coding: utf-8

# In[1]:

using Gadfly
using HDF5
using DataFrames

include("../pseudogp_student.jl"); # load in inference & plotting methods;


# In[2]:

X = h5read("../data/embeddings.h5", "monocle/embedding")
t_gt = h5read("../data/embeddings.h5", "monocle/pseudotime");



# In[4]:

# remove cells less than 0 on x
to_keep = X[:,1] .> 0
X = X[to_keep,:]
t_gt = t_gt[to_keep]

x = X[:,1] ; y = X[:,2]
xs = (x - mean(x)) / sqrt(var(x))
ys = (y - mean(y)) / sqrt(var(y))

#plot(x = X[:,1], y = X[:,2], colour = t_gt)
X_old = X
X = [xs ys];


#### Inference

srand(123)

n = size(X)[1]

n_iter = int(10e6)
burn = n_iter / 2

thin = int(n_iter / 3000)

# pseudotime parameters
eps = 1e-6
t =  rand(Uniform(.5 - eps, .5 + eps), n) # t_gt #  
tvar = 6.5e-3

# GP function
guess_at_variance = 5e-2
tau = fill(1 / guess_at_variance, size(X))
tauvar = [1, 1] * 2e-1

# kernel parameters
lambda_init = [10, 10]
sigma_init = [1 / guess_at_variance, 1 / guess_at_variance]

lvar = [1, 1] * 1e-1
svar = tauvar # [1, 1] * 1

r = 1e-3 # repulsion parameter

return_burn = true # should the burn period be returned in the traces?
cell_swap_probability = 0 # randomly swap two cells at each stage?

gamma = 50000.0

mh1 = pseudogp_student(X, n_iter, burn, thin, tau, tauvar, t, tvar, lambda_init, lvar, 
sigma_init, svar, r = r, gamma = gamma)


# In[31]:

## write traces to file for plotting with R
trace_file = "../data/10m_run_with_tau_traces.h5"

h5write(trace_file,  "tchain", mh1["tchain"])
h5write(trace_file, "tauchain", mh1["tau_chain"])
h5write(trace_file, "lambda_chain", mh1["lambda_chain"])
h5write(trace_file, "sigma_chain", mh1["sigma_chain"])
h5write(trace_file, "n_iter",  n_iter,)
h5write(trace_file, "prior", mh1["prior_chain"])
h5write(trace_file, "lambda_prior", mh1["lambda_prior_chain"])
h5write(trace_file, "t_gt", t_gt)
h5write(trace_file, "X", X)
h5write(trace_file, "to_keep", int(to_keep))


h5write(trace_file, "tvar", tvar)
h5write(trace_file,  "lvar", lvar)
h5write(trace_file, "svar", svar)


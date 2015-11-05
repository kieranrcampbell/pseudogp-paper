
# coding: utf-8

# In[29]:

using Gadfly
using HDF5
using DataFrames

include("../bgplvm.jl"); # load in inference & plotting methods;


# In[2]:

X = h5read("../data/embeddings.h5", "monocle/embedding")
t_gt = h5read("../data/embeddings.h5", "monocle/pseudotime");



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

# In[35]:


srand(123)

n = size(X)[1]

n_iter = int(4e6)
burn = n_iter / 2

thin = int(n_iter / 2000)

# pseudotime parameters
eps = 1e-6
t =  rand(Uniform(.5 - eps, .5 + eps), n) # t_gt #  
tvar = 6.5e-3

# kernel parameters
lambda_init = [1, 1] * 50
s_init = [1, 1] * 10

lvar = [1, 1] * 5e-2
svar = [1, 1] * 5e-1

r = 1e-3 # repulsion parameter

return_burn = true # should the burn period be returned in the traces?
gamma = 5e4

mh1 = B_GPLVM_MH(X, n_iter, burn, thin, t, tvar, lambda_init, lvar, 
s_init, svar, r = r, gamma = gamma)



trace_file = "../data/monocle_gausslik_traces_4m_gamma_5e4.h5"

h5write(trace_file,  "tchain", mh1["tchain"])
h5write(trace_file, "lambda_chain", mh1["lambda_chain"])
h5write(trace_file, "s_chain", mh1["s_chain"])
h5write(trace_file, "n_iter",  n_iter,)
h5write(trace_file, "log_lik", mh1["loglik_chain"])
h5write(trace_file, "prior", mh1["prior_chain"])
h5write(trace_file, "t_gt", t_gt)
h5write(trace_file, "X", X)
h5write(trace_file, "to_keep", int(to_keep))


h5write(trace_file, "tvar", tvar)
h5write(trace_file,  "lvar", lvar)
h5write(trace_file, "svar", svar)

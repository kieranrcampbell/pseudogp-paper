
using HDF5
using DataFrames
using Distributions


X = h5read("../data/embeddings.h5", "monocle/embedding")
t_gt = h5read("../data/embeddings.h5", "monocle/pseudotime");



# remove cells less than 0 on x
to_keep = X[:,1] .> 0
X = X[to_keep,:]
t_gt = t_gt[to_keep]

x = X[:,1] ; y = X[:,2]
xs = (x - mean(x)) / sqrt(var(x))
ys = (y - mean(y)) / sqrt(var(y))



X_old = X
X = [xs ys];

srand(123)
n_multisample = 100
n = size(X)[1]
sample_size = int(n/2)
output_file = "../data/multisample.h5"

for i in 1:n_multisample
	sample_name = string("sample_", i)
	to_sample = sample(1:n, sample_size)
	X_s = X[to_sample,:]

	h5write(output_file, string(sample_name, "/to_sample"), to_sample)
	h5write(output_file, string(sample_name, "X"), X_s)
end

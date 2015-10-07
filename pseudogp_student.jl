#=
Bayesian Gaussian Process latent variable models for
pseudotime inference in single-cell RNA-seq data


kieran.campbell@sjc.ox.ac.uk=#

using Distributions
using Gadfly
using DataFrames
using StatsBase

## log-pdf for non-standardized student-t distribution
function nStudent(x::Real, mu::Real, sigma::Real, df::Real)
    larg = 1 + 1 / df * ( (x - mu)/sigma )^2

    ll = lgamma( (df+1)/2 ) - lgamma(df / 2)
    ll -= 0.5 * log(pi * df) - log(sigma)
    ll -= (df + 1) / 2 * log(larg)

    larg
end

function pairwise_distance(t)
    n = length(t)
    T = zeros((n,n))
    for i in 1:n
        for j in (i+1):n
            T[i,j] = (t[i] - t[j])^2
        end
    end
    return T + transpose(T)
end

function cross_pairwise_distance(t1, t2)
    n1 = length(t1)
    n2 = length(t2)
    T = zeros((n1, n2))
    for i in 1:n1
        for j in 1:n2
            T[i,j] = (t1[i] - t2[j])^2
        end
    end
    return T
end

# Covariance matrix with normal noise (sigma)
function covariance_matrix(t, lambda, sigma)
    T = pairwise_distance(t)
    Sigma = exp(-lambda * T) + sigma * eye(length(t))
    return Sigma
end

# Covariance matrix no normal noise 
function K(t, lambda, jitter = 1e-6)
    T = pairwise_distance(t)
    Sigma = exp(-lambda * T)  + jitter * eye(length(t))
    return Sigma
end


function covariance_matrix_hnoise(t, lambda, sigma)
    T = pairwise_distance(t)
    Sigma = exp(-lambda * T) + diagm(sigma)
    return Sigma
end

function cross_covariance_matrix(t1, t2, lambda)
    T = cross_pairwise_distance(t1, t2)
    return exp(-lambda * T)
end

## Likleihood functions

# Student's log likelihood for data given GP
function student_likelihood(X, mu, sigma, df)
    N, P = size(X)
    @assert length(sigma) == P
    @assert size(X) == size(mu)

    ll = 0
    for i in 1:N
        for j in 1:P
            ll += nStudent(X[i,j], mu[i,j], sigma[j], df)
        end
    end
    return ll
end

# Probability of GP function given kernel parameters
function GP_density(mu, t, lambda)
    N, P = size(mu)

    @assert length(lambda) == P
    @assert length(t) == N

    ll = 0

    for j in 1:P
        ll += sum(logpdf(MultivariateNormal(zeros(length(t)),K(t, lambda[j])), mu[:,j]))
    end

    return ll
end


## electroGP
function corp_prior(t, r = 1)
    if r == 0 # uniform prior
        return 0
    end
    
    ll = 0
    n = length(t)
    for j in 1:n
        for i in (j+1):n
            ll += log(sin(pi * abs(t[i] - t[j])))
        end
    end
    return 2 * r * ll
end

function lambda_prior(lambda, rate = 1.0)
    lp = sum(logpdf(Exponential(rate), lambda))
    # lp = 0
    return lp
end

# function sigma_prior(sigma, rate = 1.0)
#     # sp = sum(logpdf(InverseGamma(alpha, beta), sigma))
#     sp = sum(logpdf(Exponential(rate), sigma))
#     # sp = 0
#     return sp
# end

function sigma_prior(sigma, alpha = 1.0, beta = 1.0)
    sp = sum(logpdf(InverseGamma(alpha, beta), sigma))
    return sp
end


function acceptance_ratio(X, tp, t, mu_p, mu, lambda_prop, lambda, 
    sigma_prop, sigma, r, s, gamma, df)
    """ 
    Compute the acceptance ratio for 
    @param X N-by-D data array for N points in D dimensions
    @param tp Proposed pseudotime of length N
    @param t Previous pseudotime length N
    @param thetap Propose theta = [lambda, sigma]
    @param theta Previous theta
    @param r > 0 Corp parameter
    @param s Tempering parameter: (log likelihood * prior) is raised to this value
    @param gamma Rate for exponential prior on lambda
    """
    s_likelihood = student_likelihood(X, mu_p, sigma_prop, df) - student_likelihood(X, mu, sigma, df)
    gp_likelihood = GP_density(mu_p, tp, lambda_prop) - GP_density(mu, t, lambda)
    t_prior = corp_prior(tp, r) - corp_prior(t, r)
    l_prior = lambda_prior(lambda_prop, gamma) - lambda_prior(lambda, gamma)
    s_prior = sigma_prior(sigma_prop, delta) - sigma_prior(sigma, delta)
    return s_likelihood + gp_likelihood + t_prior + l_prior + s_prior
end

function couple_update_acceptance_ratio(X, t1, t2, theta1, theta2, r, s1, s2)
    h(X, t, theta, r) = log_likelihood(X, t, theta) + corp_prior(t, r) 
    return  ( (s1 - s2) * ( h(X, t2, theta2, r) - h(X, t1, theta1, r) ) )
end

function couple_update_acceptance_ratio_likelihood_only(X, t1, t2, theta1, theta2, r, s1, s2)
    """ Same as couple_acceptance_ratio except only the likelihood is raised to s """
    h(X, t, theta, r) = log_likelihood(X, t, theta) 
    return  ( (s1 - s2) * ( h(X, t2, theta2, r) - h(X, t1, theta1, r) ) )
end
   
# TODO: propose and propose_t now same function
function propose(mean, var, lower = 0, upper = Inf)
    # sample from truncated normal of (mean, real)
    n = length(mean)
    if length(var) == 1
        var = fill(var, n)
    end

    return [rand(Truncated(Normal(mean[i], var[i]), lower, upper)) for i in 1:n]
end

function propose_t(t, var, lower = 0, upper = 1)
    #= random truncated normal for mean vector t and scalar sigma var =#
    n = length(t)
    tp = [rand(Truncated(Normal(t[i], var), lower, upper)) for i in 1:n]
    return tp
end;


function pseudogp_student(X, n_iter, burn, thin, 
    mu, muvar, t, tvar, lambda, lvar, sigma, svar, 
    r = 1, return_burn = false, cell_swap_probability = 0,
    gamma = 1.0, df = 1.0)
    
    chain_size = int(floor(n_iter / thin)) + 1 # size of the thinned chain
    burn_thin = int(floor(burn / thin)) # size of the burn region of the thinned chain
    
    N, P = size(X)

    @assert P == 2 # for now
    @assert cell_swap_probability >= 0
    @assert cell_swap_probability <= 1
    @assert length(lambda) == length(sigma) == P
    @assert burn < n_iter
    @assert length(t) == N
    @assert length(lvar) == length(svar) == P
    @assert size(X) == size(mu)
    @assert length(muvar) == P

    ## chains
    tchain = zeros((chain_size, N))
    tchain[1,:] = t

    muchain = zeros((chain_size, N, P))
    muchain[1,:,:] = mu

    lambda_chain = zeros(chain_size, P)
    lambda_chain[1,:] = lambda

    sigma_chain = zeros(chain_size, P)
    sigma_chain[1,:] = sigma
    
    accepted = zeros(n_iter)

    student_chain = zeros(chain_size)
    GP_chain = zeros(chain_size)
    prior_chain = zeros(chain_size)
    lambda_prior_chain = zeros(chain_size)
    
    student_chain[1] = student_likelihood(X, mu, sigma, df)
    GP_chain[1] = GP_density(mu, t, lambda)

    prior_chain[1] = corp_prior(t, r)
    lambda_prior_chain[1] = lambda_prior(lambda, gamma)
    #println(gamma)

    # alpha_chain = zeros(chain_size - 1)
    # rnd_chain = zeros(chain_size - 1)
    
    ## MH
    for i in 1:n_iter
        # proposals
        t_prop = propose_t(t, tvar)
        mu_prop = fill(NaN, size(mu))
        for j in 1:P
            mu_prop[:,j] = propose(mu[:,j], muvar[j], -Inf, Inf)
        end

        lambda_prop = propose(lambda, lvar)
        sigma_prop = propose(sigma, svar)

        if cell_swap_probability > 0
            if rand() < cell_swap_probability
                # swap two cells at random
                to_swap = sample(1:length(t), 2, replace = false)
                t_prop[to_swap] = t_prop[reverse(to_swap)]
            end
        end

        # calculate acceptance ratio
        alpha = acceptance_ratio(X, t_prop, t, 
                                mu_prop, mu,
                                lambda_prop, lambda, 
                                sigma_prop, sigma, 
                                r, 1, gamma, df)

        rnd = log(rand())

        # accept - reject
        if alpha > rnd
            # accept
            accepted[i] = 1
            t = t_prop
            mu = mu_prop
            lambda = lambda_prop
            sigma = sigma_prop
        end
        
    
        if i % thin == 0
            # update traces
            j = int(i / thin) + 1
            tchain[j,:] = t
            muchain[j,:,:] = mu
            lambda_chain[j,:] = lambda
            sigma_chain[j,:] = sigma
            student_chain[j] = student_likelihood(X, mu, sigma, df)
            GP_chain[j] = GP_density(mu, t, lambda)
            prior_chain[j] = corp_prior(t, r)
            lambda_prior_chain[j] = lambda_prior(lambda, gamma)
        end
    end
    
    burnt = burn_thin
    if return_burn
        burnt = 1
    end
    
    rdict = {"tchain" => tchain[burnt:end,:],
        "lambda_chain" => lambda_chain[burnt:end,:],
        "sigma_chain" => sigma_chain[burnt:end,:],
        "acceptance_rate" => (sum(accepted) / length(accepted)),
        "burn_acceptance_rate" => (sum(accepted[burnt:end]) / length(accepted[burnt:end])),
        "r" => r,
        "muchain" => muchain,
        "student_chain" => student_chain,
        "GP_chain" => GP_chain,
        "prior_chain" => prior_chain,
        "lambda_prior_chain" => lambda_prior_chain,
        "params" => {"n_iter" => n_iter,
                    "burn" => burn,
                    "thin" => thin,
                    "burn_thin" => burn_thin
            }
        }
    return rdict
end




#### Posterior predictive mean

# Posterior predictive mean given by
# $$ \mathbf{\mu_x} = K^T_* K^{-1} \mathbf{x} $$
# $$ \mathbf{\mu_y} = K^T_* K^{-1} \mathbf{y} $$
# where $K_*$ is the covariance matrix between the latent $\mathbf{t}$ and the predictive $\mathbf{t'}$, where typically $\mathbf{t'}$ is sampled as $m$ equally spaced points on the interval [0, 1].


function predict(tp, t_map, lambda_map, sigma_map, X)
    #= Returns MAP prediction of mean function given:
    @param tp Values of t at which to predict function
    @param t_map Map estimate of latent pseudotimes
    @param lambda_map Map estimate of lambda
    @param sigma_map Map estimate of sigma
    =#
    @assert length(lambda_map) == length(sigma_map) == size(X)[2]
    @assert length(t_map) == size(X)[1]
    ndim = size(X)[2]
    np = length(tp)

    Xp = fill(0.0, (np, ndim))
    for i in 1:ndim
        K_map = covariance_matrix(t_map, lambda_map[i], sigma_map[i])
        K_star_transpose = cross_covariance_matrix(tp, t_map, lambda_map[i])

        matrix_prefactor = K_star_transpose * inv(K_map)

        mu = matrix_prefactor * X[:,i]
        Xp[:,i] = vec(mu)   
    end

    return Xp
end;


#----------- Plotting functions

function plot_pseudotime_trace(mh)
    nchoose = 4
    chosen = sample(1:n, nchoose)

    df = convert(DataFrame, mh["tchain"][:, chosen])
    df[:iter] = 1:(size(df)[1])
    df = convert(DataFrame, df)
    df_melted = stack(df, [1:nchoose])
    names!(df_melted, [symbol(x) for x in ["variable", "value", "iter"]])

    return Gadfly.plot(df_melted, x = "iter", y = "value", color = "variable", Geom.line)  
end

function plot_kernel_parameter(mh, param)
    chain_name = string(param, "_chain")

    df = convert(DataFrame, mh[chain_name])
    ndim = size(df)[2]

    names!(df, [symbol(string(param, string(i))) for i in 1:ndim])
    df[:iter] = 1:(size(df)[1]) # (burn + 2)
    df_melted = stack(df, [1:2])
    return Gadfly.plot(df_melted, x = "iter", y = "value", colour = "variable", Geom.line)
end

function plot_posterior_mean(mh, tp, X)
    burn = mh["params"]["burn_thin"]
    lambda_map = mean(mh["lambda_chain"][burn:end,:], 1)
    sigma_map = mean(mh["sigma_chain"][burn:end,:], 1)

    t_map = mean(mh["tchain"][burn:end,:], 1)
    mu_p = predict(tp, t_map, lambda_map, sigma_map, X)

    return Gadfly.plot(layer(x = X[:,1], y = X[:,2], color = t_gt, Geom.point) ,
    layer(x = mu_p[:,1], y = mu_p[:,2], Geom.line(preserve_order = 1), 
    Theme(default_color=color("red"))))
end

function plot_likelihood(mh; term = "student_chain")
    df = DataFrame()
    df[:value] = mh[term]
    df[:iter] = 1:(size(df)[1]) # (burn + 2)
    Gadfly.plot(df, x = "iter", y = "value", Geom.line)
end

function plot_prior(mh)
    df = DataFrame()
    df[:value] = mh["prior_chain"]
    df[:iter] = 1:(size(df)[1]) # (burn + 2)
    Gadfly.plot(df, x = "iter", y = "value", Geom.line)
end

function plot_lambda_prior(mh)
    df = DataFrame()
    df[:value] = mh["lambda_prior_chain"]
    df[:iter] = 1:(size(df)[1]) # (burn + 2)
    Gadfly.plot(df, x = "iter", y = "value", Geom.line)
end


#----------------------- Variance measures


function kendall_tau(mh)
    #= Calculates the kendall tau non-parametric correlation
    measure along the chain between consecutive pseudotimes
    returning a vector of length (N - 1) =#
    bt = mh["params"]["burn_thin"]
    tchain = mh["tchain"][bt:end,:]
    kt = zeros(size(tchain)[1] - 1)
    for i in 1:length(kt)
        kt[i] = corkendall(vec(tchain[i+1,:]), vec(tchain[i,:]))
    end
    return kt
end

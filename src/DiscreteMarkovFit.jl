module DiscreteMarkovFit

using Distributions, LinearAlgebra
export ObservedMarkovChain, ObservedBirthDeathChain

# Data: just a vector of integers (that's fun!).

# Prior: we're particular interested in a prior for a birth-death like Markov
# Chain that is rather sparse, i.e. the only non-zero elements in P are the
# elements on the main diagonal and the diagonals above and below it.

# Note we take the transition matrix to be of the form p_ij = P(j → i), since
# in that case we respect julia's column-major ordering.

"""
    ObservedMarkovChain

An observed discrete space, discrete time Markov Chain with associated
Dirichlet prior on the columns of the transition probability matrix.

!!! note
    The transition (probability) matrix follows the convention where
    `P[i,j]`` is the probability of a `j->i` transition. The prior should
    be defined similarly.
"""
struct ObservedMarkovChain{T}
    prior   ::Matrix{T}
    observed::Matrix{Int}
    nstates ::Int
    minstate::Int
end

struct PosteriorSample{T}
    ess::T
    pps::Vector
    πs ::Matrix{T}
    Ps ::Array{T,3}
end

Base.show(io::IO, ps::PosteriorSample) = write(io,
    "ESS = $(ps.ess)\n  π = \n$(join([" ⋅ $k: $v" for (k,v) in ps.pps], "\n"))")

function ObservedMarkovChain(x::Vector{I}, ϵ=:auto) where I<:Integer
    N, nstates, minstate = observed_transitions(x)
    ϵ = ϵ == :auto ? 1.0/nstates : ϵ
    prior = zeros(nstates, nstates) .+ ϵ
    ObservedMarkovChain(prior, N, nstates, minstate)
end

function ObservedBirthDeathChain(x::Vector{I}, ϵ=1/3) where I<:Integer
    N, nstates, minstate = observed_transitions(x)
    prior = birthdeath_transition_prior(nstates, ϵ)
    ObservedMarkovChain(prior, N, nstates, minstate)
end

function observed_transitions(x::Vector{Int})
    mn, mx = extrema(x)
    z = x .- (mn - 1)  # 'normalize' the sequence
    N = zeros(Int, mx-mn+1, mx-mn+1)
    @inbounds for i=2:length(z)
        N[z[i],z[i-1]] += 1
    end
    (N=N, nstates=mx-mn+1, minstate=mn)
end

function birthdeath_transition_prior(i::Int, ϵ::T, small=1e-50) where T
    Π = diagm(0=>repeat(ϵ:ϵ, i), 1=>repeat(ϵ:ϵ, i-1), -1=>repeat(ϵ:ϵ, i-1))
    Π[1,1] = Π[2,1] = Π[i,i] = Π[i-1,i] = 3ϵ/2
    Π .+ small
end

"""
    sample(d::ObservedMarkovChain, n)

Draw `n` samples from the posterior distribution of the transition
probabilities and stationary probabilities. This will also estimate
the effective (iid) sample size associated with the Markov chain.
"""
function sample(d::ObservedMarkovChain, n)
    πs = zeros(d.nstates, n)
    Ps = zeros(d.nstates, d.nstates, n)
    @inbounds for i=1:n
        π, P = sample(d)
        πs[:,i] = π
        Ps[:,:,i] = P
    end
    PosteriorSample(ess(d, πs), pps(d, πs), πs, Ps)
end

"""
    sample(d::ObservedMarkovChain)

Draw a sample from the posterior distribution of the transition
probabilities and stationary probabilities.
"""
function sample(d::ObservedMarkovChain)
    P = zeros(d.nstates, d.nstates)
    for i=1:d.nstates
        P[:,i] = rand(Dirichlet(d.observed[:,i] .+ d.prior[:,i]))
    end
    eig = eigen(P)
    π_unnorm = real.(eig.vectors[:,end])
    return (π=π_unnorm ./ sum(π_unnorm), P=P)
end

# This is not always numerically stable apparently... Maybe implement
# the algorithm of Minka (as in Heck et al.)?
function ess(d::ObservedMarkovChain, πs::Matrix)
    fitted_dir = fit_mle(Dirichlet, πs)
    sum(fitted_dir.alpha) - sum(d.prior)
end

function pps(d::ObservedMarkovChain, πs)
    means = vec(mean(πs, dims=2))
    [i+d.minstate => means[i] for i=1:length(means)]
end

end

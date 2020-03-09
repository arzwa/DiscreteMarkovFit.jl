using DiscreteMarkovFit
using Test, Distributions, Random, DataFrames

Random.seed!(19081994)

@testset "Heck et al. (2019) toy example" begin
    πstat = [0.85, 0.13, 0.02]
    transition(z, β, d) = rand() < β ? z : rand(d)
    t = (z)->transition(z, rand(), Categorical(πstat))
    for j=1:5
        nrep = 500; n=1000; R=500; notcovered = 0
        for i=1:nrep
            z = foldl((x,y)->vcat(x, t(x[end])),
                1:n, init=rand(Categorical(πstat)))
            d = ObservedMarkovChain(z)
            out = DiscreteMarkovFit.sample(d, R) # heck used 5000
            for i=1:length(πstat)
                q1, q2 = quantile(out.πs[i,:], [0.05, 0.95])
                if !(q1 < πstat[i] < q2)
                    notcovered += 1
                end
            end
        end
        @test isapprox(notcovered/(nrep*length(πstat)), 0.1, atol=0.03)
    end
end

@testset "Beluga rjMCMC example" begin
    df = CSV.read(joinpath(@__DIR__, "example-trace.csv"))
    ks = [0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0]
    for i=2:17
        k = Array(df[7501:end,Symbol("k$i")])
        d = ObservedBirthDeathChain(k)
        out = DiscreteMarkovFit.sample(d, 1000)
        @show (i, out.ess)
        mn = vec(mean(out.πs, dims=2))
        @test map_k = findmax(mn)[2] + d.minstate - 1 == ks[i]
    end
    k = Array(df[7501:end,:k])
    d = ObservedBirthDeathChain(k)
    out = DiscreteMarkovFit.sample(d, 10000)
    mn = vec(mean(out.πs, dims=2))
    @test map_k = findmax(mn)[2] + d.minstate - 1 == 3
end

Implements the method of [Heck et al.](https://link.springer.com/article/10.1007/s11222-018-9828-0) for estimating the precision of posterior model probabilities from rjMCMC output

Maybe should find a better name for this module.

## Example

```julia
using DiscreteMarkovFit, DataFrames, CSV
df = CSV.read("test/example-trace.csv")   # reversible-jump MCMC output
k = Array(df[7501:end,:k])   # vector of model indicator variable
d = ObservedBirthDeathChain(k)
out = DiscreteMarkovFit.sample(d, 10000)
DiscreteMarkovFit.describe(d, out)
```

Example results: 

```
julia> DiscreteMarkovFit.describe(d, out)
┌ Info: ESS and posterior probabilities
│   smpl.ess = 1323.3540345266515
│   pps =
│    9-element Array{Pair{Int64,Float64},1}:
│      3 => 0.1671723937018068   
│      4 => 0.32590748729711294  
│      5 => 0.27636169998012283  
│      6 => 0.14395779909238296  
│      7 => 0.06589943580051748  
│      8 => 0.015367630793354906
│      9 => 0.0036001811224925413
│     10 => 0.0007300913023413403
└     11 => 0.0010032809098662373
```

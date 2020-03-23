Implements the method of [Heck et al.](https://link.springer.com/article/10.1007/s11222-018-9828-0) for estimating the precision of posterior model probabilities from rjMCMC output

Maybe should find a better name for this module.

## Example

```julia
using DiscreteMarkovFit, DataFrames, CSV
df = CSV.read("test/example-trace.csv")   # reversible-jump MCMC output
k = Array(df[7501:end,:k])   # vector of model indicator variable
d = ObservedBirthDeathChain(k)
out = DiscreteMarkovFit.sample(d, 10000)
```

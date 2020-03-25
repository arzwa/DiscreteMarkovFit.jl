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

```
ESS = 1367.0964411004118
  π =
 ⋅3 => (mean = 0.167, std = 0.013, q025 = 0.142, q0975 = 0.194)
 ⋅4 => (mean = 0.326, std = 0.012, q025 = 0.302, q0975 = 0.35)
 ⋅5 => (mean = 0.276, std = 0.01, q025 = 0.257, q0975 = 0.297)
 ⋅6 => (mean = 0.144, std = 0.009, q025 = 0.127, q0975 = 0.161)
 ⋅7 => (mean = 0.066, std = 0.007, q025 = 0.053, q0975 = 0.08)
 ⋅8 => (mean = 0.015, std = 0.003, q025 = 0.011, q0975 = 0.022)
 ⋅9 => (mean = 0.004, std = 0.001, q025 = 0.002, q0975 = 0.007)
 ⋅10 => (mean = 0.001, std = 0.001, q025 = 0.0, q0975 = 0.002)
 ⋅11 => (mean = 0.001, std = 0.003, q025 = 0.0, q0975 = 0.006)
```

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
ESS = 1017.1734442883951
  π =
 ⋅ 3: 0.16699099359632893
 ⋅ 4: 0.3255515069769609
 ⋅ 5: 0.2765178577364221
 ⋅ 6: 0.14403391967423124
 ⋅ 7: 0.06597299434358542
 ⋅ 8: 0.015418541657368614
 ⋅ 9: 0.003624478717210618
 ⋅ 10: 0.000737382587628031
 ⋅ 11: 0.0011523247102647948
```

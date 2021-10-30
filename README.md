# Walking Arm Trebuchet

Nice article here: https://makezine.com/projects/walking-arm-trebuchet/

## How to run it

Run `beam_optim.jl`, uncomment these lines:

```julia
res = optimize((params) -> rangeObjective(params)[1],lower,upper,x0,
               Fminbox(inner_optimizer), Optim.Options(iterations=60))
print(res)

params = Optim.minimizer(res)
```

And comment out the `params` that are hardcoded after it. This will optimize
the trebuchet and plot the results. Other useful parameters are `FPS` on line 433
and `liveview` on 534.

This was tested with Julia 1.6.

`beam_optim_small.jl` has modified upper and lower bounds and intended for throwing
a tennis ball, not a pumpkin.
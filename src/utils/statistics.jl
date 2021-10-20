#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2021-Aug-04: add this function to avoid NaN calculations
#
#######################################################################################################################################################################################################
"""
    nanmean(x::Array)

Return the mean of array by ommiting the NaN, given
- `x` Array of numbers, can be NaN

---
Example
```julia
xs = [1, 2, 4, NaN];
nmean = nanmean(xs);
```
"""
function nanmean(x::Array)
    _x = filter(!isnan, x);

    if length(_x) == 0 return NaN end;

    return mean( _x )
end

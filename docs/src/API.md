# Emerald
```@meta
CurrentModule = Emerald
```


## Types

### Soil
```@docs
AbstractSoilVC
BrooksCorey
VanGenuchten
```

### Leaf
```@docs
Leaf
LeafBiophysics
```

### Radiation
```@docs
WaveLengthSet
```


## Folium


## Soli
```@docs
Soli.soil_ψ_25
Soli.soil_ψ_25(vc::Emerald.BrooksCorey{FT}, θ::FT) where {FT<:AbstractFloat}
```


## Utilities
```@docs
calctav
nanmean
transmittance
```

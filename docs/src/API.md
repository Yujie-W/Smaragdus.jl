# Smaragdus
```@meta
CurrentModule = Smaragdus
```


## Types

### Environment
```@docs
AbstractSoilVC
BrooksCorey
VanGenuchten
SoilAir
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
Soli.soil_ψ_25(vc::Smaragdus.BrooksCorey{FT}, θ::FT) where {FT<:AbstractFloat}
```


## Utilities
```@docs
calctav
nanmean
transmittance
```

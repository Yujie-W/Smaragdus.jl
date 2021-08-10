module LeafLevel

using SpecialFunctions: expint
using UnPack: @unpack

# using external types and functions
using ..Emerald: Leaf, WaveLengthSet


function fluspect!(leaf::Leaf{FT}, wls::WaveLengthSet{FT}; APAR_car::Bool = true) where {FT<:AbstractFloat}
    BIO = leaf.BIO;
    @unpack K_ANT, K_CAB, K_CAR_V, K_CAR_Z, K_CBC, K_H₂O, K_LMA, K_PRO, K_SENES, MESOPHYLL_N, NR = BIO;

    # calculate the average absorption feature and relative Cab and Car partitions
    BIO._k_all   .= (K_CAB   .* BIO.cab .+                      # chlorophyll absorption
                     K_CAR_V .* BIO.car .* (1 - BIO.f_zeax) .+  # violaxanthin carotenoid absorption
                     K_CAR_Z .* BIO.car .* BIO.f_zeax .+        # zeaxanthin carotenoid absorption
                     K_ANT   .* BIO.ant .+                      # anthocynanin absorption absorption
                     K_SENES .* BIO.f_senes .+                  # TODO cannot be a fraction
                     K_H₂O   .* BIO.l_H₂O .+                    # water absorption
                     K_CBC   .* BIO.CBC .+                      # carbon-based constituents absorption
                     K_PRO   .* BIO.PRO .+                      # protein absorption
                     K_LMA   .* (BIO.LMA - BIO.CBC - BIO.PRO)   # dry mass absorption (if some remained)
                    ) ./ MESOPHYLL_N;
    BIO.α_cab    .= (K_CAB .* BIO.cab) ./ BIO._k_all ./ MESOPHYLL_N;
    BIO.α_cabcar .= (K_CAB .* BIO.cab .+ K_CAR_V .* BIO.car .* (1 - BIO.f_zeax) .+ K_CAR_Z .* BIO.car .* BIO.f_zeax) ./ BIO._k_all ./ MESOPHYLL_N;

    _τ = (1 .- BIO._k_all) .* exp.(-BIO._k_all) .+ BIO._k_all .^ 2 .* expint.(BIO._k_all .+ eps(FT));

    talf = calctav.(FT(40), NR)

    return nothing
end;


end

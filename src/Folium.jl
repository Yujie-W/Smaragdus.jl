module Folium

using SpecialFunctions: expint
using UnPack: @unpack

# using external types and functions
using ..Emerald: Leaf, WaveLengthSet, calctav


# export the public functions
export leaf_spectra!


#######################################################################################################################################################################################################
#
# Changes made to this function
# General
#     2021-Oct-01: rename the function to leaf_spectra! as the function updates not only fluorescence but also reflectance, transmittance, and absorption spectra
# Bug fix
# To do
#     TODO: add another method to prescribe leaf spectra such as transmittance and reflectance from broadband method
#
#######################################################################################################################################################################################################
function leaf_spectra! end


#######################################################################################################################################################################################################
#
# Changes made to this method
# General
#     2020-Mar-30: account for carotenoid absorption as APAR as well as chlorophyll
#     2021-Aug-07: replace function `expint` with that from SpecialFunctions
# Bug fix
#     2021-Aug-06: If BIO.CBC and BIO.PRO are not zero, they are accounted for twice in BIO.LMA, thus the spectrum from LMA need to subtract the contribution from CBC and PRO
# To do
#     TODO: is f_sense a ratio to scale cab down, if so BIO.cab need to be multiplied by (1 - BIO.f_senes) and K_SENES need to be multiplied by BIO.cab
#
#######################################################################################################################################################################################################
leaf_spectra!(leaf::Leaf{FT}, wls::WaveLengthSet{FT}; APAR_car::Bool = true) where {FT<:AbstractFloat} = (
    BIO = leaf.BIO_PHYSICS;
    @unpack K_ANT, K_CAB, K_CAR_V, K_CAR_Z, K_CBC, K_H₂O, K_LMA, K_PRO, K_SENES, MESOPHYLL_N, NR = BIO;

    # calculate the average absorption feature and relative Cab and Car partitions
    BIO._k_all   .= (K_CAB   .* BIO.cab .+                      # chlorophyll absorption (TODO: downscale this with f_sens?)
                     K_CAR_V .* BIO.car .* (1 - BIO.f_zeax) .+  # violaxanthin carotenoid absorption
                     K_CAR_Z .* BIO.car .* BIO.f_zeax .+        # zeaxanthin carotenoid absorption
                     K_ANT   .* BIO.ant .+                      # anthocynanin absorption absorption
                     K_SENES .* BIO.f_senes .+                  # TODO: cannot be a fraction, needs to be a concentration
                     K_H₂O   .* BIO.l_H₂O .+                    # water absorption
                     K_CBC   .* BIO.CBC .+                      # carbon-based constituents absorption
                     K_PRO   .* BIO.PRO .+                      # protein absorption
                     K_LMA   .* (BIO.LMA - BIO.CBC - BIO.PRO)   # dry mass absorption (if some remained)
                    ) ./ MESOPHYLL_N;
    BIO.α_cab    .= (K_CAB .* BIO.cab) ./ BIO._k_all ./ MESOPHYLL_N;
    BIO.α_cabcar .= (K_CAB .* BIO.cab .+ K_CAR_V .* BIO.car .* (1 - BIO.f_zeax) .+ K_CAR_Z .* BIO.car .* BIO.f_zeax) ./ BIO._k_all ./ MESOPHYLL_N;

    _τ = (1 .- BIO._k_all) .* exp.(-BIO._k_all) .+ BIO._k_all .^ 2 .* expint.(BIO._k_all .+ eps(FT));

    # TODO: continue here later
    talf = calctav.(FT(40), NR);

    return nothing
);


end

# TODO: needs more work later
"""
    transmittance(θ_in::FT, n::FT) where {FT <:AbstractFloat}

Transmittance of S and P polarizations for a plane dielectric surface for light incident, given
- `θ_i` Incoming light angle
- `n` Relative index of refraction
"""
function transmittance(θ_i::FT, n::FT) where {FT <:AbstractFloat}
    _cos_θ_i = cosd(θ_i);
    _cos_θ_t = sqrt(1 - sind(θ_i)^2 / n^2);

    # reflection of s and p polarizations
    _r_s = ((_cos_θ_i - n * _cos_θ_t) / (_cos_θ_i + n * _cos_θ_t)) ^ 2;
    _r_p = ((n * _cos_θ_i - _cos_θ_t) / (n * _cos_θ_i + _cos_θ_t)) ^ 2;

    return 1 - _r_s, 1 - _r_p
end

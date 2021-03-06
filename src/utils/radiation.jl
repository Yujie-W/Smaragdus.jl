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








# TODO: needs more work later
###############################################################################
#
# Compute transmission of isotropic radiation
#
###############################################################################
"""
    calctav(α::FT, nr::FT) where {FT<:AbstractFloat}

Computes transmission of isotropic radiation across an interface between two
    dielectrics (Stern F., 1964; Allen W.A., 1973)). From calctav.m in
    PROSPECT-D
- `α` angle of incidence
- `nr` Index of refraction
"""
function calctav(α::FT, nr::FT) where {FT<:AbstractFloat}
    a   = ((nr+1) ^ 2) / 2;
    a3  = a  ^ 3;
    n2  = nr ^ 2;
    n4  = nr ^ 4;
    n6  = nr ^ 6;
    np  = n2 + 1;
    np2 = np ^ 2;
    np3 = np ^ 3;
    nm2 = (n2 - 1) ^2;
    k   = ((n2-1) ^ 2) / -4;
    k2  = k  ^ 2;
    sa2 = sind(α) ^ 2;

    _b1 = (α==90 ? 0 : sqrt( (sa2 - np/2)^2 + k ));
    _b2 = sa2 - np/2;
    b   = _b1 - _b2;
    b3  = b ^ 3;
    ts  = ( k2 / (6*b3) + k/b - b/2 ) - ( k2 / (6*a3) + k/a - a/2 );
    tp1 = -2 * n2 * (b-a) / (np2);
    tp2 = -2 * n2 * np * log(b/a) / (nm2);
    tp3 = n2 * (1/b - 1/a) / 2;
    tp4 = 16 * n4 * (n4+1) * log((2*np*b - nm2) / (2*np*a - nm2)) / (np3*nm2);
    tp5 = 16 * n6 * (1 / (2*np*b-nm2) - 1 / (2*np*a-nm2)) / (np3);
    tp  = tp1 + tp2 + tp3 + tp4 + tp5;
    tav = (ts + tp) / (2 * sa2);

    return tav
end

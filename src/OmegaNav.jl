module OmegaNav

import Geodesy: LLA, distance

export AbstractStation, OmegaStation, VLFStation,
    AbstractReceiver, OmegaReceiver,
    isAvailable, getSNR, SNRbyDistance, distance

function omega_time()
    time() % 30.0
end

abstract type AbstractStation end

struct OmegaStation <: AbstractStation
    designation
    location_name
    position::LLA
    unique_frequency
end

struct VLFStation <: AbstractStation
    designation
    location_name
    position::LLA
    fequencies::AbstractArray
    power
end

VLFStation(designation, location_name, position, fequency::Number, power) = VLFStation(designation, location_name, position, [fequency], power)

abstract type AbstractReceiver end

struct OmegaReceiver <: AbstractReceiver
    position::LLA
    clock_offset
    clock_drift
end

function isAvailable(t::OmegaStation, r::OmegaReceiver)
    d = distance(r, t)
    if d > 500e3 && d < (π*6378137 - 2000e3)
        snr = getSNR(t, r)
        if snr > -20
            return true
        end
    end
    false
end

function getSNR(t::OmegaStation, r::OmegaReceiver)
    d = distance(r, t)
    SNRbyDistance(d)
end

function SNRbyDistance(d)
    70 + (-20 - 70) / (deg2rad(110)*6378137 - 500e3) * d
end

const f = 1/298.257222100882711

function distance(λr, Λr, λt, Λt)
    Lr = atan((1-f) * tan(λr))
    Lt = atan((1-f) * tan(λt))
    cos_d = sin(Lr) * sin(Lt) + cos(Lr) * cos(Lt) * cos(Λr - Λt)
    d = acos(cos_d)
    ϵi = f/4 * ((sin(d) - d) * (sin(Lt) + sin(Lr))^2 / (1 + cos_d) - (sin(d) + d) * (sin(Lt) - sin(Lr))^2 / (1 - cos_d))
    d + ϵi
end

# distance(r::LLA, t::LLA) = distance(deg2rad(r.lat), deg2rad(r.lon), deg2rad(t.lat), deg2rad(t.lon))
distance(x) = distance(x[1,1], x[1,2], x[2,1], x[2,2]) * 6378137

distance(r::AbstractReceiver, t::AbstractStation) = distance(r.position, t.position)
distance(t::AbstractStation, r::AbstractReceiver) = distance(r, t)

function ∇distance(λr, Λr, λt, Λt)
    Lr = atan((1-f) * tan(λr))
    Lt = atan((1-f) * tan(λt))
    cos_d = sin(Lr) * sin(Lt) + cos(Lr) * cos(Lt) * cos(Λr - Λt)
    d = acos(cos_d)

    [
        (-cos(Lr) * sin(Lt) + sin(Lr) * cos(Lt) * cos(Λr - Λt)) / sin(d) cos(Lr) * cos(Lt) * sin(Λr - Λt) / sin(d)
        (-cos(Lt) * sin(Lr) + sin(Lt) * cos(Lr) * cos(Λt - Λr)) / sin(d) cos(Lt) * cos(Lr) * sin(Λt - Λr) / sin(d)
    ]
end

function ∇distance!(out, λr, Λr, λt, Λt)
    Lr = atan((1-f) * tan(λr))
    Lt = atan((1-f) * tan(λt))
    cos_d = sin(Lr) * sin(Lt) + cos(Lr) * cos(Lt) * cos(Λr - Λt)
    d = acos(cos_d)

    out[1,1] = (-cos(Lr) * sin(Lt) + sin(Lr) * cos(Lt) * cos(Λr - Λt)) / sin(d)
    out[1,2] = cos(Lr) * cos(Lt) * sin(Λr - Λt) / sin(d)
    out[2,1] = (-cos(Lt) * sin(Lr) + sin(Lt) * cos(Lr) * cos(Λt - Λr)) / sin(d)
    out[2,2] = cos(Lt) * cos(Lr) * sin(Λt - Λr) / sin(d)
end

const sequence = [
    10.2  13.6  11.33 12.1  12.1  11.05 12.1  12.1
    12.0  10.2  13.6  11.33 12.0  12.0  11.05 12.0
    11.8  11.8  10.2  13.6  11.33 11.8  11.8  11.05
    11.05 13.1  13.1  10.2  13.6  11.33 13.1  13.1
    12.3  11.05 12.3  12.3  10.2  13.6  11.33 12.3
    12.9  12.9  11.05 12.9  12.9  10.2  13.6  11.33
    11.33 13.0  13.0  11.05 13.0  13.0  10.2  13.6
    13.6  11.33 12.8  12.8  11.05 12.8  12.8  10.2
]

const timing = [0.9, 1.0, 1.1, 1.2, 1.1, 0.9, 1.2, 1.0]

const c0 = 299792458

const k0 = 1.0

const α = 1e-7

function amplitude(lon, lat, lon0, lat0, t, f, E0 = 1.0)
    ω = 2π * f
    λ = c0 / f
    k = k0 * 2π / λ

    d = distance(lat0, lon0, lat, lon) * 6378137

    E0 * exp(-α * d) * cos(ω * t - k * d)
end

function phase(lon, lat, lon0, lat0, t, f, E0 = 1.0)
    ω = 2π * f
    λ = c0 / f
    k = k0 * 2π / λ

    d = distance(lat0, lon0, lat, lon) * 6378137

    ω * t - k * d
end

function signal(lon, lat, lon0, lat0, t, f, E0 = 1.0)
    ω = 2π * f
    λ = c0 / f
    k = k0 * 2π / λ

    d = distance(lat0, lon0, lat, lon) * 6378137

    E0 * exp(-α * d) * exp(im * (ω * t - k * d))
end

function ionosphere_conductivity(z, h′ = 70e3, β = 0.3e-3)
    2.5e5 * exp(β * (z - h′))
end

end # module

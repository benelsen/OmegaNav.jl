module OmegaNav

export AbstractStation, OmegaStation, VLFStation,
    AbstractReceiver, OmegaReceiver,
    isAvailable, getSNR, SNRbyDistance

abstract type AbstractStation end

struct OmegaStation <: AbstractStation
    designation
    location_name
    position::LLA
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
    d = distance(r, t) * 6378137
    if d > 500e3 && d < (π*6378137 - 2000e3)
        snr = getSNR(t, r)
        if snr > -20
            return true
        end
    end
    false
end

function getSNR(t::OmegaStation, r::OmegaReceiver)
    d = distance(r, t) * 6378137
    SNRbyDistance(d)
end

function SNRbyDistance(d)
    70 + (-20 - 70) / (deg2rad(110)*6378137 - 500e3) * d
end

end # module

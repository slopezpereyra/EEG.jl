using EDF
using Plots


mutable struct EEG
    signals::Dict{String, Vector{<:AbstractFloat}}
    fs::Integer
    N::Integer
    epoch_length::Integer
    epoch_count::Union{<:AbstractFloat, <:Integer}
    staging::Vector{String}

    function EEG(file, fs, epoch_length, staging)

        eeg = EDF.read(file)
        S = Dict()
        
        for signal in eeg.signals
            S[signal.header.label] = EDF.decode(signal)
        end

        lengths = length.(values(S))
        N = lengths[1]

        if !all(y -> y == lengths[1], lengths)
            @warn "Not all signals are of equal length. Setting `N` field to length of the first signal."
        end

        epoch_count = N / (fs * epoch_length)
        
        new(S, fs, N, epoch_length, epoch_count, staging)
    end

end


function epoch(eeg::EEG, n::Integer)
    """An epoch e is a function ℕ ↦ ℕ₀² s.t. if e(n) = (x, y) if and only if 
       [x, y] is the space indexes of all values in an EEG signal in the nth epoch.

        Extending this, when e : ℕ² ↦ ℕ₀², we have e(n, m) = (x, y) if and only if 
        [x, y] is the space of all indexes corresponding to values in an EEG signal 
        in epochs n, n +1, …, m.
       """
    return [( (n - 1) * eeg.fs * eeg.epoch_length ) + 1, n * eeg.fs * eeg.epoch_length]
end

function epoch(eeg::EEG, n::Integer, m::Integer)
    """An epoch e is a function ℕ ↦ ℕ₀² s.t. if e(n) = (x, y) if and only if 
       [x, y] is the space indexes of all values in an EEG signal in the nth epoch.

        Extending this, when e : ℕ² ↦ ℕ₀², we have e(n, m) = (x, y) if and only if 
        [x, y] is the space of all indexes corresponding to values in an EEG signal 
        in epochs n, n +1, …, m.
       """
    if (n == m)
        return epoch(eeg, n)
    end
    if (n > m)
        throw(ArgumentError("The second epoch should be greater than the first."))
    end
    return [( (n - 1) * eeg.fs * eeg.epoch_length ) + 1, m * eeg.fs * eeg.epoch_length]
end

function epoch(eeg::EEG, n::Integer, channel::String)
    signal = eeg.signals[channel]
    bounds = epoch(eeg, n)
    signal[bounds[1]:bounds[2]]
end

function epoch(eeg::EEG, n::Integer, m ::Integer, channel::String)
    signal = eeg.signals[channel]
    bounds = epoch(eeg, n, m)
    signal[bounds[1]:bounds[2]]
end

function time(eeg::EEG, s::AbstractFloat, e::AbstractFloat)
    """Generates time domain of an EEG signal from second `s` to second `e`."""
    step = 1/eeg.fs
    return [ i for i in (s+step):step:e ]
end

function plot_eeg(eeg::EEG, channels::String, n::Integer, m::Integer)
    if n > m 
        throw(ArgumentError("m should be greater than n"))
    end
    t = time(eeg, (n-1) * eeg.epoch_length, m * eeg.epoch_length)
    signal = epoch(eeg, n, m, channels)
    plot(t, signal)
end

function plot_eeg(eeg::EEG, channels::Vector{String}, n::Integer, m::Integer)
    if n > m 
        throw(ArgumentError("m should be greater than n"))
    end
    t = time(eeg, (n-1) * eeg.epoch_length, m * eeg.epoch_length)
    signals = [epoch(eeg, n, m, chan) for chan in channels]
    plot(t, signals, layout = (length(signals), 1), legend=false)
    xlabel!("Time")
    ylabel!("")
end

function get_stage_indexes(eeg::EEG, stages::Vector)
    """This function maps a `stages` to the array of all indexes whose values in an EEG 
        signal pertain to a stage in `stages`."""


    stage_indexes = findall(x -> x in stages, eeg.staging)
    result = []

    for index in stage_indexes
        epoch_range = epoch(eeg, index)
        v = collect(epoch_range[1]:epoch_range[2])
        result = vcat(result, v)
    end

    return result
end

function get_stage(eeg::EEG, channel::String, stages::Vector)
    indexes = get_stage_indexes(eeg, stages)
    eeg.signals[channel][indexes]
end


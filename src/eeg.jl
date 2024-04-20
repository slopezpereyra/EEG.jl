using EDF
using Plots


mutable struct EEG
    """
    A mutable struct representing EEG data with associated metadata.

    # Fields
    - `signals::Dict{String, Vector{<:AbstractFloat}}`: A dictionary mapping signal labels (strings) to arrays of floating-point values.
    - `fs::Integer`: Sampling frequency (in Hz) of the EEG signals.
    - `N::Integer`: Length of the first EEG signal (number of samples).
    - `epoch_length::Integer`: Length of each epoch (in seconds).
    - `epoch_count::Union{<:AbstractFloat, <:Integer}`: Number of epochs in the EEG data. This can be a floating-point value if exact division is not possible.
    - `staging::Vector{String}`: A vector of stage labels corresponding to each epoch.

    # Constructor
    - `EEG(file, fs, epoch_length, staging)`: Constructs an `EEG` object from an EDF file (`file`) containing EEG data. The function reads the signals, computes the necessary metadata (`fs`, `N`, `epoch_count`), and initializes the `EEG` struct with the provided `staging` vector.

    ## Example
    ```julia
    eeg_data = EEG("data.edf", 256, 30, ["Wake", "1", "1", …, "REM", "REM", "2", "3", …])
    """
    signals::Dict{String,Vector{<:AbstractFloat}}
    fs::Integer
    N::Integer
    epoch_length::Integer
    epoch_count::Union{<:AbstractFloat,<:Integer}
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
    return [((n - 1) * eeg.fs * eeg.epoch_length) + 1, n * eeg.fs * eeg.epoch_length]
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
    return [((n - 1) * eeg.fs * eeg.epoch_length) + 1, m * eeg.fs * eeg.epoch_length]
end

function epoch(eeg::EEG, n::Integer, channel::String)
    """An epoch e is a function ℕ ↦ ℕ₀² s.t. if e(n) = (x, y) if and only if 
       [x, y] is the space indexes of all values in an EEG signal in the nth epoch.

        Extending this, when e : ℕ² ↦ ℕ₀², we have e(n, m) = (x, y) if and only if 
        [x, y] is the space of all indexes corresponding to values in an EEG signal 
        in epochs n, n +1, …, m.
       """
    signal = eeg.signals[channel]
    bounds = epoch(eeg, n)
    signal[bounds[1]:bounds[2]]
end

function epoch(eeg::EEG, n::Integer, m::Integer, channel::String)
    """An epoch e is a function ℕ ↦ ℕ₀² s.t. if e(n) = (x, y) if and only if 
       [x, y] is the space indexes of all values in an EEG signal in the nth epoch.

        Extending this, when e : ℕ² ↦ ℕ₀², we have e(n, m) = (x, y) if and only if 
        [x, y] is the space of all indexes corresponding to values in an EEG signal 
        in epochs n, n +1, …, m.
       """
    signal = eeg.signals[channel]
    bounds = epoch(eeg, n, m)
    signal[bounds[1]:bounds[2]]
end

function time(eeg::EEG, s::Union{AbstractFloat,Integer}, e::Union{AbstractFloat,Integer})
    """Generates time domain of an EEG signal from second `s` to second `e`."""
    step = 1 / eeg.fs
    return [i for i in (s+step):step:e]
end

function plot_eeg(eeg::EEG, channels::String, n::Integer, m::Integer)
    """
    Plots with the active backend all EEG channels in `channels` in the 
    range from the `n`th to the `m`th epoch.
    """
    if n > m
        throw(ArgumentError("m should be greater than n"))
    end
    t = time(eeg, (n - 1) * eeg.epoch_length, m * eeg.epoch_length)
    signal = epoch(eeg, n, m, channels)
    plot(t, signal)
end

function plot_eeg(eeg::EEG, channels::Vector{String}, n::Integer, m::Integer)
    if n > m
        throw(ArgumentError("m should be greater than n"))
    end
    t = time(eeg, (n - 1) * eeg.epoch_length, m * eeg.epoch_length)
    signals = [epoch(eeg, n, m, chan) for chan in channels]
    plot(t, signals, layout=(length(signals), 1), legend=false)
    xlabel!("Time")
    ylabel!("")
end

function plot_eeg_overlay(eeg::EEG, channels::Vector{String}, n::Integer, m::Integer)
    if n > m
        throw(ArgumentError("m should be greater than n"))
    end
    t = time(eeg, (n - 1) * eeg.epoch_length, m * eeg.epoch_length)
    signals = [epoch(eeg, n, m, chan) for chan in channels]
    Plots.plot(t, signals, legend=false)
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

function filter!(eeg::EEG, channel::String, digfilter, cut_off)
    eeg.signals[channel] = filt(digitalfilter(digfilter(cut_off, fs=eeg.fs), Butterworth(4)), eeg.signals[channel])
end

function filter!(eeg::EEG, channels::Vector{<:String}, digfilter, cut_off)

    for chan in channels
        eeg.signals[chan] = filt(digitalfilter(digfilter(cut_off, eeg.fs), Butterworth(4)), eeg.signals[chan])
    end
end

function filter!(eeg::EEG, digfilter, cut_off)
    for chan in keys(eeg.signals)
        eeg.signals[chan] = filt(digitalfilter(digfilter(cut_off, eeg.fs), Butterworth(4)), eeg.signals[chan])
    end
end

function nrem(eeg::EEG, n::Integer, m::Integer)
    """
    nrem(eeg::EEG, n::Integer, m::Integer)

    This function finds all NREM periods in a stage vector. It returns a list [ [a₁, b₁], …, [aₖ, bₖ] ] s.t. 
    each of the k NREM periods of the EEG range between epochs [aᵢ, bᵢ], i ∈ ℕ.

    Given a sequence s₁, …, sₖ of sleep stages, s.t. sᵢ is the stage of the ith epoch in an EEG,
    we define the alphabet Σ = {s : s ∈ {s₁, …, sₙ})}. We assume Σ = {1, 2, 3, 4, 5, 6}, with 5 denoting 
    REM and 6 denoting wakefulness. If α = s₁… sₙ ∈ Σ*, we assume α can be decomposed into the form 

                                α = ψ₁β₁ …  ψₖβₖ ψₖ₊₁ 

    where ψᵢ is an arbitrary word and βᵢ is a word of the form 

                                φ(5ᵐ5* + 6),   φ ∈ { w ∈ {2, 3, 4}* : |w| ≥ n }

    Here, n is the minimum length imposed upon combinations of stages 2, 3, and 4, in order to account for a NREM period,
    while m is the minimum length imposed upon REM periods terminating NREM periods.

    Thus, βᵢ can be thought of as a NREM period, and the set φ(5ᵐ5* + 6) as the universe of NREM periods.
    """
    φ = Regex("[234]{$n,}(?:6|5{$(m+1),})")
    α = join(eeg.staging, "")
    matches = eachmatch(φ, α)
    return [ [match.offset, match.offset + length(match.match)] for match in matches ]
end














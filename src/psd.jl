using FFTW
using DSP

"""
    Splits a vectir `v` into segments of length `L` with an overlap `overlap_frac` expressed
    as a fraction of L. 

    Parameters
    ---------------------

        `v::Vector{T}` : The vector to be split
        `L::Int` : Length of each segment
        `overlap_frac::Union{Float64,Int})` : A number in [0, 1) which expresses the overlap as a fraction of L.
                                              If `overlap_frac` = x, then each segment will share L * x elements
                                              with its neighboring segments.
"""
function overlaps(v::Vector{T}, L::Int, overlap_frac::Union{Float64,Int}) where {T}
    if L > length(v)
        throw(ArgumentError("Segment length L must be less than or equal to the length of the vector."))
    end

    if overlap_frac < 0 || overlap_frac >= 1.0
        throw(ArgumentError("Overlap fraction must be in the range [0, 1)."))
    end

    D = L * overlap_frac
    M = Int(ceil((length(v) - L) / (L - D)))

    segments = Vector{Vector{T}}(undef, M)
    step = Int(floor((1 - overlap_frac) * L))  # Calculate step size

    for i in 1:M
        start_idx = 1 + (i - 1) * step
        end_idx = start_idx + L - 1

        # Ensure the last segment does not exceed the length of the vector
        if end_idx > length(v)
            break
        end

        segments[i] = v[start_idx:end_idx]
    end

    return segments
end



mutable struct AmplitudeSpectrum 

    freq::Vector{<:AbstractFloat}
    spectrum::Vector{<:AbstractFloat}
    formula::String
    
    function AmplitudeSpectrum(x::Vector{<:AbstractFloat}, sampling_rate::Integer, pad::Integer)
        N = length(x)
        hann = hanning(N) # Hanning window
        x = x .* hann
        if pad > 0 
            x = zero_pad(x, pad)
        end
        ft = abs.(fft(x))
        ft = ft[1:(div(N, 2)+1)] # Make one sided
        freq = [i for i in 0:(length(ft)-1)] .* sampling_rate / N
        normalization = 2 / (sum(hann .^ 2))
        spectrum = ft * normalization
        new(freq, spectrum,  "2|H(f)| / ∑ wᵢ²  with  w₁, …, wₗ a Hanning window")
    end

end

    """Structure for PSD estimations. Estimations are by default 
    one sided, with frequencies ranging from [0, fₛ/2].

    Fields
    ----------

       freq::Vector{<:AbstractFloat}: 
            Frequency range of the spectrum
       spectrum::Vector{<:AbstractFloat}
            Estimated spectral density
       method::String
            Estimation method used 
       formula::String
            A string representation of the formula used for the estimation.

    Constructor functions
    ---------------------
        (1) function PSD(x::Vector, sampling_rate::Integer) :
                
                Computes a direct PSD over a signal `x` with a given `sampling_rate`.

        (2) function PSD(x::Vector, fs::Int, L::Int, overlap::Union{ <:AbstractFloat, Integer }, normalization::Union{ <:AbstractFloat, Integer } = 1)
                
                Splits the signal `x` into segments of length L with an `overlap` in [0, 1). The overlap is understood to be a fraction of the segment length. 
                PSD is estimated within and averaged across all segments. If `overlap` is zero, this results in Barlett's method. If `overlap` is greater 
                than zero, this results in Welch's method.
                If `pad` is zero no zero-padding is done. If `pad` is greater than zero, each segment is zero-padded to a 
                length of `pad`. 
    """
struct PSD
    freq::Vector{<:AbstractFloat}
    spectrum::Vector{<:AbstractFloat}
    method::String
    formula::String

    function PSD(x::Vector{<:AbstractFloat}, sampling_rate::Integer, pad::Integer=0)
        N = length(x)
        hann = hanning(N) # Hanning window
        x = x .* hann
        if pad > 0 
            x = zero_pad(x, pad)
        end
        ft = abs2.(fft(x))
        ft = ft[1:(div(N, 2)+1)] # Make one sided
        freq = [i for i in 0:(length(ft)-1)] .* sampling_rate / N
        normalization = 2 / (sum(hann .^ 2))
        spectrum = ft * normalization
        new(freq, spectrum, "Direct (no segmentation)", "2|H(f)|² / ∑ wᵢ²  with  w₁, …, wₗ a Hanning window")
    end

    function PSD(x::Vector{<:AbstractFloat}, fs::Int, L::Int, overlap::Union{<:AbstractFloat,Integer}, normalization::Union{<:AbstractFloat,Integer}=1,
                 pad::Integer=0)


        if (L >= length(x))
            throw(ArgumentError("The length `L` of each segment should be inferior 
            to the length of the vector."))
        end
        if (overlap < 0 || overlap >= 1)
            throw(ArgumentError("The `overlap` should range in (0, 1)"))
        end

        method = overlap > 0 ? "Welch's method" : "Barlett's method"
        formula = "1/(M * normalization) ∑ ᵢᴹ [ 2|Hᵢ(f)|² / ∑  wᵢ² ]  where w₁, …, wₗ a Hanning window, M the number of segments, and Hᵢ(f) the FFT of the ith segment of the signal. "

        function H(signal::Vector, window::Vector, seg_length::Integer)
            signal = signal .* window
            if pad > 0 
                signal = zero_pad(signal, pad)
            end
            ft = abs2.(fft(signal))
            ft = ft[1:(div(seg_length, 2)+1)] # One sided
            return 2 .* ft ./ sum(window .^ 2)
        end

        hann = hanning(L)
        segs = overlaps(x, L, overlap)
        M = length(segs)
        w = sum(map(x -> H(x, hann, L), segs)) ./ (M * normalization)   # Hans uses denominator 2 * M * length(segs[1])
        freq = [i for i in 0:(div(L, 2))] .* fs / L
        new(freq, w, method, formula)
    end
end

struct Spectrogram

    time::Vector
    freq::Vector{<:AbstractFloat}
    spectrums::Matrix{<:AbstractFloat}
    segment_length::Integer

    function Spectrogram(signal::Vector{<:AbstractFloat}, fs::Integer, segment_length::Integer, overlap::AbstractFloat)
        segs = overlaps(signal, fs*segment_length, overlap)
        psds = map(x -> PSD(x, fs), segs)
        freq = psds[1].freq 
        spectrums = map(x -> pow2db.( x.spectrum ), psds) 
        
        spectrogram_data = zeros(length(spectrums), length(freq))
        for i in 1:length(spectrums)
            spectrogram_data[i, :] = spectrums[i]
        end

        new(1:length(spectrums), freq, spectrogram_data, segment_length)
    end
end

function plot_spectrogram(spec::Spectrogram, freq_lim::AbstractFloat = 30.0, type::Int=1, color = :hnipy_spectral)

    if type == 1
        return(heatmap(spec.time, spec.freq, spec.spectrums', ylims=(0, freq_lim), color = color))
    end

    if type == 2
        return(surface(spec.time, spec.freq, spec.spectrums', ylims=(0, freq_lim), 
        xlabel="Time", ylabel="Frequency (Hz)", zlabel="PSD (dB)", color = color))
    end 
    throw(ArgumentError("The plot `type` argument must be either 1 (for heatmap) or 2 (for a surface plot)."))
end

function next_power_of_two(n::Int)
    p = 1
    while p < n
        p <<= 1
    end
    return p
end

function zero_pad(v::Vector{T}, desired_length::Integer) where {T}
    current_length = length(v)
    if current_length == desired_length
        return v
    end 
    if desired_length < current_length
        throw(ArgumentError("Cannot zero-pad to a length inferior to the vector's length"))
    end
    
    padded_vector = zeros(T, desired_length)
    padded_vector[1:current_length] = v
    return padded_vector
    
end


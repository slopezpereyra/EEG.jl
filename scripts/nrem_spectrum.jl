using CSV
using DataFrames
using Statistics
using DelimitedFiles
using Base.Filesystem

include("../src/eeg.jl")
include("../src/nrem.jl")
include("../src/psd.jl")
include("../src/helpers.jl")

#plots = []

const STAGE_FILES = readdir("data/Staging")

if ARGS[1] == "baseline"
    const ANOM_FILES = readdir("bas")
    const ANOM_FOLDER = "bas"
    const EDF_FOLDER = "data/Baseline/"
elseif ARGS[1] == "disruption" 
    const ANOM_FILES = readdir("disr")
    const ANOM_FOLDER = "disr"
    const EDF_FOLDER = "data/Disruption/"
else 
    print("Invalid argument: Should be `baseline` or `disruption`")
end
const EDF_FILES = readdir(EDF_FOLDER)
const QUANT_THRESH = .95
const OUTPUT_FILE_NAME = ARGS[2] * ".csv"

function nrem_spectrum(eeg, clean_epochs, nrems)
    mean_delta_powers = []
    specs = []
    i = 1
    for nrem in nrems
        print(i, "\n")
        if length(nrem) == 0
            i = i + 1
            continue
        end
        if i > 4
            break
        end
        nrem_epochs = clean_epochs[nrem]
        #print(nrem_epochs)
        spec = Spectrogram(nrem_epochs, eeg.fs, eeg.fs*5 - 1, 0.5)
        push!(specs, spec)
        print("Spec computed\n")
        delta_band = freq_band(spec, 0.5, 3.9)
        print("Delta banda computed\n")
        mean_delta_through_freq = mean(delta_band, dims=2) 
        push!(mean_delta_powers, mean_delta_through_freq)
 #       print("Attempting to plot...\n")
 #       p = plot_spectrogram(spec, 30.0, 1, :inferno)
 #       print("Plotted a spectrogram\n")
 #       push!(plots, p)
        i = i+1
    end
    print("Attempting to output the file\n")
 #   pfinal = plot(plots...)
#    savefig(pfinal, eeg.id * "NREM_spectrums" * ".png")
    means = mean.(mean_delta_powers)
    res = hcat(means, collect(1:length(means)), [eeg.id for i in collect(1:length(means))])
    writedlm("results/" * ANOM_FOLDER * "/" * eeg.id * OUTPUT_FILE_NAME, res, ',')
end

function spectrum(edf_file)
    print("Reading " * edf_file, "\n")
    eeg = EEG(edf_file, 30, [], extract_identifier(edf_file))
    print("EEG read \n")
    set_staging(eeg, STAGE_FILES)

    c3key = filter(x -> occursin("C3", x), collect(keys(eeg.signals)))[1]
    eeg.fs = eeg.sampling_rates[c3key]
    unwanted_keys = [key for key in keys(eeg.signals) if !occursin("C3", key)]
    for key in unwanted_keys
        delete!(eeg.signals, key)
    end
    print("Proceeding to filter... \n")

    # Filtering 
    eeg.signals[c3key] = filt(digitalfilter(Lowpass(70, fs=eeg.fs), Butterworth(4)), eeg.signals[c3key] )
    eeg.signals[c3key] = filt(digitalfilter(Highpass(0.3, fs=eeg.fs), Butterworth(4)), eeg.signals[c3key] )
    eeg.signals[c3key] = filt(digitalfilter(Bandstop(50, 70, fs=eeg.fs), Butterworth(4)), eeg.signals[c3key] )

    print("Filtering completed. Detecting NREM periods...\n")
    nrems = nrem(eeg)
    anoms = flag_nrem_artifacts(eeg, ANOM_FOLDER, nrems, ANOM_FILES)
    print("Cleaning according to flagged artifacts\n")
    epochs = artifact_reject(eeg, anoms)
    epochs = filter(x -> !isempty(x), epochs)
    nrem_spectrum(eeg, epochs, nrems)
end

function main()
    for edf in EDF_FILES
        print("At " * edf)
        try
           print("\n\n")
           spectrum(EDF_FOLDER * edf)
        catch 
            print("\n WARNING : FAILED ")
            continue
        end
    end
end

main()








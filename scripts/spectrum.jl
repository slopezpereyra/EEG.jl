using CSV
using DataFrames
using Statistics
using DelimitedFiles
using Base.Filesystem

include("../src/eeg.jl")
include("../src/nrem.jl")
include("../src/psd.jl")
include("../src/helpers.jl")


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

function flag_artifacts(eeg)
    anom_files = read_anoms(eeg, ANOM_FILES)
    print("Succesfuly found anom files: ", anom_files,  "\n")
    canoms = CSV.read(ANOM_FOLDER * "/" * anom_files[1], DataFrame)
    panoms = CSV.read(ANOM_FOLDER * "/" * anom_files[2], DataFrame)
    rename!(canoms, [:group_id, :start, :end, :variate, :start_lage, :end_lage, :strength, :stat])
    panoms = filter(:strength => x -> x >= quantile(panoms.strength, QUANT_THRESH), panoms)
    canoms = filter(:strength => x -> x > quantile(canoms.strength, QUANT_THRESH), canoms)

    e_indexes = anom_to_index(canoms.group_id, canoms.end, eeg.fs)
    s_indexes = anom_to_index(canoms.group_id, canoms.start, eeg.fs)
    points = anom_to_index(panoms.group_id, panoms.location, eeg.fs)
    s = hcat(index_to_epoch.(s_indexes, eeg.fs) .+ 1, index_to_subepoch.(s_indexes, eeg.fs) .+ 1)
    e = hcat(index_to_epoch.(e_indexes, eeg.fs) .+ 1, index_to_subepoch.(e_indexes, eeg.fs) .+ 1)
    p = hcat(index_to_epoch.(points, eeg.fs) .+ 1, index_to_subepoch.(points, eeg.fs) .+ 1)
    U = vcat(s, e, p)
    U = unique(U; dims=1)
    return U
end

function spectrum(edf_file)
    print("Reading " * edf_file, "\n")
    eeg = EEG(edf_file, 30, [], extract_identifier(edf_file))
    print("EEG read \n")

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

    print("Filtering completed. Flaging artifacts...\n")
    anoms = flag_artifacts(eeg) 
    print("Cleaning according to flagged artifacts\n")
    epochs = artifact_reject(eeg, anoms)
    epochs = filter(x -> !isempty(x), epochs)
    print("Computing PSDs...\n")
    psds = map(x -> PSD(x, eeg.fs, eeg.fs*5-1, 0.5), epochs)
    print("Extracting delta band of each PSD...\n")
    deltas = map(x -> freq_band(x, 0.5, 3.9), psds)
    print("Computing mean power in each band...\n")
    means = map(mean, deltas)
    print("Computing final mean...\n")
    delta = mean(means)
    return ([eeg.id,  delta ])
end

function main()
    ids = []
    deltas = []
    for edf in EDF_FILES
        try
            print("\n\n")
            result = spectrum(EDF_FOLDER * edf)
            push!(ids, result[1])
            push!(deltas, result[2])
        catch 
            print("\n WARNING : FAILED ")
            continue
        end
    end
    writedlm(OUTPUT_FILE_NAME, hcat(ids, deltas), ',')
end

main()








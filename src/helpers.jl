using DataFrames
using CSV
using Statistics

include("psd.jl")
include("eeg.jl")

function extract_identifier(filename::AbstractString)
    # Define the regular expression pattern to match the filename format
    pattern = r"SWIP(\d+)(baseline|disruption)\.edf$"
    
    # Use the match function to find the pattern in the filename
    match_result = match(pattern, filename)
    
    if match_result === nothing
        return nothing  # Pattern not found in the filename
    end
    
    # Extract the identifier and type from the matched result
    identifier = match_result.captures[1]
    disruption_type = match_result.captures[2]
    
    # Form the desired output string
    result = identifier * disruption_type
    
    return result
end

function anom_to_index(id, i, fs)
    x = (fs* 30 * 5) .* (id .- 1)
    x .+ i
end

function index_to_epoch(index, fs)
    ep_length = fs * 30
    index ÷ ep_length
end

function index_to_subepoch(index, fs)
    sub_length = fs * 5
    ((index ÷ sub_length) % 6)
end


function read_anoms(eeg, anom_files)
    i = findall(x -> occursin(eeg.id[1:3], x), anom_files)
    print(":::::  ", anom_files[i], " ::::::::\n")
    anom_files[i]
end


function flag_nrem_artifacts(eeg, folder, nrems, anom_files)
    anom_files = read_anoms(eeg, anom_files)
    canoms = CSV.read(folder * "/" * anom_files[1], DataFrame)
    panoms = CSV.read(folder * "/" * anom_files[2], DataFrame)
    rename!(canoms, [:group_id, :start, :end, :variate, :start_lage, :end_lage, :strength, :stat])
    panoms = filter(:strength => x -> x >= quantile(panoms.strength, 0.1), panoms)
    canoms = filter(:strength => x -> x > quantile(canoms.strength, 0.1), canoms)
    e_indexes = anom_to_index(canoms.group_id, canoms.end, eeg.fs)
    s_indexes = anom_to_index(canoms.group_id, canoms.start, eeg.fs)
    points = anom_to_index(panoms.group_id, panoms.location, eeg.fs)
    s = hcat(index_to_epoch.(s_indexes, eeg.fs) .+ 1, index_to_subepoch.(s_indexes, eeg.fs) .+ 1)
    e = hcat(index_to_epoch.(e_indexes, eeg.fs) .+ 1, index_to_subepoch.(e_indexes, eeg.fs) .+ 1)
    p = hcat(index_to_epoch.(points, eeg.fs) .+ 1, index_to_subepoch.(points, eeg.fs) .+ 1)
    U = vcat(s, e, p)
    U = unique(U; dims=1)
    nrem_epochs = vcat(nrems...)
    nrem_anoms = findall(x -> x in nrem_epochs, U[:, 1])
    U = U[nrem_anoms, :]
    return U
end

function set_staging(eeg, stage_files)
    i = findall(x -> occursin(eeg.id[2:length(eeg.id)], lowercase(x)), stage_files)
    print(stage_files[i], "\n")
    staging_file = stage_files[i][1]
    print(staging_file)
    staging = CSV.read("data/Staging/" * "SleepStaging_70Baseline.csv", DataFrame)
    score = staging[!, 3]
    score = replace(score, "NS" => "?", "REM" => "5", "WK" => "6", "N1" => "1", "N2" => "2", "N3" => "3", "N4" => "4")
    eeg.staging = score
end

function artifact_reject(eeg, anom_matrix)
    c3key = filter(x -> occursin("C3", x), collect(keys(eeg.signals)))[1]
    epochs = overlaps(eeg.signals[c3key], eeg.fs*30, 0)
    windows = map(x -> overlaps(x, eeg.fs * 5, 1/(eeg.fs*5)), epochs)
    for epoch in unique(anom_matrix[:, 1])
        if epoch > length(windows)
            print("Anomaly was found on final epoch, which does not belong to the segmentation\n")
            continue
        end
        subeps_indexes = findall(x -> x == epoch, anom_matrix[:, 1])
        subeps = anom_matrix[:, 2][subeps_indexes]
        #print(epoch, " ~~~> ", subeps_indexes, "~~~~>",  subeps, "\n")
        deleteat!(windows[epoch], sort(subeps))
    end
    clean = map(x -> vcat(x...), windows)
    return clean
end

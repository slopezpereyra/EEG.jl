
const ANOM_FILES_BAS = readdir("bas")
const ANOM_FILES_DISR = readdir("disr")

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

function read_anoms(eeg)
    i = findall(x -> occursin(eeg.id[1:3], x), ANOM_FILES_BAS)
    ANOM_FILES_BAS[i]
end

function flag_artifacts(eeg)
    anom_files = read_anoms(eeg)
    canoms = CSV.read("bas/" * anom_files[1], DataFrame)
    panoms = CSV.read("bas/" * anom_files[2], DataFrame)
    rename!(canoms, [:group_id, :start, :end, :variate, :start_lage, :end_lage, :strength, :stat])
    panoms = filter(:strength => x -> x >= quantile(panoms.strength, .90), panoms)
    canoms = filter(:strength => x -> x > quantile(canoms.strength, .90), canoms)

    e_indexes = anom_to_index(canoms.group_id, canoms.end, eeg.fs)
    s_indexes = anom_to_index(canoms.group_id, canoms.start, eeg.fs)
    points = anom_to_index(panoms.group_id, panoms.location, eeg.fs)
    s = hcat(index_to_epoch.(s_indexes, eeg.fs) .+ 1, index_to_subepoch.(s_indexes, eeg.fs) .+ 1)
    e = hcat(index_to_epoch.(e_indexes, eeg.fs) .+ 1, index_to_subepoch.(e_indexes, eeg.fs) .+ 1)
    p = hcat(index_to_epoch.(points, eeg.fs) .+ 1, index_to_subepoch.(points, eeg.fs) .+ 1)
    U = vcat(s, e, p)
    U = unique(U; dims=1)
end

function compute_psd_without_artifacts(eeg, anom_matrix)
    c3key = filter(x -> occursin("C3", x), collect(keys(eeg.signals)))[1]
    epochs = overlaps(eeg.signals[c3key], eeg.fs*30, 0)
    windows = map(x -> overlaps(x, eeg.fs * 5, 1/(eeg.fs*5)), epochs)
    for epoch in unique(anom_matrix[:, 1])
        subeps_indexes = findall(x -> x == epoch, anom_matrix[:, 1])
        subeps = anom_matrix[:, 2][subeps_indexes]
        #print(epoch, " ~~~> ", subeps_indexes, "~~~~>",  subeps, "\n")
        deleteat!(windows[epoch], sort(subeps))
    end
    cleaned_epochs = map(x -> vcat(x...), windows)
    spec = Spectrogram(cleaned_epochs, eeg.fs, eeg.fs*5, 0.5)
    return spec
end



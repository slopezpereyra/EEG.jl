include("eeg.jl")

function nrem(vector, n::Integer = 30, m::Integer = 10)
    φ = Regex("[234]{$n,}(?:[5]{$(m+1),})") # Used to be φ = Regex("[1234]{$n,}(?:6|[51]{$(m+1),})")
    α = join(vector, "")
    matches = eachmatch(φ, α)
    result = [[match.offset + 1, match.offset + findfirst(Regex("[56]"), match.match)[1] - 1] for match in matches]
    return result
end

function nrem(eeg::EEG, n::Integer = 30, m::Integer = 10)
    """
    nrem(eeg::EEG, n::Integer, m::Integer)

    This function finds all NREM periods in a stage vector. It finds the kth underlying 
    NREM periods in a the staging vector and returns a list [ M₁, …, Mₖ ], with
    Mᵢ a n x 2 matrix whose first column are the epochs of the ith NREM period and whose second 
    column are the sleep stages of the corresponding epochs.

    *Further detail*. Given a sequence s₁, …, sₖ of sleep stages, s.t. sᵢ is the stage of the ith epoch in an EEG,
    we define the alphabet Σ = {s : s ∈ {s₁, …, sₙ})}. We assume Σ = {2, 3, 4, 5}, with 5 denoting 
    REM . If α = s₁… sₙ ∈ Σ*, we assume α can be decomposed into the form 

                                α = ψ₁β₁ …  ψₖβₖ ψₖ₊₁ 

    where ψᵢ is an arbitrary word and βᵢ is a word of the form 

                                φ(5ᵐ5*),   φ ∈ { w ∈ {2, 3, 4}* : |w| ≥ n }

    Here, n is the minimum length imposed upon combinations of stages 2, 3, and 4, in order to account for a NREM period,
    while m is the minimum length imposed upon REM periods terminating NREM periods.

    Thus, βᵢ can be thought of as a NREM period, and the set φ(5ᵐ5*) as the universe of NREM periods.
    """
    
    stages = hcat(1:length(eeg.staging), eeg.staging)
    X = stages[stages[:, 2] .!= string('1'), :] # Remove Stage 1 epochs.
    X = X[X[:, 2] .!= string('M'), :] # Remove M  epochs.
    X = X[X[:, 2] .!= string('6'), :] # Remove M  epochs.
    φ = Regex("[234]{$n,}(?:[5]{$(m+1),})")
    α = join(X[:, 2], "")
    print(α)
    matches = eachmatch(φ, α)

    nrem_periods = []
    for match in matches 
        indexes = (match.offset):(match.offset + findfirst(Regex("[56]"), match.match)[1] - 2)
        period = X[indexes, :]
        push!(nrem_periods, period)

    end
    return nrem_periods
end


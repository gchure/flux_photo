"""
    steady_state_growth_rate(φᵣ, ν⁺; γ⁺, Kₘ, φₒ)

Computes the steady-state growth rate  under a simple model of a self-replicating
system.

# Input
- `φᵣ::Union{Float64, Vector{Float64}}`: The prescribed allocation towards ribosomes.
- `ν⁺::Union{Float64, Vector{Float64}}` : The maximal meatabolic rate of the
   environment, provided as either a single value or a row vector. Should 
   be in units of inverse time. 
- `γ⁺::Float64=9.65`: The maximal translation rate in units of ribosomal mass units 
    per time. Default is 9.65, corresponding to a translation speed of 20 amino 
    acids per second with an *E. coli* ribosome protein mass. 
- `Kₘ`::Float64=0.03`: The Michaelis-Menten constant for binding of precursors 
   to an actively translating ribosome, in units of relative proteome abundance. 
   Default value is 0.03.
- `φₒ::Float64=0.55`: The allocation towards "other" proteins. Default value is 
   0.55.

# Output
- `λ::Union{Float64, Vector{Float64}}`: The steady-state precursor concentration 
   in units of inverse time.

"""
function steady_state_growth_rate(
    φᵣ::Union{Float64, Vector{Float64}},
    ν⁺::Union{Float64, Vector{Float64}};
    γ⁺::Float64=9.65,
    Kₘ::Float64=0.03,
    φₒ::Float64=0.55)

    # If a single value is provided, convert to a vector 
    if typeof(ν⁺) !== Vector{Float64}
        ν⁺ = [ν⁺]
    end

    # Compute maximal metabolic and translational outputs
    N = ν⁺ .* (1 .- φₒ .- φᵣ) 
    Γ = γ⁺ .* φᵣ

    # Compute the root and return
    λ = (N .+ Γ .- .√((N .+ Γ).^2 .- 4 .* (1 .- Kₘ) .* N .* Γ))./(2 * (1 .- Kₘ))
    if length(λ) == 1
        λ = λ[1]
    end
    return λ
end
export steady_state_growth_rate

"""
    steady_state_precursors(φᵣ, ν⁺; γ⁺, Kₘ, φₒ)

Computes the steady-state precursor concentration under a simple model of a 
self-replicating system 

# Input
- `φᵣ::Union{Float64, Vector{Float64}}`: The prescribed allocation towards 
   ribosomes provided as either a single value or a row vector. 
- `ν⁺::Union{Float64, Vector{Float64}}` : The maximal meatabolic rate of the
   environment, provided as either a single value or a row vector. Should 
   be in units of inverse time. 
- `γ⁺::Float64=9.65`: The maximal translation rate in units of ribosomal mass units 
    per time. Default is 9.65, corresponding to a translation speed of 20 amino 
    acids per second with an *E. coli* ribosome protein mass. 
- `Kₘ`::Float64=0.03`: The Michaelis-Menten constant for binding of precursors 
   to an actively translating ribosome, in units of relative proteome abundance. 
   Default value is 0.03.
- `φₒ::Float64=0.55`: The allocation towards "other" proteins. Default value is 
   0.55.

# Output
- `cₚ::Union{Float64, Vector{Float64}}`: The steady-state precursor concentration 
   in units of proteome abundance.
"""
function steady_state_precursors(
    φᵣ::Union{Float64, Vector{Float64}},
    ν⁺::Union{Float64, Vector{Float64}};
    γ⁺::Float64=9.65,
    Kₘ::Float64=0.03,
    φₒ::Float64=0.55)
       
    # If a single value is provided, convert to a vector 
    if typeof(ν⁺) !== Vector{Float64}
        ν⁺ = [ν⁺]
    end
    λ = steady_state_growth_rate(φᵣ, ν⁺; γ⁺, Kₘ, φₒ)
    cₚ = (ν⁺ .* (1 - φₒ - φᵣ) ./ λ) - 1
    if length(cₚ) == 1
        cₚ = cₚ[1]
    end
    return cₚ
end
export steady_state_precursors

"""

    steady_state_translation_rate(φᵣ, ν⁺; γ⁺, Kₘ, φₒ)

Computes the steady-state translation rate under a simple model of a 
self-replicating system 

# Input
- `φᵣ::Union{Float64, Vector{Float64}}`: The prescribed allocation towards 
   ribosomes provided as either a single value or a row vector. 
- `ν⁺::Union{Float64, Vector{Float64}}` : The maximal meatabolic rate of the
   environment, provided as either a single value or a row vector. Should 
   be in units of inverse time. 
- `γ⁺::Float64=9.65`: The maximal translation rate in units of ribosomal mass units 
    per time. Default is 9.65, corresponding to a translation speed of 20 amino 
    acids per second with an *E. coli* ribosome protein mass. 
- `Kₘ`::Float64=0.03`: The Michaelis-Menten constant for binding of precursors 
   to an actively translating ribosome, in units of relative proteome abundance. 
   Default value is 0.03.
- `φₒ::Float64=0.55`: The allocation towards "other" proteins. Default value is 
   0.55.

# Output
- `γ::Union{Float64, Vector{Float64}}`: The steady-state translation rate in 
   in units of ribosomal mass units per unit time.
"""
function steady_state_translation_rate(
    φᵣ::Union{Float64, Vector{Float64}},
    ν⁺::Union{Float64, Vector{Float64}};
    γ⁺::Float64=9.65,
    Kₘ::Float64=0.03,
    φₒ::Float64=0.55) 

    # If a single value is provided, convert to a vector 
    if typeof(ν⁺) !== Vector{Float64}
        ν⁺ = [ν⁺]
    end 
    cₚ = steady_state_precursors(φᵣ, ν⁺; γ⁺, Kₘ, φₒ)
    γ = γ⁺ .* cₚ ./ (cₚ + Kₘ)
    if length(γ) == 1
        γ = γ[1]
    end
    return γ
end
export steady_state_translation_rate

"""
    optimal_ribosomal_allocation(ν⁺; γ⁺, Kₘ, φₒ)

Computes the optimal allocation towards ribosomes which contextually maximizes 
the growth rate for a simple self-replicating system.

# Input
- `ν⁺::Union{Float64, Vector{Float64}}` : The maximal meatabolic rate of the
   environment, provided as either a single value or a row vector. Should 
   be in units of inverse time. 
- `γ⁺::Float64=9.65`: The maximal translation rate in units of ribosomal mass units 
    per time. Default is 9.65, corresponding to a translation speed of 20 amino 
    acids per second with an *E. coli* ribosome protein mass. 
- `Kₘ`::Float64=0.03`: The Michaelis-Menten constant for binding of precursors 
   to an actively translating ribosome, in units of relative proteome abundance. 
   Default value is 0.03.
- `φₒ::Float64=0.55`: The allocation towards "other" proteins. Default value is 
   0.55.

# Output
- `φᵣ::Union{Float64, Vector{Float64}}`: The optimal allocation towards ribosomes.
"""
function optimal_ribosomal_allocation(
    ν⁺::Union{Int64, Float64, Vector{Float64}};
    γ⁺::Float64=9.65,
    Kₘ::Float64=0.03,
    φₒ::Float64=0.55)

    # If a single value is provided, convert to a vector 
    if typeof(ν⁺) !== Vector{Float64}
        ν⁺ = [ν⁺]
    end 

    # Compute everything piecewise
    numer = γ⁺ .* ν⁺ .* (1 - 2 * Kₘ) .+ ν⁺.^2 .+ .√(Kₘ .* γ⁺ .* ν⁺) .* (γ⁺ .- ν⁺)
    denom = (γ⁺ .+ ν⁺).^2 -4 .* Kₘ .* γ⁺ .* ν⁺
    φᵣ = (1 - φₒ) .* numer ./ denom 

    if length(φᵣ) == 1
        φᵣ = φᵣ[1]
    end
    return φᵣ
end
export optimal_ribosomal_allocation

"""
    coord_conv(vals; out)
A helper function for converting coordinate systems between Ternary and Cartesian. 

# Arguments
- `vals::Vector{Float64}`: The coordinates to convert. If converting from 
  Ternary coordinates, provided values must be in the order of [a, b, c]. If 
  converting from Cartesian, provided values must be in the order of [x, y].
- `out::String="cartesian"`: The coordinate system to convert to. If converting 
  from Ternary coordinates, `out` must be "cartesian", "cart", "c", or "xy". 
  If converting from Cartesian coordinates, `out` must be "ternary", "tern", "t",
  or "abc". Default is "cartesian".

# Output
- `converted::Vector`: The converted coordinates. If converting to Cartesian
   coordinates, output will be a 2-vector in the order of [x, y]. If converting 
   to Ternary coordinates, output wil be a 3-vector in the order of [a, b, c].
"""
function coordinate_converter(vals::Vector{Float64};
                    out::String="cartesian")
    if lowercase(out) ∈ ["cartesian", "cart", "c", "xy"]
        _, b, c = vals 
        x = 0.5 * (2 * b + c) / sum(vals)
        y = (√(3)/2) * c / sum(vals)
        converted = [x,y]
    elseif  lowercase(out) ∈ ["ternary", "tern", "t", "abc"]
        x, y = vals
        c = 2 * y / √(3)
        b = (x - c) / 2
        a = 1 - b - c 
        converted = [a,b,c] 
    else 
        throw(ArgumentError("Return format not recognized must be in [\"cartesian\", \"cart\", \"c\"] for Cartesian output or [\"ternary\", \"tern\", \"t\", \"abc\"]."))
    end
    return converted
end
export coordinate_converter
using Base
using Dates
using DocStringExtensions
using DataFramesMeta

"""
FluxParity

A structure defining parameters of the flux-parity regulatory mechanism
for ribosomal allocation 

$(TYPEDFIELDS)
"""
Base.@kwdef struct FluxParityAllocator 
    "Allocation parameters of the metabolic sector for consuming mixed nutrients. "
    α::Vector{Vector{Float64}}
    "The concentration of nutrients in the feedstock. Can be either a time-dependent generating function or constant values."
    c₀::Vector
    "The degradation method" 
    δ_fun::Function
    "Maximum metabolic rates for consuming each nutrient in dimensions of inverse time."
    ν⁺::Vector{Float64} = 11.0 * ones(length(c₀))
end
export 


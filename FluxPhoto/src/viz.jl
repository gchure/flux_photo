using CairoMakie
using DataFramesMeta
using ProgressMeter
using TernaryPlots
"""
    get_colors(;allpalettes)

Generates a dictionary of my personally preferred color scheme and returns the 
dictionary and a vector of palettes or palette colors.

# Arguments

- `allpalettes::Bool=false` : If True, returned palette will be a vector of pallettes 
    in order of `dark`, `primary`, and `light`. If `false`, only `primary` will 
    be returned.
"""
function get_colors(;allpalettes::Bool=false)
    colors =Dict( 
            "dark_black" => "#2b2b2a",
            "black"=> "#3d3d3d",
            "primary_black"=> "#4c4b4c",
            "light_black"=> "#8c8c8c",
            "pale_black"=> "#afafaf",
            "dark_blue"=> "#154577",
            "blue"=> "#005da2",
            "primary_blue"=> "#3373ba",
            "light_blue"=> "#5fa6db",
            "pale_blue"=> "#8ec1e8",
            "dark_green"=> "#356835",
            "green"=> "#488d48",
            "primary_green"=> "#5cb75b",
            "light_green"=> "#99d097",
            "pale_green"=> "#b8ddb6",
            "dark_red"=> "#79302e",
            "red"=> "#a3433f",
            "primary_red"=> "#d8534f",
            "light_red"=> "#e89290",
            "pale_red"=> "#eeb3b0",
            "dark_gold"=> "#84622c",
            "gold"=> "#b1843e",
            "primary_gold"=> "#f0ad4d",
            "light_gold"=> "#f7cd8e",
            "pale_gold"=> "#f8dab0",
            "dark_purple"=> "#43355d",
            "purple"=> "#5d4a7e",
            "primary_purple"=> "#8066ad",
            "light_purple"=> "#a897c5",
            "pale_purple"=> "#c2b6d6" 
        )

        keys = ["black", "blue", "green","purple", "gold", "red"]
        dark_palette = [colors["dark_"*k] for k ∈ keys]
        primary_palette = [colors["primary_"*k] for k ∈ keys]
        light_palette = [colors["light_"*k] for k ∈ keys]

        if allpalettes
            pal = [primary_palette, dark_palette, light_palette]
        else
            pal = primary_palette
        end
        return [colors, pal]
end
export get_colors

function plotting_style(;colors::Bool=true,
                        palette::Bool=true)

    # Set font preferences
    if ~isfile(assetpath("fonts", "Lato-Regular.ttf"))
        @warn "Font 'Lato-Regular' not found. Defaulting to NotoSans..."
        regfont = assetpath("fonts", "NotoSans-Regular.ttf")
        boldfont = assetpath("fonts", "NotoSans-Regular.ttf")
    else
        regfont = assetpath("fonts", "Lato-Regular.ttf")
        boldfont = assetpath("fonts", "Lato-Bold.ttf")

    end

    # Load colors
    cors, pal = get_colors()

    # Set theme preferences 
    theme = Theme(
        Axis = (
            backgroundcolor="#f0f3f7",
            labelcolor     = "#5b5b5b",
            leftspinevisible = false,
            rightspinevisible = false,
            bottomspinevisible = false,
            topspinevisible = false,    
            grid           = true,
            xgridcolor = "#FFFFFF",
            ygridcolor = "#FFFFFF",
            gridlinewidth = 0.5,
            minorticks = false,
            xticksvisible = false,
            yticksvisible = false,

            # Fonts
            titlefont = boldfont,
            xlabelfont = regfont,
            ylabelfont = regfont,
            xticklabelfont = regfont,
            yticklabelfont = regfont,
        ),
        Cycle = (
            [c for c ∈ pal]
        ),
        backgroundcolor="#FFFFFF",
        linewidth=1
    )
    set_theme!(theme)
    out = []
    if colors
       push!(out, cors) 
    end
    if palette
        push!(out, pal)
    end
    if length(out) == 1
        return out[1]
    elseif length(out) > 0
        return [out...]
    else
        nothing
    end
end
export plotting_style
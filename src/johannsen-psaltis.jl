using SpectralFitting
using Plots
using Gradus
using GradusSpectralModels
using JLD2
# for JLD2 compression
using CodecZlib

include("additional-models.jl")

# emissivity profile table
a_range = collect(range(0.0, 0.998, 10))
h_range = collect(range(1.0, 10.0, 10))
eps_range = collect(range(-0.5, 3.0, 10))

prof_wrap = @time tabulate_emissivity_profile_jp(a_range, h_range, eps_range)

lamp_post_grid = map(eachindex(IndexCartesian(), prof_wrap)) do index
    if isnothing(prof_wrap[index])
        ai, hi, epsi = index.I
        _eps = if ai == 10
            2
        elseif ai == 9
            3
        elseif ai == 8
            5
        elseif ai ==7
            8
        else
            @error "Unhandled: $index"
            # @show index
        end
        prof_wrap[ai, hi, _eps]
    else
        prof_wrap[index]
    end
end;

interp = Gradus.MultilinearInterpolator{3}(lamp_post_grid)
lp_table = DiscProfileWrapper((a_range, h_range, eps_range), lamp_post_grid, interp)

# write to file
jldopen("johannsen_lamp_post_extended.jld2", "w"; compress = true) do f
    f["lamp_post"] = lp_table
end


table = read_or_make_transfer_table("data/johannsen-tf-grid-more-eps.jld2")

_updated_table = map(eachindex(IndexCartesian(), table)) do index
    if isnothing(table[index])
        ai, hi, epsi = index.I
        _eps = if ai == 10
            2
        elseif ai == 9
            3
        elseif ai == 8
            5
        elseif ai ==7
            8
        else
            @error "Unhandled: $index"
        end
        table[ai, hi, _eps]
    else
        table[index]
    end
end;

jldopen("johannsen_transfer_table_extended.jld2", "w"; compress = true) do f
    f["table"] = _updated_table
end

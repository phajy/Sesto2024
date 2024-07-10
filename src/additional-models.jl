using Gradus, SpectralFitting

struct DiscProfileGrid{V<:AbstractVector}
    radii::V
    emissivity::V
end

function Gradus.MultiLinearInterpolations.restructure(
    grid::DiscProfileGrid,
    vs::AbstractVector,
)
    @views begin
        DiscProfileGrid(
            vs[1:length(grid.radii)],
            vs[length(grid.radii)+1:length(grid.radii)+length(grid.emissivity)],
        )
    end
end

struct DiscProfileWrapper{N,T,I}
    # a and h
    params::NTuple{N,Vector{T}}
    grids::Array{DiscProfileGrid{Vector{T}},N}
    cache::I
end
function (d::DiscProfileWrapper)(args...)
    dpg = Gradus.interpolate!(d.cache, d.params, d.grids, promote(args...))
    Gradus.RadialDiscProfile(dpg.radii, dpg.radii, dpg.emissivity)
end

struct ThickDiscLineProfile{D,T} <: AbstractTableModel{T,Additive}
    table::D
    K::T
    "Spin"
    a::T
    "Observer inclination (degrees off of the spin axis)."
    θ::T
    "Eddington ratio"
    η::T
    "Inner radius of the accretion disc."
    rin::T
    "Outer radius of the accretion disc."
    rout::T
    "Central emission line energy (keV)."
    E₀::T
end

function ThickDiscLineProfile(
    profile,
    table::Gradus.CunninghamTransferTable;
    K = FitParam(1.0),
    a = FitParam(0.998),
    θ = FitParam(30.0),
    η = FitParam(0.2, lower_limit = 0.0, upper_limit = 0.3),
    rin = FitParam(1.0),
    rout = FitParam(100.0, upper_limit = 100.0),
    E₀ = FitParam(1.0),
    kwargs...,
)
    setup = integration_setup(
        profile,
        table((get_value(a), get_value(θ), get_value(η)));
        kwargs...,
    )
    ThickDiscLineProfile((; setup = setup, table = table), K, a, θ, η, rin, rout, E₀)
end

function SpectralFitting.invoke!(output, domain, model::ThickDiscLineProfile)
    grid = model.table.table((model.a, model.θ, model.η))
    rmin = if model.rin < grid.r_grid[1]
        grid.r_grid[1]
    else
        model.rin
    end
    Gradus.integrate_lineprofile!(
        output,
        model.table.setup,
        grid,
        domain;
        rmin = rmin,
        rmax = model.rout,
        g_scale = model.E₀,
    )
    output
end


struct LampPostLineProfile{D,T} <: AbstractTableModel{T,Additive}
    table::D
    K::T
    "Spin"
    a::T
    "Observer inclination (degrees off of the spin axis)."
    θ::T
    "Lamp post corona height"
    h::T
    "Inner radius of the accretion disc."
    rin::T
    "Outer radius of the accretion disc."
    rout::T
    "Central emission line energy (keV)."
    E₀::T
end

function LampPostLineProfile(
    profile,
    table::Gradus.CunninghamTransferTable;
    K = FitParam(1.0),
    a = FitParam(0.998),
    θ = FitParam(45.0),
    h = FitParam(5.0, lower_limit = 0.0, upper_limit = 20.0),
    rin = FitParam(1.0),
    rout = FitParam(100.0, upper_limit = 100.0),
    E₀ = FitParam(1.0),
    kwargs...,
)
    # just use anything for the emissivity profile, since we are going to override it during integration anyway
    setup = integration_setup(
        profile(get_value(a), get_value(h)),
        table((get_value(a), get_value(θ)));
        kwargs...,
    )
    LampPostLineProfile(
        (; setup = setup, profile = profile, table = table),
        K,
        a,
        θ,
        h,
        rin,
        rout,
        E₀,
    )
end

function Base.copy(m::LampPostThickDisc)
    table = Gradus.CunninghamTransferTable(m.table.table.params, m.table.table.grids)
    setup = Gradus.IntegrationSetup(
        m.table.setup.h,
        m.table.setup.time,
        m.table.setup.integrand,
        m.table.setup.pure_radial,
        m.table.setup.quadrature_rule,
        deepcopy(m.table.setup.index_cache),
        m.table.setup.g_grid_upscale,
        m.table.setup.n_radii,
    )
    profile = deepcopy(m.table.profile)
    typeof(m)(
        (; setup = setup, profile = profile, table = table),
        (copy(getproperty(m, f)) for f in fieldnames(typeof(m))[2:end])...,
    )
end

function SpectralFitting.invoke!(output, domain, model::LampPostLineProfile)
    grid = model.table.table((model.a, model.θ))
    rmin = if model.rin < grid.r_grid[1]
        grid.r_grid[1]
    else
        model.rin
    end
    Gradus.integrate_lineprofile!(
        output,
        model.table.setup,
        grid,
        domain;
        rmin = rmin,
        rmax = model.rout,
        pure_radial = r -> emissivity_at(model.table.profile(model.a, model.h), r),
        g_scale = model.E₀,
    )
    output
end

function tabulate_emissivity_profile(a_range, h_range; r_max = 150.0, n_radii = 100)
    total = length(a_range) * length(h_range)
    i = 0
    function _wrapper(a, h)
        i += 1
        @info "$i of $total"
        m = KerrMetric(1.0, a)
        d = ThinDisc(0.0, Inf)
        radii = Gradus.Grids._inverse_grid(Gradus.isco(m) + 1e-2, r_max, n_radii) |> collect
        if h < 0.5 + Gradus.inner_radius(m)
            DiscProfileGrid(radii, zeros(eltype(radii), size(radii)))
        else
            prof = @time emissivity_profile(m, d, LampPostModel(h = h), n_samples = 2000)
            ems = emissivity_at.((prof,), radii)
            DiscProfileGrid(radii, ems)
        end
    end
    vals = [_wrapper(a, h) for a in a_range, h in h_range]
    interp = Gradus.MultilinearInterpolator{2}(vals)
    DiscProfileWrapper((a_range, h_range), vals, interp)
end

function tabulate_emissivity_profile(
    a_range,
    h_range,
    η_range;
    r_max = 150.0,
    n_radii = 100,
)
    total = length(a_range) * length(h_range) * length(η_range)
    i = 0
    function _wrapper(a, h, η)
        i += 1
        @info "$i of $total"
        m = KerrMetric(1.0, a)
        d = ShakuraSunyaev(m; eddington_ratio = η)
        radii = Gradus.Grids._inverse_grid(Gradus.isco(m) + 1e-2, r_max, n_radii) |> collect
        if h < 0.5 + Gradus.inner_radius(m)
            DiscProfileGrid(radii, zeros(eltype(radii), size(radii)))
        else
            prof = @time emissivity_profile(m, d, LampPostModel(h = h), n_samples = 3000)
            ems = emissivity_at.((prof,), radii)
            DiscProfileGrid(radii, ems)
        end
    end
    vals = [_wrapper(a, h, η) for a in a_range, h in h_range, η in η_range]
    interp = Gradus.MultilinearInterpolator{3}(vals)
    DiscProfileWrapper((a_range, h_range, η_range), vals, interp)
end

struct LampPostThickDisc{D,T} <: AbstractTableModel{T,Additive}
    table::D
    K::T
    "Spin"
    a::T
    "Observer inclination (degrees off of the spin axis)."
    θ::T
    "Lamp post corona height"
    h::T
    "Eddington ratio"
    η::T
    "Inner radius of the accretion disc."
    rin::T
    "Outer radius of the accretion disc."
    rout::T
    "Central emission line energy (keV)."
    E₀::T
end

function LampPostThickDisc(
    profile,
    table::Gradus.CunninghamTransferTable;
    K = FitParam(1.0),
    a = FitParam(0.998, lower_limit = 0.0, upper_limit = 0.998),
    θ = FitParam(45.0, lower_limit = 5.0, upper_limit = 85.0),
    h = FitParam(5.0, lower_limit = 0.0, upper_limit = 20.0),
    η = FitParam(0.01, lower_limit = 0.01, upper_limit = 0.3),
    rin = FitParam(1.0),
    rout = FitParam(100.0, upper_limit = 100.0),
    E₀ = FitParam(1.0),
    kwargs...,
)
    # just use anything for the emissivity profile, since we are going to override it during integration anyway
    setup = integration_setup(
        profile(get_value(a), get_value(h), get_value(η)),
        table((get_value(a), get_value(θ), get_value(η)));
        kwargs...,
    )
    LampPostThickDisc(
        (; setup = setup, profile = profile, table = table),
        K,
        a,
        θ,
        h,
        η,
        rin,
        rout,
        E₀,
    )
end

function SpectralFitting.invoke!(output, domain, model::LampPostThickDisc)
    grid = model.table.table((model.a, model.θ, model.η))
    rmin = if model.rin < grid.r_grid[1]
        grid.r_grid[1]
    else
        model.rin
    end
    Gradus.integrate_lineprofile!(
        output,
        model.table.setup,
        grid,
        domain;
        rmin = rmin,
        rmax = model.rout,
        pure_radial = r -> emissivity_at(model.table.profile(model.a, model.h, model.η), r),
        g_scale = model.E₀,
    )
    output
end

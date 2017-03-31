"""
SACPlot.jl provides routines for plotting SAC traces.
"""
module SACPlot
# Module for plotting SAC traces

using SAC
import Plots

import SAC.SACtr

export
    plot1,
    p1,
    plot2,
    p2,
    plotpm,
    ppm,
    plotrs,
    prs,
    plotsp,
    psp

const TIME_PICKS = [:a, :t0, :t1, :t2, :t3, :t4, :t5, :t6, :t7, :t8, :t9, :f]
const NAME_PICKS = [:ka, :kt0, :kt1, :kt2, :kt3, :kt4, :kt5, :kt6, :kt7, :kt8, :kt9, :kf]

"Maximum number of samples to show by default"
const sacplot_qdp_thresh = 10_000

"""
    plot1(s::Array{SACtr}; kwargs...) -> ::Plots.Plot
    plot1(s::SACtr; kwargs...) -> ::Plots.Plot

Create a plot of the SAC trace(s) `s` and return a `Plots.Plot` object.

## Keyword arguments

| Keyword   | Example values | Description |
|:----------|:---------------|:------------|
|`label`    |`[:kcmpnm, :user0]`| A single or array of symbols listing the header values to show.|
|`line`     |`(:black, 2)`   | An argument to Plots defining the line property to use.  Arrays will apply in turn to each trace.|
|`over`     |`true`          | If `true`, plot on top of the active Plots plot object|
|`qdp`      |`false`         | Plot every single sample if `false`.  (Default produces 'quick-and-dirty plot'.)|
|`relative` |`true`, `:user0`| Plot all times relative to beginning of trace if `true` or a header value if a `Symbol`|
|`xlim`     |`10:12`         | A range, tuple or array giving the limits in time to plot.  `NaN` as a limit uses the limits of the data.|
|`ylim`     |`(NaN, 2e-9)`   | As for `xlim`, but controlling the y-axis.  `:all` will force all traces to have the same scale.|
|`xlabel`   |`"Time / s"`    | Set the independent axis label|
|`ylabel`   |`"Amplitude / nm"`| Set the dependent axis label|
"""
function plot1(a::Array{SACtr};
               label=:default, line=(:black,), over=true, qdp=true,
               relative=false, title="", xlabel="Time / s", xlim=[NaN, NaN],
               ylabel=nothing, ylim=[NaN, NaN])
    # Check arguments
    check_lims(xlim, "SACPlot.plot1")
    typeof(ylim) <: Union{String, Symbol} || check_lims(ylim, "SACPlot.plot1")
    # Make sure any single label is in an array
    typeof(label) == Symbol && (label = [label])
    for l in label
        if l != :default
            l in [:kzdate, :kztime] && continue # Special labels
            l in fieldnames(a[1]) || error("Label header '$l' is not a valid header variable")
        end
    end
    # Make plots
    n = length(a)
    b, e, depmin, depmax = lims(a, relative)
    if ! isnan(xlim[1]); b = xlim[1]; end
    if ! isnan(xlim[end]); e = xlim[end]; end
    if ! (typeof(ylim) <: Union{String,Symbol})
        if ! isnan(ylim[1]); depmin = ylim[1]; end
        if ! isnan(ylim[end]); depmax = ylim[end]; end
    end
    # Downsample plot
    iskip = qdp_skip(a, qdp)
    # Turn off x labels for all but bottom trace
    p = over ? Plots.current() : Plots.plot(layout=(n,1), grid=false, legend=false)
    for i = 1:n
        # Which samples to plot
        inds = qdp_inds(a[i], iskip)
        # Traces
        offset = lims_offset(a[i], relative)
        Plots.plot!(p[i], SAC.time(a[i])[inds] .- offset, a[i].t[inds],
            xlim=(b, e), xticks=(i==n ? true : false), l=line)
        if typeof(ylim) <: Union{String, Symbol}
            if Symbol(lowercase(string(ylim))) == :all
                Plots.ylims!(p[i], (depmin, depmax))
            end
        elseif ! all([isnan(ylim[1]), isnan(ylim[end])])
            Plots.ylims!(p[i], (depmin, depmax))
        end
        y1, y2 = Plots.ylims(p[i])
        # Add picks
        for (tp, kp) in zip(TIME_PICKS, NAME_PICKS)
            t, k = getfield(a[i], tp), strip(getfield(a[i], kp))
            if t != SAC.sac_rnull
                t -= offset
                Plots.plot!(p[i], [t, t], [y1, y2], l=:blue)
                if tp in [:a, :f] && k == SAC.sac_cnull
                    k = uppercase("$tp")
                end
                k != SAC.sac_cnull &&
                    Plots.plot!(p[i],
                        ann=(t, y1+0.05*(y2-y1), Plots.text(k, 8, :bottom, :left, :blue)))
            end
        end
        # Add text annotation
        if label != :nothing
            if :default in label
                if all(Bool[getfield(a[i], f) != SAC.sac_inull for f in [:nzyear, :nzhour, :nzmin, :nzsec, :nzmsec]])
                    date = Date(Date(0) + Dates.Year(a[i].nzyear) + Dates.Day(a[i].nzjday))
                    day = string(Dates.dayofmonth(date))
                    month = Dates.monthname(date)[1:3]
                else
                    day = month = "?"
                end
                label_string = strip(a[i].kevnm)*"\n"*strip(a[i].kstnm)*"."*
                               strip(a[i].knetwk)*"."*strip(a[i].kcmpnm) * "\n" *
                         @sprintf("%04d.%03dT%02d:%02d:%02d:%03d (%s %s)", a[i].nzyear,
                                  a[i].nzjday, a[i].nzhour, a[i].nzmin, a[i].nzsec,
                                  a[i].nzmsec, day, month)
            else
                label_string = ""
                for l in label
                    label_string *= string(l) * ": " * strip(string(getfield(a[i], l))) * "\n"
                end
            end
            Plots.plot!(p[i], ann=(b+0.97*(e-b), y1+0.95*(y2-y1),
                Plots.text(label_string, 8, :top, :right, :black)), annotationguide=:below)
        end
    end
    title != nothing && Plots.title!(p[1], title)
    xlabel != nothing && Plots.xlabel!(p[end], xlabel)
    ylabel != nothing && Plots.ylabel!(p[end÷2], ylabel)
    p
end

# Single-trace version
plot1(s::SACtr; kwargs...) = plot1([s]; kwargs...)

p1 = plot1

"""
    plot2(::Array{SACtr})

Plot all traces in array of SAC traces `a` on the same plot
"""
function plot2(a::Array{SACtr}; legend=false)
    p = Plots.plot(legend=legend, show=false)
    for i = 1:length(a)
        Plots.plot!(p, SAC.time(a[i]), a[i].t)
    end
    Plots.plot!(p, xlim=lims(a)[1:2], show=true)
    p
end

# Single-trace version
plot2(s::SACtr) = plot2([s])

p2 = plot2

"""
    plotpm(::Array{SACtr}; xlim=[NaN, NaN])

Plot the particle motion for a pair of orthogonal components, within
the time window `xlim[1]` to `xlim[2]` if provided
"""
function plotpm(a::Array{SACtr}; xlim=[NaN, NaN])
    const angle_tol = 0.1
    length(a) == 2 || error("plotpm: Can only plot two components")
    angle = (a[1].cmpaz - a[2].cmpaz)%360.
    abs(angle) - 90. > eps(angle) && abs(angle) - 270. > eps(angle) &&
        error("plotpm: Components must be orthogonal")
    # Find out whether trace 1 is clockwise of 2, or vice versa
    if abs(angle) - 90. <= angle_tol
        t1, t2 = a[1], a[2]
    else
        t1, t2 = a[2], a[1]
    end
    _, _, amin, amax = lims([t1, t2])
    p = Plots.plot(t1.t, t2.t, xlim=(amin, amax), ylim=(amin, amax),
        aspect_ratio=:equal, legend=false)
    Plots.xlabel!(p, strip(t1.kcmpnm) * " (" * string(t1.cmpaz) * ")")
    Plots.ylabel!(p, strip(t2.kcmpnm) * " (" * string(t2.cmpaz) * ")")
    p
end

plotpm(s1::SACtr, s2::SACtr; xlim=[NaN, NaN]) = plotpm([s1, s2]; xlim=xlim)
ppm = plotpm

"""
    plotsp(f::Array, S::Array, kind=\"amp\")

Plot the Fourier-transformed trace `S`, with frequencies `f`.

`kind` may be one of:\n
`amp`  : Plot amplitude\n
`phase`: Plot phase\n
`real` : Plot real part\n
`imag` : Plot imaginary part\n
"""
function plotsp(f::Array{Array, 1}, S::Array{Array, 1}, kind="amp";
                xlim=[NaN, NaN], ylim=[NaN, NaN])
    error("`plotsp` is not implemented yet")
    check_lims(xlim, "SACPlot.plotsp")
    check_lims(ylim, "SACPlot.plotsp")
end
psp = plotsp

#TODO: Make plotting limits and automatic scaling more sensible when some traces
#      have very different amplitudes
"""
    plotrs(s, align=:o; tw, dw, style="-r", size=1, y=:gcarc, over=false)

Plot a record section of the traces in `s`, aligned on the time given by
`align`.  This may be a header name, given as a Symbol, or an array of times.

By default, the y-axis is distance (`:gcarc`), but this can by any header name
passed as a symbol to the `y` keywords argument, or an arbitrary array of values
(for instance, to plot equally spaced in some order).

## Keyword arguments

|Name   |Type          |Description|
|:------|:-------------|:----------|
|`dw`   |Range or array|Set *distance window* (or other y-axis variable) for plotting|
|`fill` |`Tuple`       |Set colour fill.  Pass a 2-tuple of (positive_colour, negative_colour)|
|`label`|`Symbol`      |Label traces with header value|
|`over` |`Bool`        |If `true`, overplot this record section over the previous plot|
|`qdp`  |`Bool`        |If `true`, reduce the number of points plotted for speed. (Default is automatic)|
|`reverse`|`Bool`      |If true, reverse the direction of the y-axis.  (Default for `y=:gcarc`)|
|`tw`   |Range or array|Set *time window* for plotting|
|`size` |Real          |Scaling factor for traces|
|`line` |`Any`         |Argument passed to Plots specifying style for lines|
|`y`    |`Symbol` or array|Header value or array of values to use for y-axis|
"""
function plotrs(s::Array{SACtr}, align=0.;
        tw=[nothing, nothing], dw=[nothing, nothing], y=:gcarc, line=(:black,),
        size::Real=1., over::Bool=false, reverse=nothing, fill=(nothing, nothing),
        qdp=10_000, label=nothing,
        kwargs...)
    maxamp = maxabs([s[:depmax]; s[:depmin]])
    y_shift = _y_shifts(s, y)
    d = _x_shifts(s, align)
    reverse = if reverse == nothing
        if y == :gcarc true else false end
    else
        reverse
    end
    scale = size*abs(minimum(y_shift) - maximum(y_shift))/10
    reverse && (scale *= -1) # Set positive values to still plot 'up' the page if reversed axis
    # Determine a limit on total number of points to plot
    iskip = qdp_skip(s, qdp)
    # Create plot here
    over || Plots.plot(legend=false, grid=false, display=false)
    # Traces
    for i in 1:length(s)
        t = SAC.time(s[i]) + d[i]
        inds = qdp_inds(s[i], iskip)
        # FIXME: Filled wiggles not working with Plots yet
        any([fill[1], fill[2]] .!= nothing) && i == 1 &&
            warn("Filled wiggles not implemented yet")
        #=
        if fill[1] != nothing
            Plots.plot!(t[inds], (s[i].t[inds] + y_shift[i]).*[s[i].t[j] >= 0. ?
                scale/maxamp : NaN for j in inds], fillrange=y_shift[i], c=fill[1])
        end
        if fill[end] != nothing
            Plots.plot!(t[inds], (s[i].t[inds] + y_shift[i]).*[s[i].t[j] <= 0. ?
                scale/maxamp : NaN for j in inds], fillrange=y_shift[i], c=fill[end])
        end
        =#
        Plots.plot!(t[inds], y_shift[i] + s[i].t[inds]*scale/maxamp, l=line)
    end
    # x limits
    t1 = minimum(s[:b] + d)
    t2 = maximum(s[:e] + d)
    if !(tw[1] == tw[end] == nothing)
        tw[1] != nothing && (t1 = tw[1])
        tw[end] != nothing && (t2 = tw[end])
    end
    Plots.xlims!(t1, t2)
    # y limits
    if !(dw[1] == dw[end] == nothing)
        d1, d2 = Plots.ylims()
        dw[1] != nothing && (d1 = dw[1])
        dw[end] != nothing && (d2 = dw[end])
        Plots.ylims!(d1, d2)
    end
    reverse && Plots.yflip!()
    # Labels
    if label != nothing
        t1, t2 = Plots.xlims()
        for i in 1:length(s)
            Plots.annotate!(t2, y_shift[i], Plots.text("$(s[i][label])", 8, :left, :black),
                right_margin=10Plots.mm)
        end
    end
    # Ensure the border looks nice (not true with all backends)
    draw_borders!(;display=true, kwargs...)
    return
end
prs = plotrs

# Routines to calculate offsets in x and y for record section plots
_y_shifts(s::Array{SACtr}, y::Symbol) = s[y]
_y_shifts(s::Array{SACtr}, y::AbstractArray) =
    length(y) == length(s) ? y : error("Length of `y` not the same as number of traces")
_x_shifts(s::Array{SACtr}, t::Symbol) = -s[t]
_x_shifts(s::Array{SACtr}, t::AbstractArray) =
    length(s) == length(t) ? -t : error("Length of `align` not the same as number of traces")
_x_shifts(s::Array{SACtr}, t::Real) = -t*ones(length(s))

"""
    lims(::Array{SACtr}, relative=false) -> b, e, depmin, depmax

Get the minimum and maximum times, `b` and `e`, for an array of SAC traces.
If `relative` is true, these times are all relative to the beginning of the trace,
not the zero time.
"""
function lims(a::Array{SACtr}, relative=false)
    n = length(a)
    offset = lims_offset(a[1], relative)
    b = a[1].b - offset
    e = a[1].e - offset
    depmin = a[1].depmin
    depmax = a[1].depmax
    for i = 2:n
        offset = lims_offset(a[i], relative)
        if a[i].b - offset < b; b = a[i].b - offset; end
        if a[i].e - offset > e; e = a[i].e - offset; end
        if a[i].depmin < depmin; depmin = a[i].depmin; end
        if a[i].depmax > depmax; depmax = a[i].depmax; end
    end
    return b, e, depmin, depmax
end
lims(s::SACtr) = s.b, s.e, s.depmin, s.depmax
lims_offset(s::SACtr, relative::Bool) = relative ? s.b : zero(SAC.SACFloat)
lims_offset(s::SACtr, relative::Symbol) = getfield(s, relative)

"""
    check_lims(a, routine=\"SACPlot.check_lims\")

Throw an error if `a` is not a length-2 array with the first value lower
or equal to the second, if both are not NaN.
"""
function check_lims(a, routine="SACPlot.check_lims")
    a[1] > a[2] && # NB: (NaN {>,<,==} NaN) == false
        error(routine * ": Lower plot limit value must be less than upper")
end

"""
    draw_borders(line=(:black,1.5))

Draw a border around the current plot.
"""
function draw_borders!(p::Union{Plots.Plot,Plots.Subplot}, line=(:black,1.5); kwargs...)
    (x1, x2) = Plots.xlims(p)
    (y1, y2) = Plots.ylims(p)
    Plots.plot!(p, [x1, x2, x2, x1, x1], [y1, y1, y2, y2, y1], l=line; kwargs...)
end
draw_borders!(line=(:black, 1.5); kwargs...) =
    draw_borders!(Plots.current(), line; kwargs...)

"""
    qdp_inds(s::SACtr, iskip) -> inds

Return `inds`, a range containing which indices to use for 'quick-and-dirty
plotting' when only each `iskip`th point is needed.
"""
qdp_inds(s::SACtr, iskip::Integer) = 1:min(iskip, s.npts-1):s.npts

"""
    qdp_skip(a::Union{SACtr,Array{SACtr}}, thresh::Integer=$(sacplot_qdp_thresh))

Return `iskip`, such that only the `iskip`th sample need be plotted when using
'quick-and-dirty plotting'.
"""
qdp_skip(a::Union{SACtr,Array{SACtr}}, thresh::Integer=sacplot_qdp_thresh) =
    max(1, sum(a[:npts])÷thresh)
qdp_skip(a::Union{SACtr,Array{SACtr}}, thresh::Void) = qdp_skip(a)
qdp_skip(a::Union{SACtr,Array{SACtr}}, thresh::Bool) = thresh ? qdp_skip(a) : 1


end # module SACPlot

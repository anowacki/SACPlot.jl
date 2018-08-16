"""
SACPlot.jl provides routines for plotting SAC traces.
"""
module SACPlot

using Compat.Dates
using Compat.Printf
import Compat: Nothing, @warn
import DSP

import Plots

using SAC
import SAC: SACtr

const SACArray = AbstractArray{SACtr}

export
    plot1,
    p1,
    plot2,
    p2,
    plotpm,
    ppm,
    plotrs,
    plotrs!,
    plotspec,
    prs,
    prs!

const TIME_PICKS = [:a, :t0, :t1, :t2, :t3, :t4, :t5, :t6, :t7, :t8, :t9, :f]
const NAME_PICKS = Symbol.("k" .* String.(SACPlot.TIME_PICKS))

"Maximum number of samples to show by default"
const sacplot_qdp_thresh = 100_000
"Default plotting parameters for traces"
const sacplot_p1_defaults = Dict{Symbol,Any}()
"Default plotting parameters for record section"
const sacplot_prs_defaults = Dict{Symbol,Any}()

function __init__()
    props = Dict(:legend => false,
                 :grid => false,
                 :display => false,
                 :framestyle => :box)
    for (k, v) in props
        sacplot_p1_defaults[k] = v
    end
    props = Dict(:legend => false,
                 :grid => false,
                 :display => false,
                 :framestyle => :box)
    for (k, v) in props
        sacplot_prs_defaults[k] = v
    end
end

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
|`picks`    |`(["SKS"], [1610])`| Add arrival time picks to plot from a tuple of two arrays, one giving the phase name and the second the arrival time.|
|`qdp`      |`false`         | Plot every single sample if `false`.  (Default produces 'quick-and-dirty plot'.)|
|`relative` |`true`, `:user0`| Plot all times relative to beginning of trace if `true` or a header value if a `Symbol`|
|`xlim`     |`10:12`         | A range, tuple or array giving the limits in time to plot.  `NaN` as a limit uses the limits of the data.|
|`ylim`     |`(NaN, 2e-9)`   | As for `xlim`, but controlling the y-axis.  `:all` will force all traces to have the same scale.|
|`xlabel`   |`"Time / s"`    | Set the independent axis label|
|`ylabel`   |`"Amplitude / nm"`| Set the dependent axis label|
"""
function plot1(a::SACArray;
               label=:default, line=(:black,), over=false, picks=([],[]), qdp=true,
               relative=false, title="", xlabel="Time / s", xlim=[NaN, NaN],
               ylabel=nothing, ylim=[NaN, NaN], plots_kwargs...)
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
    p = over ? Plots.current() : Plots.plot(layout=(n,1), grid=false, legend=false,
        frame_style=:box; plots_kwargs...)
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
        # Add custom picks
        phase_names, times = picks
        length(times) == length(phase_names) ||
            error("`picks` must be a tuple of phase name and time, with both arrays the same length")
        if length(times) > 0
            ilabel = 0
            for (k, t) in zip(phase_names, times)
                if b <= t <= e
                    Plots.plot!(p[i], [t, t], [y1, y2], l=:red)
                    Plots.plot!(p[i], ann=(t, y1+(1-0.05*ilabel)*(y2-y1), Plots.text(k, 8, :top, :left, :red)))
                    ilabel = (ilabel + 1)%4
                end
            end
        end
        # Add text annotation
        if label != :nothing
            if :default in label
                if all(Bool[getfield(a[i], f) != SAC.sac_inull for f in [:nzyear, :nzhour, :nzmin, :nzsec, :nzmsec]])
                    date = Date(Date(0) + Dates.Year(a[i].nzyear) + Dates.Day(a[i].nzjday - 1))
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
# Iterable which contains SAC traces
plot1(args...; kwargs...) = plot1(SACtr[args...]; kwargs...)

p1 = plot1

"""
    plot2(::Array{SACtr}; relative=false)

Plot all traces in array of SAC traces `a` on the same plot
"""
function plot2(a::SACArray; relative=false, legend=false)
    p = Plots.plot(legend=legend, show=false)
    offset = relative ? a[:b] : zeros(length(a))
    for i = 1:length(a)
        Plots.plot!(p, SAC.time(a[i]).-offset[i], a[i].t)
    end
    Plots.plot!(p, xlim=lims(a, relative)[1:2], show=true)
    p
end

# Single-trace version
plot2(s::SACtr) = plot2([s])

p2 = plot2

"""
    plotpm(::Array{SACtr}; xlim=[NaN, NaN], geog=false, kwargs...) -> ::Plots.Plot

Plot the particle motion for a pair of orthogonal components, within
the time window `xlim[1]` to `xlim[2]` if provided.

If `geog` is `true`, then plot the particle motion as relative to north,
with north up the page, and mark on the component directions.

Pass any further plotting commands to `Plots` with `kwargs...`.
"""
function plotpm(a::SACArray; xlim=[NaN, NaN], geog=false, kwargs...)
    angle_tol = 0.1
    length(a) == 2 || error("plotpm: Can only plot two components")
    angle = SAC.angle_difference(a[1][:cmpaz], a[2][:cmpaz])
    mod(angle - 90., 180.) <= angle_tol ||
        error("plotpm: Components must be orthogonal")
    # Set so that t1 is clockwise of t2 (e.g., t1 is E and t2 is N)
    if angle - 90. <= angle_tol
        t1, t2 = a[1], a[2]
    else
        t1, t2 = a[2], a[1]
    end
    amax = sqrt(maximum(a[1].t.^2 + a[2].t.^2))
    p = Plots.plot(xlim=(-amax,amax), ylim=(-amax,amax), aspect_ratio=:equal,
        legend=false, kwargs...)
    if geog
        az1, az2 = deg2rad(t1[:cmpaz]), deg2rad(t2[:cmpaz])
        Plots.plot!(p, amax.*[sin(az1), 0, sin(az2)], amax.*[cos(az1), 0, cos(az2)],
            line=(:darkgray, 2, :dash))
        for t in (t1, t2)
            rotation = mod(180-t[:cmpaz], 180) - 90
            font = Plots.Font("sans-serif", 8, :hcenter, :bottom,
                rotation, Plots.RGB(0.0,0.0,0.0))
            Plots.annotate!(p, 0.5amax*sind(t[:cmpaz]), 0.5amax*cosd(t[:cmpaz]),
                Plots.PlotText(t[:kcmpnm], font))
        end
        amp = sqrt.(t1.t.^2 + t2.t.^2)
        phi = atan2.(t1.t, t2.t) # Angle from t2->t1 (e.g., N->E)
        Plots.plot!(p, amp.*sin.(phi + az2), amp.*cos.(phi + az2),
            xlabel="East", ylabel="North")
    else
        Plots.plot!(p, t1.t, t2.t,
            xlabel=strip(t1.kcmpnm) * " (" * string(t1.cmpaz) * ")",
            ylabel=strip(t2.kcmpnm) * " (" * string(t2.cmpaz) * ")")
    end
    p
end

plotpm(s1::SACtr, s2::SACtr; kwargs...) = plotpm([s1, s2]; kwargs...)
ppm = plotpm

const _plotrs_first_run = Ref(true)

#TODO: Make plotting limits and automatic scaling more sensible when some traces
#      have very different amplitudes
"""
    plotrs!(p::Plots.Plot, s, align=:o; plots_kwargs...) -> ::Plots.Plot
    plotrs(s, align=:o; prs_kwargs..., plots_kwargs...) -> ::Plots.Plot

Plot a record section of the traces in `s`, aligned on the time given by
`align`.  This may be a header name, given as a Symbol, or an array of times.

By default, the y-axis is distance (`:gcarc`), but this can by any header name
passed as a symbol to the `y` keywords argument, or an arbitrary array of values
(for instance, to plot equally spaced in some order).

Other plotting options are controlled by keywords arguments `prs_kwargs`:

## Keyword arguments

|Name   |Type          |Description|
|:------|:-------------|:----------|
|`dw`   |Range or array|Set *distance window* (or other y-axis variable) for plotting|
|`fill` |`Tuple`       |Set colour fill.  Pass a 2-tuple of (positive_colour, negative_colour) or (pos, neg, level) to fill from/to a certain level between -1 and 1|
|`label`|`Symbol`      |Label traces with header value or array of values|
|`qdp`  |`Bool`        |If `true`, reduce the number of points plotted for speed. (Default is automatic)|
|`reverse`|`Bool`      |If `true`, reverse the direction of the y-axis.  (Default for `y=:gcarc`)|
|`tw`   |Range or array|Set *time window* for plotting|
|`scale`|Real          |Scaling factor for traces (can be negative to reverse polarity)|
|`line` |`Any`         |Argument passed to Plots specifying style for lines.  Use `false` to not plot lines.|
|`y`    |`Symbol` or array|Header value or array of values to use for y-axis|

Keyword arguments to be passed to Plots are those appended after any of the above;
e.g., to plot a record section aligned on the `A` header, and set the size of the
plot window using the Plots `size` argument, you can do:

```
julia> prs(sample(:array) |> cut(:a, -30, :a, 30), :a, scale=0.5, size=(500,1000))
```
"""
function plotrs!(p::Union{Plots.Plot,Plots.Subplot}, s::SACArray, align=0.;
        tw=[nothing, nothing], dw=[nothing, nothing], y=:gcarc, line=(:black,),
        scale::Real=1.0, reverse=nothing, fill=(nothing, nothing),
        qdp=sacplot_qdp_thresh, label=nothing,
        kwargs...)
    maxamp = maximum(abs, (s[:depmax]; s[:depmin]))
    y_shift = _y_shifts(s, y)
    d = _x_shifts(s, align)
    reverse = if reverse == nothing
        if y == :gcarc true else false end
    else
        reverse
    end
    scale = scale*abs(minimum(y_shift) - maximum(y_shift))/10
    reverse && (scale *= -1) # Set positive values to still plot 'up' the page if reversed axis
    # Determine a limit on total number of points to plot
    iskip = qdp_skip(s, qdp)
    # Traces
    for i in 1:length(s)
        inds = qdp_inds(s[i], iskip)
        t = (SAC.time(s[i]) .+ d[i])[inds]
        trace = y_shift[i] .+ s[i].t[inds].*scale/maxamp
        # FIXME: Filled wiggles not working with Plots yet
        if _plotrs_first_run[] && any([fill[1], fill[2]] .!= nothing)
            _plotrs_first_run[] = false
            @warn("Filled wiggles do not display perfectly yet")
        end
            
        fill_level = length(fill) >= 3 ? y_shift[i]+fill[3]*scale : y_shift[i]
        if fill[1] != nothing
            Plots.plot!(p, t, [tt  < fill_level ? tt : NaN for tt in trace],
                fillrange=fill_level, c=fill[1], line=Plots.stroke(0))
        end
        if fill[2] != nothing
            Plots.plot!(p, t, [tt >= fill_level ? tt : NaN for tt in trace],
                fillrange=fill_level, c=fill[2], line=Plots.stroke(0))
        end
        
        line != false && Plots.plot!(p, t, trace, l=line)
    end
    # x limits
    t1 = minimum(s[:b] + d)
    t2 = maximum(s[:e] + d)
    if !(tw[1] == tw[end] == nothing)
        tw[1] != nothing && (t1 = tw[1])
        tw[end] != nothing && (t2 = tw[end])
    end
    Plots.xlims!(p, t1, t2)
    # y limits
    if !(dw[1] == dw[end] == nothing)
        d1, d2 = Plots.ylims(p)
        dw[1] != nothing && (d1 = dw[1])
        dw[end] != nothing && (d2 = dw[end])
        Plots.ylims!(p, d1, d2)
    end
    reverse && Plots.yflip!(p)
    # Labels
    if label != nothing
        label_text = if label isa AbstractArray
            length(label) == length(s) ||
                throw(ArgumentError("`label` vector must be same length as number of traces"))
            string.(label)
        elseif label isa Symbol
            strip.(s[label])
        else
            throw(ArgumentError("Don't know how to plot labels for input: $label"))
        end
        t1, t2 = Plots.xlims(p)
        for i in 1:length(s)
            Plots.annotate!(p, t2, y_shift[i], Plots.text(label_text[i], 8, :left, :black),
                right_margin=10Plots.mm)
        end
    end
    Plots.plot!(p; kwargs...)
    p
end
plotrs(s::SACArray, args...; kwargs...) =
    plotrs!(Plots.plot(; sacplot_prs_defaults...), s, args...; kwargs...)
prs! = plotrs!
prs = plotrs
@doc (@doc plotrs!) plotrs

# Routines to calculate offsets in x and y for record section plots
_y_shifts(s::SACArray, y::Symbol) = s[y]
_y_shifts(s::SACArray, y::AbstractArray) =
    length(y) == length(s) ? y : error("Length of `y` not the same as number of traces")
_x_shifts(s::SACArray, t::Symbol) = -s[t]
_x_shifts(s::SACArray, t::AbstractArray) =
    length(s) == length(t) ? -t : error("Length of `align` not the same as number of traces")
_x_shifts(s::SACArray, t::Real) = -t*ones(length(s))

"""
    plotspec(s::SACtr; len=trace_length/20, overlap=4len/5, trace=true, trace_frac, f=identity, plots_kwargs...) -> ::Plots.Plot

Plot a spectrogram of the SAC trace `s` with windows of length `len` s which overlap by
`overlap` s.

By default, show the trace on a separate panel above the spectrogram; set `trace`
to `false` to only return the plot of the spectrogram.  `trace_frac` sets the proportion
of the plot to be taken up by the trace.  The power is trasformed by the function `f`
elementwise.  (For example, to show the natural log, pass `f=log`.)

`plots_kwargs...` are passed to the `Plots.heatmap` function which shows the spectrogram.
"""
function plotspec(s::SACtr;
        len=(s[:e] - s[:b])/20, overlap=0.8len,
        trace=true, trace_frac=0.3, f=identity, plot_kwargs...)
    n = round(Int, len/s[:delta])
    noverlap = clamp(round(Int, overlap/s[:delta]), 0, n-1)
    spec = DSP.spectrogram(s.t, n, noverlap, fs=1/s.delta)
    dspec = step(spec.time)/2
    colours = Plots.cgrad([:white,:blue,:green,:yellow,:orange,:red])
    pspec = Plots.heatmap(spec.time .+ s[:b], spec.freq, f.(spec.power),
        xlabel="Time / s", ylabel="Frequency / Hz", cbar=false, c=colours,
        framestyle=:box;
        plot_kwargs...)
    if trace
        p = Plots.plot(
            p1(cut(s, s[:b]+spec.time[1]-dspec, s[:b]+spec.time[end]+dspec), xlabel=""),
            pspec,
            layout=Plots.grid(2, 1, heights=(trace_frac, 1 - trace_frac)),
            xticks=nothing # Remove time from all x-axes
        )
        Plots.plot!(p[end], xticks=:auto) # Add time only to bottom subplot x-axis
    else
        pspec
    end
end

"""
    lims(::Array{SACtr}, relative=false) -> b, e, depmin, depmax

Get the minimum and maximum times, `b` and `e`, for an array of SAC traces.
If `relative` is true, these times are all relative to the beginning of the trace,
not the zero time.
"""
function lims(a::SACArray, relative=false)
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
    check_lims(a, routine="SACPlot.check_lims")

Throw an error if `a` is not a length-2 array with the first value lower
or equal to the second, if both are not NaN.
"""
function check_lims(a, routine="SACPlot.check_lims")
    a[1] > a[2] && # NB: (NaN {>,<,==} NaN) == false
        error(routine * ": Lower plot limit value must be less than upper")
end

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
qdp_skip(a::Union{SACtr,SACArray}, thresh::Integer=sacplot_qdp_thresh) =
    max(1, sum(a[:npts])÷thresh)
qdp_skip(a::Union{SACtr,SACArray}, thresh::Nothing) = qdp_skip(a)
qdp_skip(a::Union{SACtr,SACArray}, thresh::Bool) = thresh ? qdp_skip(a) : 1


end # module SACPlot

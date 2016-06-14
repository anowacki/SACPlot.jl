"""
SACPlot.jl provides routines for plotting SAC traces.
"""
module SACPlot
# Module for plotting SAC traces

using SAC
using PyPlot

import SAC.SACtr
import PyPlot.plot
# For @doc doc strings
if VERSION < v"0.4"
	import Docile.@doc
end

export
	plot1,
	p1,
	plot2,
	p2,
	plotpm,
	ppm,
	plotsp,
	psp

const TIME_PICKS = [:t0, :t1, :t2, :t3, :t4, :t5, :t6, :t7, :t8, :t9]
const NAME_PICKS = [:kt0, :kt1, :kt2, :kt3, :kt4, :kt5, :kt6, :kt7, :kt8, :kt9]

@doc """
`plot1(s::Array{SACtr}; xlim=[NaN, NaN], ylim=[NaN, NaN], label=:default, title="")`
`plot1(s::SACtr; args...)

Create a plot of the SAC trace(s) `s`.

Define limits in time with `xlim`

Define dependent variable axis limits with `ylim`, which can be a 2-array
of values, or \"all\" to set all axes to have the same automatic limits.

Define the text labels with an array of sumbols, which correspond to the names
of SAC headers.
""" ->
function plot1(a::Array{SACtr}; xlim=[NaN, NaN], ylim=[NaN, NaN], label=:default,
               title="")
	# Check arguments
	check_lims(xlim, "SACPlot.plot1")
	typeof(ylim) <: AbstractString || check_lims(ylim, "SACPlot.plot1")
    # Make sure any single label is in an array
    typeof(label) == Symbol && (label = [label])
    for l in label
        if l != :default
            l in [:kzdate, :kztime] && continue # Special labels
            l in fieldnames(a[1]) || error("Label header '$l' is not a valid header variable")
        end
    end
	# Make plots
	PyPlot.clf()
	n = length(a)
	b, e, depmin, depmax = lims(a)
	if ! isnan(xlim[1]); b = xlim[1]; end
	if ! isnan(xlim[end]); e = xlim[end]; end
	if ! (typeof(ylim) <: AbstractString)
		if ! isnan(ylim[1]); depmin = ylim[1]; end
		if ! isnan(ylim[end]); depmax = ylim[end]; end
	end
	# Turn off x labels for all but bottom trace
	for i = 1:n
		PyPlot.subplot(n, 1, i)
		PyPlot.plot(SAC.time(a[i]), a[i].t)
		PyPlot.xlim([b, e])
		if typeof(ylim) <: AbstractString
			if lowercase(ylim) == "all" PyPlot.ylim([depmin, depmax]) end
		elseif ! all(isnan(ylim))
			PyPlot.ylim([depmin, depmax])
		end
        y1, y2 = PyPlot.ylim()
        # Add picks
        for (tp, kp) in zip(TIME_PICKS, NAME_PICKS)
            t, k = getfield(a[i], tp), strip(getfield(a[i], kp))
            if t != SAC.sac_rnull
                PyPlot.plot([t, t], [y1, y2], "k-")
                k != SAC.sac_cnull &&
                    PyPlot.text(t, y1, k, horizontalalignment="left",
                               verticalalignment="bottom")
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
            PyPlot.text(b+0.97*(e-b), y1+0.95*(y2-y1), label_string,
                        horizontalalignment="right", verticalalignment="top",
                        fontsize=10)
        end
        if i == 1; PyPlot.title(title); end
	end
	ticks = PyPlot.xticks()
    PyPlot.subplots_adjust(hspace=0.)
	return
end

# Single-trace version
plot1(s::SACtr) = plot1([s])

p1 = plot1

@doc """
`plot2(::Array{SACtr})`

Plot all traces in array of SAC traces `a` on the same plot
""" ->
function plot2(a::Array{SACtr}; legend=false)
	PyPlot.clf()
	for i = 1:length(a)
		PyPlot.plot(SAC.time(a[i]), a[i].t, "")
	end
	b, e = lims(a)
	PyPlot.xlim([b, e])
	PyPlot.tight_layout()
	return
end

# Single-trace version
plot2(s::SACtr) = plot2([s])

p2 = plot2

@doc """
`plotpm(::Array{SACtr}; xlim=[NaN, NaN])`

Plot the particle motion for a pair of orthogonal components, within
the time window `xlim[1]` to `xlim[2]` if provided
""" ->
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
	PyPlot.clf()
	_, _, min, max = lims([t1, t2])
	PyPlot.plot(t1.t, t2.t)
	PyPlot.xlim(min, max)
	PyPlot.ylim(min, max)
	PyPlot.xlabel(strip(t1.kcmpnm) * " (" * string(t1.cmpaz) * ")")
	PyPlot.ylabel(strip(t2.kcmpnm) * " (" * string(t2.cmpaz) * ")")
	PyPlot.axis("square")
	PyPlot.tight_layout()
	return
end

plotpm(s1::SACtr, s2::SACtr; xlim=[NaN, NaN]) = plotpm([s1, s2]; xlim=xlim)
ppm = plotpm

@doc """
`plotsp(f::Array, S::Array, kind=\"amp\")`

Plot the Fourier-transformed trace `S`, with frequencies `f`.

`kind` may be one of:\n
`amp`  : Plot amplitude\n
`phase`: Plot phase\n
`real` : Plot real part\n
`imag` : Plot imaginary part\n
""" ->
function plotsp(f::Array{Array, 1}, S::Array{Array, 1}, kind="amp";
				xlim=[NaN, NaN], ylim=[NaN, NaN])
    error("`plotsp` is not implemented yet")
	check_lims(xlim, "SACPlot.plotsp")
	check_lims(ylim, "SACPlot.plotsp")
end
psp = plotsp

@doc """
`lims(::Array{SACtr}) -> b, e`

Get the minimum and maximum times, `b` and `e`, for an array of SAC traces
""" ->
function lims(a::Array{SACtr})
	n = length(a)
	b = a[1].b
	e = a[1].e
	depmin = a[1].depmin
	depmax = a[1].depmax
	for i = 2:n
		if a[i].b < b; b = a[i].b; end
		if a[i].e > e; e = a[i].e; end
		if a[i].depmin < depmin; depmin = a[i].depmin; end
		if a[i].depmax > depmax; depmax = a[i].depmax; end
	end
	return b, e, depmin, depmax
end

lims(s::SACtr) = lims([s])

"""
`check_lims(a, routine=\"SACPlot.check_lims\")`

Throw an error if `a` is not a length-2 array with the first value lower
or equal to the second, if both are not NaN.
"""
function check_lims(a, routine="SACPlot.check_lims")
	length(a) == 2 ||
		error(routine * ": Plot limits must be given as length-2 array.")
	a[1] > a[2] && # NB: (NaN {>,<,==} NaN) == false
		error(routine * ": Lower plot limit value must be less than upper")
end


end # module SACPlot

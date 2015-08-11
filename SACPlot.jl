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
	p2

@doc """
plot1(::Array{SACtr})

   Create a plot of the SAC trace(s) `s`
""" ->
function plot1(a::Array{SACtr}; xlim=[None, None], ylim=None)
	# Check arguments
	length(xlim) == 2 || error("xlim must be length-2 array of numbers")
	if xlim[1] != None && xlim[2] != None
		xlim[2] > xlim[1] || error("Lower xlim value must be less than upper")
	end
	# Make plots
	PyPlot.clf()
	n = length(a)
	b, e = lims(a)
	if xlim[1] != None; b = xlim[1]; end
	if xlim[2] != None; e = xlim[2]; end
	# Turn off x labels for all but bottom trace
	for i = 1:n
		PyPlot.subplot(n, 1, i)
		PyPlot.plot(SAC.time(a[i]), a[i].t)
		PyPlot.xlim([b, e])
	end
	ticks = PyPlot.xticks()
	# PyPlot.xticks(ticks[1])
	PyPlot.subplots_adjust(hspace=0.)
	return
end

# Single-trace version
plot1(s::SACtr) = plot1([s])

p1 = plot1

@doc """
plot2(::Array{SACtr})

	Plot all traces in array of SAC traces `a` on the same plot
""" ->
function plot2(a::Array{SACtr})
	PyPlot.clf()
	for i = 1:length(a)
		PyPlot.plot(SAC.time(a[i]), a[i].t)
	end
	b, e = lims(a)
	PyPlot.xlim([b, e])
	return
end

# Single-trace version
plot2(s::SACtr) = plot2([s])

p2 = plot2

@doc """
lims(::Array{SACtr}) -> b, e

	Get the minimum and maximum times, b and e, for an array of SAC traces
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

end # module SACPlot

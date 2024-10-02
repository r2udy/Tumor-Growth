cd(@__DIR__)
using Pkg; Pkg.activate("."); Pkg.instantiate()

# User inputs
datapath = "data/IF Files/Surface_growth_verifalgo_(x,t).mat"
savepath = "data/IF Files/"

xspace = :fourier #can be :cheb or :fourier
tspace = :cheb
fldnames = ["heights"] # ["rho"]
spline_degree = 5
## Load packages - all automatic from here :)
include("src/data2coefficients_v0.jl")
#include("plottingstartup.jl") #plotting I use
using .Transforms
using MAT

# Read in data
println("Reading in data")
datadict = matread(datapath)

# Grids
xs = datadict["x_heights"][:]
ts = datadict["t_heights"][:]

# Variables
heights = datadict[fldnames[1]]

# Set up the grids and smoothing parameters
println("Setting up grids")
if xspace === :cheb
    xdomain = (xs[1], xs[end])
    xspacetyp = ChebyshevSpace(1, xdomain)
elseif xspace === :fourier
    dx = xs[2] - xs[1]
    xdomain = (xs[1], xs[end]+dx)
    xspacetyp = FourierSpace(1, xdomain)
end

if tspace === :cheb
    tdomain = (ts[1], ts[end])
    tspacetyp = ChebyshevSpace(2, tdomain)
elseif tspace === :fourier
    dt = ts[2] - ts[1]
    tdomain = (ts[1], ts[end]+dt)
    tspacetyp = FourierSpace(2, tdomain)
end

spaces = (xspacetyp, tspacetyp)
inputgrid = (xs, ts)

# Construct representations
println("Calculating representations")
heightsrepresentation = Representation(heights, inputgrid, spaces, k = spline_degree)

# save
println("Saving representations")
saverepresentation(savepath*"_representation_$(xspace)_$(tspace)_$(fldnames[1])_Verif_algo_h.mat", heightsrepresentation)
println("Done :)")
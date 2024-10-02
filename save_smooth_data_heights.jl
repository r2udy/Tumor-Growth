cd(@__DIR__)
using Pkg; Pkg.activate("."); Pkg.instantiate()

# Inputs
datafolderpath = "data/IF Files/"
savefolderpath = "data/IF Files/"

xspace = :fourier #can be :cheb or :fourier
tspace = :cheb
fldnames = ["heights"] # ["px", "py", "rho"]
spline_degree = 5


derivdesc_id = ["", "t", "x"]
derivdesc = ((0, 0), (1, 0), (0, 1))

## Automatic from here
# Load relevant packages
using MAT
using Printf
# Include packages and decompression functions
include("src/data2coefficients_v0.jl")
using .Transforms

#matwrite saves a v7.3 .mat file which uses a HDF5 file system. This, when loaded in Python using h5py tranposes the matrix
#Hence, here, we first transpose the matrix before saving
function transposeMats(ls)
    #ls: list of arrays
    return [permutedims(m, ndims(m):-1:1) for m in ls]
end

# Load representations
println("Loading representations")
heightsrepresentation = readrepresentation(datafolderpath*"_representation_$(xspace)_$(tspace)_$(fldnames[1])_Verif_algo_h.mat")
xs, ts = heightsrepresentation.inputgrid

#coefficient thresholds
space_thr_list = [length(xs)]
time_thr_list = [length(ts)]

# Set up the threshold tups
println("Setting up thresholds")
if xspace == :cheb
    spacethrtups = map((xthr) -> (MaxThreshold(1, xthr)), space_thr_list)
elseif xspace == :fourier
    spacethrtups = map(thr -> (IntegerRadialThreshold(1, thr)), space_thr_list)
else
    error("What are you doing big dawg?")
end

if tspace == :cheb
    timethrtups = map(thr -> MaxThreshold(2, thr), time_thr_list)
else
    error("What are you doing big dawg?")
end

println("Beginning sweep ...")

@printf("space_thr = %i, time_thr = %i\n", space_thr_list[1], time_thr_list[1])
savepath = savefolderpath * "_reconstruction_$(xspace)_$(tspace)_space$(space_thr_list[1])_time$(time_thr_list[1])_Verif_algo_h.mat"

thrtup = (spacethrtups[1], timethrtups[1])
thrheightsrepresentation = threshold(heightsrepresentation, thrtup)

##
println("$(fldnames[1]) derivatives")
heights_derivdict = map(x -> reconstruct(thrheightsrepresentation, x, k = spline_degree), derivdesc)

##
heights_derivdict_tranpose = transposeMats(heights_derivdict)

#tranpose the grid vectors
x, t = transposeMats([xs, ts])

## save as .mat
println("saving data..")
savedict = Dict(Pair.("$(fldnames[1])_".*derivdesc_id, heights_derivdict_tranpose)..., "x"=>x, "t"=>t)
matwrite(savepath, savedict)
println("Saved !")
# Changing working directory

cd("/Users/victoralexandrino/Google Drive/PhD Insper/Thesis/Paper 1 - Sovereign Default Holdings/Quant/Learning Julia/")

# Verifying directory

pwd()

# Clear working space

clearconsole()

# Installing packages

# First, type ] add InstantiateFromURL in the REPL.

using Pkg
Pkg.add("Grisu")
Pkg.add("IJulia")
Pkg.add("GR")
Pkg.add("Gadfly")

using InstantiateFromURL
using IJulia
using Grisu
using GR
using Plots
using LinearAlgebra
using Statistics

randn()

using Gadfly
plot([sin, cos], 0, 25)

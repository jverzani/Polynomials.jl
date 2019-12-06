# assert file to test polynomial implementation
using Test
using LinearAlgebra
using Polynomials
using SpecialFunctions
using RecipesBase: apply_recipe

import SparseArrays: sparse, nnz


@testset "Polynomial" begin include("Polynomial.jl") end
@testset "ChebyshevT" begin include("ChebyshevT.jl") end
#@testset "Deprecations" begin include("deprecated.jl") end
@testset "Poly (deprecated)" begin include("Poly.jl") end

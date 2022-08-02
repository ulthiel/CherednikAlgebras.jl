# CherednikAlgebras.jl

A Julia package for things around (rational) Cherednik algebras (draft, under development).

By [Ulrich Thiel](https://ulthiel.com/math) (University of Kaiserslautern).

This package uses:

* [Gapjm.jl](https://github.com/jmichel7/Gapjm.jl) for everything related to complex reflection groups
* [OSCAR](https://oscar.computeralgebra.de) for everything else

## Motivation

I already have a larger Magma package [CHAMP](https://github.com/ulthiel/Champ) for rational Cherednik algebras. The motivation for this package is my observation that the arithmetic in Julia/OSCAR is faster than in Magma and I am able to solve some problems I have not been able to solve in Magma. This concerns especially the computation of Calogero–Moser cellular characters. So, a first milestone is to re-implement in Julia the algorithm for the computation is Calogero–Moser cellular characters. This is still under development.

## Installation

```julia
julia> using Pkg
julia> Pkg.add("url=https://github.com/jmichel7/Gapjm.jl")
julia> Pkg.add("https://github.com/ulthiel/CherednikAlgebras.jl")
```

## Usage

```julia
julia> W = complex_reflection_group(4)
G₄

julia> 𝒞 = ParameterSpace(W)
Parameter space of type EG for G₄ over Abelian closure of Q

julia> c = generic_point(𝒞)
Parameter of type EG for G₄ with values [c(1), c(2)]
        
# Gaudin operators about to come...
```




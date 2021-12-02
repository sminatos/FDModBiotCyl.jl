# FDModBiotCyl.jl
A Julia implementation of 2.5-D seismic wave simulation using a finite-difference method in a cylindrical coordinate system.

## Overview
FDModBiotCyl.jl simulates seismic wavefield by solving Biot dynamic poroelasticity equations using staggered-grid finite-difference method. It assumes azimuthally symmetric wavefields in a cylindrical coordinate system.

Particular emphasis is on wave simulation involving a test hole in rocks and soils (_Borehole Geophysics_). Critical remarks are:
 - Hydrophone vertical seismic profile: Waves in a borehole due to a plane wave incidence
 - Borehole acoustics: Waves in a borehole due to a seismic source in the borehole
 - Useful to simulate waves in any materials with azimuthal symmetry (boreholes, cylinders, and pipes)
 - Acoustic, Elastic, and Poroelastic wavefields

![](img/demo_vsp.gif)
![](img/demo_acoustic_log.gif)

The package is developed and distributed hoping that Julia's performance and MATLAB-like intuitive syntax will promote this exciting field of study in both academia and industry.

## Requirements
This package requires Julia version 1.4 or later.

## Installation
First, clone the repository and move to the top-level directory:
```bash
git clone https://github.com/sminatos/FDModBiotCyl.jl.git --depth 1
cd FDModBiotCyl/
```
Then, in the Julia REPL:
```julia
julia> using Pkg
julia> Pkg.activate("./")
julia> Pkg.test()
```
If the test successfully passes, it is good to go.

## Quick Start
A sample code for simulation of a simple medium is available in [/examples/homogeneous_Elastic.jl](/examples/homogeneous_Elastic.jl). The code requires additional packages for plotting outputs. If your system does not have them:
```julia
julia> using Pkg
julia> Pkg.add("Plots")
julia> Pkg.add("Printf")
julia> Pkg.add("CPUTime")
julia> Pkg.add("ProgressMeter")
```
To run the code,
```julia
julia> include("./samples/homogeneous_Elastic.jl")
```

You can change medium properties and a source-receiver geometry by modifying the code. See [manual for /samples](/doc/manual.md) for more details.

## Examples

In `./samples` several examples for a simple medium are available ([manual for /samples](/doc/manual.md)). In `./samples_borehole` several examples for a rather complex borehole model are available ([manual for /samples_borehole](/doc/manual_samples_borehole.md)).

- [/examples/homogeneous_Elastic.jl](/examples/homogeneous_Elastic.jl) simulates elastic wavefield in a homogeneous medium due to a point source.

- [/examples/homogeneous_PoroElastic.jl](/examples/homogeneous_PoroElastic.jl) simulates poroelastic wavefield in a homogeneous medium due to a point source.

- [/samples_borehole/acoustic_logging.jl](/samples_borehole/acoustic_logging.jl) simulates pressure response in a water-filled borehole. The borehole is embedded in a homogeneous elastic medium, and a point source is located in the borehole water (monopole-source acoustic logging).

- [/samples_borehole/VSP_Elastic_2L.jl](/samples_borehole/VSP_Elastic_2L.jl) simulates pressure response in a water-filled borehole. The borehole is embedded in a layered elastic medium, and a plane P wave incident on the borehole from above (zero-offset vertical seismic profiling).

- [/samples_borehole/VSP_PoroElastic_3L.jl](/samples_borehole/VSP_PoroElastic_3L.jl) simulates pressure response in a water-filled borehole. The borehole is embedded in a layered elastic-poroelastic medium, and a plane P wave is incident from above (zero-offset vertical seismic profiling).

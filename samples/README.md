# /samples
This file explains sample files in `/samples`. They solve Biot poroelastic equations using the finite-difference method in the cylindrical coordinate system with azimuthal symmetry. The solid phase assumes isotropic elasticity. Please note that the sample files will use additional packages for visualization.

- [/samples/homogeneous_Elastic.jl](/samples/homogeneous_Elastic.jl) simulates elastic wavefield in a homogeneous medium due to a point source.

- [/samples/homogeneous_PoroElastic.jl](/samples/homogeneous_PoroElastic.jl) simulates poroelastic wavefield in a homogeneous medium due to a point source.

## Structure of the sample files
The sample files have the following structure. You can change values according to your simulation configuration.

1. Setting `dt`, `T`, and PML thicknesses (`LPML_r` and `LPML_z`)
2. Creating a model (material parameter distribution):
  - See [Material parameters](#material-parameters) below for the necessary material parameters.
  - Look at an example function, e.g., `makemodel_homogeneous_Elastic` in [/samples/homogeneous_Elastic.jl](/samples/homogeneous_Elastic.jl) for more details.
3. Creating a source wavelet `src_func`
4. Defining a source location and its amplitude-scaling factor:
  - `src_index` and `src_dn`
5. Initializing field variables by `init_fields_Por`
6. Creating a PML profile by `init_PML_profile`
7. Defining a receiver geometry
8. Setting snapshots parameters by `init_snap`
9. Running the FD main loop `main_loop!`:
  - See [FD Main Loop](#fd-main-loop) below for more details.

## Material parameters
Following parameters are matrices of the size `(nz, nr)`, where `nr` is the number of grid points in a radial direction `r`, and `nz` is that in a vertical direction `z`.
 - `Rho` : Bulk density
 - `Rhof` : Fluid density
 - `H`, `C`, `M` : Poroelastic moduli (see, e.g., [Guan and Hu, 2011; Minato et al., 2021](#references))
 - `G` : Shear modulus
 - `D1`, `D2` : Material parameters relevant to fluid flow properties (see, e.g., [Guan and Hu, 2011; Minato et al., 2021](#references))
 - `Flag_AC`, `Flag_E` : These flags are used to indicate an acoustic medium or an elastic medium at each FD cell.

## Field variables
Following parameters are matrices of the size `(nz, nr)`.
- `vr`, `vz` : solid-phase particle velocities.
- `trr`, `tzz`, `tpp`, `trz` : solid-phase stress
- `pf`, `vfr`, `vfz` : fluid pressure and fluid relative velocities.

## Grid geometry
- `(z,r)`: z is a vertical direction, and `r` is a radial direction
- Material parameters are defined at `(z,r)`. They are constant in a cell with the four corners at `(z+dz/2,r+dr/2)`, `(z-dz/2,r+dr/2)`, `(z-dz/2,r-dr/2)`, and `(z+dz/2,r-dr/2)`.
- `(z,r)` is the center of a cell, where `tzz`, `trr`, `tpp`, and `pf` are defined.
- Index `(k,j)` corresponds to `tzz(z,r)`, `tpp(z,r)`, `trr(z,r)`, `pf(z,r)`, `vr(z,r+dr/2)`, `vfr(z,r+dr/2)`, `trz(z+dz/2,r+dr/2)`, `vz(z+dz/2,r)`, `vfz(z+dz/2,r)`
- Origin `(z,r)=(0,0)` or `(k,j)=(1,1)` is at the location of `tzz`, `trr`, `tpp`, and `pf`.
- The field parameters located at `r=0` are `tzz`, `trr`, `tpp`, `pf`, `vz`, and `vfz`.
- The field parameters located at `z=0` are `tzz`, `trr`, `tpp`, `pf`, `vr`, and `vfr`.

## FD Main Loop
The function `main_loop!` in the sample files in `/samples` calculates FD using given material parameters. It then populates a receiver response matrix (e.g., `rec_vz`). It may require modifying the main loop according to a source type (e.g., body force in `z` or `r` direction) and a receiver type (e.g., `pf` or `vz`) desired in your simulation. The calculation sequence in the loop is as follows.
1. Updating velocity field (main computation region and PML region)
2. Applying boundary conditions at the left (`r=0`) and the right (`r=R`) boundaries for the velocity field
3. Injecting source at the velocity field (optional)
4. Updating stress field (main computation region and PML region)
5. Applying boundary conditions at the left (`r=0`) and the right (`r=R`) boundaries for the stress field
6. Injecting source at the stress field (optional)
7. Saving field variables at receiver locations (optional)
8. Saving field variables for snapshots (optional)

At an acoustic domain (`Flag_AC`) and an elastic domain (`Flag_E`), it does not evaluate `vfr` and `vfz`, which are zeros. When the entire modeling domain is a poroelastic formation, use `Flag_AC=zeros(nz,nr)` and `Flag_E=zeros(nz,nr)`. When the entire modeling domain is an elastic formation, use `Flag_AC=zeros(nz,nr)` and `Flag_E=ones(nz,nr)`.

## References
- Randall et al. (1991), Geophysics, 56, 1757-1769
- Peng (1994), Ph.D. Thesis, Massachusetts Institute of Technology, http://hdl.handle.net/1721.1/12218
- Mittet and Renlie (1996), Geophysics, 61, 21-33
- Guan and Hu (2011), Comun. Comput. Phys., doi: 10.4208/cicp.020810.161210a
- Ou and Wang (2019), Geophys. J. Int., doi: 10.1093/gji/ggz144
- Minato et al. (2021), arXiv:2112.03410 [physics.geo-ph] (available at http://arxiv.org/abs/2112.03410).

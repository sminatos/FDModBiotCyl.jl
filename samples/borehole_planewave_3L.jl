#testing FDTD in cylindrical coordinate
#Randall et al (1991)


#01: start. homogeneous
#02-04: module. CylFDMOD.jl
#05:testing large model before ABC
#06B:  speedup. functioned.
#07: receiver trace test.
#08: testing Interpolations
#  : testing 3dFD and DWI
#09: testing Hicks interpolated Source -> good.
#10: functionizing to speedup
#11: testing Liao's ABC
#12: testing Cerjan's sponge ABC
#13: Heterogeneity test.
#14: Snapshots test
#14B: Figures use Plots.jl+pyplot backend
#15: testing borehole model
#16: testing m=1 (assuming water @r=0, or trp=0 @r=0)
#16B: explicit imput of "m" (CylFDMod03.jl)
#16C: careful implementation of dipole source on r=0 --> not yet
#17: test functionizing main loop
#-------------Major change-----
#18: implementing 2nd order FD (Graves)
#18B: better convention
#19: minor modification (either Liao OR Cerjan ABC is better?)
#Mittet01: field variable locations are from Mittet (1996)
#02: testing all l'Hopital rule for LB (if reducing noise on axis?->not yet)-> no effect?
#  : setting back to original (same as Mittet01).
#  : "noise" was disappeared when dipole source is a single point on axis (j=1). why...?
#Nojima01: loading model created by makemodel_Nojima01.jl
#Nojima02: minor bug fix
#Nojima03: model building using module -> becomes very slow. Checking which varibales should be "const"
#        : "m","nt","dt","T"
#Nojima03: minor bug fixes due to Global/Local scopes (minor changes in modules too)

#planewave01: testing plane wave incidence
#03: B.C. of plane wave solution from Schoenberg. Test.

#OWS01: modeling OWS source (monopole)

#PML01: test implementing NPML of Wang and Tang (2003)
#PML03: CylFDMod_Mittet05tmp3.jl. bug fix about LEft BCs.
#PML04: Right PML
#PML05: checking why bottom and top is non synmetric...Good.
# temporary good version.
#PML06: implementing TopRight...Done.
#PML07: checking furhter r direction..
#PML08: make it work with the latest module (CylFDMod_Mittet06tmp6.jl)...done.
# --------------
# major change
# --------------
# Based from CylFD_PML08.jl and CylFDMod_Mittet06tmp6.jl
# --PcylFD -> Biot's poroelastic media.
# --Main -> Ou 2019, GJI, doi: 10.1093/gji/ggz144. Sub->Sidler, 2014, GJI, doi: 10.1093/gji/ggt447
# Parameters definitions are better in Sidler 2014
# --radially invariant source only (m=0)
# PcylFD_test01: no poroelastic, homogeneous + no PML, done.
# PcylFD_test02: start implementing poroelasticitiy (PCylFDMod_Ou01.jl) w/o PML looks ok.
# PcylFD_test02: start PML
# PcylFD_test04: ./PCylFDMod_Ou02.jl
# Homogeneous test looks OK!
# PcylFD_test05: start implementing borehole environment: ./PCylFDMod_Ou03.jl
# PcylFD_test06: switch to 1st order FD: ./PCylFDMod_Ou04.jl
#              : Guan's Fig3 looks OK! (copymat is important!!)
# PcylFD_planewave01: plane wave incidence (taken from CylFD_planewave_Nojima01_tmpsimple06.jl)
# PcylFD_planewave02: additional receivers along permeable layer (vfr,pf)
#                  : Checking if pf~-B*sii/3 ... done. Note: this relation required A. undrained condition, and B. zero-initial condition.
#                  : B is skempton coefficient. pf is pore pressure. sii/3 is mean "total" stress
#                  : Important-> In this FD implementation, tau are all "total" stress. vf are all "relative fluid velocity"
#                  : I keep this code so as to keep the setting for testing pf~-B*sii/3
# PcylFD_planewave03: coming back to test thin porous layer

#--Geometry convention---
# See Randall et al (1991), Geophysics, Multipole borehole acoustic waveform: Synethtic logs with beds and borehole washouts
# Mittet and Renlie (1996), Geophysics, High-order, finite-difference modeling of multipole logging in formations with anisotropic attenuation and elasticity, 61, 21-33.
# (z,r): z vertical downward
#
#
#
# Vp,Vs,Rho are defined at (z,r), and constant within the cell of (z,r),(z+dz,r),(z,r+dr),(z+dz,r+dr)
# --> (z,r) is the upper left corner of the grid
#
#--Redifining material parameters for Biot poroelasticity--
# H, C, M, mu, rhof, rho, k, eta, phi : See Ou and Wang (2019, doi: 10.1093/gji/ggz144
# will be dependent on --> D1(=0 when k(w)=k0), D2, rho, rhof, M, C, H, mu
#
#--Additional field variables for Biot poroelasticitiy--
# pf, vwr, vwz : fluid pressure, vr and vz






#--
# Threads.nthreads() displays number of available threads
# if ==1, start julia with JULIA_NUM_THREADS=4 julia
#--

#--
# If no graphics (X-windows) required, then ENV["GKSwstype"]="nul"
#--

#
using CPUTime
using ProgressMeter
#using Interpolations
using MAT
#using Base64
#using DSP
using Plots
using JLD
using DelimitedFiles
using Printf
#
using FDModBiotCyl

#---Loading Functions----
#include("./PCylFDMod_Ou1st03.jl") #main modules. With acoustic, not good, but ok.
#include("./Peng06Mod.jl") #Peng's exact solution to B.C.s of plane wave (MOD version. Use this.)

#using Plots
#gr()

#--Suppressing graphics
#ENV["GKSwstype"]="nul"

@__DIR__
@__FILE__


#================
Drawing model
================#
function drawmodel(Field_Var,nz,dz,nr,dr,src_index,index_allrec,nrec,LPML_r,LPML_z)
   println("Model drawing...")
   zvec=[0:dz:(nz-1)*dz]
   rvec=[0:dr:(nr-1)*dr]
   R=(nr-1)*dr
#   plt1=heatmap(rvec,zvec,Vp,yflip=true,xlabel="R (m)",ylabel="Z (m)",ratio=1)
   plt1=heatmap(rvec,zvec,Field_Var,yflip=true,xlabel="R (m)",ylabel="Z (m)")
   plot!([(src_index[2]-1)*dr],[(src_index[1]-1)*dz],label="",markershape=:circle,markersize=5,markerstrokewidth=0,seriescolor=:red)
   plot!((index_allrec[:,2]-ones(nrec))*dr,(index_allrec[:,1]-ones(nrec))*dz,label="",markershape=:cross,markersize=5,markerstrokewidth=0,seriescolor=:red)
   plot!([0;R],ones(2)*(LPML_z-1)*dz,label="",color=:black)
   plot!([0;R],ones(2)*((nz-LPML_z+1)*dz),label="",color=:black)
   plot!(ones(2)*((nr-LPML_r+1)*dr),[0;(nz-1)*dz],label="",color=:black)
   display(plot(plt1))
   tmp_filename=@sprintf "tmp_FDmodel.png"
   tmp_filename_full=string(@__DIR__,"/",tmp_filename)
   png(tmp_filename_full)
   println("--------")
   println("Model config figure is saved in:", tmp_filename_full)
   println("--------")
end


function tmp_makemodel()
#poroelasticity version

   nr=301
   nz=701 #for back2



   dr=0.01 #A3-A5
   dz=0.2 #A3

   #--poroelastic parameters (Sidler's)
   #solid phase
   Km=zeros(nz,nr) #frame
   Ks=zeros(nz,nr) #grain
   G=zeros(nz,nr) #bulk
   Rho=zeros(nz,nr) #bulk (weighting average of Rhos and Rhof)
   Rhos=zeros(nz,nr) #grain

   #fluid phase
   Kf=zeros(nz,nr)
   Rhof=zeros(nz,nr)
   Kappa0=zeros(nz,nr) #static permeability (m^2, 1D=9.869E-13m^2)
   Eta=zeros(nz,nr) #dynamic viscosity Pa.s=kg/(m.s)

   #other parameters
   Phi=zeros(nz,nr)
   Tot=zeros(nz,nr) #Tortuosity factor (see Sidler 2014)
   #---end poroelastic parameters

   #-Fluid all equivalent
   Kf[:,:]=ones(nz,nr)*2.25*10^9
   Rhof[:,:]=ones(nz,nr)*1000
   Kappa0[:,:]=ones(nz,nr)*1*9.869*10^(-13) #m^2
   Eta[:,:]=ones(nz,nr)*1.0*10^(-3) #water 1E-3 Pa.s


   #--homogeneous media
   Phi[:,:]=ones(nz,nr)*0.3
   #--Tortuosity factor from Sidler
   #Tot=0.5*(ones(nz,nr)+1 ./Phi) #Empirical formula (see Sidler)
   #--Tortuosity factor from Ou
   m_Ou=8.0 #
   alpha_inf_Ou=3.0 #tortuosity
   Tot=ones(nz,nr)*(1+2/m_Ou)*alpha_inf_Ou

   #-Solid parameters from dry_Vp,dry_Vs,grain_density, and grain_bulkmodulus (Table 1 of Ou2019)
   #
   #== -Frame 1: Vpsat=4000,Vssat=2000,RhoB=2500
   Vp_dry=4165.0
   Vs_dry=2135.0
   Rhos[:,:]=ones(nz,nr)*3143.0
   Ks[:,:]=ones(nz,nr)*5*1E10
   ==#
   #== -Frame 2: Vpsat=5011,Vssat=3002,RhoB=2500
   Vp_dry=5200.0
   Vs_dry=3200.0
   Rhos[:,:]=ones(nz,nr)*3143.0
   Ks[:,:]=ones(nz,nr)*9*1E10
   ==#
   # -Frame 2B: Vpsat=5000,Vssat=3000,RhoB=2500
   Vp_dry=5170.0
   Vs_dry=3198.0
   Rhos[:,:]=ones(nz,nr)*3143.0
   Ks[:,:]=ones(nz,nr)*10*1E10
   #


   Rho_dry=(ones(nz,nr)-Phi).*Rhos
   G[:,:]=Vs_dry^2*Rho_dry[:,:]
   Km[:,:]=Vp_dry^2*Rho_dry[:,:]-4/3*G[:,:]
   Rho=(ones(nz,nr)-Phi).*Rhos+Phi.*Rhof
   #





   #--additional borehole (acoustic media)
   Flag_AC=zeros(nz,nr)
   #
   #ir_wall=Int(round(0.1/dr))
   ir_wall=6 #A3-A5
   for iz=1:nz
       Flag_AC[iz,1:ir_wall]=ones(ir_wall)
       Rho[iz,1:ir_wall]=ones(ir_wall)*Rhof[1,1] #Rho=Rhof
       #--
       G[iz,1:ir_wall]=ones(ir_wall)*0.0 #G=0
       Phi[iz,1:ir_wall]=ones(ir_wall) #Phi=1
       Km[iz,1:ir_wall]=ones(ir_wall)*0.0 #Km=0 (G=0,phi=0->M=C=H=Kf, Ou's)
       #--
#       Tot[iz,1:ir]=Phi[iz,1:ir] #Tot=Phi (-> D1=rhof, Ou's)
#       Tot[iz,1:ir]=Phi[iz,1:ir] #Tot=Phi (-> D1=1, Guan's)
       Eta[iz,1:ir_wall]=ones(ir_wall)*0.0 #eta=0 (-> D2=0, Ou's)
   end
   #




   #3Layer model (Elastic-PoroElastic-Elastic)
   #boundary @iz=271, L~1m or k=5samples?
   #== --additional elastic media
   # background 1
   Vp_elastic=4000.0 #Or, Specify K_elastic
#   K_elastic=1.12*1E10 #Or, Specify Vp_elastic
   Vs_elastic=2000.0
   Rho_elastic=2500.0
   ==#

   # background 2
   Flag_E=zeros(nz,nr)
   Vp_elastic=5000.0 #Or, Specify K_elastic
   Vs_elastic=3000.0
   Rho_elastic=2500.0
   #

#
   for iz=1:270 #Upper Elastic media

       G_elastic=Vs_elastic^2*Rho_elastic
       K_elastic=Vp_elastic^2*Rho_elastic-4/3*G_elastic #Activate here if you have specified Vp_elastic

       Flag_E[iz,ir_wall+1:end]=ones(nr-ir_wall)
       Rho[iz,ir_wall+1:end]=ones(nr-ir_wall)*Rho_elastic #Rho=Rho_elastic
       Rhof[iz,ir_wall+1:end]=ones(nr-ir_wall)*Rho_elastic #Rhof=Rho_elastic
       G[iz,ir_wall+1:end]=ones(nr-ir_wall)*G_elastic #G=Gs (Ou's)
       #--
       Phi[iz,ir_wall+1:end]=ones(nr-ir_wall)*0.0 #Phi=0
       Km[iz,ir_wall+1:end]=ones(nr-ir_wall)*K_elastic #Km
       Ks[iz,ir_wall+1:end]=ones(nr-ir_wall)*K_elastic #Ks=Km (-> M=Inf, Ou's)
       #--
       Kappa0[iz,ir_wall+1:end]=ones(nr-ir_wall)*0.0 #k0=0 (-> D1=D2=Inf, Ou's)
   end

#==
   for iz=276:nz #Lower Elastic media

       G_elastic=Vs_elastic^2*Rho_elastic
       K_elastic=Vp_elastic^2*Rho_elastic-4/3*G_elastic #Activate here if you have specified Vp_elastic

       Flag_E[iz,ir_wall+1:end]=ones(nr-ir_wall)
       Rho[iz,ir_wall+1:end]=ones(nr-ir_wall)*Rho_elastic #Rho=Rho_elastic
       Rhof[iz,ir_wall+1:end]=ones(nr-ir_wall)*Rho_elastic #Rhof=Rho_elastic
       G[iz,ir_wall+1:end]=ones(nr-ir_wall)*G_elastic #G=Gs (Ou's)
       #--
       Phi[iz,ir_wall+1:end]=ones(nr-ir_wall)*0.0 #Phi=0
       Km[iz,ir_wall+1:end]=ones(nr-ir_wall)*K_elastic #Km
       Ks[iz,ir_wall+1:end]=ones(nr-ir_wall)*K_elastic #Ks=Km (-> M=Inf, Ou's)
       #--
       Kappa0[iz,ir_wall+1:end]=ones(nr-ir_wall)*0.0 #k0=0 (-> D1=D2=Inf, Ou's)
   end
==#
#

#==TMPTMP everywhere elastic TMPTMP
   K_elastic=1.12*1E10 #Or, Specify Vp_elastic
   Vs_elastic=2000.0
#   Vs_elastic=1500.0
   Rho_elastic=2800.0

   for iz=501:nz
   #   for iz=1:nz

       G_elastic=Vs_elastic^2*Rho_elastic
   #       K_elastic=Vp_elastic^2*Rho_elastic-4/3*G_elastic #Activate here if you have specified Vp_elastic

       Flag_E[iz,ir+1:end]=ones(nr-ir)
       Rho[iz,ir+1:end]=ones(nr-ir)*Rho_elastic #Rho=Rho_elastic
       Rhof[iz,ir+1:end]=ones(nr-ir)*Rho_elastic #Rhof=Rho_elastic
       G[iz,ir+1:end]=ones(nr-ir)*G_elastic #G=Gs (Ou's)
       #--
       Phi[iz,ir+1:end]=ones(nr-ir)*0.0 #Phi=0
       Km[iz,ir+1:end]=ones(nr-ir)*K_elastic #Km
       Ks[iz,ir+1:end]=ones(nr-ir)*K_elastic #Ks=Km (-> M=Inf, Ou's)
       #--
       Kappa0[iz,ir+1:end]=ones(nr-ir)*0.0 #k0=0 (-> D1=D2=Inf, Ou's)
   end
==# #TMPTMP everywhere elastic TMPTMP


   #--Converting Sidler's param into Ou's param
   #Rho, Rhof, M, C, H, G, D1, D2
   #Rho as is
   #Rhof as is
   Alpha=ones(nz,nr)-Km./Ks
   M=( (Alpha-Phi)./Ks+Phi./Kf ).^(-1)
   C=M.*Alpha
   H=Km+4/3*G+M.*Alpha.^2
   #G as is
   D1=Tot.*Rhof./Phi #(1+2/m)*alpha_inf*rhof/phi
   D2=Eta./Kappa0
#   D1=ones(nz,nr)*(1+2/8)*3*1000/0.2
#   D2=ones(nz,nr)*1E-3/(1*1E-12)

   #correction of H and C for elastic media
   for ir=1:nr
      for iz=1:nz
         if(Flag_E[iz,ir]==1)
            H[iz,ir]=Ks[iz,ir]+4/3*G[iz,ir]
            C[iz,ir]=Ks[iz,ir] #Ou's
#            C[iz,ir]=0.0 #C does not matter!
         end

         if(Flag_AC[iz,ir]==1)
#            D1[iz,ir]=1.0 #Guan's
            D1[iz,ir]=1000. #Ou's D1=rhof
            D2[iz,ir]=0. #Ou's D2=0
         end
      end
   end


   return nr,nz,dr,dz,Rho,Rhof,M,C,H,G,D1,D2,Flag_AC,Flag_E,ir_wall
end



function drawmodel_tmp(Data1,Data2,nz,dz,nr,dr,src_index,index_allrec_tii,nrec,LPML_r,LPML_z)
   println("Model drawing...")
   zvec=[0:dz:(nz-1)*dz]
   rvec=[0:dr:(nr-1)*dr]
   R=(nr-1)*dr
#   plt1=heatmap(rvec,zvec,Vp,yflip=true,xlabel="R (m)",ylabel="Z (m)",ratio=1)
   plt1=heatmap(rvec,zvec,Data1,yflip=true,xlabel="R (m)",ylabel="Z (m)",legend=:none)
#   plot!([(src_index[2]-1)*dr],[(src_index[1]-1)*dz],label="",markershape=:circle,markersize=5,markerstrokewidth=0,seriescolor=:red)
#   plot!((index_allrec_tii[:,2]-ones(nrec))*dr,(index_allrec_tii[:,1]-ones(nrec))*dz,label="",markershape=:cross,markersize=5,markerstrokewidth=0,seriescolor=:red)
   plot!([0;R],ones(2)*(LPML_z-1)*dz,label="",color=:black)
   plot!([0;R],ones(2)*((nz-LPML_z+1)*dz),label="",color=:black)
   plot!(ones(2)*((nr-LPML_r+1)*dr),[0;(nz-1)*dz],label="",color=:black)
   plt2=heatmap(rvec,zvec,Data2,yflip=true,xlabel="R (m)",ylabel="Z (m)",ylims=(50,60),legend=:none)
#   plot!([(src_index[2]-1)*dr],[(src_index[1]-1)*dz],label="",markershape=:circle,markersize=5,markerstrokewidth=0,seriescolor=:red)
#   plot!((index_allrec_tii[:,2]-ones(nrec))*dr,(index_allrec_tii[:,1]-ones(nrec))*dz,label="",markershape=:cross,markersize=5,markerstrokewidth=0,seriescolor=:red)
   plot!([0;R],ones(2)*(LPML_z-1)*dz,label="",color=:black)
   plot!([0;R],ones(2)*((nz-LPML_z+1)*dz),label="",color=:black)
   plot!(ones(2)*((nr-LPML_r+1)*dr),[0;(nz-1)*dz],label="",color=:black)


   display(plot(plt1,plt2,layout=(1,2)))

end

function drawmodel_tmp2(Data1,Data2,Data3,Data4,nz,dz,nr,dr,src_index,index_allrec_tii,nrec,LPML_r,LPML_z)
   println("Model drawing...")
   zvec=[0:dz:(nz-1)*dz]
   rvec=[0:dr:(nr-1)*dr]
   R=(nr-1)*dr
#   plt1=heatmap(rvec,zvec,Vp,yflip=true,xlabel="R (m)",ylabel="Z (m)",ratio=1)
   plt1=heatmap(rvec,zvec,Data1,yflip=true,xlabel="R (m)",ylabel="Z (m)",legend=:none)
   plot!([0;R],ones(2)*(LPML_z-1)*dz,label="",color=:black)
   plot!([0;R],ones(2)*((nz-LPML_z+1)*dz),label="",color=:black)
   plot!(ones(2)*((nr-LPML_r+1)*dr),[0;(nz-1)*dz],label="",color=:black)

   plt2=heatmap(rvec,zvec,Data2,yflip=true,xlabel="R (m)",ylabel="Z (m)",ylims=(50,60),legend=:none)
   plot!([0;R],ones(2)*(LPML_z-1)*dz,label="",color=:black)
   plot!([0;R],ones(2)*((nz-LPML_z+1)*dz),label="",color=:black)
   plot!(ones(2)*((nr-LPML_r+1)*dr),[0;(nz-1)*dz],label="",color=:black)

   plt3=heatmap(rvec,zvec,Data3,yflip=true,xlabel="R (m)",ylabel="Z (m)",ylims=(50,60),legend=:none)
   plot!([0;R],ones(2)*(LPML_z-1)*dz,label="",color=:black)
   plot!([0;R],ones(2)*((nz-LPML_z+1)*dz),label="",color=:black)
   plot!(ones(2)*((nr-LPML_r+1)*dr),[0;(nz-1)*dz],label="",color=:black)

   plt4=heatmap(rvec,zvec,Data4,yflip=true,xlabel="R (m)",ylabel="Z (m)",ylims=(50,60),legend=:none)
   plot!([0;R],ones(2)*(LPML_z-1)*dz,label="",color=:black)
   plot!([0;R],ones(2)*((nz-LPML_z+1)*dz),label="",color=:black)
   plot!(ones(2)*((nr-LPML_r+1)*dr),[0;(nz-1)*dz],label="",color=:black)

   display(plot(plt1,plt2,plt3,plt4,layout=(2,2)))

end


#---function to create receiver geometry for plane wave experiment ()temporary
function make_receivers()
   #--Receiver positions----
   #==
      nrec=8
      recgeom=zeros(nrec,2) #(z,r)
      for irec=1:nrec
         recgeom[irec,1]=srcgeom[1,1]-2.7432-(irec-1)*0.1524
         recgeom[irec,2]=0.0
      end
   ==#
   #---borehole center
   nrec1=nz
   recgeom1=zeros(nrec1,2) #(z,r)
   for irec=1:nrec1
      recgeom1[irec,1]=dz*(irec-1)
      recgeom1[irec,2]=0.0
   end

   #--Receivers without interpolation
   rec_vr1,rec_vz1,rec_tii1,index_allrec_vr1,index_allrec_vz1,index_allrec_tii1=init_receiver(recgeom1,nrec1,dr,dz,nr,nz,nt)

   #---corresponding borehole wall
   # exact boundrary is at (ir_wall-1)*dr+dr/2
   # coincident location: vr and trz
   # dr/2 inside elastci media: vz,tii,vphi,tpz
   recgeom2=zeros(nrec1,2) #(z,r)
   #--collect indexes for vr
   #---homogeneous R version
   for irec=1:nrec1
      recgeom2[irec,1]=dz*(irec-1)
      recgeom2[irec,2]=(ir_wall-1)*dr+dr/2
   end
   #---inhomogeneous R version (assuming Vp[1:ir_wall@z]=1500)
#   for irec=1:nrec1
#      ir_wall_now=maximum(findall(x->x==1500,Vp[irec,:])) #assuming nrec=nz
#      recgeom2[irec,1]=dz*(irec-1)
#      recgeom2[irec,2]=(ir_wall_now-1)*dr+dr/2 #this is a position of wall (good for vr)
#   end

   #--Receivers without interpolation
   rec_vr,rec_vz,rec_tii,index_allrec_vr2,index_allrec_vz_dummy,index_allrec_tii_dummy=init_receiver(recgeom2,nrec1,dr,dz,nr,nz,nt)
   # repeat for tii
   for irec=1:nrec1
      recgeom2[irec,1]=dz*(irec-1)
      recgeom2[irec,2]=(ir_wall-1)*dr+dr
   end
#   for irec=1:nrec1
#      ir_wall_now=maximum(findall(x->x==1500,Vp[irec,:])) #assuming nrec=nz
#      recgeom2[irec,1]=dz*(irec-1) #this is a position good for tii (same as Vp)
#      recgeom2[irec,2]=(ir_wall_now-1)*dr+dr #this is a position of wall+dr/2 (good for vz and tii)
#   end
   #--Receivers without interpolation
#   rec_vr,rec_vz,rec_tii,index_allrec_vr_dummy,index_allrec_vz2,index_allrec_tii2=init_receiver(recgeom2,nrec1,dr,dz,nr,nz,nt)
   rec_vr,rec_vz,rec_tii,index_allrec_vr_dummy,index_allrec_vz_dummy,index_allrec_tii2=init_receiver(recgeom2,nrec1,dr,dz,nr,nz,nt)

   # repeat for vz (be careful Vp is defined at tii and vr, and vz is located in between. Where to define "solid phase")
   Vp_E=(H./Rho).^0.5
   Vs_E=(G./Rho).^0.5
   for irec=1:nrec1-1
      ir_wall_up=maximum(findall(x->x==1500,Vp_E[irec,:])) #assuming nrec=nz
      ir_wall_down=maximum(findall(x->x==1500,Vp_E[irec+1,:])) #assuming nrec=nz
      if (ir_wall_up==ir_wall_down)
         recgeom2[irec,1]=dz*(irec-1)+dz/2 #this is a position good for vz
         recgeom2[irec,2]=(ir_wall_up-1)*dr+dr #this is a position of wall+dr/2 (good for vz and tii)
      elseif (ir_wall_up>ir_wall_down) #radius gets smaller
         recgeom2[irec,1]=dz*(irec-1)+dz/2 #this is a position good for vz
         recgeom2[irec,2]=(ir_wall_up-1)*dr+dr #
      elseif (ir_wall_up<ir_wall_down) #radius gets larger
         recgeom2[irec,1]=dz*(irec-1)+dz/2 #this is a position good for vz
         recgeom2[irec,2]=(ir_wall_down-1)*dr+dr #
      else
         error("ERROR in RECGEOM!")
      end

   end
   #--Receivers without interpolation
   rec_vr,rec_vz,rec_tii,index_allrec_vr_dummy,index_allrec_vz2,index_allrec_tii_dummy=init_receiver(recgeom2,nrec1,dr,dz,nr,nz,nt)


   #--initializing arrays with correct nrec
   nrec=nrec1*2
   rec_vr=zeros(nt,nrec)
   rec_vz=zeros(nt,nrec)
   rec_tii=zeros(nt,nrec)
   #---merge indexes
   index_allrec_vr=zeros(nrec,2) #z,r
   index_allrec_vz=zeros(nrec,2) #z,r
   index_allrec_tii=zeros(nrec,2) #z,r
   index_allrec_vr[1:nrec1,:]=index_allrec_vr1[1:nrec1,:]
   index_allrec_vr[nrec1+1:2*nrec1,:]=index_allrec_vr2[1:nrec1,:]
   index_allrec_vz[1:nrec1,:]=index_allrec_vz1[1:nrec1,:]
   index_allrec_vz[nrec1+1:2*nrec1,:]=index_allrec_vz2[1:nrec1,:]
   index_allrec_tii[1:nrec1,:]=index_allrec_tii1[1:nrec1,:]
   index_allrec_tii[nrec1+1:2*nrec1,:]=index_allrec_tii2[1:nrec1,:]

   index_allrec_vr=map(Int,index_allrec_vr)
   index_allrec_vz=map(Int,index_allrec_vz)
   index_allrec_tii=map(Int,index_allrec_tii)
   #--Receivers with Interpolations.jl (not in main loop yet)
   #rec_vr,rec_vz,rec_tzz,
   #wis_allrec_vr,wis_allrec_vz,wis_allrec_tzz=init_receiver_interp(recgeom,nrec,dr,dz,nr,nz,nt)



   #--Receivers vfr and pf----
   nrec_PE=nr
   recgeom_PE=zeros(nrec_PE,2) #(z,r)
   for irec=1:nrec_PE
#         recgeom_PE[irec,1]=(273-1)*dz
      recgeom_PE[irec,1]=(320-1)*dz
      recgeom_PE[irec,2]=dr*(irec-1)
   end
   #--Receivers without interpolation
   rec_vfr_PE,rec_vfz_PE,rec_pf_PE,index_allrec_vfr_PE,index_allrec_vfz_PE,index_allrec_pf_PE=init_receiver(recgeom_PE,nrec_PE,dr,dz,nr,nz,nt)
   rec_tii_PE=zeros(size(rec_pf_PE)) #for mean normal stress

return  nrec,rec_vr,rec_vz,rec_tii,index_allrec_vr,index_allrec_vz,index_allrec_tii,
      nrec_PE,rec_vfr_PE,rec_pf_PE,rec_tii_PE,index_allrec_vfr_PE,index_allrec_pf_PE

end





#----MAIN loop function----
function main_loop!(nr,nz,dr,dz,Rho,Rhof,M,C,H,G,D1,D2,dt,nt,T,
   vr,vz,trr,tpp,tzz,trz,
   vfr,vfz,pf,
   ir_wall,
   Flag_AC,Flag_E,
   nrec,index_allrec_vr,index_allrec_vz,index_allrec_tii,rec_vr,rec_vz,rec_tii,
   nrec_PE,index_allrec_vfr_PE,index_allrec_pf_PE,rec_vfr_PE,rec_pf_PE,rec_tii_PE,
   LPML_r,LPML_z,PML_Wr,PML_Wz,PML_IWr,PML_Wr2,PML_Wz2,PML_IWr2,
   snapshots_vr,snapshots_vz,snapshots_trr,nsnap,itvec_snap,nskip)

#println("start")

#--
Flag_vf_zero=get_Flag_vf_zero(Flag_AC,Flag_E,nr,nz)
divV=zeros(nz,nr)
curV=zeros(nz,nr)

#initilize matrices of PML
Prz_T,Pzz_T,Rz_T,Srz_T,PzzPE_T,RzPE_T,
Prz_B,Pzz_B,Rz_B,Srz_B,PzzPE_B,RzPE_B,
Prr_R,Qrp_R,Pzr_R,Qzp_R,Rr_R,Rp_R,Rrz_R,PrrPE_R,RrPE_R,RpPE_R,
Prr_TR,Qrp_TR,Prz_TR,Pzr_TR,Qzp_TR,Pzz_TR,Rr_TR,Rp_TR,Rz_TR,Rrz_TR,Srz_TR,PrrPE_TR,PzzPE_TR,RrPE_TR,RpPE_TR,RzPE_TR,
Prr_BR,Qrp_BR,Prz_BR,Pzr_BR,Qzp_BR,Pzz_BR,Rr_BR,Rp_BR,Rz_BR,Rrz_BR,Srz_BR,PrrPE_BR,PzzPE_BR,RrPE_BR,RpPE_BR,RzPE_BR=
init_PML_Por(LPML_r,LPML_z,nr,nz)

memT_vr,memT_vz,memT_trr,memT_tpp,memT_tzz,memT_trz,memT_vfr,memT_vfz,memT_pf,
memB_vr,memB_vz,memB_trr,memB_tpp,memB_tzz,memB_trz,memB_vfr,memB_vfz,memB_pf,
memR_vr,memR_vz,memR_trr,memR_tpp,memR_tzz,memR_trz,memR_vfr,memR_vfz,memR_pf=
init_memory_variables_Por(LPML_r,LPML_z,nr,nz)
#return

#required for BC
vr_old=zeros(nz,nr)
vz_old=zeros(nz,nr)
vfr_old=zeros(nz,nr)
vfz_old=zeros(nz,nr)

mean_tii=zeros(nz,nr)

# Making sure RightBC for stress (zero values are trp and trz)
ApplyBCRight_stress!(vr,vz, #Use with flag_zero=1 (see ApplyBCRight_stress1D01)
    trz,
    G,nr,nz,dr,dz,dt)

#error()
#@showprogress for ii=1:nt
@showprogress for ii=1:1001
#for ii=1:1



#----Time at (ii-1)*dt---(updating velocities)-----

#--PML: save velocity at previous step
PML_save_vel!(memT_vr,memT_vz,memT_vfr,memT_vfz,
              memB_vr,memB_vz,memB_vfr,memB_vfz,
              memR_vr,memR_vz,memR_vfr,memR_vfz,
              vr,vz,vfr,vfz,nr,nz,LPML_z,LPML_r)
#println("update vel")

#keep current values before updating velocity for Acoustic-Poroelastic BC later
mycopy_mat(vr,vr_old,nz,nr)
mycopy_mat(vz,vz_old,nz,nr)
mycopy_mat(vfr,vfr_old,nz,nr)
mycopy_mat(vfz,vfz_old,nz,nr)

# Main velocity update!
update_velocity_1st_Por!(vr,vz,trr,tpp,tzz,trz,vfr,vfz,pf,Rho,Rhof,D1,D2,Flag_vf_zero,nr,nz,dr,dz,dt,LPML_z,LPML_r)

#error()

#PML: update velocity
PML_update_vel!(vr,vz,trr,tpp,tzz,trz,vfr,vfz,pf,
              Rho,Rhof,D1,D2,Flag_vf_zero,nr,nz,dr,dz,dt,
              Prz_T,Pzz_T,PzzPE_T,
              Prz_B,Pzz_B,PzzPE_B,
              Prr_R,Qrp_R,Pzr_R,Qzp_R,Rr_R,Rp_R,Rrz_R,
              PrrPE_R,RrPE_R,RpPE_R,
              Prz_TR,Pzz_TR,PzzPE_TR,
              Prr_TR,Qrp_TR,Pzr_TR,Qzp_TR,Rr_TR,Rp_TR,Rrz_TR,
              PrrPE_TR,RrPE_TR,RpPE_TR,
              Prz_BR,Pzz_BR,PzzPE_BR,
              Prr_BR,Qrp_BR,Pzr_BR,Qzp_BR,Rr_BR,Rp_BR,Rrz_BR,
              PrrPE_BR,RrPE_BR,RpPE_BR,
              LPML_z,LPML_r)
#println("update vel Left")
#if(ii==1)
#   error()
#end

# Making sure RightBC for velocity (zero values are vr)
ApplyBCRight_vel!(vr,vz,
      trr,tpp,tzz,trz,
      vfr,vfz,pf,
      Rho,Rhof,D1,D2,Flag_vf_zero,nr,nz,dr,dz,dt)

#---velocity B.Cs-----
#Periodic Left Edge
ApplyBCLeft_vel!(vr,vz,trr,tpp,tzz,trz,vfr,vfz,pf,
                          Rho,Rhof,D1,D2,Flag_vf_zero,
                          nr,nz,dr,dz,dt,
                          Prz_T,Pzz_T,PzzPE_T,
                          Prz_B,Pzz_B,PzzPE_B,
                          LPML_z)

#---Fluid-PE BC @borehole wall (test): replacing velocity values
#ir_wall=10
update_vr_vfr_1st_vertical(vr,trr,tpp,tzz,trz,vfr,pf,
    vr_old,vfr_old,
    Rho,Rhof,D1,D2,nr,nz,dr,dz,dt,LPML_z,
    Prz_T,Prz_B,
    ir_wall,Flag_AC,Flag_E)


#PML: update memory variables for stress (Rx and Sxx) using velocity at two time steps
PML_update_memRS!(Rz_T,Srz_T,RzPE_T,
                  Rz_B,Srz_B,RzPE_B,
                  Rr_R,Rp_R,Rrz_R,RrPE_R,RpPE_R,
                  Rr_TR,Rp_TR,Rz_TR,Rrz_TR,Srz_TR,RrPE_TR,RpPE_TR,RzPE_TR,
                  Rr_BR,Rp_BR,Rz_BR,Rrz_BR,Srz_BR,RrPE_BR,RpPE_BR,RzPE_BR,
                  memT_vr,memT_vz,memT_vfr,memT_vfz,
                  memB_vr,memB_vz,memB_vfr,memB_vfz,
                  memR_vr,memR_vz,memR_vfr,memR_vfz,
                  vr,vz,vfr,vfz,nr,nz,dr,dz,dt,LPML_z,LPML_r,PML_Wz,PML_Wr,PML_IWr,PML_Wz2,PML_Wr2,PML_IWr2)

#--src injection (velocity)
#srcamp=src_func[ii]
#srcapply!(vz,src_index,src_dn,srcamp)
#srcapply!(vr,src_index,src_dn,srcamp)

#----Time at (ii-1)*dt+dt/2----(updating stress)---


#--PML: save stress at previous step
PML_save_stress!(memT_trr,memT_tpp,memT_tzz,memT_trz,memT_pf,
                 memB_trr,memB_tpp,memB_tzz,memB_trz,memB_pf,
                 memR_trr,memR_tpp,memR_tzz,memR_trz,memR_pf,
                 trr,tpp,tzz,trz,pf,nr,nz,LPML_z,LPML_r)

#println("update stress")

# Main stress update!
update_stress_1st_Por!(vr,vz,trr,tpp,tzz,trz,vfr,vfz,pf,
     M,C,H,G,nr,nz,dr,dz,dt,LPML_z,LPML_r)
#error()

#PML: stress update
PML_update_stress!(vr,vz,trr,tpp,tzz,trz,vfr,vfz,pf,
                  M,C,H,G,nr,nz,dr,dz,dt,
                  Rz_T,Srz_T,RzPE_T,
                  Rz_B,Srz_B,RzPE_B,
                  Rr_R,Rp_R,Rrz_R,RrPE_R,RpPE_R,
                  Rr_TR,Rp_TR,Rz_TR,Rrz_TR,Srz_TR,RrPE_TR,RpPE_TR,RzPE_TR,
                  Rr_BR,Rp_BR,Rz_BR,Rrz_BR,Srz_BR,RrPE_BR,RpPE_BR,RzPE_BR,
                  LPML_z,LPML_r)

# Making sure RightBC for stress (zero values are trp and trz)
ApplyBCRight_stress!(vr,vz, #Use with flag_zero=1 (see ApplyBCRight_stress1D01)
     trz,
     G,nr,nz,dr,dz,dt)

#---stress B.Cs
#println("update stress Left")
#ApplyBCLeft_stress_2nd!(vr,vphi,vz,trr,tpp,tzz,trp,trz,tpz,lmat,mmat,m,nr,nz,dr,dz)
ApplyBCLeft_stress!(vr,vz,trr,tpp,tzz,trz,vfr,vfz,pf,
                    M,C,H,G,nr,nz,dr,dz,dt,
                    Rz_T,Srz_T,RzPE_T,
                    Rz_B,Srz_B,RzPE_B,
                    LPML_z)
#return

#--src injection (stress:monopole/dipole src)
#srcamp=-src_func[ii]
#srcapply!(pf,src_index,src_dn,-srcamp)
#srcapply!(tzz,src_index,src_dn,srcamp)
#srcapply!(trr,src_index,src_dn,srcamp)
#srcapply!(tpp,src_index,src_dn,srcamp)


#---Additional BC when Flag_Acoustic/Flag_Elastic
ApplyBC_stress_AcousticMedia_TEST!(trr,tpp,tzz,trz,pf,Flag_AC,nr,nz)
ApplyBC_stress_ElasticMedia_Ou!(pf,Flag_E,nr,nz)


#--PML: update memory variables for velocity (Pxx and Qxx) using stress at two time steps
PML_update_memPQ!(Prz_T,Pzz_T,PzzPE_T,
                  Prz_B,Pzz_B,PzzPE_B,
                  Prr_R,Qrp_R,Pzr_R,Qzp_R,PrrPE_R,
                  Prr_TR,Qrp_TR,Prz_TR,Pzr_TR,Qzp_TR,Pzz_TR,PrrPE_TR,PzzPE_TR,
                  Prr_BR,Qrp_BR,Prz_BR,Pzr_BR,Qzp_BR,Pzz_BR,PrrPE_BR,PzzPE_BR,
                  memT_trr,memT_tpp,memT_tzz,memT_trz,memT_pf,
                  memB_trr,memB_tpp,memB_tzz,memB_trz,memB_pf,
                  memR_trr,memR_tpp,memR_tzz,memR_trz,memR_pf,
                  trr,tpp,tzz,trz,pf,nr,nz,dr,dz,dt,LPML_z,LPML_r,PML_Wz,PML_Wr,PML_IWr,PML_Wz2,PML_Wr2,PML_IWr2)

#-----End of updating field----

#println("receiver")
#-Receiver field (extraction)
getRecData_from_index!(vr,rec_vr,index_allrec_vr,nrec,ii)
getRecData_from_index!(vz,rec_vz,index_allrec_vz,nrec,ii)
getRecData_from_index!(trr,rec_tii,index_allrec_tii,nrec,ii)

#Along porous layer
getRecData_from_index!(vfr,rec_vfr_PE,index_allrec_vfr_PE,nrec_PE,ii)
getRecData_from_index!(pf,rec_pf_PE,index_allrec_pf_PE,nrec_PE,ii)
mean_tii=(trr+tpp+tzz)/3
getRecData_from_index!(mean_tii,rec_tii_PE,index_allrec_pf_PE,nrec_PE,ii)



#--Snapshots
#println("check snap")
check_snap=findall(x -> x==ii,itvec_snap)
if (length(check_snap)!=0)
   println("ii=",ii)
   cnt_snap=Int(check_snap[1])
   get_snapshots!(snapshots_vr,snapshots_vz,cnt_snap,vfr,vz,nr,nz)
   get_snapshots_t!(snapshots_trr,cnt_snap,pf,nr,nz)
#   display(heatmap(snapshots_vz[1:5:end,1:5:end,cnt_snap]))
#   display(heatmap(snapshots_vr[100:150,1:20,cnt_snap]))

#   display(heatmap(snapshots_trr[:,:,cnt_snap]))

#   drawmodel02(snapshots_vz[:,:,cnt_snap],nz,dz,nr,dr,
#            src_index,index_allrec_tii,nrec,LPML_r,LPML_z)

#   drawmodel02(snapshots_trr[:,:,cnt_snap],nz,dz,nr,dr,
#            src_index,index_allrec_tii,nrec,LPML_r,LPML_z)

#   drawmodel_tmp(snapshots_vz[:,:,cnt_snap],snapshots_trr[:,:,cnt_snap],nz,dz,nr,dr,
#            src_index,index_allrec_tii,nrec,LPML_r,LPML_z)
   drawmodel_tmp2(snapshots_vz[:,:,cnt_snap],snapshots_trr[:,:,cnt_snap],snapshots_vr[:,:,cnt_snap],snapshots_vz[:,:,cnt_snap],nz,dz,nr,dr,
            src_index,index_allrec_tii,nrec,LPML_r,LPML_z)

end

end #i Time

end #function
#------------------------------

#-----ENTRY POINT HERE----
CPUtic()
start=time()

#======================
General variables
======================#

#==time samples
const dt=0.125E-5
const nt=16001
const T=(nt-1)*dt
==#
dt=0.125E-5
nt=16001
T=(nt-1)*dt
tvec=range(0.0,T,length=nt) # range object (no memory allocation)
tvec=collect(tvec) # a vector

#PML thickness in samples
LPML_r=51
LPML_z=61

#======================
Creating model
======================#
nr,nz,dr,dz,Rho,Rhof,M,C,H,G,D1,D2,Flag_AC,Flag_E,ir_wall=tmp_makemodel()
#error()

#======================
Stability check
======================#
Vmax=maximum((H./Rho).^(0.5)) #for PML
check_stability01(dt,dr,Vmax,0)
check_stability01(dt,dz,Vmax,0)

#======================
Creating src wavelet
======================#
f0=200 #src Freq Ou
delay=1/f0*1.0
src_func=myricker2(tvec,f0,delay,2) #when using 2nd derivative Gaussian (good for DWI)
#src_func=myricker2(tvec,f0,delay,1) #when using 1st derivative Gaussian (Randall?)
#src_func=myricker2(tvec,f0,delay,3) #when using 3rd derivative Gaussian
tmp_maxamp=maximum(map(abs,src_func))
src_func=src_func/tmp_maxamp
display(plot(tvec,src_func[:],title="src function"))
tmpfilename_src=string(@__DIR__,"/src_function.png")
png(tmpfilename_src)
println("Source signature figure saved in ",tmpfilename_src)

#================================
Src depth for initial plane wave
================================#
srcdepth=28. #meter


#==============================
Initializig field variables
==============================#
vr,vz,trr,tpp,tzz,trz,
vfr,vfz,pf,
vz_srcBC,vr_srcBC,
trr_srcBC,tpp_srcBC,tzz_srcBC,trz_srcBC,
vz_init,vr_init,trr_init,tpp_init,tzz_init,
trz_init,
src_iz,noffset,
f0=initialize_planewave(nr,nz,dr,dz,ir_wall,Rho,H,G,nt,tvec,dt,src_func,srcdepth,f0,delay)

LPML_r,LPML_z,PML_Wr,PML_Wz,PML_IWr,PML_Wr2,PML_Wz2,PML_IWr2=init_PML_profile(LPML_r,LPML_z,Vmax,dr,dz,nr)
#PML_check(LPML_r,LPML_z,Vmax,dr,dz,f0)


#error()

#==============================
Receiver geometry (plane wave experiments)
==============================#
nrec,rec_vr,rec_vz,rec_tii,index_allrec_vr,index_allrec_vz,index_allrec_tii,
nrec_PE,rec_vfr_PE,rec_pf_PE,rec_tii_PE,index_allrec_vfr_PE,index_allrec_pf_PE=make_receivers()


#==============================
Check model
==============================#
src_index=[1 1]; #dummy
drawmodel(vz_init,nz,dz,nr,dr,src_index,index_allrec_tii,nrec,LPML_r,LPML_z)


#==============================
Snapshot settings
==============================#
nskip,snapshots_trr,snapshots_vr,snapshots_vz,
nsnap,itvec_snap=init_snap(nt,nz,nr,100)
#error()

#==============================
Start main FD Loop
==============================#

main_loop!(nr,nz,dr,dz,Rho,Rhof,M,C,H,G,D1,D2,dt,nt,T,
   vr_init,vz_init,trr_init,tpp_init,tzz_init,trz_init,
   vfr,vfz,pf,
   ir_wall,
   Flag_AC,Flag_E,
   nrec,index_allrec_vr,index_allrec_vz,index_allrec_tii,rec_vr,rec_vz,rec_tii,
   nrec_PE,index_allrec_vfr_PE,index_allrec_pf_PE,rec_vfr_PE,rec_pf_PE,rec_tii_PE,
   LPML_r,LPML_z,PML_Wr,PML_Wz,PML_IWr,PML_Wr2,PML_Wz2,PML_IWr2,
   snapshots_vr,snapshots_vz,snapshots_trr,nsnap,itvec_snap,nskip)


#finalizing
#Vp1D,Vs1D,Rho1D,Caliper1D,Caliper_ir1D=model1D_load_Nojima(tmpfilename_JLD)
#save_data(zsrc,tvec,nr,nz,dr,dz,m,dt,nt,T,rec_tii,src_func,index_allrec_tii,src_index,snapshots_trr,itvec_snap,Vp1D,Vs1D,Rho1D,Caliper1D,Caliper_ir1D)
#end

CPUtoc()
println("elapsed real time: ", round(time() - start;digits=3)," seconds")

error("Stopped w/o problem.")



   println("Model drawing...")
   zvec=[0:dz:(nz-1)*dz]
   rvec=[0:dr:(nr-1)*dr]
   R=(nr-1)*dr
   plt1=heatmap(rvec,zvec,vz,yflip=true,xlabel="R (m)",ylabel="Z (m)",ratio=1)
   plot!([(src_index[2]-1)*dr],[(src_index[1]-1)*dz],label="",markershape=:circle,markersize=5,markerstrokewidth=0,seriescolor=:red)
   plot!((index_allrec_tii[:,2]-ones(nrec))*dr,(index_allrec_tii[:,1]-ones(nrec))*dz,label="",markershape=:cross,markersize=5,markerstrokewidth=0,seriescolor=:red)
   plot!([0;R],ones(2)*(LPML_z-1)*dz,label="",color=:black)
   plot!([0;R],ones(2)*((nz-LPML_z+1)*dz),label="",color=:black)
   plot!(ones(2)*((nr-LPML_r+1)*dr),[0;(nz-1)*dz],label="",color=:black)
   display(plot(plt1))



tmp_filename="out.jld"
tmp_filename_full=string(@__DIR__,"/",tmp_filename)
save(tmp_filename_full,
"tvec",tvec,"nr",nr,"nz",nz,"dr",dr,"dz",dz,"m",m,"dt",dt,"nt",nt,"T",T,"rec_tii",rec_tii,"src_func",src_func,
"index_allrec_tii",index_allrec_tii,"src_index",src_index,"snapshots_trr",snapshots_trr,
"itvec_snap",itvec_snap)


#tmp_filename="out_reference_back2.mat"
tmp_filename="out_frame2B_back2_thick.mat"

iskip_rec=10
tvec_rec=tvec[1:iskip_rec:end]

tmp_filename_full=string(@__DIR__,"/",tmp_filename)
matwrite(tmp_filename_full, Dict(
"tvec"=>tvec,"nr"=>nr,"nz"=>nz,"dr"=>dr,"dz"=>dz,
"dt"=>dt,"nt"=>nt,"T"=>T,
"tvec_rec"=>tvec_rec,
"rec_tii"=>rec_tii[1:iskip_rec:end,:],
"rec_vz"=>rec_vz[1:iskip_rec:end,:],
#"rec_vr"=>rec_vr,
"rec_tii_PE"=>rec_tii_PE[1:iskip_rec:end,:],
"rec_pf_PE"=>rec_pf_PE[1:iskip_rec:end,:],
"src_func"=>src_func,
"index_allrec_tii"=>index_allrec_tii,
"index_allrec_vz"=>index_allrec_vz,
"index_allrec_vr"=>index_allrec_vr,
"index_allrec_pf_PE"=>index_allrec_pf_PE,
"src_index"=>src_index,
#"snapshots_trr"=>snapshots_trr,
"snapshots_pf"=>snapshots_trr,
#"snapshots_vr"=>snapshots_vr,
"snapshots_vz"=>snapshots_vz,
"itvec_snap"=>itvec_snap); compress = true)

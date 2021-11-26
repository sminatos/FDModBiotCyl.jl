#Poroelastic FDTD in the cylindrical coordinate system
#Guan and Hu (2011), Coomun. Comput. Phys., doi: 10.4208/cicp.020810.161210a
#Ou and Wang (2019), Geophys. J. Int., doi: 10.1093/gji/ggz144
#--Geometry convention---
# (z,r): z vertical downward
# Material parameters are defined at (z,r), and constant within the cell of (z+dz/2,r+dr/2),(z-dz/2,r+dr/2),(z-dz/2,r-dr/2),(z+dz/2,r-dr/2)
# --> (z,r) is the center of the cell, where tzz,trr,tpp,pf are defined.
#
# index (k,j) corresponds to:
# tzz(z,r),tpp(z,r),trr(z,r),pf(z,r),
# vr(z,r+dr/2),vfr(z,r+dr/2),
# trz(z+dz/2,r+dr/2),
# vz(z+dz/2,r),vfz(z+dz/2,r)
#
# Origin (z=0,r=0) or (k=1,j=1) is at the location of tzz, trr, tpp, and pf.
# The field parameters located at r=0 are tzz,trr,tpp,pf,vz,vfz.
# The field parameters located at z=0 are tzz,trr,tpp,pf,vr,vfr.
#
#--Redefining material parameters for Biot poroelasticity--
# H, C, M, mu, rhof, rho, k, eta, phi : See Ou and Wang (2019, doi: 10.1093/gji/ggz144
# will be dependent on --> D1(=0 when k(w)=k0), D2, rho, rhof, M, C, H, mu
#
#--Additional field variables for Biot poroelasticity--
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
using Dierckx

#--Suppressing graphics
#ENV["GKSwstype"]="nul"

@__DIR__
@__FILE__


#================
Drawing model
================#
function drawmodel(Field_Var,nz,dz,nr,dr,LPML_r,LPML_z)
   println("Model drawing...")
   zvec=[0:dz:(nz-1)*dz]
   rvec=[0:dr:(nr-1)*dr]
   R=(nr-1)*dr
#   plt1=heatmap(rvec,zvec,Vp,yflip=true,xlabel="R (m)",ylabel="Z (m)",ratio=1)
   plt1=heatmap(rvec,zvec,Field_Var,yflip=true,xlabel="R (m)",ylabel="Z (m)")
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

#=========================
Creating homogeneous model
=========================#
function makemodel_homogeneous_PoroElastic()
#Nessesary input matrices to main loop function are:
#Rho: Bulk density
#Rhof: Fluid density
#M,C,H: Poroelastic moduli
#G: Formation shear moduli
#D1,D2: Poroelastic moduli in Ou's formulation
#Flag_AC,Flag_E: flag specifying acoustic region or elastic region
#Other input parameters to the main loop function are:
#ir_wall: grid number in r direction where a borehole wall starts (single value)
#
# This function creates the above input parameters assuming a
# homogeneous poroelastic media

#Model size
   nr=251 #samples
   nz=351 #samples
   dr=0.5 #meter
   dz=0.5 #meter

   #==============================================
   Initializing poroelastic parameters (Sidler's)
   ==============================================#
   #First creating Sidler's poroelastic parameters and then converting to
   #nessesary parameters for our FD
   #----solid phase-----
   Km=zeros(nz,nr) #frame
   Ks=zeros(nz,nr) #grain
   G=zeros(nz,nr) #bulk
   Rho=zeros(nz,nr) #bulk (weighted average of Rhos and Rhof)
   Rhos=zeros(nz,nr) #grain
   #---fluid phase----
   Kf=zeros(nz,nr)
   Rhof=zeros(nz,nr)
   Kappa0=zeros(nz,nr) #static permeability (m^2, 1D=9.869E-13m^2)
   Eta=zeros(nz,nr) #dynamic viscosity Pa.s=kg/(m.s)
   #---other parameters---
   Phi=zeros(nz,nr)
   Tot=zeros(nz,nr) #Tortuosity factor (see Sidler 2014)
   #End initializing poroelastic parameters

   #==============================================
   Fluid parameters
   ==============================================#
   #Assuming that fluid properties are identical everywhere
   Kf[:,:]=ones(nz,nr)*2.25*10^9
   Rhof[:,:]=ones(nz,nr)*1000
   Eta[:,:]=ones(nz,nr)*1.0*10^(-3) #water 1E-3 Pa.s

   Flag_AC=zeros(nz,nr)
   ir_wall=0

   #==========================
   Homogeneous PoroElastic model
   ==========================#
   Kappa0[:,:]=ones(nz,nr)*1*9.869*10^(-13) #m^2
   Phi[:,:]=ones(nz,nr)*0.3
   #--Tortuosity factor from Ou
   m_Ou=8.0 #
   alpha_inf_Ou=3.0 #tortuosity
   Tot=ones(nz,nr)*(1+2/m_Ou)*alpha_inf_Ou
   #-Solid parameters from dry_Vp,dry_Vs,grain_density, and grain_bulkmodulus (Table 1 of Ou2019)
   Vp_dry=2170.0
   Vs_dry=1360.0
   Rhos[:,:]=ones(nz,nr)*3143.0
   Ks[:,:]=ones(nz,nr)*5*1E10
   #
   Rho_dry=(ones(nz,nr)-Phi).*Rhos
   G[:,:]=Vs_dry^2*Rho_dry[:,:]
   Km[:,:]=Vp_dry^2*Rho_dry[:,:]-4/3*G[:,:]
   Rho=(ones(nz,nr)-Phi).*Rhos+Phi.*Rhof


   Flag_E=zeros(nz,nr)


   #=====================================================
   Converting Sidler's poroelastic parameters into Ou's
   =====================================================#
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

   #Correction of H and C for elastic media
   for ir=1:nr
      for iz=1:nz
         if(Flag_E[iz,ir]==1)
            H[iz,ir]=Ks[iz,ir]+4/3*G[iz,ir]
            C[iz,ir]=Ks[iz,ir] #Ou's
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


function drawsnap(Data1,nz,dz,nr,dr,LPML_r,LPML_z)
   println("Snapshot drawing...")
   zvec=[0:dz:(nz-1)*dz]
   rvec=[0:dr:(nr-1)*dr]
   R=(nr-1)*dr
   plt1=heatmap(rvec,zvec,Data1,yflip=true,xlabel="R (m)",ylabel="Z (m)",legend=:none)
   plot!([0;R],ones(2)*(LPML_z-1)*dz,label="",color=:black)
   plot!([0;R],ones(2)*((nz-LPML_z+1)*dz),label="",color=:black)
   plot!(ones(2)*((nr-LPML_r+1)*dr),[0;(nz-1)*dz],label="",color=:black)
   display(plot(plt1))

end

#=========================
Computing Vp_sat, Vs_sat
from poroelastic properties
=========================#
function Gassman(vp_dry,vs_dry,rhos,phi,Ks)
#   vp_dry=5170; #dry vp
#   vs_dry=3198; #dry vs
#   rhos=3143; #grain density
#   phi=0.3; #porosity
#   Ks=10E10; #grain bulk modulus
   Kf=1000*1500^2; #fluid bulk modulus
   Rhof=1000;

   rho_dry=(1-phi)*rhos;
   rho_bulk=(1-phi)*rhos+phi*Rhof;
   G=vs_dry^2*rho_dry;
   K_dry=vp_dry^2*rho_dry-4/3*G;
   alpha=1-K_dry./Ks;
   M=( (alpha-phi)./Ks+phi./Kf ).^(-1);
   C=M.*alpha;
   H=K_dry+4/3*G+M.*alpha.^2; #H=K_undrained+4/3*G (see 2.6c, Guan2011)


   #Gassman K_undrained=alpha^2*M+K_dry (see Guan2011)
   #Gassman K_sat=K_undrained
   Ksat=K_dry+(1-K_dry/Ks)^2/(
       phi/Kf+(1-phi)/Ks-K_dry/Ks^2);
   K_undrained=K_dry+alpha^2*M;
   Vp0=sqrt(H/rho_bulk); #->sqrt((K_undrained+4/3G)/rho_bulk)
   Vs0=sqrt(G/rho_bulk); #->sqrt((K_undrained+4/3G)/rho_bulk)
   #Skempton Coefficient
   B=(1/K_dry-1/Ks)/(1/K_dry-1/Ks+phi*(1/Kf-1/Ks));

   return Vp0,Vs0,rho_bulk
end



#----MAIN loop function with src----
function main_loop!(nr,nz,dr,dz,Rho,Rhof,M,C,H,G,D1,D2,dt,nt,T,
   vr,vz,trr,tpp,tzz,trz,
   vfr,vfz,pf,
   ir_wall,
   src_index,src_dn,
   Flag_AC,Flag_E,
   nrec,rec_vr,rec_vz,rec_pf,index_allrec_vr,index_allrec_vz,index_allrec_pf,
   LPML_r,LPML_z,PML_Wr,PML_Wz,PML_IWr,PML_Wr2,PML_Wz2,PML_IWr2,
   snapshots_vr,snapshots_vz,snapshots_trr,nsnap,itvec_snap,nskip)

#println("start")

#--
Flag_vf_zero=get_Flag_vf_zero(Flag_AC,Flag_E,nr,nz)

#initializing matrices of PML
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


@showprogress for ii=1:nt

#----Time at (ii-1)*dt---(updating velocities)-----

#--PML: save velocity at previous step
PML_save_vel!(memT_vr,memT_vz,memT_vfr,memT_vfz,
              memB_vr,memB_vz,memB_vfr,memB_vfz,
              memR_vr,memR_vz,memR_vfr,memR_vfz,
              vr,vz,vfr,vfz,nr,nz,LPML_z,LPML_r)

# Main velocity update!
update_velocity_1st_Por!(vr,vz,trr,tpp,tzz,trz,vfr,vfz,pf,Rho,Rhof,D1,D2,Flag_vf_zero,nr,nz,dr,dz,dt,LPML_z,LPML_r)

# PML: updating velocity
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

# Making sure RightBC for velocity (zero values -> vr)
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


#PML: updating memory variables for stress (Rx and Sxx) using velocity at two-time steps
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
#srcapply!(vfz,src_index,src_dn,srcamp)
#srcapply!(vfr,src_index,src_dn,srcamp)

#----Time at (ii-1)*dt+dt/2----(updating stress)---


#--PML: saving stress at the previous step
PML_save_stress!(memT_trr,memT_tpp,memT_tzz,memT_trz,memT_pf,
                 memB_trr,memB_tpp,memB_tzz,memB_trz,memB_pf,
                 memR_trr,memR_tpp,memR_tzz,memR_trz,memR_pf,
                 trr,tpp,tzz,trz,pf,nr,nz,LPML_z,LPML_r)


# Main stress update!
update_stress_1st_Por!(vr,vz,trr,tpp,tzz,trz,vfr,vfz,pf,
     M,C,H,G,nr,nz,dr,dz,dt,LPML_z,LPML_r)

# PML: updating stress
PML_update_stress!(vr,vz,trr,tpp,tzz,trz,vfr,vfz,pf,
                  M,C,H,G,nr,nz,dr,dz,dt,
                  Rz_T,Srz_T,RzPE_T,
                  Rz_B,Srz_B,RzPE_B,
                  Rr_R,Rp_R,Rrz_R,RrPE_R,RpPE_R,
                  Rr_TR,Rp_TR,Rz_TR,Rrz_TR,Srz_TR,RrPE_TR,RpPE_TR,RzPE_TR,
                  Rr_BR,Rp_BR,Rz_BR,Rrz_BR,Srz_BR,RrPE_BR,RpPE_BR,RzPE_BR,
                  LPML_z,LPML_r)


# Making sure RightBC for stress (zero values -> trp and trz)
ApplyBCRight_stress!(vr,vz, #Use with flag_zero=1 (see ApplyBCRight_stress1D01)
     trz,
     G,nr,nz,dr,dz,dt)

#---stress B.Cs
ApplyBCLeft_stress!(vr,vz,trr,tpp,tzz,trz,vfr,vfz,pf,
                    M,C,H,G,nr,nz,dr,dz,dt,
                    Rz_T,Srz_T,RzPE_T,
                    Rz_B,Srz_B,RzPE_B,
                    LPML_z)

#--src injection (stress:monopole src)
#srcamp=-src_func[ii]
#srcapply!(tzz,src_index,src_dn,srcamp)
#srcapply!(trr,src_index,src_dn,srcamp)
#srcapply!(tpp,src_index,src_dn,srcamp)
srcamp=src_func[ii]
srcamp=srcamp*dt/Rhof[1,1]
srcapply!(pf,src_index,src_dn,srcamp)


#--PML: updating memory variables for velocity (Pxx and Qxx) using stress at two-time steps
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

#Receiver field (extraction)
getRecData_from_index!(vr,rec_vr,index_allrec_vr,nrec,ii)
getRecData_from_index!(vz,rec_vz,index_allrec_vz,nrec,ii)
#getRecData_from_index!(trr,rec_tii,index_allrec_tii,nrec,ii)
getRecData_from_index!(pf,rec_pf,index_allrec_pf,nrec,ii)



#--Snapshots
check_snap=findall(x -> x==ii,itvec_snap)
if (length(check_snap)!=0)
   println("ii=",ii)
   cnt_snap=Int(check_snap[1])
   get_snapshots!(snapshots_vr,snapshots_vz,cnt_snap,vr,vz,nr,nz)
   get_snapshots_t!(snapshots_trr,cnt_snap,pf,nr,nz)

   drawsnap(snapshots_vz[:,2:end,cnt_snap],nz,dz,nr-1,dr,
            LPML_r,LPML_z)

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

#time samples
dt=0.5E-4
nt=2001

T=(nt-1)*dt
tvec=range(0.0,T,length=nt) # range object (no memory allocation)
tvec=collect(tvec) # a vector

#PML thickness in samples
LPML_r=30
LPML_z=30

#======================
Creating model
======================#
nr,nz,dr,dz,Rho,Rhof,M,C,H,G,D1,D2,Flag_AC,Flag_E,ir_wall=makemodel_homogeneous_PoroElastic()
drawmodel(H,nz,dz,nr,dr,LPML_r,LPML_z)

#======================
Stability check
======================#
Vmax=maximum((H./Rho).^(0.5))
check_stability01(dt,dr,Vmax,0)
check_stability01(dt,dz,Vmax,0)

#======================
Creating src wavelet
======================#
f0=100 #src Freq
delay=1/f0*1.0
#src_func=myricker2(tvec,f0,delay,2) #when using 2nd derivative Gaussian (good for DWI)
#src_func=myricker2(tvec,f0,delay,1) #when using 1st derivative Gaussian (Randall?)
src_func=myricker2(tvec,f0,delay,3) #when using 3rd derivative Gaussian
tmp_maxamp=maximum(map(abs,src_func))
src_func=src_func/tmp_maxamp
display(plot(tvec,src_func[:],title="src function"))
tmpfilename_src=string(@__DIR__,"/src_function.png")
png(tmpfilename_src)
println("Source signature figure saved in ",tmpfilename_src)

#================================
Point src geometry
================================#
srcgeom=zeros(1,2) #(z,r)
srcgeom[1,1]=(nz-1)*dz/2 #z meter
#srcgeom[1,2]=(nr-LPML_r+1)*dr #r meter
srcgeom[1,2]=0*dr #r meter
#Gaussian point src
wsize=7 #window size (odd number)
wsigma=1 #std
src_index,src_dn=get_srcindex_pGauss(srcgeom,dr,dz,wsize,wsigma)
plot(src_index[1,:],src_index[2,:],src_dn[:],marker=:circle,seriestype=:scatter,zcolor=src_dn[:],camera=(90,90))

#==============================
Initializing field variables
==============================#
vr,vz,trr,tpp,tzz,trz,vfr,vfz,pf=init_fields_Por(nz,nr) #

# Initializing PML field variables
LPML_r,LPML_z,PML_Wr,PML_Wz,
PML_IWr,PML_Wr2,PML_Wz2,PML_IWr2=init_PML_profile(LPML_r,LPML_z,Vmax,dr,dz,nr,f0)
#PML_check(LPML_r,LPML_z,Vmax,dr,dz,f0)

#===================
Receiver geometry
====================#
#--Fluid pressure at borehole center
nrec=nz
recgeom=zeros(nrec,2) #(z,r)
for irec=1:nrec
   recgeom[irec,1]=dz*(irec-1)
   recgeom[irec,2]=25.0
end
rec_vr,rec_vz,index_allrec_vr,index_allrec_vz=init_receiver_geophone(recgeom,nrec,dr,dz,nr,nz,nt)
rec_pf,index_allrec_pf=init_receiver_hydrophone(recgeom,nrec,dr,dz,nr,nz,nt)



#==============================
Snapshot settings
==============================#
nskip,snapshots_trr,snapshots_vr,snapshots_vz,
nsnap,itvec_snap=init_snap(nt,nz,nr,30)
#error()

#==============================
Start main FD Loop
==============================#

main_loop!(nr,nz,dr,dz,Rho,Rhof,M,C,H,G,D1,D2,dt,nt,T,
   vr,vz,trr,tpp,tzz,trz,
   vfr,vfz,pf,
   ir_wall,
   src_index,src_dn,
   Flag_AC,Flag_E,
   nrec,rec_vr,rec_vz,rec_pf,index_allrec_vr,index_allrec_vz,index_allrec_pf,
   LPML_r,LPML_z,PML_Wr,PML_Wz,PML_IWr,PML_Wr2,PML_Wz2,PML_IWr2,
   snapshots_vr,snapshots_vz,snapshots_trr,nsnap,itvec_snap,nskip)


CPUtoc()
println("elapsed real time: ", round(time() - start;digits=3)," seconds")

error("Stopped w/o problem.")


offset=recgeom[:,1]-srcgeom[1]*ones(size(recgeom[:,1]))
heatmap(offset[:],tvec[:],rec_vz,xlims=(-50,50),ylims=(0,0.1),clims=(-1E-15,1E-15))
heatmap!(flip=true)

heatmap(offset[:],tvec[:],rec_vr,xlims=(-50,50),ylims=(0,0.1),clims=(-1E-15,1E-15))
heatmap!(flip=true)

#Check
rec_vz_Aki=Green_Aki(2525.7,1275.8,2500.1,tvec,src_func,srcgeom,recgeom)
rec_vz_Aki=rec_vz_Aki/maximum(rec_vz_Aki)
rec_vz_norm=rec_vz/maximum(rec_vz)
plt1=heatmap(recgeom[:,1],tvec[:],rec_vz_norm,xlabel="Z (m)",ylabel="time (s)",title="FD")
plt2=heatmap(recgeom[:,1],tvec[:],rec_vz_Aki,xlabel="Z (m)",ylabel="time (s)",title="Theory")
plt3=plot(tvec[:],rec_vz_norm[:,151],xlabel="time (s)",label="FD")
plot!(tvec[:],rec_vz_Aki[:,151],xlabel="time (s)",label="Theory")
display(plot(plt1,plt2,plt3,layout=(3,1)))

tol=0.05
sqrt(sum((rec_vz_norm[:,151] - rec_vz_Aki[:,151]) .^ 2)) / sqrt(sum(rec_vz_norm[:,151] .^ 2)) < tol

error()

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
tmp_filename="out_2L.mat"

iskip_rec=10
tvec_rec=tvec[1:iskip_rec:end]

tmp_filename_full=string(@__DIR__,"/",tmp_filename)
matwrite(tmp_filename_full, Dict(
"tvec"=>tvec,"nr"=>nr,"nz"=>nz,"dr"=>dr,"dz"=>dz,
"dt"=>dt,"nt"=>nt,"T"=>T,
"tvec_rec"=>tvec_rec,
"rec_tii"=>rec_tii[1:iskip_rec:end,:],
#"rec_vz"=>rec_vz[1:iskip_rec:end,:],
#"rec_vr"=>rec_vr,
#"rec_tii_PE"=>rec_tii_PE[1:iskip_rec:end,:],
#"rec_pf_PE"=>rec_pf_PE[1:iskip_rec:end,:],
"src_func"=>src_func,
"index_allrec_tii"=>index_allrec_tii,
#"index_allrec_vz"=>index_allrec_vz,
#"index_allrec_vr"=>index_allrec_vr,
#"index_allrec_pf_PE"=>index_allrec_pf_PE,
#"src_index"=>src_index,
#"snapshots_trr"=>snapshots_trr,
#"snapshots_pf"=>snapshots_trr,
#"snapshots_vr"=>snapshots_vr,
"snapshots_vz"=>snapshots_vz,
"itvec_snap"=>itvec_snap); compress = true)

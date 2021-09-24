#testing FDTD in cylindrical coordinate
#Randall et al (1991)

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
Creating two layer model
=========================#
function makemodel_homogeneous_Elastic()
#Nessesary input matrices to main loop function are:
#Rho: Bulk density
#Rhof: Fluid density
#M,C,H: Poroelastic moduli
#G: Formation shear moduli
#D1,D2: Poroelastic moduli in Ou's formulation
#Flag_AC,Flag_E: flag specifying acoustic region or elastic region
#Other input parameters to main loop function are:
#ir_wall: grid number in r direction where a borehole wall starts (single value)
#
# This function creates above input paramers assuming a borehole embedded in
# homogeneous elastic media

#Model size
   nr=351 #samples
   nz=1401 #samples
   dr=0.005 #meter
   dz=0.005 #meter


   #==============================================
   Initializing poroelastic parameters (Sidler's)
   ==============================================#
   #First creating Sidler's poroelastic parameters and then converting to
   #nessesary parameters for our FD
   #----solid phase-----
   Km=zeros(nz,nr) #frame
   Ks=zeros(nz,nr) #grain
   G=zeros(nz,nr) #bulk
   Rho=zeros(nz,nr) #bulk (weighting average of Rhos and Rhof)
   Rhos=zeros(nz,nr) #grain
   #---fluid phase----
   Kf=zeros(nz,nr)
   Rhof=zeros(nz,nr)
   Kappa0=zeros(nz,nr) #static permeability (m^2, 1D=9.869E-13m^2)
   Eta=zeros(nz,nr) #dynamic viscosity Pa.s=kg/(m.s)
   #---other parameters---
   Phi=zeros(nz,nr)
   Tot=zeros(nz,nr) #Tortuosity factor (see Sidler 2014)
   #End initilizing poroelastic parameters

   #==============================================
   Fluid parameters
   ==============================================#
   #Assuming that fluid properties are identical everywhere
   Kf[:,:]=ones(nz,nr)*2.25*10^9
   Rhof[:,:]=ones(nz,nr)*1000
   Eta[:,:]=ones(nz,nr)*1.0*10^(-3) #water 1E-3 Pa.s

   #========================================
   Inclusion of a borehole (acoustic media)
   ========================================#
   Flag_AC=zeros(nz,nr)
   #Borehole wall (homogeneous radius)
   ir_wall=Int(round(0.1/dr))

   for iz=1:nz
       Flag_AC[iz,1:ir_wall]=ones(ir_wall)
       Rho[iz,1:ir_wall]=ones(ir_wall)*Rhof[1,1] #Rho=Rhof
       G[iz,1:ir_wall]=ones(ir_wall)*0.0 #G=0
       Phi[iz,1:ir_wall]=ones(ir_wall) #Phi=1
       Km[iz,1:ir_wall]=ones(ir_wall)*0.0 #Km=0 (G=0,phi=0->M=C=H=Kf, Ou's)
       Ks[iz,1:ir_wall]=ones(ir_wall)*Kf[1,1] #
       Eta[iz,1:ir_wall]=ones(ir_wall)*0.0 #eta=0 (-> D2=0, Ou's)
   end

   #==========================
   Homogeneous Elastic model
   ==========================#
   Flag_E=zeros(nz,nr)
   Vp1_elastic=2500.0 #Or, Specify K_elastic
   Vs1_elastic=2000.0
   Rho1_elastic=2500.0

#
   for iz=1:nz #Elastic media

       G_elastic=Vs1_elastic^2*Rho1_elastic
       K_elastic=Vp1_elastic^2*Rho1_elastic-4/3*G_elastic #Activate here if you have specified Vp_elastic

       Flag_E[iz,ir_wall+1:end]=ones(nr-ir_wall)
       Rho[iz,ir_wall+1:end]=ones(nr-ir_wall)*Rho1_elastic #Rho=Rho_elastic
       Rhof[iz,ir_wall+1:end]=ones(nr-ir_wall)*Rho1_elastic #Rhof=Rho_elastic
       G[iz,ir_wall+1:end]=ones(nr-ir_wall)*G_elastic #G=Gs (Ou's)
       #--
       Phi[iz,ir_wall+1:end]=ones(nr-ir_wall)*0.0 #Phi=0
       Km[iz,ir_wall+1:end]=ones(nr-ir_wall)*K_elastic #Km
       Ks[iz,ir_wall+1:end]=ones(nr-ir_wall)*K_elastic #Ks=Km (-> M=Inf, Ou's)
       #--
       Kappa0[iz,ir_wall+1:end]=ones(nr-ir_wall)*0.0 #k0=0 (-> D1=D2=Inf, Ou's)
   end

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




#----MAIN loop function with src----
function main_loop!(nr,nz,dr,dz,Rho,Rhof,M,C,H,G,D1,D2,dt,nt,T,
   vr,vz,trr,tpp,tzz,trz,
   vfr,vfz,pf,
   ir_wall,
   src_index,src_dn,
   Flag_AC,Flag_E,
   nrec,index_allrec_tii,rec_tii,
   LPML_r,LPML_z,PML_Wr,PML_Wz,PML_IWr,PML_Wr2,PML_Wz2,PML_IWr2,
   snapshots_vr,snapshots_vz,snapshots_trr,nsnap,itvec_snap,nskip)

#println("start")

#--
Flag_vf_zero=get_Flag_vf_zero(Flag_AC,Flag_E,nr,nz)

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


# Making sure RightBC for stress (zero values are trp and trz)
ApplyBCRight_stress!(vr,vz, #Use with flag_zero=1 (see ApplyBCRight_stress1D01)
    trz,
    G,nr,nz,dr,dz,dt)

#error()
@showprogress for ii=1:nt
#@showprogress for ii=1:201
#for ii=1:1



#----Time at (ii-1)*dt---(updating velocities)-----

#--PML: save velocity at previous step
PML_save_vel!(memT_vr,memT_vz,memT_vfr,memT_vfz,
              memB_vr,memB_vz,memB_vfr,memB_vfz,
              memR_vr,memR_vz,memR_vfr,memR_vfz,
              vr,vz,vfr,vfz,nr,nz,LPML_z,LPML_r)
#println("update vel")

#keep current values before updating velocity for Acoustic-Poroelastic BC later
#mycopy_mat(vr,vr_old,nz,nr)
#mycopy_mat(vz,vz_old,nz,nr)
#mycopy_mat(vfr,vfr_old,nz,nr)
#mycopy_mat(vfz,vfz_old,nz,nr)

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
#update_vr_vfr_1st_vertical(vr,trr,tpp,tzz,trz,vfr,pf,
#    vr_old,vfr_old,
#    Rho,Rhof,D1,D2,nr,nz,dr,dz,dt,LPML_z,
#    Prz_T,Prz_B,
#    ir_wall,Flag_AC,Flag_E)


#PML: update memory variables for stress (Rx and Sxx) using velocity at two time steps
PML_update_memRS!(Rz_T,Srz_T,RzPE_T,
                  Rz_B,Srz_B,RzPE_B,
                  Rr_R,Rp_R,Rrz_R,RrPE_R,RpPE_R,
                  Rr_TR,Rp_TR,Rz_TR,Rrz_TR,Srz_TR,RrPE_TR,RpPE_TR,RzPE_TR,
                  Rr_BR,Rp_BR,Rz_BR,Rrz_BR,Srz_BR,RrPE_BR,RpPE_BR,RzPE_BR,
                  memT_vr,memT_vz,memT_vfr,memT_vfz,
                  memB_vr,memB_vz,memT_vfr,memT_vfz,
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

#--src injection (stress:monopole src)
srcamp=-src_func[ii]
srcapply!(tzz,src_index,src_dn,srcamp)
srcapply!(trr,src_index,src_dn,srcamp)
srcapply!(tpp,src_index,src_dn,srcamp)


#---Additional BC when Flag_Acoustic/Flag_Elastic
#ApplyBC_stress_AcousticMedia_TEST!(trr,tpp,tzz,trz,pf,Flag_AC,nr,nz) #pf=-1/3tii
#ApplyBC_stress_ElasticMedia_Ou!(pf,Flag_E,nr,nz) #when Flag_E==1, then pf=0


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
#Receiver field (extraction)
#getRecData_from_index!(vr,rec_vr,index_allrec_vr,nrec,ii)
#getRecData_from_index!(vz,rec_vz,index_allrec_vz,nrec,ii)
getRecData_from_index!(trr,rec_tii,index_allrec_tii,nrec,ii)



#--Snapshots
#println("check snap")
check_snap=findall(x -> x==ii,itvec_snap)
if (length(check_snap)!=0)
   println("ii=",ii)
   cnt_snap=Int(check_snap[1])
   get_snapshots!(snapshots_vr,snapshots_vz,cnt_snap,vfr,vz,nr,nz)
   get_snapshots_t!(snapshots_trr,cnt_snap,pf,nr,nz)

   drawsnap(snapshots_vz[:,:,cnt_snap],nz,dz,nr,dr,
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

#==time samples
const dt=0.125E-5
const nt=16001
const T=(nt-1)*dt
==#
dt=1E-6
nt=2501
T=(nt-1)*dt
tvec=range(0.0,T,length=nt) # range object (no memory allocation)
tvec=collect(tvec) # a vector

#PML thickness in samples
LPML_r=60
LPML_z=60

#======================
Creating model
======================#
nr,nz,dr,dz,Rho,Rhof,M,C,H,G,D1,D2,Flag_AC,Flag_E,ir_wall=makemodel_homogeneous_Elastic()
drawmodel(H,nz,dz,nr,dr,LPML_r,LPML_z)
#error()

#======================
Stability check
======================#
Vmax=maximum((H./Rho).^(0.5))
check_stability01(dt,dr,Vmax,0)
check_stability01(dt,dz,Vmax,0)
#error()
#======================
Creating src wavelet
======================#
f0=5000 #src Freq
delay=1/f0*1.0
#src_func=myricker2(tvec,f0,delay,2) #when using 2nd derivative Gaussian (good for DWI)
src_func=myricker2(tvec,f0,delay,1) #when using 1st derivative Gaussian (Randall?)
#src_func=myricker2(tvec,f0,delay,3) #when using 3rd derivative Gaussian
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
srcgeom[1,1]=(nz-1)*dz/3 #z meter
srcgeom[1,2]=0.1*dr #r meter
src_index,src_dn=get_srcindex_monopole(srcgeom,dr,dz)

#error()
#==============================
Initializig field variables
==============================#
vr,vz,trr,tpp,tzz,trz,vfr,vfz,pf=init_fields_Por(nz,nr) #

# Initializing PML field variables
LPML_r,LPML_z,PML_Wr,PML_Wz,
PML_IWr,PML_Wr2,PML_Wz2,PML_IWr2=init_PML_profile(LPML_r,LPML_z,Vmax,dr,dz,nr)
#PML_check(LPML_r,LPML_z,Vmax,dr,dz,f0)

#===================
Receiver geometry
====================#
#--Fluid pressure at borehole center
nrec=nz
recgeom=zeros(nrec,2) #(z,r)
for irec=1:nrec
   recgeom[irec,1]=dz*(irec-1)
   recgeom[irec,2]=0.0
end
rec_tii,index_allrec_tii=init_receiver_hydrophone(recgeom,nrec,dr,dz,nr,nz,nt)



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
   vr,vz,trr,tpp,tzz,trz,
   vfr,vfz,pf,
   ir_wall,
   src_index,src_dn,
   Flag_AC,Flag_E,
   nrec,index_allrec_tii,rec_tii,
   LPML_r,LPML_z,PML_Wr,PML_Wz,PML_IWr,PML_Wr2,PML_Wz2,PML_IWr2,
   snapshots_vr,snapshots_vz,snapshots_trr,nsnap,itvec_snap,nskip)


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

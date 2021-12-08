#--
# Threads.nthreads() displays number of available threads
# if ==1, start julia with JULIA_NUM_THREADS=4 julia
#--

#--
# If no graphics (X-windows) required, then ENV["GKSwstype"]="nul"
#--

using FDModBiotCyl
#
using CPUTime
using ProgressMeter
using Plots
using Printf
using Dierckx # for Aki-Richards Green's function

#--Suppressing graphics
#ENV["GKSwstype"]="nul"

@__DIR__
@__FILE__


#================
Drawing a model
================#
function drawmodel(Field_Var,nz,dz,nr,dr,LPML_r,LPML_z)
   println("Drawing a model...")
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

#================
Drawing a snapshot
================#
function drawsnap(Data1,nz,dz,nr,dr,LPML_r,LPML_z,title_str)
   println("Drawing a snapshot...")
   zvec=[0:dz:(nz-1)*dz]
   rvec=[0:dr:(nr-1)*dr]
   R=(nr-1)*dr
   plt1=heatmap(rvec,zvec,Data1,yflip=true,xlabel="R (m)",ylabel="Z (m)",legend=:none)
   heatmap!(title=title_str)
   plot!([0;R],ones(2)*(LPML_z-1)*dz,label="",color=:black)
   plot!([0;R],ones(2)*((nz-LPML_z+1)*dz),label="",color=:black)
   plot!(ones(2)*((nr-LPML_r+1)*dr),[0;(nz-1)*dz],label="",color=:black)
   display(plot(plt1))

end

#=========================
Creating a homogeneous model
=========================#
function makemodel_homogeneous_Elastic()
   #Necessary input matrices to the main loop function are:
   #Rho: Bulk density
   #Rhof: Fluid density
   #M,C,H: Poroelastic moduli
   #G: Formation shear modulus
   #D1,D2: Poroelastic moduli (Guan and Hu, 2011; Ou and Wang, 2019)
   #Flag_AC,Flag_E: flag specifying acoustic region or elastic region
   # This function creates the above input parameters assuming a
   # homogeneous elastic media

   # Model size
   nr=251 #samples
   nz=351 #samples
   dr=0.5 #meter
   dz=0.5 #meter


   #==============================================
   Initializing grain/frame parameters
   ==============================================#
   #First creating grain/frame parameters and then converting them
   #to necessary poroelastic parameters for our FD
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

   #==========================
   Homogeneous Elastic model
   ==========================#
   Flag_E=zeros(nz,nr)
   Vp1_elastic=2500.0 #Or, Specify K_elastic
   Vs1_elastic=2000.0
   Rho1_elastic=2000.0

   G_elastic=Vs1_elastic^2*Rho1_elastic
   K_elastic=Vp1_elastic^2*Rho1_elastic-4/3*G_elastic

   Flag_E=ones(nz,nr)
   Rho=Rho1_elastic*ones(nz,nr) #Rho=Rho_elastic
   G=ones(nz,nr)*G_elastic #G=Gs (Ou's)
   #--
   Km=ones(nz,nr)*K_elastic #Km
   Ks=ones(nz,nr)*K_elastic #Ks=Km (-> M=Inf, Ou's)
   #--

   #=====================================================
   Converting grain/frame parameters into poroelastic parameters
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


   return nr,nz,dr,dz,Rho,Rhof,M,C,H,G,D1,D2,Flag_AC,Flag_E
end





#=========================
Green's function
=========================#
function Green_Aki(vp,vs,rho,tvec,src_func,srcgeom,recgeom)
#Aki and Rechards eq 4.23
#vp=2500
#vs=2000
#rho=2500

tmp_R=recgeom-repeat(srcgeom,nrec,1) #(z,r)
R=map(x->sqrt(x),sum(tmp_R.^2,dims=2))
gammax=tmp_R[:,2]./R
gammaz=tmp_R[:,1]./R
src_spl=Spline1D(tvec[:],src_func[:],k=1)
    #================
          Gzz
    ================#
     amp_p=1/(4pi*rho*vp^2)*gammaz.*gammaz./R
     amp_s=-1/(4pi*rho*vs^2)*(gammaz.*gammaz-ones(size(gammax)))./R
     Gzz_far=zeros(nt,nrec)
     for irec=1:nrec
        tmp=amp_p[irec]*src_spl(tvec[:]-R[irec]/vp*ones(size(tvec)))+
            amp_s[irec]*src_spl(tvec[:]-R[irec]/vs*ones(size(tvec)))
        Gzz_far[:,irec]=tmp #far-field terms only (displacement-body force)
     end
     #inclusion of a near-field term
     Gzz_near=zeros(nt,nrec)
     amp_nf=1/(4pi*rho)*(3*gammaz.*gammaz-ones(size(gammax)))./(R.^3)

     for irec=1:nrec
        tau=range(R[irec]/vp,R[irec]/vs,length=nt) # range object (no memory allocation)
        tau=collect(tau) # a vector
        dt_tau=tau[2]-tau[1]
        tmp_nf=zeros(nt,1)
        for it=1:nt
           src_intp=src_spl(tvec[it]*ones(size(tau))-tau[:])
           tmp_nf[it]=sum(src_intp.*tau)*dt_tau
        end
        Gzz_near[:,irec]=tmp_nf*amp_nf[irec]
     end

     Gzz=Gzz_far+Gzz_near

#converting displacement to particle velocity
  tmp=(Gzz[2:end,:]-Gzz[1:end-1,:])/dt
  vz_an=zeros(nt,nrec)
  for irec=1:nrec
    tmp_spl=Spline1D(tvec[1:end-1]+dt/2*ones(size(tmp[:,irec])),tmp[:,irec],k=1)
    vz_an[:,irec]=tmp_spl(tvec)
  end

  return vz_an
end



#----MAIN loop function with src----
function main_loop!(nr,nz,dr,dz,Rho,Rhof,M,C,H,G,D1,D2,dt,nt,T,
   vr,vz,trr,tpp,tzz,trz,
   vfr,vfz,pf,
   src_index,src_dn,src_func,
   Flag_AC,Flag_E,
   nrec,rec_vr,rec_vz,index_allrec_vr,index_allrec_vz,
   LPML_r,LPML_z,PML_Wr,PML_Wz,PML_IWr,PML_Wr2,PML_Wz2,PML_IWr2,
   snapshots_vr,snapshots_vz,snapshots_trr,nsnap,itvec_snap,nskip)

#println("start")

#--
Flag_vf_zero=get_Flag_vf_zero(Flag_AC,Flag_E,nr,nz)

#initialize matrices of PML
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
#println("update vel")

# Main velocity update!
update_velocity_1st_Por!(vr,vz,trr,tpp,tzz,trz,vfr,vfz,pf,Rho,Rhof,D1,D2,Flag_vf_zero,nr,nz,dr,dz,dt,LPML_z,LPML_r)

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
srcamp=src_func[ii]
srcamp=srcamp*dt/Rho[1,1]
srcapply!(vz,src_index,src_dn,srcamp)
#srcapply!(vr,src_index,src_dn,srcamp)

#----Time at (ii-1)*dt+dt/2----(updating stress)---


#--PML: save stress at previous step
PML_save_stress!(memT_trr,memT_tpp,memT_tzz,memT_trz,memT_pf,
                 memB_trr,memB_tpp,memB_tzz,memB_trz,memB_pf,
                 memR_trr,memR_tpp,memR_tzz,memR_trz,memR_pf,
                 trr,tpp,tzz,trz,pf,nr,nz,LPML_z,LPML_r)


# Main stress update!
update_stress_1st_Por!(vr,vz,trr,tpp,tzz,trz,vfr,vfz,pf,
     M,C,H,G,nr,nz,dr,dz,dt,LPML_z,LPML_r)

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
ApplyBCLeft_stress!(vr,vz,trr,tpp,tzz,trz,vfr,vfz,pf,
                    M,C,H,G,nr,nz,dr,dz,dt,
                    Rz_T,Srz_T,RzPE_T,
                    Rz_B,Srz_B,RzPE_B,
                    LPML_z)

#--src injection (stress:monopole src)
#srcamp=-src_func[ii]
#srcapply!(trz,src_index,src_dn,srcamp)
#srcapply!(tzz,src_index,src_dn,srcamp)
#srcapply!(trr,src_index,src_dn,srcamp)
#srcapply!(tpp,src_index,src_dn,srcamp)



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
#Receiver field
getRecData_from_index!(vr,rec_vr,index_allrec_vr,nrec,ii)
getRecData_from_index!(vz,rec_vz,index_allrec_vz,nrec,ii)



#--Snapshots
check_snap=findall(x -> x==ii,itvec_snap)
if (length(check_snap)!=0)
   println("ii=",ii)
   cnt_snap=Int(check_snap[1])
   get_snapshots!(snapshots_vz,cnt_snap,vz,nr,nz)
   get_snapshots!(snapshots_vr,cnt_snap,vr,nr,nz)
   get_snapshots!(snapshots_trr,cnt_snap,trr,nr,nz)


   drawsnap(snapshots_vz[:,:,cnt_snap],nz,dz,nr,dr,
            LPML_r,LPML_z,"Snapshots: Vz")

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
dt=1E-4
nt=751
#nt=101
T=(nt-1)*dt
tvec=range(0.0,T,length=nt) # range object (no memory allocation)
tvec=collect(tvec) # a vector

#PML thickness in samples
LPML_r=20
LPML_z=20

#======================
Creating a model
======================#
nr,nz,dr,dz,Rho,Rhof,M,C,H,G,D1,D2,Flag_AC,Flag_E=makemodel_homogeneous_Elastic()
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
Creating a src wavelet
======================#
f0=100 #src Freq
delay=1/f0*1.0
#src_func=myricker2(tvec,f0,delay,2) #when using 2nd derivative Gaussian (Ricker with a negative sign)
src_func=myricker2(tvec,f0,delay,1) #when using 1st derivative Gaussian
#src_func=myricker2(tvec,f0,delay,3) #when using 3rd derivative Gaussian
tmp_maxamp=maximum(map(abs,src_func))
src_func=src_func/tmp_maxamp*1E+9 #this scaling is just for visualization

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
srcgeom[1,2]=0*dr #r meter.
#Gaussian point src
wsize=7 #window size (odd number)
wsigma=1 #std
src_index,src_dn=get_srcindex_pGauss(srcgeom,dr,dz,wsize,wsigma)
# if you want to see the amplitude distribution,
# plot(src_index[1,:],src_index[2,:],src_dn[:],marker=:circle,seriestype=:scatter,zcolor=src_dn[:],camera=(90,90))

#error()
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



#==============================
Snapshot settings
==============================#
nskip=30 #change here to adjust time sampling (1~nt)
snapshots_vz,nsnap,itvec_snap=init_snapshots(nt,nz,nr,nskip)
snapshots_vr,nsnap,itvec_snap=init_snapshots(nt,nz,nr,nskip)
snapshots_trr,nsnap,itvec_snap=init_snapshots(nt,nz,nr,nskip)

#error()

#==============================
Start main FD Loop
==============================#

main_loop!(nr,nz,dr,dz,Rho,Rhof,M,C,H,G,D1,D2,dt,nt,T,
   vr,vz,trr,tpp,tzz,trz,
   vfr,vfz,pf,
   src_index,src_dn,src_func,
   Flag_AC,Flag_E,
   nrec,rec_vr,rec_vz,index_allrec_vr,index_allrec_vz,
   LPML_r,LPML_z,PML_Wr,PML_Wz,PML_IWr,PML_Wr2,PML_Wz2,PML_IWr2,
   snapshots_vr,snapshots_vz,snapshots_trr,nsnap,itvec_snap,nskip)


CPUtoc()
println("elapsed real time: ", round(time() - start;digits=3)," seconds")

println("FD finished without a problem.")


#Following is to check the FD results with the analytical solution
println("Calculating analytical solution and plotting final results...")
#Analytical Green's function
rec_vz_Aki=Green_Aki(2500,2000,2000,tvec,src_func,srcgeom,recgeom) #vp,vs,rho
#Plotting results
plt1=heatmap(recgeom[:,1],tvec[:],rec_vz,xlabel="Z (m)",ylabel="time (s)",title="FD")
plt2=heatmap(recgeom[:,1],tvec[:],rec_vz_Aki,xlabel="Z (m)",ylabel="time (s)",title="Theory")
plt3=plot(tvec[:],rec_vz[:,151],xlabel="time (s)",label="FD")
plot!(tvec[:],rec_vz_Aki[:,151],xlabel="time (s)",label="Theory")
display(plot(plt1,plt2,layout=(1,2)))
display(plot(plt3))

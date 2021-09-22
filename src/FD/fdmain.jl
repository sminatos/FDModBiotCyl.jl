# Poroelastic FD in cylindrical coordinate system
# main staggered-grid finite differencing functions

#-----------------------------
# initializing field variables
#-----------------------------
function init_fields_Por(nz,nr)
    vr=zeros(Float64,nz,nr)
    vz=zeros(Float64,nz,nr)
    trr=zeros(Float64,nz,nr)
    tpp=zeros(Float64,nz,nr)
    tzz=zeros(Float64,nz,nr)
    trz=zeros(Float64,nz,nr)

    vfr=zeros(Float64,nz,nr)
    vfz=zeros(Float64,nz,nr)
    pf=zeros(Float64,nz,nr)
    return vr,vz,trr,tpp,tzz,trz,vfr,vfz,pf
end


#-----------------------------
# updating velocity field
#-----------------------------
function update_vr(___vr_now,___trr_f,___trr_b,___trz_f,___trz_b,___trr_av,___tpp_av,
    ___trp_now,___dz,___dr,___dt,___r_now,___Rho,___m)
    #first-order difference
    #X_now: current time, current position [k,j]
    #X_f: foward position (trr[k,j],trz[k+1,j])
    #X_b: backward position (trr[k,j-1],trz[k,j])
    #_av: averaged values (tpp_av,trr_av)
    #trp is at the same location
    #r_now: current radial position @vr
    #m: order (0: monopole, 1:dipole)
    return ___vr_now+___dt/___Rho*(
                +(___trr_f-___trr_b)/___dr #drtrr
                +(___trr_av-___tpp_av+___m*___trp_now)/(___r_now) #trr/r-tpp/r+m*trp/r
                +(___trz_f-___trz_b)/___dz #dztrz
                )
end
function update_vz(___vz_now,___trz_f,___trz_b,___tzz_f,___tzz_b,___trz_av,
    ___tpz_now,___dz,___dr,___dt,___r_now,___Rho,___m)
    #first-order difference
    #X_now: current time, current position [k,j]
    #X_f: foward position (trz[k,j+1],tzz[k,j])
    #X_b: backward position (trz[k,j],tzz[k-1,j])
    #_av: averaged values (trz_av)
    #tpz is at the same location
    #r_now: current radial position @vz
    #m: order (0: monopole, 1:dipole)
    return ___vz_now+___dt/___Rho*(
              +(___trz_f-___trz_b)/___dr #drtrz
              +(___trz_av+___m*___tpz_now)/(___r_now) #trz/r+m*tpz/r
              +(___tzz_f-___tzz_b)/___dz #dztzz
              )
end

function update_vr_1st_Por(vr_now,vfr_now,
    trr_f,trr_b,trz_f,trz_b,trr_av,tpp_av,
    pf_f,pf_b,
    dz,dr,dt,r_now,Rho,Rhof,D1,D2,Flag_vf_zero)
    #first-order difference, 2nd-order accuracy (Graves 1996)
    #X_now: current time, current position [k,j]
    #X_f: foward position (trr[k,j],trz[k+1,j]) 1:nearest, 2:next nearest
    #X_b: backward position (trr[k,j-1],trz[k,j]) 1:nearest, 2:next nearest
    #_av: averaged values (tpp_av,trr_av) --> will be calculated in this routine
    #trp is at the same location (not used)
    #r_now: current radial position @vr
    #m: order (0: monopole, 1:dipole). m=0 only

    drtrr=(trr_f-trr_b)/dr
    dztrz=(trz_f-trz_b)/dz

    A=drtrr+(trr_av-tpp_av)/(r_now)+dztrz #common term in solid and fluid phase

    drpf=(pf_f-pf_b)/dr

    #updating fluid phase
    if(Flag_vf_zero!=1)
        vfr_old=vfr_now
        vfr_now=(D2/2+(D1-Rhof^2/Rho)/dt)^(-1) * (
                    (-D2/2+(D1-Rhof^2/Rho)/dt)*vfr_old-drpf-Rhof/Rho*A
                    )
    else
        vfr_old=vfr_now
        vfr_now=0.0
    end


    dert_vfr=(vfr_now-vfr_old)/dt

    #check if
    #updating solid phase
    vr_old=vr_now
    vr_now=vr_old+dt/Rho*(A-Rhof*dert_vfr)

    return vr_now,vfr_now
end

function update_vz_1st_Por(vz_now,vfz_now,
    trz_f,trz_b,tzz_f,tzz_b,trz_av,
    pf_f,pf_b,
    dz,dr,dt,r_now,Rho,Rhof,D1,D2,Flag_vf_zero)
    #first-order difference, 2nd-order accuracy (Graves 1996)
    #X_now: current time, current position [k,j]
    #X_f: foward position (trz[k,j+1],tzz[k,j]) 1:nearest, 2:next nearest
    #X_b: backward position (trz[k,j],tzz[k-1,j]) 1:nearest, 2:next nearest
    #_av: averaged values (trz_av) --> will be calculated in this routine
    #tpz is at the same location
    #r_now: current radial position @vz
    #m: order (0: monopole, 1:dipole), m=0 only!

    drtrz=(trz_f-trz_b)/dr
    dztzz=(tzz_f-tzz_b)/dz

    A=drtrz+trz_av/(r_now)+dztzz #common term in solid and fluid phase

    dzpf=(pf_f-pf_b)/dz

    #updating fluid phase
    if(Flag_vf_zero!=1)
        vfz_old=vfz_now
        vfz_now=(D2/2+(D1-Rhof^2/Rho)/dt)^(-1) * (
                    (-D2/2+(D1-Rhof^2/Rho)/dt)*vfz_old-dzpf-Rhof/Rho*A
                    )
    else
        vfz_old=vfz_now
        vfz_now=0.0
    end

    dert_vfz=(vfz_now-vfz_old)/dt

    #updating solid phase
    vz_old=vz_now
    vz_now=vz_old+dt/Rho*(A-Rhof*dert_vfz)

    return vz_now,vfz_now
end

function update_velocity_1st_Por!(vr,vz,trr,tpp,tzz,trz,vfr,vfz,pf,
    rhomat,rhofmat,D1mat,D2mat,Flagmat_vf_zero,nr,nz,dr,dz,dt,LPML_z,LPML_r)

#PE-to-E horizontal boundary incorpolated by checking if Flagmat_vf_zero[k,j]==0 && D1[k+1,j]==Inf

#2nd order FD (Graves, 1996)
#assuming PML!
#@inbounds Threads.@threads for j=3:nr-2
#      for k=3:nz-2

@inbounds Threads.@threads for j=2:nr-LPML_r
    for k=LPML_z+1:nz-LPML_z
          #I am accesing [k,j], but the definition of r_now changes with each component!
         r_now=(j-1)*dr+dr/2.0 #for vr (Mittet)
         vr_now=vr[k,j]
         trr_f1,trr_b1=trr[k,j+1],trr[k,j]
         trz_f1,trz_b1=trz[k,j],trz[k-1,j]
         rho_av=(rhomat[k,j]+rhomat[k,j+1])/2.0 #be carefull

         trr_av=0.5*(trr_f1+trr_b1)
         tpp_av=0.5*(tpp[k,j+1]+tpp[k,j]) #see trr

         vfr_now=vfr[k,j]
         pf_f1,pf_b1=pf[k,j+1],pf[k,j]
         rhof_av=(rhofmat[k,j]+rhofmat[k,j+1])/2.0 #be carefull
         D1_av=(D1mat[k,j]+D1mat[k,j+1])/2.0 #be carefull
         D2_av=(D2mat[k,j]+D2mat[k,j+1])/2.0 #be carefull

         vr[k,j],vfr[k,j]=update_vr_1st_Por(vr_now,vfr_now,
             trr_f1,trr_b1,trz_f1,trz_b1,trr_av,tpp_av,
             pf_f1,pf_b1,
             dz,dr,dt,r_now,rho_av,rhof_av,D1_av,D2_av,Flagmat_vf_zero[k,j])

         r_now=(j-1)*dr #for vphi and vz (Mittet)

         vz_now=vz[k,j]
         trz_f1,trz_b1=trz[k,j],trz[k,j-1]
         tzz_f1,tzz_b1=tzz[k+1,j],tzz[k,j]
         rho_av=(rhomat[k,j]+rhomat[k+1,j])/2.0 #be careful

         trz_av=0.5*(trz_f1+trz_b1)

         vfz_now=vfz[k,j]
         pf_f1,pf_b1=pf[k+1,j],pf[k,j]
         rhof_av=(rhofmat[k,j]+rhofmat[k+1,j])/2.0 #be careful
         D1_av=(D1mat[k,j]+D1mat[k+1,j])/2.0 #be careful
         D2_av=(D2mat[k,j]+D2mat[k+1,j])/2.0 #be careful

         vz[k,j],vfz[k,j]=update_vz_1st_Por(vz_now,vfz_now,
             trz_f1,trz_b1,tzz_f1,tzz_b1,trz_av,
             pf_f1,pf_b1,
             dz,dr,dt,r_now,rho_av,rhof_av,D1_av,D2_av,Flagmat_vf_zero[k,j])

         #implementing PE to E horizontal boundary
         if(Flagmat_vf_zero[k,j]==0 && D1mat[k+1,j]==Inf)
             vz[k,j],vfz[k,j]=update_vz_1st_Por(vz_now,vfz_now,
                 trz_f1,trz_b1,tzz_f1,tzz_b1,trz_av,
                 pf_f1,pf_b1,
                 dz,dr,dt,r_now,rho_av,rhof_av,D1_av,D2_av,1)
         end

     end #k (z)
end #j (r)

 #Surrounding region (exept for edge) is evaluated using 1st order FD
 # k=2 & j=2:nr-1, k=nz-1 & j=2:nr-1, j=nr-1 & k=2:nz-1
 #--> inside PML. Skipped.
end #function


#-----------------------------
# updating stress field
#-----------------------------
function update_trr_1st_Por(trr_now,vr_f,vr_b,vz_f,vz_b,vr_av,
    vfr_f,vfr_b,vfz_f,vfz_b,vfr_av,
    dz,dr,dt,r_now,C,H,G)

    #first-order difference
    #X_now: current time, current position [k,j]
    #X_f: foward position (vr[k,j+1],vz[k+1,j])
    #X_b: backward position (vr[k,j],vz[k,j])
    #_av: averaged values (vr_av)
    #vphi is at the same location
    #r_now: current radial position @trr
    #m: order (0: monopole, 1:dipole)
    drvr=vr_f-vr_b
    drvr=drvr/dr
    dzvz=vz_f-vz_b
    dzvz=dzvz/dz

    drvfr=vfr_f-vfr_b
    drvfr=drvfr/dr
    dzvfz=vfz_f-vfz_b
    dzvfz=dzvfz/dz

    return trr_now+dt*(
                      +(H-2.0*G)*(vr_av/r_now+dzvz)
                      +H*drvr
                      +C*(drvfr+vfr_av/r_now+dzvfz)
                      )
end

function update_tpp_1st_Por(tpp_now,vr_f,vr_b,vz_f,vz_b,vr_av,
    vfr_f,vfr_b,vfz_f,vfz_b,vfr_av,
    dz,dr,dt,r_now,C,H,G)

    #first-order difference
    #X_now: current time, current position [k,j]
    #X_f: foward position (vr[k,j+1],vz[k+1,j])
    #X_b: backward position (vr[k,j],vz[k,j])
    #_av: averaged values (vr_av)
    #vphi is at the same location
    #r_now: current radial position @trr
    #m: order (0: monopole, 1:dipole)
    drvr=vr_f-vr_b
    drvr=drvr/dr
    dzvz=vz_f-vz_b
    dzvz=dzvz/dz

    drvfr=vfr_f-vfr_b
    drvfr=drvfr/dr
    dzvfz=vfz_f-vfz_b
    dzvfz=dzvfz/dz

    return tpp_now+dt*(
                      +(H-2.0*G)*(drvr+dzvz)
                      +H*vr_av/r_now
                      +C*(drvfr+vfr_av/r_now+dzvfz)
                      )
end

function update_tzz_1st_Por(tzz_now,vr_f,vr_b,vz_f,vz_b,vr_av,
    vfr_f,vfr_b,vfz_f,vfz_b,vfr_av,
    dz,dr,dt,r_now,C,H,G)

    #first-order difference
    #X_now: current time, current position [k,j]
    #X_f: foward position (vr[k,j+1],vz[k+1,j])
    #X_b: backward position (vr[k,j],vz[k,j])
    #_av: averaged values (vr_av)
    #vphi is at the same location
    #r_now: current radial position @trr
    #m: order (0: monopole, 1:dipole)
    drvr=vr_f-vr_b
    drvr=drvr/dr
    dzvz=vz_f-vz_b
    dzvz=dzvz/dz

    drvfr=vfr_f-vfr_b
    drvfr=drvfr/dr
    dzvfz=vfz_f-vfz_b
    dzvfz=dzvfz/dz

    return tzz_now+dt*(
                      +(H-2.0*G)*(vr_av/r_now+drvr)
                      +H*dzvz
                      +C*(drvfr+vfr_av/r_now+dzvfz)
                      )
end

function update_trz(___trz_now,___vr_f,___vr_b,___vz_f,___vz_b,
    ___dz,___dr,___dt,___r_now,___mu,___m)
    #first-order difference
    #X_now: current time, current position [k,j]
    #X_f: foward position (vr[k,j],vz[k,j])
    #X_b: backward position (vr[k-1,j],vz[k,j-1])
    #_av: averaged values ()
    #-- is at the same location
    #r_now: current radial position @trz
    #m: order (0: monopole, 1:dipole)
    return ___trz_now+___dt*(
               +___mu*(___vr_f-___vr_b)/___dz #mu*dzvr
               +___mu*(___vz_f-___vz_b)/___dr #mu*drvz
               )
end

function update_pf_1st_Por(pf_now,vr_f,vr_b,vz_f,vz_b,vr_av,
    vfr_f,vfr_b,vfz_f,vfz_b,vfr_av,
    dz,dr,dt,r_now,C,M)
    #first-order difference, 2nd-order accuracy (Graves 1996)
    #X_now: current time, current position [k,j]
    #X_f: foward position (vr[k,j+1],vz[k+1,j]) 1:nearest, 2:next nearest
    #X_b: backward position (vr[k,j],vz[k,j]) 1:nearest, 2:next nearest
    #_av: averaged values (vr_av) --> will be calculated in this routine
    #vphi is at the same location
    #r_now: current radial position @trr
    #m: order (0: monopole, 1:dipole), m=0 only!

    drvr=vr_f-vr_b
    drvr=drvr/dr
    dzvz=vz_f-vz_b
    dzvz=dzvz/dz

    drvfr=vfr_f-vfr_b
    drvfr=drvfr/dr
    dzvfz=vfz_f-vfz_b
    dzvfz=dzvfz/dz


    return pf_now+dt*(
                      -C*(drvr+vr_av/r_now+dzvz)
                      -M*(drvfr+vfr_av/r_now+dzvfz)
                      )

end

function update_stress_1st_Por!(vr,vz,trr,tpp,tzz,trz,vfr,vfz,pf,
    Mmat,Cmat,Hmat,Gmat,nr,nz,dr,dz,dt,LPML_z,LPML_r)
#2nd order FD (Graves, 1996)
#assuming PML!
#@inbounds Threads.@threads
#println("--in")
#@inbounds Threads.@threads for j=3:nr-2
#      for k=3:nz-2
@inbounds Threads.@threads for j=2:nr-LPML_r
      for k=LPML_z+1:nz-LPML_z
         #I am accesing [k,j], but definition of r_now changes with each component!
         r_now=(j-1)*dr #for trr,tpp,tzz,tpz Mittet

         vr_f1,vr_b1=vr[k,j],vr[k,j-1]
         vz_f1,vz_b1=vz[k,j],vz[k-1,j]

         vfr_f1,vfr_b1=vfr[k,j],vfr[k,j-1]
         vfz_f1,vfz_b1=vfz[k,j],vfz[k-1,j]

         H=Hmat[k,j] #be careful
         C=Cmat[k,j] #be careful
         G=Gmat[k,j] #be careful
         M=Mmat[k,j] #be careful

         vr_av=0.5*(vr_f1+vr_b1)
         vfr_av=0.5*(vfr_f1+vfr_b1)

         trr_now=trr[k,j]
         trr[k,j]=update_trr_1st_Por(trr_now,vr_f1,vr_b1,vz_f1,vz_b1,vr_av,
                vfr_f1,vfr_b1,vfz_f1,vfz_b1,vfr_av,
                dz,dr,dt,r_now,C,H,G)

         tpp_now=tpp[k,j]
         tpp[k,j]=update_tpp_1st_Por(tpp_now,vr_f1,vr_b1,vz_f1,vz_b1,vr_av,
                vfr_f1,vfr_b1,vfz_f1,vfz_b1,vfr_av,
                dz,dr,dt,r_now,C,H,G)

         tzz_now=tzz[k,j]
         tzz[k,j]=update_tzz_1st_Por(tzz_now,vr_f1,vr_b1,vz_f1,vz_b1,vr_av,
                vfr_f1,vfr_b1,vfz_f1,vfz_b1,vfr_av,
                dz,dr,dt,r_now,C,H,G)

         pf_now=pf[k,j]
         pf[k,j]=update_pf_1st_Por(pf_now,vr_f1,vr_b1,vz_f1,vz_b1,vr_av,
                vfr_f1,vfr_b1,vfz_f1,vfz_b1,vfr_av,
                dz,dr,dt,r_now,C,M)

         r_now=(j-1)*dr+dr/2.0 #for trp,trz Mittet
         trz_now=trz[k,j]
         #---Special averaging---
#         if (mmat[k,j]+mmat[k+1,j+1]==0.0)
#             mu_av=0.0
#         else
#             mu_av=2.0*mmat[k,j]*mmat[k+1,j+1]/(mmat[k,j]+mmat[k+1,j+1]) #be careful
#         end
         if (Gmat[k,j]*Gmat[k,j+1]*Gmat[k+1,j]*Gmat[k+1,j+1]==0.0)
              G_av=0.0
         else
             gj1=2.0*Gmat[k,j]*Gmat[k+1,j]/(Gmat[k,j]+Gmat[k+1,j])
             gj2=2.0*Gmat[k,j+1]*Gmat[k+1,j+1]/(Gmat[k,j+1]+Gmat[k+1,j+1])
             G_av=2.0*gj1*gj2/(gj1+gj2)
         end
#         println(mu_av)
         #---End special averaging---

         vr_f1,vr_b1=vr[k+1,j],vr[k,j]
         vz_f1,vz_b1=vz[k,j+1],vz[k,j]
         trz[k,j]=update_trz(trz_now,vr_f1,vr_b1,vz_f1,vz_b1,
             dz,dr,dt,r_now,G_av,0)

   end #k (z)
end #j (r)

#Surrounding region (exept for edge) is evaluated using 1st order FD
# k=2 & j=2:nr-1, k=nz-1 & j=2:nr-1, j=nr-1 & k=2:nz-1
#--> inside PML. Skipped.

end #function


#-----------------------------
# Boundary conditions at r=0 (velocity)
#-----------------------------
#Entry function
function ApplyBCLeft_vel!(vr,vz,trr,tpp,tzz,trz,vfr,vfz,pf,
                          rhomat,rhofmat,D1mat,D2mat,Flagmat_vf_zero,
                          nr,nz,dr,dz,dt,
                          Prz_T,Pzz_T,PzzPE_T,
                          Prz_B,Pzz_B,PzzPE_B,
                          LPML_z)

    ApplyBCLeft_velocity_1st_Por!(vr,vz,trr,tpp,tzz,trz,vfr,vfz,pf,
                                  rhomat,rhofmat,D1mat,D2mat,Flagmat_vf_zero,nr,nz,dr,dz,dt,LPML_z)
    ApplyBCLeft_velocity_1st_atPML_Top_Por!(vr,vz,trr,tpp,tzz,trz,vfr,vfz,pf,
                                            rhomat,rhofmat,D1mat,D2mat,Flagmat_vf_zero,nr,nz,dr,dz,dt,
                                            Prz_T,Pzz_T,PzzPE_T,
                                            LPML_z)
    ApplyBCLeft_velocity_1st_atPML_Bottom_Por!(vr,vz,trr,tpp,tzz,trz,vfr,vfz,pf,
                                               rhomat,rhofmat,D1mat,D2mat,Flagmat_vf_zero,nr,nz,dr,dz,dt,
                                               Prz_B,Pzz_B,PzzPE_B,
                                               LPML_z)

end


#-------------
function ApplyBCLeft_velocity_1st_Por!(vr,vz,trr,tpp,tzz,trz,vfr,vfz,pf,
    rhomat,rhofmat,D1mat,D2mat,Flagmat_vf_zero,nr,nz,dr,dz,dt,LPML_z)
#Filling field values at j==1 (r) & k=2:nz-1
#j=1 corresponds to r=0 for variables of vz[], vphi[], tii[], tpz[]. Other values are @dr/2, not r=0
#2nd order FD (Graves, 1996)
#Supplement discussion from Mittet, 1996, Geophysics,61,1,21-33.
# m=0, monopole: vr(r)=-vr(-r),vphi=-vphi(-r),vz(r)=vz(-r), tij(r)=tij(-r)
# m=1, dipole: vr(r)=vr(-r),vphi=vphi(-r),vz(r)=-vz(-r), tij(r)=-tij(-r)
# nonzero@r=0 --> only "symmetric" components, otherwize 0.
#m=0 only!
#new version to ignore PML region


#THIS FUNCTION MIGHT NEED UPDATE IN ORDER TO SEPARATELY APPLY BC AT PML-Z REGION!

#Borehole axis (Left Edge): d/dr=0 (Symmetric BC)
j=1
#Variables@left-edge are vz[], vphi[], tii[], tpz[] Mittet
#Case m=0: vr=0, trz=0 (Randall, Mittet->fluid), trp=0 (Randall, Mittet->fluid)
for k=LPML_z+1:nz-LPML_z #caution !

        r_now=(j-1)*dr+dr/2.0 #for vr (Mittet)

        vr_now=vr[k,j]
        trr_f1,trr_b1=trr[k,j+1],trr[k,j]
        trz_f1,trz_b1=trz[k,j],trz[k-1,j]
        rho_av=(rhomat[k,j]+rhomat[k,j+1])/2.0 #be carefull

        trr_av=0.5*(trr_f1+trr_b1)
        tpp_av=0.5*(tpp[k,j+1]+tpp[k,j]) #see trr

        vfr_now=vfr[k,j]
        pf_f1,pf_b1=pf[k,j+1],pf[k,j]
        rhof_av=(rhofmat[k,j]+rhofmat[k,j+1])/2.0 #be carefull
        D1_av=(D1mat[k,j]+D1mat[k,j+1])/2.0 #be carefull
        D2_av=(D2mat[k,j]+D2mat[k,j+1])/2.0 #be carefull

        vr[k,j],vfr[k,j]=update_vr_1st_Por(vr_now,vfr_now,
            trr_f1,trr_b1,trz_f1,trz_b1,trr_av,tpp_av,
            pf_f1,pf_b1,
            dz,dr,dt,r_now,rho_av,rhof_av,D1_av,D2_av,Flagmat_vf_zero[k,j])

        r_now=(j-1)*dr #for vphi and vz (Mittet)

        #----THIS PART NEEDS CONSIDERATION FOR Poroelasticity
        #--vz requires special attention (see Mittet)
        # D2*vfz+(D1-rhof^2/rho)*dvfz/dt=-dp/dz-rhof/rho*A, A=2*dtrz/dr+dtzz/dz (l'Hopital)
        # Mittet says trz(r)=trz(-r) -> dtrz/dr=0
        # Randall says trz(r=0)=0 -> dtrz/dr¥=0
        # when fluid is around r=0, they are same

        vz_now=vz[k,j]
        tzz_f1,tzz_b1=tzz[k+1,j],tzz[k,j]
        rho_av=(rhomat[k,j]+rhomat[k+1,j])/2.0 #be careful

        dztzz=(tzz_f1-tzz_b1)
        dztzz=dztzz/dz

        #Experimental (Randall type)
        trz_f1,trz_b1=trz[k,j],-trz[k,j] #symtr (assume Randall's trz=0. not Mittet sym)

        drtrz=(trz_f1-trz_b1)
        drtrz=drtrz/dr

        A=2.0*drtrz+dztzz

        vfz_now=vfz[k,j]
        pf_f1,pf_b1=pf[k+1,j],pf[k,j]
        rhof_av=(rhofmat[k,j]+rhofmat[k+1,j])/2.0 #be careful
        D1_av=(D1mat[k,j]+D1mat[k+1,j])/2.0 #be careful
        D2_av=(D2mat[k,j]+D2mat[k+1,j])/2.0 #be careful

        dzpf=(pf_f1-pf_b1)
        dzpf=dzpf/dz

        vfz_old=vfz_now
        if(Flagmat_vf_zero[k,j]==1)
            vfz_now=0.0
        else
        vfz_now=(D2_av/2+(D1_av-rhof_av^2/rho_av)/dt)^(-1) * (
                    (-D2_av/2+(D1_av-rhof_av^2/rho_av)/dt)*vfz_old-dzpf-rhof_av/rho_av*A
                    )
        end
        dert_vfz=(vfz_now-vfz_old)/dt
        vfz[k,j]=vfz_now

        vz[k,j]=vz_now+dt/rho_av*(A-rhof_av*dert_vfz)
        #----END THIS PART NEEDS CONSIDERATION FOR Poroelasticity

end #k (z)

end #function


#-----------------------------
# Boundary conditions at r=0 (stress)
#-----------------------------
#Entry function
function ApplyBCLeft_stress!(vr,vz,trr,tpp,tzz,trz,vfr,vfz,pf,
                             M,C,H,G,nr,nz,dr,dz,dt,
                             Rz_T,Srz_T,RzPE_T,
                             Rz_B,Srz_B,RzPE_B,
                             LPML_z)


    ApplyBCLeft_stress_1st_Por!(vr,vz,trr,tpp,tzz,trz,vfr,vfz,pf,
                                M,C,H,G,nr,nz,dr,dz,dt,LPML_z)
    ApplyBCLeft_stress_1st_atPML_Top_Por!(vr,vz,trr,tpp,tzz,trz,vfr,vfz,pf,
                                          M,C,H,G,nr,nz,dr,dz,dt,
                                          Rz_T,Srz_T,RzPE_T,
                                          LPML_z)
    ApplyBCLeft_stress_1st_atPML_Bottom_Por!(vr,vz,trr,tpp,tzz,trz,vfr,vfz,pf,
                                             M,C,H,G,nr,nz,dr,dz,dt,
                                             Rz_B,Srz_B,RzPE_B,
                                             LPML_z)
    
end


#---
function ApplyBCLeft_stress_1st_Por!(vr,vz,trr,tpp,tzz,trz,vfr,vfz,pf,
    Mmat,Cmat,Hmat,Gmat,nr,nz,dr,dz,dt,LPML_z)
#Filling field values at j==1 (r) & k=2:nz-1
#j=1 corresponds to r=0 for variables of vz[], vphi[], tii[], tpz[]. Other values are @dr/2, not r=0
#2nd order FD (Graves, 1996)
#Supplement discussion from Mittet, 1996, Geophysics,61,1,21-33.
# m=0, monopole: vr(r)=-vr(-r),vphi=-vphi(-r),vz(r)=vz(-r), tij(r)=tij(-r)
# m=1, dipole: vr(r)=vr(-r),vphi=vphi(-r),vz(r)=-vz(-r), tij(r)=-tij(-r)
# nonzero@r=0 --> only "symmetric" components, otherwize 0.
#ignore PML region

#Borehole axis (Left Edge): d/dr=0 (Symmetric BC)
j=1
#Variables@left-edge are vz[], vphi[], tii[], tpz[] Mittet
#Case m=0: vr=0, trz=0 (Randall, Mittet->fluid), trp=0 (Randall, Mittet->fluid)
for k=LPML_z+1:nz-LPML_z

    #trz[k,j]=0.0 #Not on Left edge
    #trp[k,j]=0.0 #Not on Left edge

        r_now=(j-1)*dr #for trr,tpp,tzz,tpz Mittet

        trr_now=trr[k,j]
#            vr_f1,vr_b1=vr[k,j],vr[k,j-1]
        vr_f1,vr_b1=vr[k,j],-vr[k,j]

        vz_f1,vz_b1=vz[k,j],vz[k-1,j]

        vfr_f1,vfr_b1=vfr[k,j],-vfr[k,j] #symmtr
        vfz_f1,vfz_b1=vfz[k,j],vfz[k-1,j]

        H=Hmat[k,j] #be careful
        C=Cmat[k,j] #be careful
        G=Gmat[k,j] #be careful
        M=Mmat[k,j] #be careful

        #Special treatment @r=0
        #trr=tzz=tpp=lambda*(2*dr/dr+dvz/dz) FLUID

        drvr=(vr_f1-vr_b1)
        drvr=drvr/dr
        dzvz=(vz_f1-vz_b1)
        dzvz=dzvz/dz

        drvfr=(vfr_f1-vfr_b1)
        drvfr=drvfr/dr
        dzvfz=(vfz_f1-vfz_b1)
        dzvfz=dzvfz/dz


        trr[k,j]=trr[k,j]+dt*(
                             +(H-2.0*G)*(drvr+dzvz)
                             +H*drvr
                             +C*(2.0*drvfr+dzvfz)
                             )
        tpp[k,j]=tpp[k,j]+dt*(
                            +(H-2.0*G)*(drvr+dzvz)
                            +H*drvr
                            +C*(2.0*drvfr+dzvfz)
                            )
        tzz[k,j]=tzz[k,j]+dt*(
                            +(H-2.0*G)*2.0*drvr
                            +H*dzvz
                            +C*(2.0*drvfr+dzvfz)
                            )
        pf[k,j]=pf[k,j]+dt*(
                            -C*(2.0*drvr+dzvz)
                            -M*(2.0*drvfr+dzvfz)
                            )

        #---trp and trz not on left edge!
        r_now=(j-1)*dr+dr/2.0 #for trp,trz Mittet
        trz_now=trz[k,j]
        #---Special averaging---
#         if (mmat[k,j]+mmat[k+1,j+1]==0.0)
#             mu_av=0.0
#         else
#             mu_av=2.0*mmat[k,j]*mmat[k+1,j+1]/(mmat[k,j]+mmat[k+1,j+1]) #be careful
#         end
        if (Gmat[k,j]*Gmat[k,j+1]*Gmat[k+1,j]*Gmat[k+1,j+1]==0.0)
             G_av=0.0
        else
            gj1=2.0*Gmat[k,j]*Gmat[k+1,j]/(Gmat[k,j]+Gmat[k+1,j])
            gj2=2.0*Gmat[k,j+1]*Gmat[k+1,j+1]/(Gmat[k,j+1]+Gmat[k+1,j+1])
            G_av=2.0*gj1*gj2/(gj1+gj2)
        end
#         println(mu_av)
        #---End special averaging---

        vr_f1,vr_b1=vr[k+1,j],vr[k,j]
        vz_f1,vz_b1=vz[k,j+1],vz[k,j]
        trz[k,j]=update_trz(trz_now,vr_f1,vr_b1,vz_f1,vz_b1,
            dz,dr,dt,r_now,G_av,0)


 end #k (z)

end #function

#-----------------------------
# Boundary conditions in acoustic media
#-----------------------------
function ApplyBC_stress_AcousticMedia_TEST!(trr,tpp,tzz,trz,pf,Flag_AC,nr,nz)
    #when FlagAC==1, then set pf=-mean(tii) and tau=-pf*eye
    @inbounds Threads.@threads for j=1:nr
        for k=1:nz
            if(Flag_AC[k,j]==1)
                pf[k,j]=-(trr[k,j]+tpp[k,j]+tzz[k,j])/3
                trr[k,j]=-pf[k,j]
                tpp[k,j]=-pf[k,j]
                tzz[k,j]=-pf[k,j]
            end
        end
    end
end

#-----------------------------
# Boundary conditions in elastic media
#-----------------------------
function ApplyBC_stress_ElasticMedia_Ou!(pf,Flag_E,nr,nz)
    #when FlagE==1, then set pf=0
    @inbounds Threads.@threads for j=1:nr
        for k=1:nz
            if(Flag_E[k,j]==1)
                pf[k,j]=0.0 #Ou's
            end
        end
    end
end


#-----------------------------
# Boundary conditions at r=R (assuming plane-wave propagation)
#-----------------------------
#Entry function
function ApplyBCRight_vel!(vr,vz,
                           trr,tpp,tzz,trz,
                           vfr,vfz,pf,
                           Rho,Rhof,D1,D2,Flag_vf_zero,nr,nz,dr,dz,dt)


    ApplyBCRight_velocity1D_Por01!(0.0,vr,vz,
                                   trr,tpp,tzz,trz,
                                   vfr,vfz,pf,
                                   Rho,Rhof,D1,D2,Flag_vf_zero,nr,nz,dr,dz,dt)

end



function ApplyBCRight_stress!(vr,vz,
                              trz,
                              Gmat,nr,nz,dr,dz,dt)


    ApplyBCRight_stress1D_Por01!(0.0,vr,vz,
                                 trz,
                                 Gmat,nr,nz,dr,dz,dt)
end

#-
function ApplyBCRight_velocity1D_Por01!(origin_r,vr,vz,
    trr,tpp,tzz,trz,
    vfr,vfz,pf,
    Rho,Rhof,D1,D2,Flag_vf_zero,nr,nz,dr,dz,dt)
    #updating velocity on Right most edge (it was Neumann BC), assuming only 1D wave propagation
    #testing for plane wave inciden ce
    #assuming LVTS, region 2 (may be this assumption is not nessesary)

    #Por01: Poroelatic version (Note: 1st order)

    #nonzero values at j=nr are
    #j=nr->vphi,vz (1st), vr (zero)
    #Make sure regular [update] function has been filled till the right edge!


    #Following flag is 1 when no 1D assumption, but simply zeroing blank component (i.e., vr)
    flag_zero=1

if (flag_zero==1)
    j=nr
@inbounds for k=1:nz
    vr[k,j]=0.0
    vfr[k,j]=0.0
end

else

    j=nr
@inbounds for k=2:nz-1
    r_now=(j-1)*dr+dr/2.0 #for vr (Mittet)
    r_now=r_now+origin_r #for LVTS, region 2


    #---option A: 1D elastic wave, i.e., vr is same as elastic, and vfr is zero
    vr_now=vr[k,j]
    trz_f1,trz_b1=trz[k,j],trz[k-1,j]
    dztrz=(trz_f1-trz_b1)
    dztrz=dztrz/dz

    vr[k,j]=vr_now+dt/Rho[k,j]*(
                +(trr[k,j]-tpp[k,j])/r_now
                +dztrz
                )
    vfr[k,j]=0.0

  end
end #if
end #function


function ApplyBCRight_stress1D_Por01!(origin_r,vr,vz,
    trz,
    Gmat,nr,nz,dr,dz,dt)
    #updating stress on Right most edge (it was Neumann BC), assuming only 1D wave propagation
    #testing for plane wave inciden ce
    #assuming LVTS, region 2 (may be this assumption is not nessesary)

    #Por01: Poroelasic version (same as elastic, but reduced number of variables)

    #nonzero values @j=nr are
    #j=nr->tii,tpz (1st), trp,trz (zero)
    #Make sure regular [update] function has been filled till the right edge!

    #Following flag is 1 when no 1D assumption, but simply zeroing blank component (i.e., trp,trz)
    flag_zero=1

if (flag_zero==1)
    j=nr
@inbounds for k=3:nz-2
    trz[k,j]=0.0
end

else

    j=nr
@inbounds for k=1:nz

    r_now=(j-1)*dr+dr/2.0 #for trp,trz Mittet
    r_now=r_now+origin_r #for LVTS, region 2

    trz_now=trz[k,j]

    vr_f1,vr_b1=vr[k+1,j],vr[k,j]
    dzvr=(vr_f1-vr_b1)
    dzvr=dzvr/dz
    trz[k,j]=trz_now+dt*(
               +G[k,j]*dzvr #mu*dzvr
               )
  end

 end #if

end #function


#------------------------------------------------------------
# main entry functions for PML calculations
#------------------------------------------------------------
#PML: save velocity at previous step
function PML_save_vel!(memT_vr,memT_vz,memT_vfr,memT_vfz,
                       memB_vr,memB_vz,memB_vfr,memB_vfz,
                       memR_vr,memR_vz,memR_vfr,memR_vfz,
                       vr,vz,vfr,vfz,nr,nz,LPML_z,LPML_r)
    PML_save_vel_Top_Por!(memT_vr,memT_vz,memT_vfr,memT_vfz,vr,vz,vfr,vfz,nr,nz,LPML_z,LPML_r)
    PML_save_vel_Bottom_Por!(memB_vr,memB_vz,memB_vfr,memB_vfz,vr,vz,vfr,vfz,nr,nz,LPML_z,LPML_r)
    PML_save_vel_Right_Por!(memR_vr,memR_vz,memR_vfr,memR_vfz,vr,vz,vfr,vfz,nr,nz,LPML_z,LPML_r)
end 

#PML: save stress at previous step
function PML_save_stress!(memT_trr,memT_tpp,memT_tzz,memT_trz,memT_pf,
                          memB_trr,memB_tpp,memB_tzz,memB_trz,memB_pf,
                          memR_trr,memR_tpp,memR_tzz,memR_trz,memR_pf,
                          trr,tpp,tzz,trz,pf,nr,nz,LPML_z,LPML_r)

    PML_save_stress_Top_Por!(memT_trr,memT_tpp,memT_tzz,memT_trz,memT_pf,
                             trr,tpp,tzz,trz,pf,nr,nz,LPML_z,LPML_r)
    PML_save_stress_Bottom_Por!(memB_trr,memB_tpp,memB_tzz,memB_trz,memB_pf,
                                trr,tpp,tzz,trz,pf,nr,nz,LPML_z,LPML_r)
    PML_save_stress_Right_Por!(memR_trr,memR_tpp,memR_tzz,memR_trz,memR_pf,
                               trr,tpp,tzz,trz,pf,nr,nz,LPML_z,LPML_r)
end

#PML: update velocity
function PML_update_vel!(vr,vz,trr,tpp,tzz,trz,vfr,vfz,pf,
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
                       
    PML_update_velocity_1st_Top_Por!(vr,vz,trr,tpp,tzz,trz,vfr,vfz,pf,
                                     Rho,Rhof,D1,D2,Flag_vf_zero,nr,nz,dr,dz,dt,Prz_T,Pzz_T,PzzPE_T,LPML_z,LPML_r)
    PML_update_velocity_1st_Bottom_Por!(vr,vz,trr,tpp,tzz,trz,vfr,vfz,pf,
                                        Rho,Rhof,D1,D2,Flag_vf_zero,nr,nz,dr,dz,dt,Prz_B,Pzz_B,PzzPE_B,LPML_z,LPML_r)
    PML_update_velocity_1st_Right_Por!(vr,vz,trr,tpp,tzz,trz,vfr,vfz,pf,
                                       Rho,Rhof,D1,D2,Flag_vf_zero,nr,nz,dr,dz,dt,
                                       Prr_R,Qrp_R,Pzr_R,Qzp_R,Rr_R,Rp_R,Rrz_R,
                                       PrrPE_R,RrPE_R,RpPE_R,
                                       LPML_z,LPML_r)
    PML_update_velocity_1st_TopRight_Por!(vr,vz,trr,tpp,tzz,trz,vfr,vfz,pf,
                                          Rho,Rhof,D1,D2,Flag_vf_zero,nr,nz,dr,dz,dt,
                                          Prz_TR,Pzz_TR,PzzPE_TR,
                                          Prr_TR,Qrp_TR,Pzr_TR,Qzp_TR,Rr_TR,Rp_TR,Rrz_TR,
                                          PrrPE_TR,RrPE_TR,RpPE_TR,
                                          LPML_z,LPML_r)
    PML_update_velocity_1st_BottomRight_Por!(vr,vz,trr,tpp,tzz,trz,vfr,vfz,pf,
                                             Rho,Rhof,D1,D2,Flag_vf_zero,nr,nz,dr,dz,dt,
                                             Prz_BR,Pzz_BR,PzzPE_BR,
                                             Prr_BR,Qrp_BR,Pzr_BR,Qzp_BR,Rr_BR,Rp_BR,Rrz_BR,
                                             PrrPE_BR,RrPE_BR,RpPE_BR,
                                             LPML_z,LPML_r)
end

#PML: update memory variables for stress (Rx and Sxx) using velocity at two time steps
function PML_update_memRS!(Rz_T,Srz_T,RzPE_T,
                           Rz_B,Srz_B,RzPE_B,
                           Rr_R,Rp_R,Rrz_R,RrPE_R,RpPE_R,
                           Rr_TR,Rp_TR,Rz_TR,Rrz_TR,Srz_TR,RrPE_TR,RpPE_TR,RzPE_TR,
                           Rr_BR,Rp_BR,Rz_BR,Rrz_BR,Srz_BR,RrPE_BR,RpPE_BR,RzPE_BR,
                           memT_vr,memT_vz,memT_vfr,memT_vfz,
                           memB_vr,memB_vz,memB_vfr,memB_vfz,
                           memR_vr,memR_vz,memR_vfr,memR_vfz,
                           vr,vz,vfr,vfz,nr,nz,dr,dz,dt,LPML_z,LPML_r,PML_Wz,PML_Wr,PML_IWr,PML_Wz2,PML_Wr2,PML_IWr2)
    
    PML_update_memRS_1st_Top_Por!(Rz_T,Srz_T,RzPE_T,
                                  memT_vr,memT_vz,memT_vfr,memT_vfz,
                                  vr,vz,vfr,vfz,nr,nz,dr,dz,dt,LPML_z,LPML_r,PML_Wz,PML_Wz2)
    PML_update_memRS_1st_Bottom_Por!(Rz_B,Srz_B,RzPE_B,
                                     memB_vr,memB_vz,memB_vfr,memB_vfz,
                                     vr,vz,vfr,vfz,nr,nz,dr,dz,dt,LPML_z,LPML_r,PML_Wz,PML_Wz2)
    PML_update_memRS_1st_Right_Por!(Rr_R,Rp_R,Rrz_R,RrPE_R,RpPE_R,
                                    memR_vr,memR_vz,memR_vfr,memR_vfz,
                                    vr,vz,vfr,vfz,nr,nz,dr,dz,dt,LPML_z,LPML_r,PML_Wr,PML_IWr,PML_Wr2,PML_IWr2)
    PML_update_memRS_1st_TopRight_Por!(Rr_TR,Rp_TR,Rz_TR,Rrz_TR,Srz_TR,RrPE_TR,RpPE_TR,RzPE_TR,
                                       memT_vr,memT_vz,memT_vfr,memT_vfz,
                                       memR_vr,memR_vz,memR_vfr,memR_vfz,
                                       vr,vz,vfr,vfz,nr,nz,dr,dz,dt,LPML_z,LPML_r,PML_Wz,PML_Wr,PML_IWr,PML_Wz2,PML_Wr2,PML_IWr2)
    PML_update_memRS_1st_BottomRight_Por!(Rr_BR,Rp_BR,Rz_BR,Rrz_BR,Srz_BR,RrPE_BR,RpPE_BR,RzPE_BR,
                                          memB_vr,memB_vz,memB_vfr,memB_vfz,
                                          memR_vr,memR_vz,memR_vfr,memR_vfz,
                                          vr,vz,vfr,vfz,nr,nz,dr,dz,dt,LPML_z,LPML_r,PML_Wz,PML_Wr,PML_IWr,PML_Wz2,PML_Wr2,PML_IWr2)
end


#PML: stress update
function PML_update_stress!(vr,vz,trr,tpp,tzz,trz,vfr,vfz,pf,
                            M,C,H,G,nr,nz,dr,dz,dt,
                            Rz_T,Srz_T,RzPE_T,
                            Rz_B,Srz_B,RzPE_B,
                            Rr_R,Rp_R,Rrz_R,RrPE_R,RpPE_R,
                            Rr_TR,Rp_TR,Rz_TR,Rrz_TR,Srz_TR,RrPE_TR,RpPE_TR,RzPE_TR,
                            Rr_BR,Rp_BR,Rz_BR,Rrz_BR,Srz_BR,RrPE_BR,RpPE_BR,RzPE_BR,
                            LPML_z,LPML_r)


    PML_update_stress_1st_Top_Por!(vr,vz,trr,tpp,tzz,trz,vfr,vfz,pf,
                                   M,C,H,G,nr,nz,dr,dz,dt,Rz_T,Srz_T,RzPE_T,LPML_z,LPML_r)
    PML_update_stress_1st_Bottom_Por!(vr,vz,trr,tpp,tzz,trz,vfr,vfz,pf,
                                      M,C,H,G,nr,nz,dr,dz,dt,Rz_B,Srz_B,RzPE_B,LPML_z,LPML_r)
    PML_update_stress_1st_Right_Por!(vr,vz,trr,tpp,tzz,trz,vfr,vfz,pf,
                                     M,C,H,G,nr,nz,dr,dz,dt,
                                     Rr_R,Rp_R,Rrz_R,RrPE_R,RpPE_R,
                                     LPML_z,LPML_r)
    PML_update_stress_1st_TopRight_Por!(vr,vz,trr,tpp,tzz,trz,vfr,vfz,pf,
                                        M,C,H,G,nr,nz,dr,dz,dt,
                                        Rr_TR,Rp_TR,Rz_TR,Rrz_TR,Srz_TR,RrPE_TR,RpPE_TR,RzPE_TR,
                                        LPML_z,LPML_r)
    PML_update_stress_1st_BottomRight_Por!(vr,vz,trr,tpp,tzz,trz,vfr,vfz,pf,
                                           M,C,H,G,nr,nz,dr,dz,dt,
                                           Rr_BR,Rp_BR,Rz_BR,Rrz_BR,Srz_BR,RrPE_BR,RpPE_BR,RzPE_BR,
                                           LPML_z,LPML_r)
end


#PML: update memory variables for velocity (Pxx and Qxx) using stress at two time steps
function PML_update_memPQ!(Prz_T,Pzz_T,PzzPE_T,
                           Prz_B,Pzz_B,PzzPE_B,
                           Prr_R,Qrp_R,Pzr_R,Qzp_R,PrrPE_R,
                           Prr_TR,Qrp_TR,Prz_TR,Pzr_TR,Qzp_TR,Pzz_TR,PrrPE_TR,PzzPE_TR,
                           Prr_BR,Qrp_BR,Prz_BR,Pzr_BR,Qzp_BR,Pzz_BR,PrrPE_BR,PzzPE_BR,
                           memT_trr,memT_tpp,memT_tzz,memT_trz,memT_pf,
                           memB_trr,memB_tpp,memB_tzz,memB_trz,memB_pf,
                           memR_trr,memR_tpp,memR_tzz,memR_trz,memR_pf,
                           trr,tpp,tzz,trz,pf,nr,nz,dr,dz,dt,LPML_z,LPML_r,PML_Wz,PML_Wr,PML_IWr,PML_Wz2,PML_Wr2,PML_IWr2)


    PML_update_memPQ_1st_Top_Por!(Prz_T,Pzz_T,PzzPE_T,
                                  memT_trr,memT_tpp,memT_tzz,memT_trz,memT_pf,
                                  trr,tpp,tzz,trz,pf,nr,nz,dr,dz,dt,LPML_z,LPML_r,PML_Wz,PML_Wz2)
    PML_update_memPQ_1st_Bottom_Por!(Prz_B,Pzz_B,PzzPE_B,
                                     memB_trr,memB_tpp,memB_tzz,memB_trz,memB_pf,
                                     trr,tpp,tzz,trz,pf,nr,nz,dr,dz,dt,LPML_z,LPML_r,PML_Wz,PML_Wz2)
    PML_update_memPQ_1st_Right_Por!(Prr_R,Qrp_R,Pzr_R,Qzp_R,PrrPE_R,
                                    memR_trr,memR_tpp,memR_tzz,memR_trz,memR_pf,
                                    trr,tpp,tzz,trz,pf,nr,nz,dr,dz,dt,LPML_z,LPML_r,PML_Wr,PML_IWr,PML_Wr2,PML_IWr2)
    PML_update_memPQ_1st_TopRight_Por!(Prr_TR,Qrp_TR,Prz_TR,Pzr_TR,Qzp_TR,Pzz_TR,PrrPE_TR,PzzPE_TR,
                                       memT_trr,memT_tpp,memT_tzz,memT_trz,memT_pf,
                                       memR_trr,memR_tpp,memR_tzz,memR_trz,memR_pf,
                                       trr,tpp,tzz,trz,pf,nr,nz,dr,dz,dt,LPML_z,LPML_r,PML_Wz,PML_Wr,PML_IWr,PML_Wz2,PML_Wr2,PML_IWr2)
    PML_update_memPQ_1st_BottomRight_Por!(Prr_BR,Qrp_BR,Prz_BR,Pzr_BR,Qzp_BR,Pzz_BR,PrrPE_BR,PzzPE_BR,
                                          memB_trr,memB_tpp,memB_tzz,memB_trz,memB_pf,
                                          memR_trr,memR_tpp,memR_tzz,memR_trz,memR_pf,
                                          trr,tpp,tzz,trz,pf,nr,nz,dr,dz,dt,LPML_z,LPML_r,PML_Wz,PML_Wr,PML_IWr,PML_Wz2,PML_Wr2,PML_IWr2)
end




#------------------------------
# initializing PML variables
# (16 additional variables)
#------------------------------
function init_PML_Por(LPML_r,LPML_z,nr,nz)
 #LPML_ : PML thickness (points)

 #Top and Bottom (z streching only)
 Prz_T=zeros(LPML_z,nr-LPML_r)
 Pzz_T=zeros(LPML_z,nr-LPML_r)
 Rz_T=zeros(LPML_z,nr-LPML_r)
 Srz_T=zeros(LPML_z,nr-LPML_r)
 PzzPE_T=zeros(LPML_z,nr-LPML_r)
 RzPE_T=zeros(LPML_z,nr-LPML_r)

 Prz_B=zeros(LPML_z,nr-LPML_r)
 Pzz_B=zeros(LPML_z,nr-LPML_r)
 Rz_B=zeros(LPML_z,nr-LPML_r)
 Srz_B=zeros(LPML_z,nr-LPML_r)
 PzzPE_B=zeros(LPML_z,nr-LPML_r)
 RzPE_B=zeros(LPML_z,nr-LPML_r)

 #Right (r streching only)
 Prr_R=zeros(nz-2*LPML_z,LPML_r)
 Qrp_R=zeros(nz-2*LPML_z,LPML_r)
 Pzr_R=zeros(nz-2*LPML_z,LPML_r)
 Qzp_R=zeros(nz-2*LPML_z,LPML_r)
 Rr_R=zeros(nz-2*LPML_z,LPML_r)
 Rp_R=zeros(nz-2*LPML_z,LPML_r)
 Rrz_R=zeros(nz-2*LPML_z,LPML_r)
 PrrPE_R=zeros(nz-2*LPML_z,LPML_r)
 RrPE_R=zeros(nz-2*LPML_z,LPML_r)
 RpPE_R=zeros(nz-2*LPML_z,LPML_r)

 #TopRight and BottomRight (both z and r streching)
 Prr_TR=zeros(LPML_z,LPML_r)
 Qrp_TR=zeros(LPML_z,LPML_r)
 Prz_TR=zeros(LPML_z,LPML_r)
 Pzr_TR=zeros(LPML_z,LPML_r)
 Qzp_TR=zeros(LPML_z,LPML_r)
 Pzz_TR=zeros(LPML_z,LPML_r)
 Rr_TR=zeros(LPML_z,LPML_r)
 Rp_TR=zeros(LPML_z,LPML_r)
 Rz_TR=zeros(LPML_z,LPML_r)
 Rrz_TR=zeros(LPML_z,LPML_r)
 Srz_TR=zeros(LPML_z,LPML_r)
 PrrPE_TR=zeros(LPML_z,LPML_r)
 PzzPE_TR=zeros(LPML_z,LPML_r)
 RrPE_TR=zeros(LPML_z,LPML_r)
 RpPE_TR=zeros(LPML_z,LPML_r)
 RzPE_TR=zeros(LPML_z,LPML_r)

 Prr_BR=zeros(LPML_z,LPML_r)
 Qrp_BR=zeros(LPML_z,LPML_r)
 Prz_BR=zeros(LPML_z,LPML_r)
 Pzr_BR=zeros(LPML_z,LPML_r)
 Qzp_BR=zeros(LPML_z,LPML_r)
 Pzz_BR=zeros(LPML_z,LPML_r)
 Rr_BR=zeros(LPML_z,LPML_r)
 Rp_BR=zeros(LPML_z,LPML_r)
 Rz_BR=zeros(LPML_z,LPML_r)
 Rrz_BR=zeros(LPML_z,LPML_r)
 Srz_BR=zeros(LPML_z,LPML_r)
 PrrPE_BR=zeros(LPML_z,LPML_r)
 PzzPE_BR=zeros(LPML_z,LPML_r)
 RrPE_BR=zeros(LPML_z,LPML_r)
 RpPE_BR=zeros(LPML_z,LPML_r)
 RzPE_BR=zeros(LPML_z,LPML_r)

#----
return Prz_T,Pzz_T,Rz_T,Srz_T,PzzPE_T,RzPE_T,
Prz_B,Pzz_B,Rz_B,Srz_B,PzzPE_B,RzPE_B,
Prr_R,Qrp_R,Pzr_R,Qzp_R,Rr_R,Rp_R,Rrz_R,PrrPE_R,RrPE_R,RpPE_R,
Prr_TR,Qrp_TR,Prz_TR,Pzr_TR,Qzp_TR,Pzz_TR,Rr_TR,Rp_TR,Rz_TR,Rrz_TR,Srz_TR,PrrPE_TR,PzzPE_TR,RrPE_TR,RpPE_TR,RzPE_TR,
Prr_BR,Qrp_BR,Prz_BR,Pzr_BR,Qzp_BR,Pzz_BR,Rr_BR,Rp_BR,Rz_BR,Rrz_BR,Srz_BR,PrrPE_BR,PzzPE_BR,RrPE_BR,RpPE_BR,RzPE_BR

end


#------------------------------
# creating PML profile function
#------------------------------
function init_PML_profile(LPML_r,LPML_z,Vmax,dr,dz,nr,f0)
 #LPML_ : PML thickness (points)
 #Vmax: maximum velocity

 #Try to look at PWL_Wz: when slope is very sharp reflection can occur.
 #                     : The highest point is controled by Ralpha,
 #                     : slope continuing from non-PML to PML is controlled by a, b
 #                     : large Ralpha -> too sharp. large a -> continuation may be discontinous

 #---Wang and Tang 2003 eqs (42-43)
# a=0.25
# b=0.75

 #==--default---
 a=0.0
 b=1.0
 c=0.0
# Ralpha=1E-8 #required reflection magnitude order: ~log()=3pi(m+1)? m is order of OMEGA function

 Ralpha=1E-3
    
 dvec_r=collect(1:1:LPML_r)
 dvec_z=collect(1:1:LPML_z)
 dvec_r=dvec_r*dr-dr*ones(length(dvec_r)) #0,dr,2dr,3dr... do not change!
 dvec_z=dvec_z*dz-dz*ones(length(dvec_z))

 LPML_r_meter=(LPML_r-1)*dr
 LPML_z_meter=(LPML_z-1)*dz

#0,dr,2dr,3dr...
 PML_Wr=-Vmax*log(Ralpha)/LPML_r_meter*(a*dvec_r/LPML_r_meter+b*dvec_r.^2/LPML_r_meter^2+c*dvec_r.^3/LPML_r_meter^3)
 PML_Wz=-Vmax*log(Ralpha)/LPML_z_meter*(a*dvec_z/LPML_z_meter+b*dvec_z.^2/LPML_z_meter^2+c*dvec_z.^3/LPML_z_meter^3)
#dr/2,dr/2+dr,dr/2+2dr...
 dvec_r2=collect(1:1:LPML_r)
 dvec_z2=collect(1:1:LPML_z)
 dvec_r2=dvec_r2*dr-dr/2.0*ones(length(dvec_r2))
 dvec_z2=dvec_z2*dz-dz/2.0*ones(length(dvec_z2))
 PML_Wr2=-Vmax*log(Ralpha)/LPML_r_meter*(a*dvec_r2/LPML_r_meter+b*dvec_r2.^2/LPML_r_meter^2+c*dvec_r2.^3/LPML_r_meter^3)
 PML_Wz2=-Vmax*log(Ralpha)/LPML_z_meter*(a*dvec_z2/LPML_z_meter+b*dvec_z2.^2/LPML_z_meter^2+c*dvec_z2.^3/LPML_z_meter^3)

#Analytical integral assuming specific form PML_Wr=P0(ar/L+br^2/L^2)
#0,dr,2dr,3dr... from r0
r0=(nr-LPML_r)*dr #PML just starts at this R (Wr=0): location of tii
dvec_r_global=dvec_r+r0*ones(length(dvec_r)) #r0+0, r0+dr, r0+2dr, ...
P0=-Vmax*log(Ralpha)/LPML_r_meter
PML_IWr=P0*(1.0/2.0*a/LPML_r_meter*(dvec_r).^2+
            1.0/3.0*b/LPML_r_meter^2*(dvec_r).^3+
            1.0/4.0*c/LPML_r_meter^3*(dvec_r).^4
            )./dvec_r_global

#dr/2,dr/2+dr,dr/2+2dr... from r0
#  PML_IWr2=-Vmax*log(Ralpha)/LPML_r_meter*(1.0/2.0*a*dvec_r2/LPML_r_meter+1.0/3.0*b*dvec_r2.^2/LPML_r_meter^2)
dvec_r2_global=dvec_r2+r0*ones(length(dvec_r))
PML_IWr2=P0*(1.0/2.0*a/LPML_r_meter*(dvec_r2).^2+
            1.0/3.0*b/LPML_r_meter^2*(dvec_r2).^3+
            1.0/4.0*b/LPML_r_meter^3*(dvec_r2).^4
            )./dvec_r2_global

  ==#

    #--test---Komatisch profile, N=2
    #PML_Wr=d0(r^2/L^2)-amax(r/L-1)
    LPML_r_meter=(LPML_r-1)*dr
    LPML_z_meter=(LPML_z-1)*dz
    Ralpha=1E-3 #Komatitsch
#    Ralpha=1E-6 #
    d0_r=-3*Vmax*log(Ralpha)/(2*LPML_r_meter)
    d0_z=-3*Vmax*log(Ralpha)/(2*LPML_z_meter)
    amax=pi*f0
        
    dvec_r=collect(1:1:LPML_r)
    dvec_z=collect(1:1:LPML_z)
    dvec_r=dvec_r*dr-dr*ones(length(dvec_r)) #0,dr,2dr,3dr... do not change!
    dvec_z=dvec_z*dz-dz*ones(length(dvec_z))


    #0,dr,2dr,3dr...
    PML_Wr=d0_r*(dvec_r.^2/LPML_r_meter^2)-amax*(dvec_r/LPML_r_meter-ones(size(dvec_r)))
    PML_Wz=d0_r*(dvec_z.^2/LPML_z_meter^2)-amax*(dvec_z/LPML_z_meter-ones(size(dvec_z)))
    #dr/2,dr/2+dr,dr/2+2dr...
    dvec_r2=collect(1:1:LPML_r)
    dvec_z2=collect(1:1:LPML_z)
    dvec_r2=dvec_r2*dr-dr/2.0*ones(length(dvec_r2))
    dvec_z2=dvec_z2*dz-dz/2.0*ones(length(dvec_z2))
    PML_Wr2=d0_r*(dvec_r2.^2/LPML_r_meter^2)-amax*(dvec_r2/LPML_r_meter-ones(size(dvec_r2)))
    PML_Wz2=d0_r*(dvec_z2.^2/LPML_z_meter^2)-amax*(dvec_z2/LPML_z_meter-ones(size(dvec_z2)))

    #Analytical integral (1/r*int) assuming specific form PML_Wr=d0(r^2/L^2)-amax(r/L-1)
    #0,dr,2dr,3dr... from r0
    r0=(nr-LPML_r)*dr #PML just starts at this R (Wr=0): location of tii
    dvec_r_global=dvec_r+r0*ones(length(dvec_r)) #r0+0, r0+dr, r0+2dr, ...
    PML_IWr=(
             d0_r/3*1/LPML_r_meter^2*(dvec_r).^3-
             amax*(1/2*1/LPML_r_meter*(dvec_r).^2-dvec_r)
                )./dvec_r_global

    #dr/2,dr/2+dr,dr/2+2dr... from r0
    #  PML_IWr2=d0*(1.0/2.0*a*dvec_r2/LPML_r_meter+1.0/3.0*b*dvec_r2.^2/LPML_r_meter^2)
    dvec_r2_global=dvec_r2+r0*ones(length(dvec_r))
    PML_IWr2=(
             d0_r/3*1/LPML_r_meter^2*(dvec_r2).^3-
             amax*(1/2*1/LPML_r_meter*(dvec_r2).^2-dvec_r2)
                )./dvec_r2_global

  #


 return  LPML_r,LPML_z,PML_Wr,PML_Wz,PML_IWr,PML_Wr2,PML_Wz2,PML_IWr2
end


#------------------------------
# initializing memory variables
#------------------------------
function init_memory_variables_Por(LPML_r,LPML_z,nr,nz)
    #--memory variables to calculate PML
    #for 2nd order FD, storing +2 sample layer larger than PML layer
    #m=0 only
    #Top and Bottom
    memT_vr=zeros(LPML_z+2,nr)
    memT_vz=zeros(LPML_z+2,nr)
    memT_trr=zeros(LPML_z+2,nr)
    memT_tpp=zeros(LPML_z+2,nr)
    memT_tzz=zeros(LPML_z+2,nr)
    memT_trz=zeros(LPML_z+2,nr)
    memT_vfr=zeros(LPML_z+2,nr)
    memT_vfz=zeros(LPML_z+2,nr)
    memT_pf=zeros(LPML_z+2,nr)

    memB_vr=zeros(LPML_z+2,nr)
    memB_vz=zeros(LPML_z+2,nr)
    memB_trr=zeros(LPML_z+2,nr)
    memB_tpp=zeros(LPML_z+2,nr)
    memB_tzz=zeros(LPML_z+2,nr)
    memB_trz=zeros(LPML_z+2,nr)
    memB_vfr=zeros(LPML_z+2,nr)
    memB_vfz=zeros(LPML_z+2,nr)
    memB_pf=zeros(LPML_z+2,nr)

    #Right
    memR_vr=zeros(nz,LPML_r+2)
    memR_vz=zeros(nz,LPML_r+2)
    memR_trr=zeros(nz,LPML_r+2)
    memR_tpp=zeros(nz,LPML_r+2)
    memR_tzz=zeros(nz,LPML_r+2)
    memR_trz=zeros(nz,LPML_r+2)
    memR_vfr=zeros(nz,LPML_r+2)
    memR_vfz=zeros(nz,LPML_r+2)
    memR_pf=zeros(nz,LPML_r+2)

return    memT_vr,memT_vz,memT_trr,memT_tpp,memT_tzz,memT_trz,memT_vfr,memT_vfz,memT_pf,
    memB_vr,memB_vz,memB_trr,memB_tpp,memB_tzz,memB_trz,memB_vfr,memB_vfz,memB_pf,
    memR_vr,memR_vz,memR_trr,memR_tpp,memR_tzz,memR_trz,memR_vfr,memR_vfz,memR_pf

end

#-----------------------------------
# Functions to save memory variables
#-----------------------------------
function save_memory_Top!(save_var,input_var,nr,nz,LPML_z,LPML_r)
    # storing field variables at top PML layer [index 1 = Top boundary]
    #for 2nd order FD, storing +2 sample layer larger than PML layer
@inbounds Threads.@threads  for ir=1:nr
        for iz=1:LPML_z+2
            save_var[iz,ir]=input_var[iz,ir]
        end
    end
end

function save_memory_Bottom!(save_var,input_var,nr,nz,LPML_z,LPML_r)
    # storing field variables at bottom PML layer [index 1 = PML-2 ]
    #for 2nd order FD, storing +2 sample layer larger than PML layer
@inbounds Threads.@threads for ir=1:nr
        for iz=1:LPML_z+2
            save_var[iz,ir]=input_var[nz-LPML_z-2+iz,ir]
        end
    end
end

function save_memory_Right!(save_var,input_var,nr,nz,LPML_z,LPML_r)
    # storing field variables at right PML layer [index 1 = PML-2]
    #for 2nd order FD, storing +2 sample layer larger than PML layer
@inbounds Threads.@threads for ir=1:LPML_r+2
        for iz=1:nz
            save_var[iz,ir]=input_var[iz,nr-LPML_r-2+ir]
        end
    end
end

function PML_save_vel_Top_Por!(memT_vr,memT_vz,memT_vfr,memT_vfz,vr,vz,vfr,vfz,nr,nz,LPML_z,LPML_r)
    save_memory_Top!(memT_vr,vr,nr,nz,LPML_z,LPML_r)
    save_memory_Top!(memT_vz,vz,nr,nz,LPML_z,LPML_r)
    save_memory_Top!(memT_vfr,vfr,nr,nz,LPML_z,LPML_r)
    save_memory_Top!(memT_vfz,vfz,nr,nz,LPML_z,LPML_r)
end

function PML_save_vel_Bottom_Por!(memB_vr,memB_vz,memB_vfr,memB_vfz,vr,vz,vfr,vfz,nr,nz,LPML_z,LPML_r)
    save_memory_Bottom!(memB_vr,vr,nr,nz,LPML_z,LPML_r)
    save_memory_Bottom!(memB_vz,vz,nr,nz,LPML_z,LPML_r)
    save_memory_Bottom!(memB_vfr,vfr,nr,nz,LPML_z,LPML_r)
    save_memory_Bottom!(memB_vfz,vfz,nr,nz,LPML_z,LPML_r)
end

function PML_save_vel_Right_Por!(memR_vr,memR_vz,memR_vfr,memR_vfz,vr,vz,vfr,vfz,nr,nz,LPML_z,LPML_r)
    save_memory_Right!(memR_vr,vr,nr,nz,LPML_z,LPML_r)
    save_memory_Right!(memR_vz,vz,nr,nz,LPML_z,LPML_r)
    save_memory_Right!(memR_vfr,vfr,nr,nz,LPML_z,LPML_r)
    save_memory_Right!(memR_vfz,vfz,nr,nz,LPML_z,LPML_r)
end

function PML_save_stress_Top_Por!(memT_trr,memT_tpp,memT_tzz,memT_trz,memT_pf,
    trr,tpp,tzz,trz,pf,nr,nz,LPML_z,LPML_r)
    save_memory_Top!(memT_trr,trr,nr,nz,LPML_z,LPML_r)
    save_memory_Top!(memT_tpp,tpp,nr,nz,LPML_z,LPML_r)
    save_memory_Top!(memT_tzz,tzz,nr,nz,LPML_z,LPML_r)
    save_memory_Top!(memT_trz,trz,nr,nz,LPML_z,LPML_r)
    save_memory_Top!(memT_pf,pf,nr,nz,LPML_z,LPML_r)
end

function PML_save_stress_Bottom_Por!(memB_trr,memB_tpp,memB_tzz,memB_trz,memB_pf,
    trr,tpp,tzz,trz,pf,nr,nz,LPML_z,LPML_r)
    save_memory_Bottom!(memB_trr,trr,nr,nz,LPML_z,LPML_r)
    save_memory_Bottom!(memB_tpp,tpp,nr,nz,LPML_z,LPML_r)
    save_memory_Bottom!(memB_tzz,tzz,nr,nz,LPML_z,LPML_r)
    save_memory_Bottom!(memB_trz,trz,nr,nz,LPML_z,LPML_r)
    save_memory_Bottom!(memB_pf,pf,nr,nz,LPML_z,LPML_r)
end

function PML_save_stress_Right_Por!(memR_trr,memR_tpp,memR_tzz,memR_trz,memR_pf,
    trr,tpp,tzz,trz,pf,nr,nz,LPML_z,LPML_r)
    save_memory_Right!(memR_trr,trr,nr,nz,LPML_z,LPML_r)
    save_memory_Right!(memR_tpp,tpp,nr,nz,LPML_z,LPML_r)
    save_memory_Right!(memR_tzz,tzz,nr,nz,LPML_z,LPML_r)
    save_memory_Right!(memR_trz,trz,nr,nz,LPML_z,LPML_r)
    save_memory_Right!(memR_pf,pf,nr,nz,LPML_z,LPML_r)
end

#-----------------------------------
# Updating velocity in PML
#-----------------------------------
function PML_update_velocity_1st_Top_Por!(vr,vz,trr,tpp,tzz,trz,vfr,vfz,pf,
    rhomat,rhofmat,D1mat,D2mat,Flagmat_vf_zero,nr,nz,dr,dz,dt,Prz_T,Pzz_T,PzzPE_T,LPML_z,LPML_r)
    #PML_variables (PQRS) index starts from first layer in PML (maybe reversed order)
    #z streching only
    #1st order
    #new in 06tmp5: k=2, vr,vphi->1st, vz->2nd
@inbounds Threads.@threads for j=2:nr-LPML_r #assumig LPML_r>2
          for k=1:LPML_z

              #I am  accesing [k,j], but definition of r_now changes with each component!
             if(k!=1)
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

                 vfr_old=vfr_now
                 if(Flagmat_vf_zero[k,j]==1)
                     vfr[k,j]=0.0
                 else
                     vr_dummy,vfr[k,j]=update_vr_1st_Por(vr_now,vfr_now,
                            trr_f1,trr_b1,trz_f1,trz_b1,trr_av,tpp_av,
                            pf_f1,pf_b1,
                            dz,dr,dt,r_now,rho_av,rhof_av,D1_av,D2_av,Flagmat_vf_zero[k,j])

                     #PML addition (vfr)
                     vfr[k,j]=vfr[k,j]+
                              (D2_av/2+(D1_av-rhof_av^2/rho_av)/dt)^(-1)*(
                                                -rhof_av/rho_av*Prz_T[LPML_z-k+1,j]
                                                )
                 end

                 #PML (vr) calculate elastic one then vfr addition
                 vr[k,j]=update_vr(vr_now,trr_f1,trr_b1,trz_f1,trz_b1,trr_av,tpp_av,
                        0,dz,dr,dt,r_now,rho_av,0)
                 vr[k,j]=vr[k,j]+dt/rho_av*Prz_T[LPML_z-k+1,j] #elastic+PML
                 dert_vfr=(vfr[k,j]-vfr_old)/dt
                 vr[k,j]=vr[k,j]-dt*rhof_av/rho_av*dert_vfr #poroelastic PML

            end #if k==1

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

             vfz_old=vfz_now
             if(Flagmat_vf_zero[k,j]==1)
                 vfz[k,j]=0.0
             else
                 vz_dummy,vfz[k,j]=update_vz_1st_Por(vz_now,vfz_now,
                     trz_f1,trz_b1,tzz_f1,tzz_b1,trz_av,
                     pf_f1,pf_b1,
                     dz,dr,dt,r_now,rho_av,rhof_av,D1_av,D2_av,Flagmat_vf_zero[k,j])

                 #PML addition (vfz)
                 vfz[k,j]=vfz[k,j]+
                          (D2_av/2+(D1_av-rhof_av^2/rho_av)/dt)^(-1)*(
                                            -PzzPE_T[LPML_z-k+1,j]
                                            -rhof_av/rho_av*Pzz_T[LPML_z-k+1,j]
                                            )
             end

             #PML (vz) calculate elastic one then vfz addition
             vz[k,j]=update_vz(vz_now,trz_f1,trz_b1,tzz_f1,tzz_b1,trz_av,
                 0,dz,dr,dt,r_now,rho_av,0)
             vz[k,j]=vz[k,j]+dt/rho_av*Pzz_T[LPML_z-k+1,j] #elastic+PML
             dert_vfz=(vfz[k,j]-vfz_old)/dt
             vz[k,j]=vz[k,j]-dt*rhof_av/rho_av*dert_vfz #poroelatic PML

             #tmptmptmp
#             if(j==5)
#             if(k==LPML_z)
#                 println("TopPML vz(LPML_z,5)=",vz[k,j])
#             end
#             end


         end #k (z)
    end #j (r)

end

#Poroelastic
function PML_update_velocity_1st_Bottom_Por!(vr,vz,trr,tpp,tzz,trz,vfr,vfz,pf,
    rhomat,rhofmat,D1mat,D2mat,Flagmat_vf_zero,nr,nz,dr,dz,dt,Prz_B,Pzz_B,PzzPE_B,LPML_z,LPML_r)
    #PML_variables (PQRS) index starts from first layer in PML (maybe reversed order)
    #z streching only
    #1st order
    #new in 06tmp5:
    #    : 2nd (in z)->vphi,vr(k=>nz-1), vz(k=>nz-2)
    #    : 1st (in z)->vphi,vr(k=nz), vz(k=nz-1)
    #    : vz(k=nz) is blank (Dirichlet BC)

@inbounds Threads.@threads for j=2:nr-LPML_r #assumig LPML_r>2
#          for k=3:LPML_z
          for k=nz-LPML_z+1:nz
              #I am  accesing [k,j], but the definition of r_now changes with each component!
             r_now=(j-1)*dr+dr/2.0 #for vr (Mittet)

             vr_now=vr[k,j]
             trr_f1,trr_b1=trr[k,j+1],trr[k,j]
             trz_f1,trz_b1=trz[k,j],trz[k-1,j]
             rho_av=(rhomat[k,j]+rhomat[k,j+1])/2.0 #be carefull

             tpp_av=0.5*(tpp[k,j+1]+tpp[k,j]) #see trr
             trr_av=0.5*(trr_f1+trr_b1)

             vfr_now=vfr[k,j]
             pf_f1,pf_b1=pf[k,j+1],pf[k,j]
             rhof_av=(rhofmat[k,j]+rhofmat[k,j+1])/2.0 #be carefull
             D1_av=(D1mat[k,j]+D1mat[k,j+1])/2.0 #be carefull
             D2_av=(D2mat[k,j]+D2mat[k,j+1])/2.0 #be carefull

             vfr_old=vfr_now
             if(Flagmat_vf_zero[k,j]==1)
                 vfr[k,j]=0.0
             else
                 vr_dummy,vfr[k,j]=update_vr_1st_Por(vr_now,vfr_now,
                        trr_f1,trr_b1,trz_f1,trz_b1,trr_av,tpp_av,
                        pf_f1,pf_b1,
                        dz,dr,dt,r_now,rho_av,rhof_av,D1_av,D2_av,Flagmat_vf_zero[k,j])

                 #PML addition (vfr)
                 vfr[k,j]=vfr[k,j]+
                          (D2_av/2+(D1_av-rhof_av^2/rho_av)/dt)^(-1)*(
                                            -rhof_av/rho_av*Prz_B[k-nz+LPML_z,j]
                                            )
             end

             #PML (vr) calculate elastic one then vfr addition
             vr[k,j]=update_vr(vr_now,trr_f1,trr_b1,trz_f1,trz_b1,trr_av,tpp_av,
                    0,dz,dr,dt,r_now,rho_av,0)
             vr[k,j]=vr[k,j]+dt/rho_av*Prz_B[k-nz+LPML_z,j] #elastic PML
             dert_vfr=(vfr[k,j]-vfr_old)/dt
             vr[k,j]=vr[k,j]-dt*rhof_av/rho_av*dert_vfr #poroelastic PML


             if(k!=nz)
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

                 vfz_old=vfz_now
                 if(Flagmat_vf_zero[k,j]==1)
                     vfz[k,j]=0.0
                 else
                     vz_dummy,vfz[k,j]=update_vz_1st_Por(vz_now,vfz_now,
                         trz_f1,trz_b1,tzz_f1,tzz_b1,trz_av,
                         pf_f1,pf_b1,
                         dz,dr,dt,r_now,rho_av,rhof_av,D1_av,D2_av,Flagmat_vf_zero[k,j])

                     #PML addition (vfz)
                     vfz[k,j]=vfz[k,j]+
                              (D2_av/2+(D1_av-rhof_av^2/rho_av)/dt)^(-1)*(
                                                -PzzPE_B[k-nz+LPML_z,j]
                                                -rhof_av/rho_av*Pzz_B[k-nz+LPML_z,j]
                                                )

                 end

                 #PML (vz) calculate elastic one then vfz addition
                 vz[k,j]=update_vz(vz_now,trz_f1,trz_b1,tzz_f1,tzz_b1,trz_av,
                     0,dz,dr,dt,r_now,rho_av,0)
                 vz[k,j]=vz[k,j]+dt/rho_av*Pzz_B[k-nz+LPML_z,j] #elastic PML
                 dert_vfz=(vfz[k,j]-vfz_old)/dt
                 vz[k,j]=vz[k,j]-dt*rhof_av/rho_av*dert_vfz #poroelastic PML
             end #k==nz

         end #k (z)
    end #j (r)

end

function PML_update_velocity_1st_Right_Por!(vr,vz,trr,tpp,tzz,trz,vfr,vfz,pf,
    rhomat,rhofmat,D1mat,D2mat,Flagmat_vf_zero,nr,nz,dr,dz,dt,
    Prr_R,Qrp_R,Pzr_R,Qzp_R,Rr_R,Rp_R,Rrz_R,
    PrrPE_R,RrPE_R,RpPE_R,
    LPML_z,LPML_r)
    #PML_variables (PQRS) index starts from first layer in PML (maybe reversed order)
    #r streching only

    #PE-to-E horizontal boundary incorpolated by checking if Flagmat_vf_zero[k,j]==0 && D1[k+1,j]==Inf

    #j=nr-1->vphi,vz (2nd), vr (1st)
    #j=nr->vphi,vz (1st), vr (zero)

    #1st order
@inbounds Threads.@threads  for j=nr-LPML_r+1:nr
          for k=LPML_z+1:nz-LPML_z #assuming LPML_r>2
              #I am alywas accesing [k,j], but definition of r_now changes with components!
             if(j!=nr)

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

                 vfr_old=vfr_now
                 if(Flagmat_vf_zero[k,j]==1)
                     vfr[k,j]=0.0
                 else
                     vr_dummy,vfr[k,j]=update_vr_1st_Por(vr_now,vfr_now,
                         trr_f1,trr_b1,trz_f1,trz_b1,trr_av,tpp_av,
                         pf_f1,pf_b1,
                         dz,dr,dt,r_now,rho_av,rhof_av,D1_av,D2_av,Flagmat_vf_zero[k,j])

                     #PML addition (vfr)
                     vfr[k,j]=vfr[k,j]+
                              (D2_av/2+(D1_av-rhof_av^2/rho_av)/dt)^(-1)*(
                                                -PrrPE_R[k-LPML_z,j-nr+LPML_r]
                                                -rhof_av/rho_av*(
                                                    Prr_R[k-LPML_z,j-nr+LPML_r]+1.0/r_now*Qrp_R[k-LPML_z,j-nr+LPML_r]
                                                    )
                                                )
                 end

                 #PML (vr) calculate elastic one then vfr addition
                 vr[k,j]=update_vr(vr_now,trr_f1,trr_b1,trz_f1,trz_b1,trr_av,tpp_av,
                        0,dz,dr,dt,r_now,rho_av,0)
                 vr[k,j]=vr[k,j]+dt/rho_av*(
                      Prr_R[k-LPML_z,j-nr+LPML_r]+1.0/r_now*Qrp_R[k-LPML_z,j-nr+LPML_r]) #elastic PML
                 dert_vfr=(vfr[k,j]-vfr_old)/dt
                 vr[k,j]=vr[k,j]-dt*rhof_av/rho_av*dert_vfr #poroelastic PML
             end #if j==nr

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

             vfz_old=vfz_now
             #implementing PE to E horizontal boundary
             if(Flagmat_vf_zero[k,j]==1 || (Flagmat_vf_zero[k,j]==0 && D1mat[k+1,j]==Inf) )
#             if(Flagmat_vf_zero[k,j]==1)
                 vfz[k,j]=0.0
             else
                 vz_dummy,vfz[k,j]=update_vz_1st_Por(vz_now,vfz_now,
                     trz_f1,trz_b1,tzz_f1,tzz_b1,trz_av,
                     pf_f1,pf_b1,
                     dz,dr,dt,r_now,rho_av,rhof_av,D1_av,D2_av,Flagmat_vf_zero[k,j])

                 #PML addition (vfz)
                 vfz[k,j]=vfz[k,j]+
                          (D2_av/2+(D1_av-rhof_av^2/rho_av)/dt)^(-1)*(
                              -rhof_av/rho_av*(
                                  Pzr_R[k-LPML_z,j-nr+LPML_r]+1.0/r_now*Qzp_R[k-LPML_z,j-nr+LPML_r]
                                  )
                            )

             end
             #PML (vz) calculate elastic one then vfz addition
             vz[k,j]=update_vz(vz_now,trz_f1,trz_b1,tzz_f1,tzz_b1,trz_av,
                 0,dz,dr,dt,r_now,rho_av,0)
             vz[k,j]=vz[k,j]+dt/rho_av*(
                Pzr_R[k-LPML_z,j-nr+LPML_r]+1.0/r_now*Qzp_R[k-LPML_z,j-nr+LPML_r]) #elastic PML
             dert_vfz=(vfz[k,j]-vfz_old)/dt
             vz[k,j]=vz[k,j]-dt*rhof_av/rho_av*dert_vfz #poroelastic PML

         end #k (z)
    end #j (r)

end

function PML_update_velocity_1st_TopRight_Por!(vr,vz,trr,tpp,tzz,trz,vfr,vfz,pf,
    rhomat,rhofmat,D1mat,D2mat,Flagmat_vf_zero,nr,nz,dr,dz,dt,
    Prz_TR,Pzz_TR,PzzPE_TR,
    Prr_TR,Qrp_TR,Pzr_TR,Qzp_TR,Rr_TR,Rp_TR,Rrz_TR,
    PrrPE_TR,RrPE_TR,RpPE_TR,
    LPML_z,LPML_r)
    #PML_variables (PQRS) index starts from first layer in PML (maybe reversed order)
    #z streching only
    #1st order
@inbounds Threads.@threads for j=nr-LPML_r+1:nr #assumig LPML_r>2
          for k=1:LPML_z
              #I am accesing [k,j], but the definition of r_now changes with each component!
             if(k!=1 && j!=nr)
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

                 vfr_old=vfr_now
                if(Flagmat_vf_zero[k,j]==1)
                    vfr[k,j]=0.0
                else
                     vr_dummy,vfr[k,j]=update_vr_1st_Por(vr_now,vfr_now,
                            trr_f1,trr_b1,trz_f1,trz_b1,trr_av,tpp_av,
                            pf_f1,pf_b1,
                            dz,dr,dt,r_now,rho_av,rhof_av,D1_av,D2_av,Flagmat_vf_zero[k,j])

                     #PML addition (vfr)
                     vfr[k,j]=vfr[k,j]+
                              (D2_av/2+(D1_av-rhof_av^2/rho_av)/dt)^(-1)*(
                                                -PrrPE_TR[LPML_z-k+1,j-nr+LPML_r]
                                                -rhof_av/rho_av*(
                                                    Prr_TR[LPML_z-k+1,j-nr+LPML_r]
                                                    +1.0/r_now*Qrp_TR[LPML_z-k+1,j-nr+LPML_r]
                                                    +Prz_TR[LPML_z-k+1,j-nr+LPML_r]
                                                    )
                                                )
                 end

                 #PML (vr) calculate elastic one then vfr addition
                 vr[k,j]=update_vr(vr_now,trr_f1,trr_b1,trz_f1,trz_b1,trr_av,tpp_av,
                        0,dz,dr,dt,r_now,rho_av,0)
                 vr[k,j]=vr[k,j]+dt/rho_av*(
                    Prz_TR[LPML_z-k+1,j-nr+LPML_r]+
                    Prr_TR[LPML_z-k+1,j-nr+LPML_r]+
                    1.0/r_now*Qrp_TR[LPML_z-k+1,j-nr+LPML_r]) #elastic PML
                 dert_vfr=(vfr[k,j]-vfr_old)/dt
                 vr[k,j]=vr[k,j]-dt*rhof_av/rho_av*dert_vfr #poroelastic PML
             end #if k==1, j==nr


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

             vfz_old=vfz_now
            if(Flagmat_vf_zero[k,j]==1)
               vfz[k,j]=0.0
            else
                 vz_dummy,vfz[k,j]=update_vz_1st_Por(vz_now,vfz_now,
                     trz_f1,trz_b1,tzz_f1,tzz_b1,trz_av,
                     pf_f1,pf_b1,
                     dz,dr,dt,r_now,rho_av,rhof_av,D1_av,D2_av,Flagmat_vf_zero[k,j])

                 #PML addition (vfz)
                 vfz[k,j]=vfz[k,j]+
                          (D2_av/2+(D1_av-rhof_av^2/rho_av)/dt)^(-1)*(
                              -PzzPE_TR[LPML_z-k+1,j-nr+LPML_r]
                              -rhof_av/rho_av*(
                                  +Pzr_TR[LPML_z-k+1,j-nr+LPML_r]
                                  +Pzz_TR[LPML_z-k+1,j-nr+LPML_r]
                                  +1.0/r_now*Qzp_TR[LPML_z-k+1,j-nr+LPML_r]
                                  )
                            )
             end
             #PML (vz) calculate elastic one then vfz addition
             vz[k,j]=update_vz(vz_now,trz_f1,trz_b1,tzz_f1,tzz_b1,trz_av,
                 0,dz,dr,dt,r_now,rho_av,0)

             vz[k,j]=vz[k,j]+dt/rho_av*(
                Pzz_TR[LPML_z-k+1,j-nr+LPML_r]+
                Pzr_TR[LPML_z-k+1,j-nr+LPML_r]+
                1.0/r_now*Qzp_TR[LPML_z-k+1,j-nr+LPML_r]) #elastic PML
             dert_vfz=(vfz[k,j]-vfz_old)/dt
             vz[k,j]=vz[k,j]-dt*rhof_av/rho_av*dert_vfz #poroelastic PML

         end #k (z)
    end #j (r)

end


function PML_update_velocity_1st_BottomRight_Por!(vr,vz,trr,tpp,tzz,trz,vfr,vfz,pf,
    rhomat,rhofmat,D1mat,D2mat,Flagmat_vf_zero,nr,nz,dr,dz,dt,
    Prz_BR,Pzz_BR,PzzPE_BR,
    Prr_BR,Qrp_BR,Pzr_BR,Qzp_BR,Rr_BR,Rp_BR,Rrz_BR,
    PrrPE_BR,RrPE_BR,RpPE_BR,
    LPML_z,LPML_r)
    #PML_variables (PQRS) index starts from first layer in PML (maybe reversed order)
    #z streching only
    #1st order
@inbounds Threads.@threads for j=nr-LPML_r+1:nr #assumig LPML_r>2
    for k=nz-LPML_z+1:nz #this points normal field index
              #I am accesing [k,j], but the definition of r_now changes with each component!
             if(j!=nr)
                 r_now=(j-1)*dr+dr/2.0 #for vr (Mittet)

                 vr_now=vr[k,j]
                 trr_f1,trr_b1=trr[k,j+1],trr[k,j]
                 trz_f1,trz_b1=trz[k,j],trz[k-1,j]
                 rho_av=(rhomat[k,j]+rhomat[k,j+1])/2.0 #be carefull

                 tpp_av=0.5*(tpp[k,j+1]+tpp[k,j]) #see trr
                 trr_av=0.5*(trr_f1+trr_b1)

                 vfr_now=vfr[k,j]
                 pf_f1,pf_b1=pf[k,j+1],pf[k,j]
                 rhof_av=(rhofmat[k,j]+rhofmat[k,j+1])/2.0 #be carefull
                 D1_av=(D1mat[k,j]+D1mat[k,j+1])/2.0 #be carefull
                 D2_av=(D2mat[k,j]+D2mat[k,j+1])/2.0 #be carefull

                 vfr_old=vfr_now
                 if(Flagmat_vf_zero[k,j]==1)
                     vfr[k,j]=0.0
                 else
                     vr_dummy,vfr[k,j]=update_vr_1st_Por(vr_now,vfr_now,
                            trr_f1,trr_b1,trz_f1,trz_b1,trr_av,tpp_av,
                            pf_f1,pf_b1,
                            dz,dr,dt,r_now,rho_av,rhof_av,D1_av,D2_av,Flagmat_vf_zero[k,j])

                     #PML addition (vfr)
                     vfr[k,j]=vfr[k,j]+
                              (D2_av/2+(D1_av-rhof_av^2/rho_av)/dt)^(-1)*(
                                                -PrrPE_BR[k-nz+LPML_z,j-nr+LPML_r]
                                                -rhof_av/rho_av*(
                                                    Prr_BR[k-nz+LPML_z,j-nr+LPML_r]
                                                    +1.0/r_now*Qrp_BR[k-nz+LPML_z,j-nr+LPML_r]
                                                    +Prz_BR[k-nz+LPML_z,j-nr+LPML_r]
                                                    )
                                                )

                 end
                 #PML (vr) calculate elastic one then vfr addition
                 vr[k,j]=update_vr(vr_now,trr_f1,trr_b1,trz_f1,trz_b1,trr_av,tpp_av,
                        0,dz,dr,dt,r_now,rho_av,0)
                 vr[k,j]=vr[k,j]+dt/rho_av*(
                    Prz_BR[k-nz+LPML_z,j-nr+LPML_r]+
                    Prr_BR[k-nz+LPML_z,j-nr+LPML_r]+
                    1.0/r_now*Qrp_BR[k-nz+LPML_z,j-nr+LPML_r]) #elastic PML
                 dert_vfr=(vfr[k,j]-vfr_old)/dt
                 vr[k,j]=vr[k,j]-dt*rhof_av/rho_av*dert_vfr #poroelastic PML
             end #j==nr

             if(k!=nz)
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

                 vfz_old=vfz_now
                 if(Flagmat_vf_zero[k,j]==1)
                     vfz[k,j]=0.0
                 else
                     vz_dummy,vfz[k,j]=update_vz_1st_Por(vz_now,vfz_now,
                         trz_f1,trz_b1,tzz_f1,tzz_b1,trz_av,
                         pf_f1,pf_b1,
                         dz,dr,dt,r_now,rho_av,rhof_av,D1_av,D2_av,Flagmat_vf_zero[k,j])

                     #PML addition (vfz)
                     vfz[k,j]=vfz[k,j]+
                              (D2_av/2+(D1_av-rhof_av^2/rho_av)/dt)^(-1)*(
                                  -PzzPE_BR[k-nz+LPML_z,j-nr+LPML_r]
                                  -rhof_av/rho_av*(
                                      +Pzr_BR[k-nz+LPML_z,j-nr+LPML_r]
                                      +Pzz_BR[k-nz+LPML_z,j-nr+LPML_r]
                                      +1.0/r_now*Qzp_BR[k-nz+LPML_z,j-nr+LPML_r]
                                      )
                                )

                 end

                 #PML (vz) calculate elastic one then vfz addition
                 vz[k,j]=update_vz(vz_now,trz_f1,trz_b1,tzz_f1,tzz_b1,trz_av,
                     0,dz,dr,dt,r_now,rho_av,0)
                 vz[k,j]=vz[k,j]+dt/rho_av*(
                    Pzz_BR[k-nz+LPML_z,j-nr+LPML_r]+
                    Pzr_BR[k-nz+LPML_z,j-nr+LPML_r]+
                    1.0/r_now*Qzp_BR[k-nz+LPML_z,j-nr+LPML_r]) #elastic PML
                 dert_vfz=(vfz[k,j]-vfz_old)/dt
                 vz[k,j]=vz[k,j]-dt*rhof_av/rho_av*dert_vfz  #poroelastic PML
             end #if k==nz

         end #k (z)
    end #j (r)

end

#-----------------------------------
# Updating stress in PML
#-----------------------------------
function PML_update_stress_1st_Top_Por!(vr,vz,trr,tpp,tzz,trz,vfr,vfz,pf,
    Mmat,Cmat,Hmat,Gmat,nr,nz,dr,dz,dt,Rz_T,Srz_T,RzPE_T,LPML_z,LPML_r)
    #PML_variables (PQRS) index starts from first layer in PML (maybe reversed order)
    #z streching only
    #1st order
    #new in 06tmp5
    #    : 2nd (in z)->tpz,trz(k=>2), tii,trp(k=>3)
    #    : 1st (in z)->tpz,trz(k=1), tii,trp(k=2)
    #    : tii,trp(k=1) is blank (Dirichlet BC)
    # modified rhomat->lmat,mmat
@inbounds Threads.@threads for j=2:nr-LPML_r
          for k=1:LPML_z
             #I am accesing [k,j], but the definition of r_now changes with each component!
             if(k!=1)
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

                 #PML addition
                 trr[k,j]=trr[k,j]
                          +dt*(H-2G)*Rz_T[LPML_z-k+1,j]
                          +dt*C*RzPE_T[LPML_z-k+1,j]


                 tpp_now=tpp[k,j]
                 tpp[k,j]=update_tpp_1st_Por(tpp_now,vr_f1,vr_b1,vz_f1,vz_b1,vr_av,
                        vfr_f1,vfr_b1,vfz_f1,vfz_b1,vfr_av,
                        dz,dr,dt,r_now,C,H,G)

                 #PML addition
                 tpp[k,j]=tpp[k,j]
                          +dt*(H-2G)*Rz_T[LPML_z-k+1,j]
                          +dt*C*RzPE_T[LPML_z-k+1,j]


                 tzz_now=tzz[k,j]
                 tzz[k,j]=update_tzz_1st_Por(tzz_now,vr_f1,vr_b1,vz_f1,vz_b1,vr_av,
                        vfr_f1,vfr_b1,vfz_f1,vfz_b1,vfr_av,
                        dz,dr,dt,r_now,C,H,G)

                 #PML addition
                 tzz[k,j]=tzz[k,j]
                          +dt*H*Rz_T[LPML_z-k+1,j]
                          +dt*C*RzPE_T[LPML_z-k+1,j]


                 pf_now=pf[k,j]
                 pf[k,j]=update_pf_1st_Por(pf_now,vr_f1,vr_b1,vz_f1,vz_b1,vr_av,
                        vfr_f1,vfr_b1,vfz_f1,vfz_b1,vfr_av,
                        dz,dr,dt,r_now,C,M)

                 #PML addition
                 pf[k,j]=pf[k,j]
                         -dt*C*Rz_T[LPML_z-k+1,j]
                         -dt*M*RzPE_T[LPML_z-k+1,j]
             end #if k==1


             r_now=(j-1)*dr+dr/2.0 #for trp,trz Mittet

             trz_now=trz[k,j]
             if (Gmat[k,j]*Gmat[k,j+1]*Gmat[k+1,j]*Gmat[k+1,j+1]==0.0)
                  G_av=0.0
             else
                 gj1=2.0*Gmat[k,j]*Gmat[k+1,j]/(Gmat[k,j]+Gmat[k+1,j])
                 gj2=2.0*Gmat[k,j+1]*Gmat[k+1,j+1]/(Gmat[k,j+1]+Gmat[k+1,j+1])
                 G_av=2.0*gj1*gj2/(gj1+gj2)
             end

             vr_f1,vr_b1=vr[k+1,j],vr[k,j]
             vz_f1,vz_b1=vz[k,j+1],vz[k,j]
             trz[k,j]=update_trz(trz_now,vr_f1,vr_b1,vz_f1,vz_b1,
                 dz,dr,dt,r_now,G_av,0)

             #PML addition
             trz[k,j]=trz[k,j]+dt*G_av*Srz_T[LPML_z-k+1,j]
             #experimental: when k==LPML_z, trz is outside PML

       end #k (z)
    end #j (r)

end


function PML_update_stress_1st_Bottom_Por!(vr,vz,trr,tpp,tzz,trz,vfr,vfz,pf,
    Mmat,Cmat,Hmat,Gmat,nr,nz,dr,dz,dt,Rz_B,Srz_B,RzPE_B,LPML_z,LPML_r)
    #PML_variables (PQRS) index starts from first layer in PML (maybe reversed order)
    #z streching only
    #1st order
    #new in 06tmp5
    #    : 2nd (in z)->tii,trp(k=>nz-1), tpz,trz(k=>nz-2)
    #    : 1st (in z)->tii,trp(k=nz), tpz,trz(k=nz-1)
    #    : tpz,trz(k=nz) is blank (Dirichlet BC)
    # modified rhomat->lmat,mmat

@inbounds Threads.@threads for j=2:nr-LPML_r
          for k=nz-LPML_z+1:nz
             #I am accesing [k,j], but the definition of r_now changes with each component!
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


             #PML addition
             trr[k,j]=trr[k,j]
                     +dt*(H-2G)*Rz_B[k-nz+LPML_z,j]
                     +dt*C*RzPE_B[k-nz+LPML_z,j]


             tpp_now=tpp[k,j]
             tpp[k,j]=update_tpp_1st_Por(tpp_now,vr_f1,vr_b1,vz_f1,vz_b1,vr_av,
                    vfr_f1,vfr_b1,vfz_f1,vfz_b1,vfr_av,
                    dz,dr,dt,r_now,C,H,G)

             #PML addition
             tpp[k,j]=tpp[k,j]
                      +dt*(H-2G)*Rz_B[k-nz+LPML_z,j]
                      +dt*C*RzPE_B[k-nz+LPML_z,j]


             tzz_now=tzz[k,j]
             tzz[k,j]=update_tzz_1st_Por(tzz_now,vr_f1,vr_b1,vz_f1,vz_b1,vr_av,
                    vfr_f1,vfr_b1,vfz_f1,vfz_b1,vfr_av,
                    dz,dr,dt,r_now,C,H,G)

             #PML addition
             tzz[k,j]=tzz[k,j]
                       +dt*H*Rz_B[k-nz+LPML_z,j]
                       +dt*C*RzPE_B[k-nz+LPML_z,j]


             pf_now=pf[k,j]
             pf[k,j]=update_pf_1st_Por(pf_now,vr_f1,vr_b1,vz_f1,vz_b1,vr_av,
                    vfr_f1,vfr_b1,vfz_f1,vfz_b1,vfr_av,
                    dz,dr,dt,r_now,C,M)

             #PML addition
             pf[k,j]=pf[k,j]
                   -dt*C*Rz_B[k-nz+LPML_z,j]
                   -dt*M*RzPE_B[k-nz+LPML_z,j]



             if(k!=nz)
                 r_now=(j-1)*dr+dr/2.0 #for trp,trz Mittet

                 trz_now=trz[k,j]
                 if (Gmat[k,j]*Gmat[k,j+1]*Gmat[k+1,j]*Gmat[k+1,j+1]==0.0)
                     G_av=0.0
                 else
                     gj1=2.0*Gmat[k,j]*Gmat[k+1,j]/(Gmat[k,j]+Gmat[k+1,j])
                     gj2=2.0*Gmat[k,j+1]*Gmat[k+1,j+1]/(Gmat[k,j+1]+Gmat[k+1,j+1])
                     G_av=2.0*gj1*gj2/(gj1+gj2)
                 end

                 vr_f1,vr_b1=vr[k+1,j],vr[k,j]
                 vz_f1,vz_b1=vz[k,j+1],vz[k,j]
                 trz[k,j]=update_trz(trz_now,vr_f1,vr_b1,vz_f1,vz_b1,
                     dz,dr,dt,r_now,G_av,0)
                 #PML addition
                 trz[k,j]=trz[k,j]+dt*G_av*Srz_B[k-nz+LPML_z,j]
             end #k==nz

       end #k (z)
    end #j (r)
end

function PML_update_stress_1st_Right_Por!(vr,vz,trr,tpp,tzz,trz,vfr,vfz,pf,
    Mmat,Cmat,Hmat,Gmat,nr,nz,dr,dz,dt,
    Rr_R,Rp_R,Rrz_R,RrPE_R,RpPE_R,
    LPML_z,LPML_r)
    #PML_variables (PQRS) index starts from first layer in PML (maybe reversed order)
    #r streching only
    #1st order
    # modified rhomat->lmat,mmat
    #j=nr-1->tii,tpz (2nd), trp,trz (1st)
    #j=nr->tii,tpz (1st), trp,trz (zero)

@inbounds Threads.@threads  for j=nr-LPML_r+1:nr
        for k=LPML_z+1:nz-LPML_z #assuming LPML_r>2
             #I am accesing [k,j], but the definition of r_now changes with each component!
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

             #PML addition
             trr[k,j]=trr[k,j]
                      +dt*(H-2G)*(1.0/r_now*Rp_R[k-LPML_z,j-nr+LPML_r])
                      +dt*H*Rr_R[k-LPML_z,j-nr+LPML_r]
                      +dt*C*(
                        RrPE_R[k-LPML_z,j-nr+LPML_r]+1.0/r_now*RpPE_R[k-LPML_z,j-nr+LPML_r])



             tpp_now=tpp[k,j]
             tpp[k,j]=update_tpp_1st_Por(tpp_now,vr_f1,vr_b1,vz_f1,vz_b1,vr_av,
                    vfr_f1,vfr_b1,vfz_f1,vfz_b1,vfr_av,
                    dz,dr,dt,r_now,C,H,G)

             #PML addition
             tpp[k,j]=tpp[k,j]
                    +dt*(H-2G)*(Rr_R[k-LPML_z,j-nr+LPML_r])
                    +dt*H*(1.0/r_now*Rp_R[k-LPML_z,j-nr+LPML_r])
                    +dt*C*(
                    RrPE_R[k-LPML_z,j-nr+LPML_r]+1.0/r_now*RpPE_R[k-LPML_z,j-nr+LPML_r])

             tzz_now=tzz[k,j]
             tzz[k,j]=update_tzz_1st_Por(tzz_now,vr_f1,vr_b1,vz_f1,vz_b1,vr_av,
                    vfr_f1,vfr_b1,vfz_f1,vfz_b1,vfr_av,
                    dz,dr,dt,r_now,C,H,G)

             #PML addition
             tzz[k,j]=tzz[k,j]
                     +dt*(H-2G)*(1.0/r_now*Rp_R[k-LPML_z,j-nr+LPML_r]+Rr_R[k-LPML_z,j-nr+LPML_r])
                     +dt*C*(
                     RrPE_R[k-LPML_z,j-nr+LPML_r]+1.0/r_now*RpPE_R[k-LPML_z,j-nr+LPML_r])


             pf_now=pf[k,j]
             pf[k,j]=update_pf_1st_Por(pf_now,vr_f1,vr_b1,vz_f1,vz_b1,vr_av,
                    vfr_f1,vfr_b1,vfz_f1,vfz_b1,vfr_av,
                    dz,dr,dt,r_now,C,M)

             #PML addition
             pf[k,j]=pf[k,j]
                     -dt*C*(Rr_R[k-LPML_z,j-nr+LPML_r]+1.0/r_now*Rp_R[k-LPML_z,j-nr+LPML_r])
                     -dt*M*(RrPE_R[k-LPML_z,j-nr+LPML_r]+1.0/r_now*RpPE_R[k-LPML_z,j-nr+LPML_r])



             if(j!=nr)
                 r_now=(j-1)*dr+dr/2.0 #for trp,trz Mittet

                 trz_now=trz[k,j]
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
                 #PML addition : caution k and j index (starts from 1)
                 trz[k,j]=trz[k,j]+dt*(
                     G_av*Rrz_R[k-LPML_z,j-nr+LPML_r])
             end #if j==nr

       end #k (z)
    end #j (r)

end


function PML_update_stress_1st_TopRight_Por!(vr,vz,trr,tpp,tzz,trz,vfr,vfz,pf,
    Mmat,Cmat,Hmat,Gmat,nr,nz,dr,dz,dt,
    Rr_TR,Rp_TR,Rz_TR,Rrz_TR,Srz_TR,RrPE_TR,RpPE_TR,RzPE_TR,
    LPML_z,LPML_r)
    #PML_variables (PQRS) index starts from first layer in PML (maybe reversed order)
    #z streching only
    #2nd order
    # modified rhomat->lmat,mmat

@inbounds Threads.@threads  for j=nr-LPML_r+1:nr #
          for k=1:LPML_z
             #I am accesing [k,j], but the definition of r_now changes with each component!
             if(k!=1)
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

                #PML addition
                trr[k,j]=trr[k,j]
                       +dt*(H-2G)*(1.0/r_now*Rp_TR[LPML_z-k+1,j-nr+LPML_r]+Rz_TR[LPML_z-k+1,j-nr+LPML_r])
                       +dt*H*Rr_TR[LPML_z-k+1,j-nr+LPML_r]
                       +dt*C*(
                         RrPE_TR[LPML_z-k+1,j-nr+LPML_r]
                         +1.0/r_now*RpPE_TR[LPML_z-k+1,j-nr+LPML_r]
                         +RzPE_TR[LPML_z-k+1,j-nr+LPML_r])


                 tpp_now=tpp[k,j]
                 tpp[k,j]=update_tpp_1st_Por(tpp_now,vr_f1,vr_b1,vz_f1,vz_b1,vr_av,
                        vfr_f1,vfr_b1,vfz_f1,vfz_b1,vfr_av,
                        dz,dr,dt,r_now,C,H,G)

                 #PML addition
                 tpp[k,j]=tpp[k,j]
                     +dt*(H-2G)*(Rz_TR[LPML_z-k+1,j-nr+LPML_r]+Rr_TR[LPML_z-k+1,j-nr+LPML_r])
                     +dt*H*(1.0/r_now*Rp_TR[LPML_z-k+1,j-nr+LPML_r])
                     +dt*C*(
                       RrPE_TR[LPML_z-k+1,j-nr+LPML_r]
                       +1.0/r_now*RpPE_TR[LPML_z-k+1,j-nr+LPML_r]
                       +RzPE_TR[LPML_z-k+1,j-nr+LPML_r])


                 tzz_now=tzz[k,j]
                 tzz[k,j]=update_tzz_1st_Por(tzz_now,vr_f1,vr_b1,vz_f1,vz_b1,vr_av,
                        vfr_f1,vfr_b1,vfz_f1,vfz_b1,vfr_av,
                        dz,dr,dt,r_now,C,H,G)


                 #PML addition
                 tzz[k,j]=tzz[k,j]
                       +dt*(H-2G)*(1.0/r_now*Rp_TR[LPML_z-k+1,j-nr+LPML_r]
                                    +Rr_TR[LPML_z-k+1,j-nr+LPML_r])
                       +dt*H*Rz_TR[LPML_z-k+1,j-nr+LPML_r]
                       +dt*C*(
                         RrPE_TR[LPML_z-k+1,j-nr+LPML_r]
                         +1.0/r_now*RpPE_TR[LPML_z-k+1,j-nr+LPML_r]
                         +RzPE_TR[LPML_z-k+1,j-nr+LPML_r])



                 pf_now=pf[k,j]
                 pf[k,j]=update_pf_1st_Por(pf_now,vr_f1,vr_b1,vz_f1,vz_b1,vr_av,
                        vfr_f1,vfr_b1,vfz_f1,vfz_b1,vfr_av,
                        dz,dr,dt,r_now,C,M)


                 #PML addition
                 pf[k,j]=pf[k,j]
                         -dt*C*(Rr_TR[LPML_z-k+1,j-nr+LPML_r]
                            +1.0/r_now*Rp_TR[LPML_z-k+1,j-nr+LPML_r]
                            +Rz_TR[LPML_z-k+1,j-nr+LPML_r])
                         -dt*M*(RrPE_TR[LPML_z-k+1,j-nr+LPML_r]
                            +1.0/r_now*RpPE_TR[LPML_z-k+1,j-nr+LPML_r]
                            +RzPE_TR[LPML_z-k+1,j-nr+LPML_r])
             end #if k==1

             if(j!=nr)
                 r_now=(j-1)*dr+dr/2.0 #for trp,trz Mittet

                 trz_now=trz[k,j]
                 if (Gmat[k,j]*Gmat[k,j+1]*Gmat[k+1,j]*Gmat[k+1,j+1]==0.0)
                      G_av=0.0
                 else
                     gj1=2.0*Gmat[k,j]*Gmat[k+1,j]/(Gmat[k,j]+Gmat[k+1,j])
                     gj2=2.0*Gmat[k,j+1]*Gmat[k+1,j+1]/(Gmat[k,j+1]+Gmat[k+1,j+1])
                     G_av=2.0*gj1*gj2/(gj1+gj2)
                 end

                 vr_f1,vr_b1=vr[k+1,j],vr[k,j]
                 vz_f1,vz_b1=vz[k,j+1],vz[k,j]
                 trz[k,j]=update_trz(trz_now,vr_f1,vr_b1,vz_f1,vz_b1,
                     dz,dr,dt,r_now,G_av,0)
                 #PML addition
                 trz[k,j]=trz[k,j]+dt*(
                    G_av*Srz_TR[LPML_z-k+1,j-nr+LPML_r]+
                    G_av*Rrz_TR[LPML_z-k+1,j-nr+LPML_r])
            end #if j==nr

       end #k (z)
    end #j (r)

end


function PML_update_stress_1st_BottomRight_Por!(vr,vz,trr,tpp,tzz,trz,vfr,vfz,pf,
    Mmat,Cmat,Hmat,Gmat,nr,nz,dr,dz,dt,
    Rr_BR,Rp_BR,Rz_BR,Rrz_BR,Srz_BR,RrPE_BR,RpPE_BR,RzPE_BR,
    LPML_z,LPML_r)
    #PML_variables (PQRS) index starts from first layer in PML (maybe reversed order)
    #z streching only
    #1st order
    # modified rhomat->lmat,mmat

@inbounds Threads.@threads  for j=nr-LPML_r+1:nr #assumig LPML_r>2
    for k=nz-LPML_z+1:nz
             #I am accesing [k,j], but the definition of r_now changes with each component!
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

            #PML addition
            trr[k,j]=trr[k,j]
                   +dt*(H-2G)*(1.0/r_now*Rp_BR[k-nz+LPML_z,j-nr+LPML_r]+Rz_BR[k-nz+LPML_z,j-nr+LPML_r])
                   +dt*H*Rr_BR[k-nz+LPML_z,j-nr+LPML_r]
                   +dt*C*(
                     RrPE_BR[k-nz+LPML_z,j-nr+LPML_r]
                     +1.0/r_now*RpPE_BR[k-nz+LPML_z,j-nr+LPML_r]
                     +RzPE_BR[k-nz+LPML_z,j-nr+LPML_r])


             tpp_now=tpp[k,j]
             tpp[k,j]=update_tpp_1st_Por(tpp_now,vr_f1,vr_b1,vz_f1,vz_b1,vr_av,
                    vfr_f1,vfr_b1,vfz_f1,vfz_b1,vfr_av,
                    dz,dr,dt,r_now,C,H,G)

             #PML addition
             tpp[k,j]=tpp[k,j]
                 +dt*(H-2G)*(Rz_BR[k-nz+LPML_z,j-nr+LPML_r]+Rr_BR[k-nz+LPML_z,j-nr+LPML_r])
                 +dt*H*(1.0/r_now*Rp_BR[k-nz+LPML_z,j-nr+LPML_r])
                 +dt*C*(
                   RrPE_BR[k-nz+LPML_z,j-nr+LPML_r]
                   +1.0/r_now*RpPE_BR[k-nz+LPML_z,j-nr+LPML_r]
                   +RzPE_BR[k-nz+LPML_z,j-nr+LPML_r])


             tzz_now=tzz[k,j]
             tzz[k,j]=update_tzz_1st_Por(tzz_now,vr_f1,vr_b1,vz_f1,vz_b1,vr_av,
                    vfr_f1,vfr_b1,vfz_f1,vfz_b1,vfr_av,
                    dz,dr,dt,r_now,C,H,G)

             #PML addition
             tzz[k,j]=tzz[k,j]
                   +dt*(H-2G)*(1.0/r_now*Rp_BR[k-nz+LPML_z,j-nr+LPML_r]
                                +Rr_BR[k-nz+LPML_z,j-nr+LPML_r])
                   +dt*H*Rz_BR[k-nz+LPML_z,j-nr+LPML_r]
                   +dt*C*(
                     RrPE_BR[k-nz+LPML_z,j-nr+LPML_r]
                     +1.0/r_now*RpPE_BR[k-nz+LPML_z,j-nr+LPML_r]
                     +RzPE_BR[k-nz+LPML_z,j-nr+LPML_r])



             pf_now=pf[k,j]
             pf[k,j]=update_pf_1st_Por(pf_now,vr_f1,vr_b1,vz_f1,vz_b1,vr_av,
                    vfr_f1,vfr_b1,vfz_f1,vfz_b1,vfr_av,
                    dz,dr,dt,r_now,C,M)

             #PML addition
             pf[k,j]=pf[k,j]
                     -dt*C*(Rr_BR[k-nz+LPML_z,j-nr+LPML_r]
                        +1.0/r_now*Rp_BR[k-nz+LPML_z,j-nr+LPML_r]
                        +Rz_BR[k-nz+LPML_z,j-nr+LPML_r])
                     -dt*M*(RrPE_BR[k-nz+LPML_z,j-nr+LPML_r]
                        +1.0/r_now*RpPE_BR[k-nz+LPML_z,j-nr+LPML_r]
                        +RzPE_BR[k-nz+LPML_z,j-nr+LPML_r])


             if(k!=nz && j!=nr)
                 r_now=(j-1)*dr+dr/2.0 #for trp,trz Mittet

                 trz_now=trz[k,j]
                 if (Gmat[k,j]*Gmat[k,j+1]*Gmat[k+1,j]*Gmat[k+1,j+1]==0.0)
                      G_av=0.0
                 else
                     gj1=2.0*Gmat[k,j]*Gmat[k+1,j]/(Gmat[k,j]+Gmat[k+1,j])
                     gj2=2.0*Gmat[k,j+1]*Gmat[k+1,j+1]/(Gmat[k,j+1]+Gmat[k+1,j+1])
                     G_av=2.0*gj1*gj2/(gj1+gj2)
                 end

                 vr_f1,vr_b1=vr[k+1,j],vr[k,j]
                 vz_f1,vz_b1=vz[k,j+1],vz[k,j]
                 trz[k,j]=update_trz(trz_now,vr_f1,vr_b1,vz_f1,vz_b1,
                     dz,dr,dt,r_now,G_av,0)
                 #PML addition
                 trz[k,j]=trz[k,j]+dt*(
                    G_av*Srz_BR[k-nz+LPML_z,j-nr+LPML_r]+
                    G_av*Rrz_BR[k-nz+LPML_z,j-nr+LPML_r])
            end #if k==nz, j==nr

       end #k (z)
    end #j (r)

end



#--------------------------------
# updating PML variables:
# stress (Rx and Sxx)
#--------------------------------
function PML_update_memRS_1st_Top_Por!(Rz,Srz,RzPE,
    memT_vr,memT_vz,memT_vfr,memT_vfz,
    vr,vz,vfr,vfz,nr,nz,dr,dz,dt,LPML_z,LPML_r,PML_Wz,PML_Wz2)
#modified 2: Srz (k=2) 2nd, Srz (k=1) 1st.
#  LPML corresponds to the location of tii, and at tii(k=LPML), OMEGA function is zero
#  Mittet grid. origin is tii
#  k=1: trz is nonzero, tii is zero
#  k=LPML_z: trz is redundant

@inbounds Threads.@threads for j=1:nr-LPML_r #assumig LPML_r>2, and no r derivative (thus j=1 included)
        for k=1:LPML_z #this points normal field index
        #PML profile and PQRS memory var: index starts from first layer
        #memT_ variables: index starts from Top boundary and until LPML+2

        #Rz: corresponds updating tii (same location as tii)
        #dvz/dz of previous and current time
        if(k!=1)

            vz_f1,vz_b1=memT_vz[k,j],memT_vz[k-1,j]
            dzvz_prev=(vz_f1-vz_b1)
            dzvz_prev=dzvz_prev/dz

            vz_f1,vz_b1=vz[k,j],vz[k-1,j]
            dzvz_now=(vz_f1-vz_b1)
            dzvz_now=dzvz_now/dz

            PML_Wz_now=PML_Wz[LPML_z-k+1] #Check
            Rz[LPML_z-k+1,j]=Rz[LPML_z-k+1,j]*exp(-PML_Wz_now*dt)-
                1.0/2.0*PML_Wz_now*dt*(
                exp(-PML_Wz_now*dt)*dzvz_prev+dzvz_now)
        end

        #Srz: corresponds updating trz (same location as trz)
        #dvr/dz of previous and current time

        if(k!=LPML_z)
            vr_f1,vr_b1=memT_vr[k+1,j],memT_vr[k,j]
            dzvr_prev=(vr_f1-vr_b1)
            dzvr_prev=dzvr_prev/dz

            vr_f1,vr_b1=vr[k+1,j],vr[k,j]
            dzvr_now=(vr_f1-vr_b1)
            dzvr_now=dzvr_now/dz

            PML_Wz_now=PML_Wz2[LPML_z-k] #CHECK2 (index is different from PQRS below)

            Srz[LPML_z-k+1,j]=Srz[LPML_z-k+1,j]*exp(-PML_Wz_now*dt)-
                  1.0/2.0*PML_Wz_now*dt*(
                  exp(-PML_Wz_now*dt)*dzvr_prev+dzvr_now)

        end

        #RzPE: Poroelastic term of Rz

        if(k!=1)
            vfz_f1,vfz_b1=memT_vfz[k,j],memT_vfz[k-1,j]
            dzvfz_prev=(vfz_f1-vfz_b1)
            dzvfz_prev=dzvfz_prev/dz

            vfz_f1,vfz_b1=vfz[k,j],vfz[k-1,j]
            dzvfz_now=(vfz_f1-vfz_b1)
            dzvfz_now=dzvfz_now/dz

            PML_Wz_now=PML_Wz[LPML_z-k+1] #Check
            RzPE[LPML_z-k+1,j]=RzPE[LPML_z-k+1,j]*exp(-PML_Wz_now*dt)-
                1.0/2.0*PML_Wz_now*dt*(
                exp(-PML_Wz_now*dt)*dzvfz_prev+dzvfz_now)
        end
    end #k
end#j

end


#Poroelastic
function PML_update_memRS_1st_Bottom_Por!(Rz,Srz,RzPE,
    memB_vr,memB_vz,memB_vfr,memB_vfz,
    vr,vz,vfr,vfz,nr,nz,dr,dz,dt,LPML_z,LPML_r,PML_Wz,PML_Wz2)
    #modefied in 06tmp5
    #mod in Ou1st02tmp:
    #  LPML corresponds to the location of tii, and at tii(k=nz-LPML+1), OMEGA function is zero
    #  Mittet grid. origin is tii
    #
    # k=nz: only vr and tii. trz and vz=0

@inbounds Threads.@threads for j=1:nr-LPML_r #assumig LPML_r>2, and no r derivative (thus j=1 included)
        for k=nz-LPML_z+1:nz #this points normal field index
        #PML profile and PQRS memory var: index starts from first layer
        #memB_ variables: index starts from nz-LPML-1 and until nz

        #Rz: corresponds updating tii (same location as tii)
        #dvz/dz of previous and current time

        k_now=k-nz+LPML_z+2 #caution (starts from 3)
        vz_f1,vz_b1=memB_vz[k_now,j],memB_vz[k_now-1,j]
        dzvz_prev=(vz_f1-vz_b1)
        dzvz_prev=dzvz_prev/dz

        vz_f1,vz_b1=vz[k,j],vz[k-1,j]
        dzvz_now=(vz_f1-vz_b1)
        dzvz_now=dzvz_now/dz

        PML_Wz_now=PML_Wz[k-nz+LPML_z] #Check
        Rz[k-nz+LPML_z,j]=Rz[k-nz+LPML_z,j]*exp(-PML_Wz_now*dt)-
            1.0/2.0*PML_Wz_now*dt*(
            exp(-PML_Wz_now*dt)*dzvz_prev+dzvz_now)

        #Srz: corresponds updating trz (same location as trz)
        #dvr/dz of previous and current time
        if(k!=nz)
            vr_f1,vr_b1=memB_vr[k_now+1,j],memB_vr[k_now,j]
            dzvr_prev=(vr_f1-vr_b1)
            dzvr_prev=dzvr_prev/dz

            vr_f1,vr_b1=vr[k+1,j],vr[k,j]
            dzvr_now=(vr_f1-vr_b1)
            dzvr_now=dzvr_now/dz

            PML_Wz_now=PML_Wz2[k-nz+LPML_z] #Check
            Srz[k-nz+LPML_z,j]=Srz[k-nz+LPML_z,j]*exp(-PML_Wz_now*dt)-
                1.0/2.0*PML_Wz_now*dt*(
                exp(-PML_Wz_now*dt)*dzvr_prev+dzvr_now)
        end #k==nz

        #RzPE: Poroelastic term of Rz

        vfz_f1,vfz_b1=memB_vfz[k_now,j],memB_vfz[k_now-1,j]
        dzvfz_prev=(vfz_f1-vfz_b1)
        dzvfz_prev=dzvfz_prev/dz

        vfz_f1,vfz_b1=vfz[k,j],vfz[k-1,j]
        dzvfz_now=(vfz_f1-vfz_b1)
        dzvfz_now=dzvfz_now/dz

        PML_Wz_now=PML_Wz[k-nz+LPML_z] #Check

        RzPE[k-nz+LPML_z,j]=RzPE[k-nz+LPML_z,j]*exp(-PML_Wz_now*dt)-
            1.0/2.0*PML_Wz_now*dt*(
            exp(-PML_Wz_now*dt)*dzvfz_prev+dzvfz_now)

        end
    end
end


function PML_update_memRS_1st_Right_Por!(Rr_R,Rp_R,Rrz_R,RrPE_R,RpPE_R,
    memR_vr,memR_vz,memR_vfr,memR_vfz,
    vr,vz,vfr,vfz,nr,nz,dr,dz,dt,LPML_z,LPML_r,PML_Wr,PML_IWr,PML_Wr2,PML_IWr2)

    #j-derivatives: Rr_R(j-2.j+1), Rrp_R(j-1.j+2), Rrz_R(j-1.j+2)
    #no j-derivatives: Rp_R,Srp_R,Spz_R

    #j=nr-1: Rp_R,Srp_R,Spz_R (as-is), Rr_R(2nd), Rrp_R,Rrz_R(1st)
    #j=nr: Rp_R,Srp_R,Spz_R (as-is), Rr_R(1st), Rrp_R,Rrz_R(zero)

    #mod in Ou1st02tmp:
    #  LPML corresponds to the location of tii, and at tii(j=nr-LPML+1), OMEGA function is zero
    #  Mittet grid. origin is tii
    # j=nr: only tii and vz. vr and trz are zero

@inbounds Threads.@threads for j=nr-LPML_r+1:nr #this points normal field index
        for k=LPML_z+1:nz-LPML_z #assumig LPML_r>2
        #PML profile and PQRS memory var: index starts from first layer
        #memR_ variables: index starts from nr-LPML-1 and until nr

        #Rr: corresponds updating tii (same location as tii)
        #dvr/dr of previous and current time
        j_now=j-nr+LPML_r+2 #caution (starts from 3)

        vr_f1,vr_b1=memR_vr[k,j_now],memR_vr[k,j_now-1]
        drvr_prev=(vr_f1-vr_b1)
        drvr_prev=drvr_prev/dr

        vr_f1,vr_b1=vr[k,j],vr[k,j-1]
        drvr_now=(vr_f1-vr_b1)
        drvr_now=drvr_now/dr

        #caution k and j index (starts from 1)
        PML_Wr_now=PML_Wr[j-nr+LPML_r] #Check
        Rr_R[k-LPML_z,j-nr+LPML_r]=Rr_R[k-LPML_z,j-nr+LPML_r]*exp(-PML_Wr_now*dt)-
            1.0/2.0*PML_Wr_now*dt*(
            exp(-PML_Wr_now*dt)*drvr_prev+drvr_now)


        #Rp: corresponds updating tii (same location as tii)
        #vr+m*vphi of previous and current time
        j_now=j-nr+LPML_r+2 #caution (starts from 3)

        vr_prev=memR_vr[k,j_now]
        vr_now=vr[k,j]

        v_prev=vr_prev
        v_now=vr_now

        #caution k and j index (starts from 1)
        PML_IWr_now=PML_IWr[j-nr+LPML_r] #Check
        Rp_R[k-LPML_z,j-nr+LPML_r]=Rp_R[k-LPML_z,j-nr+LPML_r]*exp(-PML_IWr_now*dt)-
            1.0/2.0*PML_IWr_now*dt*(
            exp(-PML_IWr_now*dt)*v_prev+v_now)

        #Rrz: corresponds updating trz (same location as trz)
        #dvz/dr of previous and current time
        if(j!=nr)
            j_now=j-nr+LPML_r+2 #caution (starts from 3)

            vz_f1,vz_b1=memR_vz[k,j_now+1],memR_vz[k,j_now]
            drvz_prev=(vz_f1-vz_b1)
            drvz_prev=drvz_prev/dr

            vz_f1,vz_b1=vz[k,j+1],vz[k,j]
            drvz_now=(vz_f1-vz_b1)
            drvz_now=drvz_now/dr

            #caution k and j index (starts from 1)
            PML_Wr_now=PML_Wr2[j-nr+LPML_r] #Check

            Rrz_R[k-LPML_z,j-nr+LPML_r]=Rrz_R[k-LPML_z,j-nr+LPML_r]*exp(-PML_Wr_now*dt)-
                1.0/2.0*PML_Wr_now*dt*(
                exp(-PML_Wr_now*dt)*drvz_prev+drvz_now)
        end


        #RrPE: Poroelasic term of Rr
        j_now=j-nr+LPML_r+2 #caution (starts from 3)
        vfr_f1,vfr_b1=memR_vfr[k,j_now],memR_vfr[k,j_now-1]
        drvfr_prev=(vfr_f1-vfr_b1)
        drvfr_prev=drvfr_prev/dr

        vfr_f1,vfr_b1=vfr[k,j],vfr[k,j-1]
        drvfr_now=(vfr_f1-vfr_b1)
        drvfr_now=drvfr_now/dr

        #caution k and j index (starts from 1)
        PML_Wr_now=PML_Wr[j-nr+LPML_r] #Check
        RrPE_R[k-LPML_z,j-nr+LPML_r]=RrPE_R[k-LPML_z,j-nr+LPML_r]*exp(-PML_Wr_now*dt)-
            1.0/2.0*PML_Wr_now*dt*(
            exp(-PML_Wr_now*dt)*drvfr_prev+drvfr_now)

        #RpPE: Poroelasic term of Rp
        j_now=j-nr+LPML_r+2 #caution (starts from 3)

        vfr_prev=memR_vfr[k,j_now]
        vfr_now=vfr[k,j]

        v_prev=vfr_prev
        v_now=vfr_now

        #caution k and j index (starts from 1)
        PML_IWr_now=PML_IWr[j-nr+LPML_r] #Check
        RpPE_R[k-LPML_z,j-nr+LPML_r]=RpPE_R[k-LPML_z,j-nr+LPML_r]*exp(-PML_IWr_now*dt)-
            1.0/2.0*PML_IWr_now*dt*(
            exp(-PML_IWr_now*dt)*v_prev+v_now)

        end
    end

end


function PML_update_memRS_1st_TopRight_Por!(Rr_TR,Rp_TR,Rz_TR,Rrz_TR,Srz_TR,RrPE_TR,RpPE_TR,RzPE_TR,
    memT_vr,memT_vz,memT_vfr,memT_vfz,
    memR_vr,memR_vz,memR_vfr,memR_vfz,
    vr,vz,vfr,vfz,nr,nz,dr,dz,dt,LPML_z,LPML_r,
    PML_Wz,PML_Wr,PML_IWr,PML_Wz2,PML_Wr2,PML_IWr2)
    #mod in Ou1st02tmp:
    #  LPML corresponds to the location of tii, and at tii(k=LPML), OMEGA function is zero
    #  LPML corresponds to the location of tii, and at tii(j=nr-LPML+1), OMEGA function is zero
    #  Mittet grid. origin is tii
    # j=nr: only tii and vz. vr and trz are zero
    # k=1: only vz and trz. tii and vr=0
    # k=LPML_z: vz and trz is redundant (outside PML)

@inbounds Threads.@threads for j=nr-LPML_r+1:nr #assumig LPML_r>2,
        for k=1:LPML_z #this points normal field index
        #PML profile and PQRS memory var: index starts from first layer
        #memT_ variables: index starts from Top boundary and until LPML+2

        #--z strechings

        #Rz: corresponds updating tii (same location as tii)
        #dvz/dz of previous and current time
        if(k!=1)

            vz_f1,vz_b1=memT_vz[k,j],memT_vz[k-1,j]
            dzvz_prev=(vz_f1-vz_b1)
            dzvz_prev=dzvz_prev/dz

            vz_f1,vz_b1=vz[k,j],vz[k-1,j]
            dzvz_now=(vz_f1-vz_b1)
            dzvz_now=dzvz_now/dz

            PML_Wz_now=PML_Wz[LPML_z-k+1] #Check
            Rz_TR[LPML_z-k+1,j-nr+LPML_r]=Rz_TR[LPML_z-k+1,j-nr+LPML_r]*exp(-PML_Wz_now*dt)-
                1.0/2.0*PML_Wz_now*dt*(
                exp(-PML_Wz_now*dt)*dzvz_prev+dzvz_now)
        end

        #Srz: corresponds updating trz (same location as trz)
        #dvr/dz of previous and current time
        if(k!=LPML_z)

            vr_f1,vr_b1=memT_vr[k+1,j],memT_vr[k,j]
            dzvr_prev=(vr_f1-vr_b1)
            dzvr_prev=dzvr_prev/dz

            vr_f1,vr_b1=vr[k+1,j],vr[k,j]
            dzvr_now=(vr_f1-vr_b1)
            dzvr_now=dzvr_now/dz

            PML_Wz_now=PML_Wz2[LPML_z-k] #CHECK2 (index is different from PQRS below)

            Srz_TR[LPML_z-k+1,j-nr+LPML_r]=Srz_TR[LPML_z-k+1,j-nr+LPML_r]*exp(-PML_Wz_now*dt)-
                  1.0/2.0*PML_Wz_now*dt*(
                  exp(-PML_Wz_now*dt)*dzvr_prev+dzvr_now)
       end
       #RzPE: Poroelastic term of Rz
       if(k!=1)

           vfz_f1,vfz_b1=memT_vfz[k,j],memT_vfz[k-1,j]
           dzvfz_prev=(vfz_f1-vfz_b1)
           dzvfz_prev=dzvfz_prev/dz

           vfz_f1,vfz_b1=vfz[k,j],vfz[k-1,j]
           dzvfz_now=(vfz_f1-vfz_b1)
           dzvfz_now=dzvfz_now/dz

           PML_Wz_now=PML_Wz[LPML_z-k+1] #Check
           RzPE_TR[LPML_z-k+1,j-nr+LPML_r]=RzPE_TR[LPML_z-k+1,j-nr+LPML_r]*exp(-PML_Wz_now*dt)-
                1.0/2.0*PML_Wz_now*dt*(
                exp(-PML_Wz_now*dt)*dzvfz_prev+dzvfz_now)

        end

        # r strechings
        #PML profile and PQRS memory var: index starts from first layer
        #memR_ variables: index starts from nr-LPML-1 and until nr
        #Rr: corresponds updating tii (same location as tii)
        #dvr/dr of previous and current time
        j_now=j-nr+LPML_r+2 #caution (starts from 3)

        vr_f1,vr_b1=memR_vr[k,j_now],memR_vr[k,j_now-1]
        drvr_prev=(vr_f1-vr_b1)
        drvr_prev=drvr_prev/dr

        vr_f1,vr_b1=vr[k,j],vr[k,j-1]
        drvr_now=(vr_f1-vr_b1)
        drvr_now=drvr_now/dr

        #caution k and j index (starts from 1)
        PML_Wr_now=PML_Wr[j-nr+LPML_r] #Check
        Rr_TR[LPML_z-k+1,j-nr+LPML_r]=Rr_TR[LPML_z-k+1,j-nr+LPML_r]*exp(-PML_Wr_now*dt)-
            1.0/2.0*PML_Wr_now*dt*(
            exp(-PML_Wr_now*dt)*drvr_prev+drvr_now)

        #Rp: corresponds updating tii (same location as tii)
        #vr+m*vphi of previous and current time
        j_now=j-nr+LPML_r+2 #caution (starts from 3)

        vr_prev=memR_vr[k,j_now]
        vr_now=vr[k,j]

        v_prev=vr_prev
        v_now=vr_now

        #caution k and j index (starts from 1)
        PML_IWr_now=PML_IWr[j-nr+LPML_r] #Check
        Rp_TR[LPML_z-k+1,j-nr+LPML_r]=Rp_TR[LPML_z-k+1,j-nr+LPML_r]*exp(-PML_IWr_now*dt)-
            1.0/2.0*PML_IWr_now*dt*(
            exp(-PML_IWr_now*dt)*v_prev+v_now)


        #Rrz: corresponds updating trz (same location as trz)
        #dvz/dr of previous and current time
        if(j!=nr)
            j_now=j-nr+LPML_r+2 #caution (starts from 3)

            vz_f1,vz_b1=memR_vz[k,j_now+1],memR_vz[k,j_now]
            drvz_prev=(vz_f1-vz_b1)
            drvz_prev=drvz_prev/dr

            vz_f1,vz_b1=vz[k,j+1],vz[k,j]
            drvz_now=(vz_f1-vz_b1)
            drvz_now=drvz_now/dr

            #caution k and j index (starts from 1)
            PML_Wr_now=PML_Wr2[j-nr+LPML_r] #Check

            Rrz_TR[LPML_z-k+1,j-nr+LPML_r]=Rrz_TR[LPML_z-k+1,j-nr+LPML_r]*exp(-PML_Wr_now*dt)-
                1.0/2.0*PML_Wr_now*dt*(
                exp(-PML_Wr_now*dt)*drvz_prev+drvz_now)
        end

        #RrPE: Poroelasic term of Rr
        j_now=j-nr+LPML_r+2 #caution (starts from 3)
        vfr_f1,vfr_b1=memR_vfr[k,j_now],memR_vfr[k,j_now-1]
        drvfr_prev=(vfr_f1-vfr_b1)
        drvfr_prev=drvfr_prev/dr

        vfr_f1,vfr_b1=vfr[k,j],vfr[k,j-1]
        drvfr_now=(vfr_f1-vfr_b1)
        drvfr_now=drvfr_now/dr

        #caution k and j index (starts from 1)
        PML_Wr_now=PML_Wr[j-nr+LPML_r] #Check

        RrPE_TR[LPML_z-k+1,j-nr+LPML_r]=RrPE_TR[LPML_z-k+1,j-nr+LPML_r]*exp(-PML_Wr_now*dt)-
            1.0/2.0*PML_Wr_now*dt*(
            exp(-PML_Wr_now*dt)*drvfr_prev+drvfr_now)

        #RpPE: Poroelasic term of Rp
        j_now=j-nr+LPML_r+2 #caution (starts from 3)

        vfr_prev=memR_vfr[k,j_now]
        vfr_now=vfr[k,j]

        v_prev=vfr_prev
        v_now=vfr_now

        #caution k and j index (starts from 1)
        PML_IWr_now=PML_IWr[j-nr+LPML_r] #Check

        RpPE_TR[LPML_z-k+1,j-nr+LPML_r]=RpPE_TR[LPML_z-k+1,j-nr+LPML_r]*exp(-PML_IWr_now*dt)-
            1.0/2.0*PML_IWr_now*dt*(
            exp(-PML_IWr_now*dt)*v_prev+v_now)

        end #k
    end #j


end


function PML_update_memRS_1st_BottomRight_Por!(Rr_BR,Rp_BR,Rz_BR,Rrz_BR,Srz_BR,RrPE_BR,RpPE_BR,RzPE_BR,
    memB_vr,memB_vz,memB_vfr,memB_vfz,
    memR_vr,memR_vz,memR_vfr,memR_vfz,
    vr,vz,vfr,vfz,nr,nz,dr,dz,dt,LPML_z,LPML_r,PML_Wz,PML_Wr,PML_IWr,PML_Wz2,PML_Wr2,PML_IWr2)
    #mod in Ou1st02tmp:
    #  LPML corresponds to the location of tii, and at tii(k=nz-LPML+1), OMEGA function is zero
    #  LPML corresponds to the location of tii, and at tii(j=nr-LPML+1), OMEGA function is zero
    #  Mittet grid. origin is tii
    # j=nr: only tii and vz. vr and trz are zero
    # k=nz: only vr and tii. trz and vz=0

@inbounds Threads.@threads for j=nr-LPML_r+1:nr #assumig LPML_r>2,
    for k=nz-LPML_z+1:nz #this points normal field index
        #PML profile and PQRS memory var: index starts from first layer
        #memT_ variables: index starts from Top boundary and until LPML+2

        #--z strechings

        #Rz: corresponds updating tii (same location as tii)
        #dvz/dz of previous and current time
        k_now=k-nz+LPML_z+2 #caution (starts from 3)
        vz_f1,vz_b1=memB_vz[k_now,j],memB_vz[k_now-1,j]
        dzvz_prev=(vz_f1-vz_b1)
        dzvz_prev=dzvz_prev/dz

        vz_f1,vz_b1=vz[k,j],vz[k-1,j]
        dzvz_now=(vz_f1-vz_b1)
        dzvz_now=dzvz_now/dz

        PML_Wz_now=PML_Wz[k-nz+LPML_z] #Check

        Rz_BR[k-nz+LPML_z,j-nr+LPML_r]=Rz_BR[k-nz+LPML_z,j-nr+LPML_r]*exp(-PML_Wz_now*dt)-
            1.0/2.0*PML_Wz_now*dt*(
            exp(-PML_Wz_now*dt)*dzvz_prev+dzvz_now)

        #Srz: corresponds updating trz (same location as trz)
        #dvr/dz of previous and current time
        if(k!=nz)
            vr_f1,vr_b1=memB_vr[k_now+1,j],memB_vr[k_now,j]
            dzvr_prev=(vr_f1-vr_b1)
            dzvr_prev=dzvr_prev/dz

            vr_f1,vr_b1=vr[k+1,j],vr[k,j]
            dzvr_now=(vr_f1-vr_b1)
            dzvr_now=dzvr_now/dz

            PML_Wz_now=PML_Wz2[k-nz+LPML_z] #Check
            Srz_BR[k-nz+LPML_z,j-nr+LPML_r]=Srz_BR[k-nz+LPML_z,j-nr+LPML_r]*exp(-PML_Wz_now*dt)-
                1.0/2.0*PML_Wz_now*dt*(
                exp(-PML_Wz_now*dt)*dzvr_prev+dzvr_now)
        end #k==nz

        #RzPE: Poroelastic term of Rz
        k_now=k-nz+LPML_z+2 #caution (starts from 3)
        vfz_f1,vfz_b1=memB_vfz[k_now,j],memB_vfz[k_now-1,j]
        dzvfz_prev=(vfz_f1-vfz_b1)
        dzvfz_prev=dzvfz_prev/dz

        vfz_f1,vfz_b1=vfz[k,j],vfz[k-1,j]
        dzvfz_now=(vfz_f1-vfz_b1)
        dzvfz_now=dzvfz_now/dz

        PML_Wz_now=PML_Wz[k-nz+LPML_z] #Check

        RzPE_BR[k-nz+LPML_z,j-nr+LPML_r]=RzPE_BR[k-nz+LPML_z,j-nr+LPML_r]*exp(-PML_Wz_now*dt)-
            1.0/2.0*PML_Wz_now*dt*(
            exp(-PML_Wz_now*dt)*dzvfz_prev+dzvfz_now)


        # r strechings
        #PML profile and PQRS memory var: index starts from first layer
        #memR_ variables: index starts from nr-LPML-1 and until nr
        #Rr: corresponds updating tii (same location as tii)
        #dvr/dr of previous and current time
        j_now=j-nr+LPML_r+2 #caution (starts from 3)

        vr_f1,vr_b1=memR_vr[k,j_now],memR_vr[k,j_now-1]
        drvr_prev=(vr_f1-vr_b1)
        drvr_prev=drvr_prev/dr

        vr_f1,vr_b1=vr[k,j],vr[k,j-1]
        drvr_now=(vr_f1-vr_b1)
        drvr_now=drvr_now/dr

        #caution k and j index (starts from 1)
        PML_Wr_now=PML_Wr[j-nr+LPML_r] #Check

        Rr_BR[k-nz+LPML_z,j-nr+LPML_r]=Rr_BR[k-nz+LPML_z,j-nr+LPML_r]*exp(-PML_Wr_now*dt)-
            1.0/2.0*PML_Wr_now*dt*(
            exp(-PML_Wr_now*dt)*drvr_prev+drvr_now)

        #Rp: corresponds updating tii (same location as tii)
        #vr+m*vphi of previous and current time
        j_now=j-nr+LPML_r+2 #caution (starts from 3)

        vr_prev=memR_vr[k,j_now]
        vr_now=vr[k,j]

        v_prev=vr_prev
        v_now=vr_now

        #caution k and j index (starts from 1)
        PML_IWr_now=PML_IWr[j-nr+LPML_r] #Check
        Rp_BR[k-nz+LPML_z,j-nr+LPML_r]=Rp_BR[k-nz+LPML_z,j-nr+LPML_r]*exp(-PML_IWr_now*dt)-
            1.0/2.0*PML_IWr_now*dt*(
            exp(-PML_IWr_now*dt)*v_prev+v_now)


        #Rrz: corresponds updating trz (same location as trz)
        #dvz/dr of previous and current time
        if(j!=nr)
            j_now=j-nr+LPML_r+2 #caution (starts from 3)

            vz_f1,vz_b1=memR_vz[k,j_now+1],memR_vz[k,j_now]
            drvz_prev=(vz_f1-vz_b1)
            drvz_prev=drvz_prev/dr

            vz_f1,vz_b1=vz[k,j+1],vz[k,j]
            drvz_now=(vz_f1-vz_b1)
            drvz_now=drvz_now/dr

            #caution k and j index (starts from 1)
            PML_Wr_now=PML_Wr2[j-nr+LPML_r] #Check

            Rrz_BR[k-nz+LPML_z,j-nr+LPML_r]=Rrz_BR[k-nz+LPML_z,j-nr+LPML_r]*exp(-PML_Wr_now*dt)-
                1.0/2.0*PML_Wr_now*dt*(
                exp(-PML_Wr_now*dt)*drvz_prev+drvz_now)
        end #j==nr


        #RrPE: Poroelasic term of Rr
        j_now=j-nr+LPML_r+2 #caution (starts from 3)

        vfr_f1,vfr_b1=memR_vfr[k,j_now],memR_vfr[k,j_now-1]
        drvfr_prev=(vfr_f1-vfr_b1)
        drvfr_prev=drvfr_prev/dr

        vfr_f1,vfr_b1=vfr[k,j],vfr[k,j-1]
        drvfr_now=(vfr_f1-vfr_b1)
        drvfr_now=drvfr_now/dr

        #caution k and j index (starts from 1)
        PML_Wr_now=PML_Wr[j-nr+LPML_r] #Check
        RrPE_BR[k-nz+LPML_z,j-nr+LPML_r]=Rr_BR[k-nz+LPML_z,j-nr+LPML_r]*exp(-PML_Wr_now*dt)-
            1.0/2.0*PML_Wr_now*dt*(
            exp(-PML_Wr_now*dt)*drvfr_prev+drvfr_now)

        #RpPE: Poroelasic term of Rp
        j_now=j-nr+LPML_r+2 #caution (starts from 3)

        vfr_prev=memR_vfr[k,j_now]
        vfr_now=vfr[k,j]

        v_prev=vfr_prev
        v_now=vfr_now

        #caution k and j index (starts from 1)
        PML_IWr_now=PML_IWr[j-nr+LPML_r] #Check
        RpPE_BR[k-nz+LPML_z,j-nr+LPML_r]=RpPE_BR[k-nz+LPML_z,j-nr+LPML_r]*exp(-PML_IWr_now*dt)-
            1.0/2.0*PML_IWr_now*dt*(
            exp(-PML_IWr_now*dt)*v_prev+v_now)


        end #k
    end #j

end


#--------------------------------
# updating PML variables:
# velocity (Pxx)
#--------------------------------
function PML_update_memPQ_1st_Top_Por!(Prz,Pzz,PzzPE,
    memT_trr,memT_tpp,memT_tzz,memT_trz,memT_pf,
    trr,tpp,tzz,trz,pf,nr,nz,dr,dz,dt,LPML_z,LPML_r,PML_Wz,PML_Wz2)
#mod in 06tmp5
#mod in Ou1st02tmp:
#  LPML corresponds to the location of tii, and at tii(k=LPML), OMEGA function is zero
#  Mittet grid. origin is tii
#
# k=1: only vz and trz. tii and vr=0
# k=LPML_z: vz and trz is redundant (outside PML)
@inbounds Threads.@threads for j=1:nr-LPML_r #assumig LPML_r>2, and no r derivative (thus j=1 included)
        for k=1:LPML_z #this points normal field index
        #PML profile and PQ memory var: index starts from first layer
        #memT_ variables: index starts from Top boundary and until LPML+2

        #Prz: corresponds updating vr (same location as vr)
        #dtrz/dz of previous and current time
        if(k!=1)
            trz_f1,trz_b1=memT_trz[k,j],memT_trz[k-1,j]
            dztrz_prev=(trz_f1-trz_b1)
            dztrz_prev=dztrz_prev/dz

            trz_f1,trz_b1=trz[k,j],trz[k-1,j]
            dztrz_now=(trz_f1-trz_b1)
            dztrz_now=dztrz_now/dz

            PML_Wz_now=PML_Wz[LPML_z-k+1] #Check
            Prz[LPML_z-k+1,j]=Prz[LPML_z-k+1,j]*exp(-PML_Wz_now*dt)-
                1.0/2.0*PML_Wz_now*dt*(
                exp(-PML_Wz_now*dt)*dztrz_prev+dztrz_now)
        end

        #Pzz: corresponds updating vz (same location as vz)
        #dtzz/dz of previous and current time
        if(k!=LPML_z)
            tzz_f1,tzz_b1=memT_tzz[k+1,j],memT_tzz[k,j]
            dztzz_prev=(tzz_f1-tzz_b1)
            dztzz_prev=dztzz_prev/dz

            tzz_f1,tzz_b1=tzz[k+1,j],tzz[k,j]
            dztzz_now=(tzz_f1-tzz_b1)
            dztzz_now=dztzz_now/dz

            PML_Wz_now=PML_Wz2[LPML_z-k] #CHECK2(index is different from PQRS below)

            Pzz[LPML_z-k+1,j]=Pzz[LPML_z-k+1,j]*exp(-PML_Wz_now*dt)-
                  1.0/2.0*PML_Wz_now*dt*(
                  exp(-PML_Wz_now*dt)*dztzz_prev+dztzz_now)
        end


        #PzzPE: Poroelastic term of Pzz
        if(k!=LPML_z)
            pf_f1,pf_b1=memT_pf[k+1,j],memT_pf[k,j]
            dzpf_prev=(pf_f1-pf_b1)
            dzpf_prev=dzpf_prev/dz

            pf_f1,pf_b1=pf[k+1,j],pf[k,j]
            dzpf_now=(pf_f1-pf_b1)
            dzpf_now=dzpf_now/dz


            PML_Wz_now=PML_Wz2[LPML_z-k] #CHECK2(index is different from PQRS below)

            PzzPE[LPML_z-k+1,j]=PzzPE[LPML_z-k+1,j]*exp(-PML_Wz_now*dt)-
                  1.0/2.0*PML_Wz_now*dt*(
                  exp(-PML_Wz_now*dt)*dzpf_prev+dzpf_now)
        end


        end

    end

end


function PML_update_memPQ_1st_Bottom_Por!(Prz,Pzz,PzzPE,
    memB_trr,memB_tpp,memB_tzz,memB_trz,memB_pf,
    trr,tpp,tzz,trz,pf,nr,nz,dr,dz,dt,LPML_z,LPML_r,PML_Wz,PML_Wz2)
#  LPML corresponds to the location of tii, and at tii(k=nz-LPML+1), OMEGA function is zero
#  Mittet grid. origin is tii
#
# k=nz: only vr and tii. trz and vz=0

@inbounds Threads.@threads for j=1:nr-LPML_r #assumig LPML_r>2, and no r derivative (thus j=1 included)
#        for k=3:LPML_z #this points normal field index
        for k=nz-LPML_z+1:nz #this points normal field index
        #PML profile and PQ memory var: index starts from first layer
        #memB_ variables: index starts from nz-LPML-1 and until nz

        #Prz: corresponds updating vr (same location as vr)
        #dtrz/dz of previous and current time
        k_now=k-nz+LPML_z+2 #caution (starts from 3)
        trz_f1,trz_b1=memB_trz[k_now,j],memB_trz[k_now-1,j]
        dztrz_prev=(trz_f1-trz_b1)
        dztrz_prev=dztrz_prev/dz

        trz_f1,trz_b1=trz[k,j],trz[k-1,j]
        dztrz_now=(trz_f1-trz_b1)
        dztrz_now=dztrz_now/dz

        PML_Wz_now=PML_Wz[k-nz+LPML_z] #Check
        Prz[k-nz+LPML_z,j]=Prz[k-nz+LPML_z,j]*exp(-PML_Wz_now*dt)-
            1.0/2.0*PML_Wz_now*dt*(
            exp(-PML_Wz_now*dt)*dztrz_prev+dztrz_now)


        #Pzz: corresponds updating vz (same location as vz)
        #dtzz/dz of previous and current time
        if(k!=nz)
            tzz_f1,tzz_b1=memB_tzz[k_now+1,j],memB_tzz[k_now,j]
            dztzz_prev=(tzz_f1-tzz_b1)
            dztzz_prev=dztzz_prev/dz

            tzz_f1,tzz_b1=tzz[k+1,j],tzz[k,j]
            dztzz_now=(tzz_f1-tzz_b1)
            dztzz_now=dztzz_now/dz

             PML_Wz_now=PML_Wz2[k-nz+LPML_z] #Check
             Pzz[k-nz+LPML_z,j]=Pzz[k-nz+LPML_z,j]*exp(-PML_Wz_now*dt)-
                     1.0/2.0*PML_Wz_now*dt*(
                     exp(-PML_Wz_now*dt)*dztzz_prev+dztzz_now)
         end

         #PzzPE: Poroelasic term of Pzz
         if(k!=nz)
             pf_f1,pf_b1=memB_pf[k_now+1,j],memB_pf[k_now,j]
             dzpf_prev=(pf_f1-pf_b1)
             dzpf_prev=dzpf_prev/dz

             pf_f1,pf_b1=pf[k+1,j],pf[k,j]
             dzpf_now=(pf_f1-pf_b1)
             dzpf_now=dzpf_now/dz

              PML_Wz_now=PML_Wz2[k-nz+LPML_z] #Check
              PzzPE[k-nz+LPML_z,j]=PzzPE[k-nz+LPML_z,j]*exp(-PML_Wz_now*dt)-
                      1.0/2.0*PML_Wz_now*dt*(
                      exp(-PML_Wz_now*dt)*dzpf_prev+dzpf_now)
          end #if k==nz

      end #k
  end #j


end


function PML_update_memPQ_1st_Right_Por!(Prr_R,Qrp_R,Pzr_R,Qzp_R,PrrPE_R,
    memR_trr,memR_tpp,memR_tzz,memR_trz,memR_pf,
    trr,tpp,tzz,trz,pf,nr,nz,dr,dz,dt,LPML_z,LPML_r,PML_Wr,PML_IWr,PML_Wr2,PML_IWr2)

    #j-derivaitves Prr(j-1.j+2), Ppr(j-2.j+1), Pzr(j-2.j+1)
    #no j-drivatives Qrp,Qpp,Qzp

    #j=nr-1->Qrp,Qpp,Qzp(2nd), Ppr,Pzr(2nd), Prr (1st)
    #j=nr->Qrp,Qpp,Qzp(2nd), Ppr,Pzr(1st), Prr (zero)

    #mod in Ou1st02tmp:
    #  LPML corresponds to the location of tii, and at tii(j=nr-LPML+1), OMEGA function is zero
    #  Mittet grid. origin is tii
    # j=nr: only tii and vz. vr and trz are zero

@inbounds Threads.@threads for j=nr-LPML_r+1:nr #this points normal field index
        for k=LPML_z+1:nz-LPML_z #assumig LPML_r>2
        #PML profile and PQ memory var: index starts from first layer
        #memR_ variables: index starts from nr-LPML-1 and until nr
        #Prr: corresponds updating vr (same location as vr)
        #dtrr/dr of previous and current time
        if(j!=nr)
            j_now=j-nr+LPML_r+2 #caution (starts from 3)
            trr_f1,trr_b1=memR_trr[k,j_now+1],memR_trr[k,j_now]
            drtrr_prev=(trr_f1-trr_b1)
            drtrr_prev=drtrr_prev/dr

            trr_f1,trr_b1=trr[k,j+1],trr[k,j]
            drtrr_now=(trr_f1-trr_b1)
            drtrr_now=drtrr_now/dr

            PML_Wr_now=PML_Wr2[j-nr+LPML_r] #CHECK

            Prr_R[k-LPML_z,j-nr+LPML_r]=Prr_R[k-LPML_z,j-nr+LPML_r]*exp(-PML_Wr_now*dt)-
                1.0/2.0*PML_Wr_now*dt*(
                exp(-PML_Wr_now*dt)*drtrr_prev+drtrr_now)
        end

        #Qrp: corresponds updating vr (same location as vr)
        #trr-tpp+mtrp of previous and current time
        if(j!=nr)
            j_now=j-nr+LPML_r+2 #caution (starts from 3)
            trr_prev=memR_trr[k,j_now]
            tpp_prev=memR_tpp[k,j_now]
            trr_now=trr[k,j]
            tpp_now=tpp[k,j]

            tii_prev=trr_prev-tpp_prev
            tii_now=trr_now-tpp_now

            PML_IWr_now=PML_IWr2[j-nr+LPML_r] #Check

            Qrp_R[k-LPML_z,j-nr+LPML_r]=Qrp_R[k-LPML_z,j-nr+LPML_r]*exp(-PML_IWr_now*dt)-
                1.0/2.0*PML_IWr_now*dt*(
                exp(-PML_IWr_now*dt)*tii_prev+tii_now)
        end

        #Pzr: corresponds updating vz (same location as vz)
        #dtrz/dr of previous and current time
        j_now=j-nr+LPML_r+2 #caution (starts from 3)
        trz_f1,trz_b1=memR_trz[k,j_now],memR_trz[k,j_now-1]
        drtrz_prev=(trz_f1-trz_b1)
        drtrz_prev=drtrz_prev/dr

        trz_f1,trz_b1=trz[k,j],trz[k,j-1]
        drtrz_now=(trz_f1-trz_b1)
        drtrz_now=drtrz_now/dr

        #caution k and j index (starts from 1)
        PML_Wr_now=PML_Wr[j-nr+LPML_r] #Check

        Pzr_R[k-LPML_z,j-nr+LPML_r]=Pzr_R[k-LPML_z,j-nr+LPML_r]*exp(-PML_Wr_now*dt)-
            1.0/2.0*PML_Wr_now*dt*(
            exp(-PML_Wr_now*dt)*drtrz_prev+drtrz_now)

        #Qzp: corresponds updating vz (same location as vz)
        #trz+m*tpz of previous and current time
        j_now=j-nr+LPML_r+2 #caution (starts from 3)
        trz_prev=memR_trz[k,j_now]
        trz_now=trz[k,j]

        tii_prev=trz_prev
        tii_now=trz_now
        #caution k and j index (starts from 1)
        PML_IWr_now=PML_IWr[j-nr+LPML_r] #Check

        Qzp_R[k-LPML_z,j-nr+LPML_r]=Qzp_R[k-LPML_z,j-nr+LPML_r]*exp(-PML_IWr_now*dt)-
            1.0/2.0*PML_IWr_now*dt*(
            exp(-PML_IWr_now*dt)*tii_prev+tii_now)


        #PrrPE: Poroelasic term of Prr
        if(j!=nr)
            j_now=j-nr+LPML_r+2 #caution (starts from 3)
            pf_f1,pf_b1=memR_pf[k,j_now+1],memR_pf[k,j_now]
            drpf_prev=(pf_f1-pf_b1)
            drpf_prev=drpf_prev/dr

            pf_f1,pf_b1=pf[k,j+1],pf[k,j]
            drpf_now=(pf_f1-pf_b1)
            drpf_now=drpf_now/dr

            PML_Wr_now=PML_Wr2[j-nr+LPML_r] #check
            PrrPE_R[k-LPML_z,j-nr+LPML_r]=PrrPE_R[k-LPML_z,j-nr+LPML_r]*exp(-PML_Wr_now*dt)-
                1.0/2.0*PML_Wr_now*dt*(
                exp(-PML_Wr_now*dt)*drpf_prev+drpf_now)
        end #if

        end
    end

end



function PML_update_memPQ_1st_TopRight_Por!(Prr_TR,Qrp_TR,Prz_TR,Pzr_TR,Qzp_TR,Pzz_TR,PrrPE_TR,PzzPE_TR,
    memT_trr,memT_tpp,memT_tzz,memT_trz,memT_pf,
    memR_trr,memR_tpp,memR_tzz,memR_trz,memR_pf,
    trr,tpp,tzz,trz,pf,nr,nz,dr,dz,dt,LPML_z,LPML_r,
    PML_Wz,PML_Wr,PML_IWr,PML_Wz2,PML_Wr2,PML_IWr2)

    #j-derivatives: Prr(j-1.j+2), Ppr(j-2.j+1), Pzr(j-2.j+1)
    #no j-detrivatives: Prz(k-2.k+1),Ppz(k-2.k+1),Pzz(k-1.k+2), Qrp,Qpp,Qzp

    #mod in Ou1st02tmp:
    #  LPML corresponds to the location of tii, and at tii(k=LPML), OMEGA function is zero
    #  LPML corresponds to the location of tii, and at tii(j=nr-LPML+1), OMEGA function is zero
    #  Mittet grid. origin is tii
    # j=nr: only tii and vz. vr and trz are zero
    # k=1: only vz and trz. tii and vr=0
    # k=LPML_z: vz and trz is redundant (outside PML)


#@inbounds Threads.@threads
for j=nr-LPML_r+1:nr #assumig LPML_r>2,
        for k=1:LPML_z #this points normal field index
        #PML profile and PQ memory var: index starts from first layer
        #memT_ variables: k-index starts from Top boundary and until LPML+2
        #               : j-index is same as original

        #--z strechings

        #Prz: corresponds updating vr (same location as vr)
        #dtrz/dz of previous and current time
        if(k!=1)
            trz_f1,trz_b1=memT_trz[k,j],memT_trz[k-1,j]
            dztrz_prev=(trz_f1-trz_b1)
            dztrz_prev=dztrz_prev/dz

            trz_f1,trz_b1=trz[k,j],trz[k-1,j]
            dztrz_now=(trz_f1-trz_b1)
            dztrz_now=dztrz_now/dz

            #caution k and j index (k starts from deeper PML, j starts from 1)
            PML_Wz_now=PML_Wz[LPML_z-k+1] #Check
            Prz_TR[LPML_z-k+1,j-nr+LPML_r]=Prz_TR[LPML_z-k+1,j-nr+LPML_r]*exp(-PML_Wz_now*dt)-
                1.0/2.0*PML_Wz_now*dt*(
                exp(-PML_Wz_now*dt)*dztrz_prev+dztrz_now)
        end


        #Pzz: corresponds updating vz (same location as vz)
        #dtzz/dz of previous and current time
        if(k!=LPML_z)
            tzz_f1,tzz_b1=memT_tzz[k+1,j],memT_tzz[k,j]
            dztzz_prev=(tzz_f1-tzz_b1)
            dztzz_prev=dztzz_prev/dz

            tzz_f1,tzz_b1=tzz[k+1,j],tzz[k,j]
            dztzz_now=(tzz_f1-tzz_b1)
            dztzz_now=dztzz_now/dz

            PML_Wz_now=PML_Wz2[LPML_z-k] #CHECK2(index is different from PQRS below)

            Pzz_TR[LPML_z-k+1,j-nr+LPML_r]=Pzz_TR[LPML_z-k+1,j-nr+LPML_r]*exp(-PML_Wz_now*dt)-
                  1.0/2.0*PML_Wz_now*dt*(
                  exp(-PML_Wz_now*dt)*dztzz_prev+dztzz_now)
        end

      #PzzPE: Poroelastic term of Pzz
      if(k!=LPML_z)
          pf_f1,pf_b1=memT_pf[k+1,j],memT_pf[k,j]
          dzpf_prev=(pf_f1-pf_b1)
          dzpf_prev=dzpf_prev/dz

          pf_f1,pf_b1=pf[k+1,j],pf[k,j]
          dzpf_now=(pf_f1-pf_b1)
          dzpf_now=dzpf_now/dz

          PML_Wz_now=PML_Wz2[LPML_z-k] #CHECK2(index is different from PQRS below)

          PzzPE_TR[LPML_z-k+1,j-nr+LPML_r]=PzzPE_TR[LPML_z-k+1,j-nr+LPML_r]*exp(-PML_Wz_now*dt)-
               1.0/2.0*PML_Wz_now*dt*(
               exp(-PML_Wz_now*dt)*dzpf_prev+dzpf_now)
      end

      #--r strechings
      #PML profile and PQ memory var: index starts from first layer
      #memR_ variables: index starts from nr-LPML-1 and until nr
      #Prr: corresponds updating vr (same location as vr)
      #dtrr/dr of previous and current time
      if(j!=nr)
          j_now=j-nr+LPML_r+2 #caution (starts from 3)
          trr_f1,trr_b1=memR_trr[k,j_now+1],memR_trr[k,j_now]
          drtrr_prev=(trr_f1-trr_b1)
          drtrr_prev=drtrr_prev/dr

          trr_f1,trr_b1=trr[k,j+1],trr[k,j]
          drtrr_now=(trr_f1-trr_b1)
          drtrr_now=drtrr_now/dr

          PML_Wr_now=PML_Wr2[j-nr+LPML_r] #CHECK

          #caution k and j index
          Prr_TR[LPML_z-k+1,j-nr+LPML_r]=Prr_TR[LPML_z-k+1,j-nr+LPML_r]*exp(-PML_Wr_now*dt)-
              1.0/2.0*PML_Wr_now*dt*(
              exp(-PML_Wr_now*dt)*drtrr_prev+drtrr_now)
     end

      #Qrp: corresponds updating vr (same location as vr)
      #trr-tpp+mtrp of previous and current time
      if(j!=nr)
          j_now=j-nr+LPML_r+2 #caution (starts from 3)
          trr_prev=memR_trr[k,j_now]
          tpp_prev=memR_tpp[k,j_now]
          trr_now=trr[k,j]
          tpp_now=tpp[k,j]

          tii_prev=trr_prev-tpp_prev
          tii_now=trr_now-tpp_now
          #caution k and j index (starts from 1)
          PML_IWr_now=PML_IWr2[j-nr+LPML_r] #Check

          Qrp_TR[LPML_z-k+1,j-nr+LPML_r]=Qrp_TR[LPML_z-k+1,j-nr+LPML_r]*exp(-PML_IWr_now*dt)-
              1.0/2.0*PML_IWr_now*dt*(
              exp(-PML_IWr_now*dt)*tii_prev+tii_now)
      end

      #Pzr: corresponds updating vz (same location as vz)
      #dtrz/dr of previous and current time
      j_now=j-nr+LPML_r+2 #caution (starts from 3)
      trz_f1,trz_b1=memR_trz[k,j_now],memR_trz[k,j_now-1]
      drtrz_prev=(trz_f1-trz_b1)
      drtrz_prev=drtrz_prev/dr

      trz_f1,trz_b1=trz[k,j],trz[k,j-1]
      drtrz_now=(trz_f1-trz_b1)
      drtrz_now=drtrz_now/dr

      #caution k and j index
      PML_Wr_now=PML_Wr[j-nr+LPML_r] #Check

      Pzr_TR[LPML_z-k+1,j-nr+LPML_r]=Pzr_TR[LPML_z-k+1,j-nr+LPML_r]*exp(-PML_Wr_now*dt)-
          1.0/2.0*PML_Wr_now*dt*(
          exp(-PML_Wr_now*dt)*drtrz_prev+drtrz_now)

      #Qzp: corresponds updating vz (same location as vz)
      #trz+m*tpz of previous and current time
      j_now=j-nr+LPML_r+2 #caution (starts from 3)
      trz_prev=memR_trz[k,j_now]
      trz_now=trz[k,j]

      tii_prev=trz_prev
      tii_now=trz_now
      #caution k and j index
      PML_IWr_now=PML_IWr[j-nr+LPML_r] #Check

      Qzp_TR[LPML_z-k+1,j-nr+LPML_r]=Qzp_TR[LPML_z-k+1,j-nr+LPML_r]*exp(-PML_IWr_now*dt)-
          1.0/2.0*PML_IWr_now*dt*(
          exp(-PML_IWr_now*dt)*tii_prev+tii_now)

      #PrrPE: Poroelasic term of Prr
      if(j!=nr)
          j_now=j-nr+LPML_r+2 #caution (starts from 3)
          pf_f1,pf_b1=memR_pf[k,j_now+1],memR_pf[k,j_now]
          drpf_prev=(pf_f1-pf_b1)
          drpf_prev=drpf_prev/dr

          pf_f1,pf_b1=pf[k,j+1],pf[k,j]
          drpf_now=(pf_f1-pf_b1)
          drpf_now=drpf_now/dr

          #caution k and j index
          PML_Wr_now=PML_Wr2[j-nr+LPML_r] #Check
          PrrPE_TR[LPML_z-k+1,j-nr+LPML_r]=PrrPE_TR[LPML_z-k+1,j-nr+LPML_r]*exp(-PML_Wr_now*dt)-
              1.0/2.0*PML_Wr_now*dt*(
              exp(-PML_Wr_now*dt)*drpf_prev+drpf_now)
      end #if

      end #k
  end #j

end


function PML_update_memPQ_1st_BottomRight_Por!(Prr_BR,Qrp_BR,Prz_BR,Pzr_BR,Qzp_BR,Pzz_BR,PrrPE_BR,PzzPE_BR,
    memB_trr,memB_tpp,memB_tzz,memB_trz,memB_pf,
    memR_trr,memR_tpp,memR_tzz,memR_trz,memR_pf,
    trr,tpp,tzz,trz,pf,nr,nz,dr,dz,dt,LPML_z,LPML_r,PML_Wz,PML_Wr,PML_IWr,PML_Wz2,PML_Wr2,PML_IWr2)
    #mod in Ou1st02tmp:
    #  LPML corresponds to the location of tii, and at tii(k=nz-LPML+1), OMEGA function is zero
    #  LPML corresponds to the location of tii, and at tii(j=nr-LPML+1), OMEGA function is zero
    #  Mittet grid. origin is tii
    # j=nr: only tii and vz. vr and trz are zero
    # k=nz: only vr and tii. trz and vz=0

@inbounds Threads.@threads for j=nr-LPML_r+1:nr #assumig LPML_r>2,
        for k=nz-LPML_z+1:nz #this points normal field index
        #PML profile and PQ memory var: index starts from first layer
        #memT_ variables: k-index starts from Top boundary and until LPML+2
        #               : j-index is same as original

        #--z strechings

        #Prz: corresponds updating vr (same location as vr)
        #dtrz/dz of previous and current time
        k_now=k-nz+LPML_z+2 #caution (starts from 3)
        trz_f1,trz_b1=memB_trz[k_now,j],memB_trz[k_now-1,j]
        dztrz_prev=(trz_f1-trz_b1)
        dztrz_prev=dztrz_prev/dz

        trz_f1,trz_b1=trz[k,j],trz[k-1,j]
        dztrz_now=(trz_f1-trz_b1)
        dztrz_now=dztrz_now/dz

        #caution k and j index (k starts from 1, j starts from 1)
        PML_Wz_now=PML_Wz[k-nz+LPML_z] #Check

        Prz_BR[k-nz+LPML_z,j-nr+LPML_r]=Prz_BR[k-nz+LPML_z,j-nr+LPML_r]*exp(-PML_Wz_now*dt)-
            1.0/2.0*PML_Wz_now*dt*(
            exp(-PML_Wz_now*dt)*dztrz_prev+dztrz_now)


        #Pzz: corresponds updating vz (same location as vz)
        #dtzz/dz of previous and current time
        if(k!=nz)
            tzz_f1,tzz_b1=memB_tzz[k_now+1,j],memB_tzz[k_now,j]
            dztzz_prev=(tzz_f1-tzz_b1)
            dztzz_prev=dztzz_prev/dz

            tzz_f1,tzz_b1=tzz[k+1,j],tzz[k,j]
            dztzz_now=(tzz_f1-tzz_b1)
            dztzz_now=dztzz_now/dz

            PML_Wz_now=PML_Wz2[k-nz+LPML_z] #Check
            Pzz_BR[k-nz+LPML_z,j-nr+LPML_r]=Pzz_BR[k-nz+LPML_z,j-nr+LPML_r]*exp(-PML_Wz_now*dt)-
                1.0/2.0*PML_Wz_now*dt*(
                exp(-PML_Wz_now*dt)*dztzz_prev+dztzz_now)
        end #k==nz

      #PzzPE: Poroelastic term of Pzz
      if(k!=nz)
          pf_f1,pf_b1=memB_pf[k_now+1,j],memB_pf[k_now,j]
          dzpf_prev=(pf_f1-pf_b1)
          dzpf_prev=dzpf_prev/dz

          pf_f1,pf_b1=pf[k+1,j],pf[k,j]
          dzpf_now=(pf_f1-pf_b1)
          dzpf_now=dzpf_now/dz

       #accounting for dz/2 difference in PML_Wz when Pzz or vz
          PML_Wz_now=PML_Wz2[k-nz+LPML_z] #Check
          PzzPE_BR[k-nz+LPML_z,j-nr+LPML_r]=PzzPE_BR[k-nz+LPML_z,j-nr+LPML_r]*exp(-PML_Wz_now*dt)-
              1.0/2.0*PML_Wz_now*dt*(
              exp(-PML_Wz_now*dt)*dzpf_prev+dzpf_now)
      end #k==nz

      #--r strechings
      #PML profile and PQ memory var: index starts from first layer
      #memR_ variables: index starts from nr-LPML-1 and until nr
      #Prr: corresponds updating vr (same location as vr)
      #dtrr/dr of previous and current time
      if(j!=nr)
          j_now=j-nr+LPML_r+2 #caution (starts from 3)
          trr_f1,trr_b1=memR_trr[k,j_now+1],memR_trr[k,j_now]
          drtrr_prev=(trr_f1-trr_b1)
          drtrr_prev=drtrr_prev/dr

          trr_f1,trr_b1=trr[k,j+1],trr[k,j]
          drtrr_now=(trr_f1-trr_b1)
          drtrr_now=drtrr_now/dr

          #caution k and j index
          PML_Wr_now=PML_Wr2[j-nr+LPML_r] #CHECK
          Prr_BR[k-nz+LPML_z,j-nr+LPML_r]=Prr_BR[k-nz+LPML_z,j-nr+LPML_r]*exp(-PML_Wr_now*dt)-
              1.0/2.0*PML_Wr_now*dt*(
              exp(-PML_Wr_now*dt)*drtrr_prev+drtrr_now)
      end #j==nr

      #Qrp: corresponds updating vr (same location as vr)
      #trr-tpp+mtrp of previous and current time
      if(j!=nr)
          j_now=j-nr+LPML_r+2 #caution (starts from 3)
          trr_prev=memR_trr[k,j_now]
          tpp_prev=memR_tpp[k,j_now]
          trr_now=trr[k,j]
          tpp_now=tpp[k,j]

          tii_prev=trr_prev-tpp_prev
          tii_now=trr_now-tpp_now

          #caution k and j index (starts from 1)
          PML_IWr_now=PML_IWr2[j-nr+LPML_r] #Check

          Qrp_BR[k-nz+LPML_z,j-nr+LPML_r]=Qrp_BR[k-nz+LPML_z,j-nr+LPML_r]*exp(-PML_IWr_now*dt)-
              1.0/2.0*PML_IWr_now*dt*(
              exp(-PML_IWr_now*dt)*tii_prev+tii_now)
      end #j==nr

      #Pzr: corresponds updating vz (same location as vz)
      #dtrz/dr of previous and current time
      j_now=j-nr+LPML_r+2 #caution (starts from 3)
      trz_f1,trz_b1=memR_trz[k,j_now],memR_trz[k,j_now-1]
      drtrz_prev=(trz_f1-trz_b1)
      drtrz_prev=drtrz_prev/dr

      trz_f1,trz_b1=trz[k,j],trz[k,j-1]
      drtrz_now=(trz_f1-trz_b1)
      drtrz_now=drtrz_now/dr

      #caution k and j index
      PML_Wr_now=PML_Wr[j-nr+LPML_r] #Check

      Pzr_BR[k-nz+LPML_z,j-nr+LPML_r]=Pzr_BR[k-nz+LPML_z,j-nr+LPML_r]*exp(-PML_Wr_now*dt)-
          1.0/2.0*PML_Wr_now*dt*(
          exp(-PML_Wr_now*dt)*drtrz_prev+drtrz_now)

      #Qzp: corresponds updating vz (same location as vz)
      #trz+m*tpz of previous and current time
      j_now=j-nr+LPML_r+2 #caution (starts from 3)
      trz_prev=memR_trz[k,j_now]
      trz_now=trz[k,j]

      tii_prev=trz_prev
      tii_now=trz_now

      #caution k and j index
      PML_IWr_now=PML_IWr[j-nr+LPML_r] #Check

      Qzp_BR[k-nz+LPML_z,j-nr+LPML_r]=Qzp_BR[k-nz+LPML_z,j-nr+LPML_r]*exp(-PML_IWr_now*dt)-
          1.0/2.0*PML_IWr_now*dt*(
          exp(-PML_IWr_now*dt)*tii_prev+tii_now)


      #PrrPE: Poroelasic term of Prr
      if(j!=nr)
          j_now=j-nr+LPML_r+2 #caution (starts from 3)
          pf_f1,pf_b1=memR_pf[k,j_now+1],memR_pf[k,j_now]
          drpf_prev=(pf_f1-pf_b1)
          drpf_prev=drpf_prev/dr

          pf_f1,pf_b1=pf[k,j+1],pf[k,j]
          drpf_now=(pf_f1-pf_b1)
          drpf_now=drpf_now/dr

          #caution k and j index
          PML_Wr_now=PML_Wr2[j-nr+LPML_r] #check
          PrrPE_BR[k-nz+LPML_z,j-nr+LPML_r]=PrrPE_BR[k-nz+LPML_z,j-nr+LPML_r]*exp(-PML_Wr_now*dt)-
              1.0/2.0*PML_Wr_now*dt*(
              exp(-PML_Wr_now*dt)*drpf_prev+drpf_now)
      end #j==nr

      end #k
  end #j

end


#----------------------------
# Boundary condition at r=0
# in PML
#----------------------------
function ApplyBCLeft_velocity_1st_atPML_Top_Por!(vr,vz,trr,tpp,tzz,trz,vfr,vfz,pf,
    rhomat,rhofmat,D1mat,D2mat,Flagmat_vf_zero,nr,nz,dr,dz,dt,
    Prz_T,Pzz_T,PzzPE_T,
    LPML_z)
#Filling field values at j==1 (r) & k=2:nz-1
#j=1 corresponds to r=0 for variables of vz[], vphi[], tii[], tpz[]. Other values are @dr/2, not r=0
#2nd order FD (Graves, 1996)
#Supplement discussion from Mittet, 1996, Geophysics,61,1,21-33.
# m=0, monopole: vr(r)=-vr(-r),vphi=-vphi(-r),vz(r)=vz(-r), tij(r)=tij(-r)
# m=1, dipole: vr(r)=vr(-r),vphi=vphi(-r),vz(r)=-vz(-r), tij(r)=-tij(-r)
# nonzero@r=0 --> only "symmetric" components, otherwize 0.
#m=0 only!
#PML-only version

#Borehole axis (Left Edge): d/dr=0 (Symmetric BC)
j=1
#Variables@left-edge are vz[], vphi[], tii[], tpz[] Mittet
#Case m=0: vr=0, trz=0 (Randall, Mittet->fluid), trp=0 (Randall, Mittet->fluid)

#PML (TOP)
for k=1:LPML_z

    #I am alywas accesing [k,j], but definition of r_now changes with components!
   if(k!=1)
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

       vfr_old=vfr_now
       if(Flagmat_vf_zero[k,j]==1)
           vfr[k,j]=0.0
       else
           vr_dummy,vfr[k,j]=update_vr_1st_Por(vr_now,vfr_now,
                  trr_f1,trr_b1,trz_f1,trz_b1,trr_av,tpp_av,
                  pf_f1,pf_b1,
                  dz,dr,dt,r_now,rho_av,rhof_av,D1_av,D2_av,Flagmat_vf_zero[k,j])

           #PML addition (vfr)
           vfr[k,j]=vfr[k,j]+
                (D2_av/2+(D1_av-rhof_av^2/rho_av)/dt)^(-1)*(
                                  -rhof_av/rho_av*Prz_T[LPML_z-k+1,j]
                                  )
        end

       #PML (vr) calculate elastic one then PML addition
       vr[k,j]=update_vr(vr_now,trr_f1,trr_b1,trz_f1,trz_b1,trr_av,tpp_av,
              0,dz,dr,dt,r_now,rho_av,0)

       vr[k,j]=vr[k,j]+dt/rho_av*Prz_T[LPML_z-k+1,j] #elastic PML
       dert_vfr=(vfr[k,j]-vfr_old)/dt
       vr[k,j]=vr[k,j]-dt*rhof_av/rho_av*dert_vfr #poroelastic PML
   end #if k==1

   r_now=(j-1)*dr #for vphi and vz (Mittet)

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
       vfz[k,j]=0.0
   else
       vfz_now=(D2_av/2+(D1_av-rhof_av^2/rho_av)/dt)^(-1) * (
                   (-D2_av/2+(D1_av-rhof_av^2/rho_av)/dt)*vfz_old-dzpf-rhof_av/rho_av*A
                   )

       #PML addition (vfz)
       vfz[k,j]=vfz_now+
                (D2_av/2+(D1_av-rhof_av^2/rho_av)/dt)^(-1)*(
                                  -PzzPE_T[LPML_z-k+1,j]
                                  -rhof_av/rho_av*Pzz_T[LPML_z-k+1,j]
                                  )
   end

   #PML (vz) calculate elastic one then vfz addition
   vz[k,j]=vz_now+dt/rho_av*A #elastic one with r=0
   vz[k,j]=vz[k,j]+dt/rho_av*Pzz_T[LPML_z-k+1,j] #elastic PML
   dert_vfz=(vfz[k,j]-vfz_old)/dt
   vz[k,j]=vz[k,j]-dt/rho_av*rhof_av*dert_vfz #poroelastic PML

end #k (z)

end #function

function ApplyBCLeft_velocity_1st_atPML_Bottom_Por!(vr,vz,trr,tpp,tzz,trz,vfr,vfz,pf,
    rhomat,rhofmat,D1mat,D2mat,Flagmat_vf_zero,nr,nz,dr,dz,dt,
    Prz_B,Pzz_B,PzzPE_B,
    LPML_z)
#Filling field values at j==1 (r) & k=2:nz-1
#j=1 corresponds to r=0 for variables of vz[], vphi[], tii[], tpz[]. Other values are @dr/2, not r=0
#2nd order FD (Graves, 1996)
#Supplement discussion from Mittet, 1996, Geophysics,61,1,21-33.
# m=0, monopole: vr(r)=-vr(-r),vphi=-vphi(-r),vz(r)=vz(-r), tij(r)=tij(-r)
# m=1, dipole: vr(r)=vr(-r),vphi=vphi(-r),vz(r)=-vz(-r), tij(r)=-tij(-r)
# nonzero@r=0 --> only "symmetric" components, otherwize 0.
#m=0 only!
#PML-only version

#Borehole axis (Left Edge): d/dr=0 (Symmetric BC)
j=1
#Variables@left-edge are vz[], vphi[], tii[], tpz[] Mittet
#Case m=0: vr=0, trz=0 (Randall, Mittet->fluid), trp=0 (Randall, Mittet->fluid)

#PML (Bottom)
for k=nz-LPML_z+1:nz

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

   vfr_old=vfr_now
   if(Flagmat_vf_zero[k,j]==1)
       vfr[k,j]=0.0
   else
       vr_dummy,vfr[k,j]=update_vr_1st_Por(vr_now,vfr_now,
              trr_f1,trr_b1,trz_f1,trz_b1,trr_av,tpp_av,
              pf_f1,pf_b1,
              dz,dr,dt,r_now,rho_av,rhof_av,D1_av,D2_av,Flagmat_vf_zero[k,j])

       #PML addition (vfr)
       vfr[k,j]=vfr[k,j]+
                (D2_av/2+(D1_av-rhof_av^2/rho_av)/dt)^(-1)*(
                                  -rhof_av/rho_av*Prz_B[k-nz+LPML_z,j]
                                  )
   end
   #PML (vr) calculate elastic one then vfr addition
   vr[k,j]=update_vr(vr_now,trr_f1,trr_b1,trz_f1,trz_b1,trr_av,tpp_av,
          0,dz,dr,dt,r_now,rho_av,0)
   vr[k,j]=vr[k,j]+dt/rho_av*Prz_B[k-nz+LPML_z,j] #elastic PML
   dert_vfr=(vfr[k,j]-vfr_old)/dt
   vr[k,j]=vr[k,j]-dt*rhof_av/rho_av*dert_vfr #poroelastic PML

   if(k!=nz)
       r_now=(j-1)*dr #for vphi and vz (Mittet)

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
           vfz[k,j]=0.0
       else
           vfz_now=(D2_av/2+(D1_av-rhof_av^2/rho_av)/dt)^(-1) * (
                       (-D2_av/2+(D1_av-rhof_av^2/rho_av)/dt)*vfz_old-dzpf-rhof_av/rho_av*A
                       )

           #PML addition (vfz)
           vfz[k,j]=vfz_now+
                    (D2_av/2+(D1_av-rhof_av^2/rho_av)/dt)^(-1)*(
                                      -PzzPE_B[k-nz+LPML_z,j]
                                      -rhof_av/rho_av*Pzz_B[k-nz+LPML_z,j]
                                      )
       end

       #PML (vz) calculate elastic one then vfz addition
       vz[k,j]=vz_now+dt/rho_av*A #elastic one with r=0
       vz[k,j]=vz[k,j]+dt/rho_av*Pzz_B[k-nz+LPML_z,j] #elastic PML
       dert_vfz=(vfz[k,j]-vfz_old)/dt
       vz[k,j]=vz[k,j]-dt/rho_av*rhof_av*dert_vfz #poroelastic PML
   end #if k==nz

end #k (z)

end #function


function ApplyBCLeft_stress_1st_atPML_Top_Por!(vr,vz,trr,tpp,tzz,trz,vfr,vfz,pf,
    Mmat,Cmat,Hmat,Gmat,nr,nz,dr,dz,dt,
    Rz_T,Srz_T,RzPE_T,
    LPML_z)
#Filling field values at j==1 (r) & k=2:nz-1
#j=1 corresponds to r=0 for variables of vz[], vphi[], tii[], tpz[]. Other values are @dr/2, not r=0
#2nd order FD (Graves, 1996)
#Supplement discussion from Mittet, 1996, Geophysics,61,1,21-33.
# m=0, monopole: vr(r)=-vr(-r),vphi=-vphi(-r),vz(r)=vz(-r), tij(r)=tij(-r)
# m=1, dipole: vr(r)=vr(-r),vphi=vphi(-r),vz(r)=-vz(-r), tij(r)=-tij(-r)
# nonzero@r=0 --> only "symmetric" components, otherwize 0.
#PML_only version


#Borehole axis (Left Edge): d/dr=0 (Symmetric BC)
j=1
#Variables@left-edge are vz[], vphi[], tii[], tpz[] Mittet
#Case m=0: vr=0, trz=0 (Randall, Mittet->fluid), trp=0 (Randall, Mittet->fluid)

 #PML(Top)
 for k=1:LPML_z #2nd order FD
    #I am alywas accesing [k,j], but definition of r_now changes with components!
    if(k!=1)
        r_now=(j-1)*dr #for trr,tpp,tzz,tpz Mittet

        trr_now=trr[k,j]
        vr_f1,vr_b1=vr[k,j],-vr[k,j] #symmtr
        vz_f1,vz_b1=vz[k,j],vz[k-1,j]

        vfr_f1,vfr_b1=vfr[k,j],-vfr[k,j] #symmtr
        vfz_f1,vfz_b1=vfz[k,j],vfz[k-1,j]

        H=Hmat[k,j] #be careful
        C=Cmat[k,j] #be careful
        G=Gmat[k,j] #be careful
        M=Mmat[k,j] #be careful

        #Special treatment @r=0
        drvr=(vr_f1-vr_b1)
        drvr=drvr/dr
        dzvz=(vz_f1-vz_b1)
        dzvz=dzvz/dz

        drvfr=(vfr_f1-vfr_b1)
        drvfr=drvfr/dr
        dzvfz=(vfz_f1-vfz_b1)
        dzvfz=dzvfz/dz

        #w/o PML then PML addition
        trr[k,j]=trr[k,j]+dt*(
                             +(H-2G)*(drvr+dzvz)
                             +H*drvr
                             +C*(2drvfr+dzvfz)
                             )
        trr[k,j]=trr[k,j]
                +dt*(H-2G)*Rz_T[LPML_z-k+1,j]
                +dt*C*RzPE_T[LPML_z-k+1,j]


        tpp[k,j]=tpp[k,j]+dt*(
                            +(H-2G)*(drvr+dzvz)
                            +H*drvr
                            +C*(2drvfr+dzvfz)
                            )
        tpp[k,j]=tpp[k,j]
                 +dt*(H-2G)*Rz_T[LPML_z-k+1,j]
                 +dt*C*RzPE_T[LPML_z-k+1,j]


        tzz[k,j]=tzz[k,j]+dt*(
                            +(H-2.0*G)*2.0*drvr
                            +H*dzvz
                            +C*(2.0*drvfr+dzvfz)
                            )
        tzz[k,j]=tzz[k,j]
                 +dt*H*Rz_T[LPML_z-k+1,j]
                 +dt*C*RzPE_T[LPML_z-k+1,j]


        pf[k,j]=pf[k,j]+dt*(
                            -C*(2.0*drvr+dzvz)
                            -M*(2.0*drvfr+dzvfz)
                            )
        pf[k,j]=pf[k,j]
                -dt*C*Rz_T[LPML_z-k+1,j]
                -dt*M*RzPE_T[LPML_z-k+1,j]
    end #if k==1

    #---trp and trz not on left edge!
    r_now=(j-1)*dr+dr/2.0 #for trp,trz Mittet
    trz_now=trz[k,j]
    if (Gmat[k,j]*Gmat[k,j+1]*Gmat[k+1,j]*Gmat[k+1,j+1]==0.0)
         G_av=0.0
    else
        gj1=2.0*Gmat[k,j]*Gmat[k+1,j]/(Gmat[k,j]+Gmat[k+1,j])
        gj2=2.0*Gmat[k,j+1]*Gmat[k+1,j+1]/(Gmat[k,j+1]+Gmat[k+1,j+1])
        G_av=2.0*gj1*gj2/(gj1+gj2)
    end

    vr_f1,vr_b1=vr[k+1,j],vr[k,j]
    vz_f1,vz_b1=vz[k,j+1],vz[k,j]
    trz[k,j]=update_trz(trz_now,vr_f1,vr_b1,vz_f1,vz_b1,dz,dr,dt,r_now,G_av,0)

    #PML addition
    trz[k,j]=trz[k,j]+dt*G_av*Srz_T[LPML_z-k+1,j]
    #experimental: when k==LPML_z, trz is outside PML

end #k (z)


end #function


function ApplyBCLeft_stress_1st_atPML_Bottom_Por!(vr,vz,trr,tpp,tzz,trz,vfr,vfz,pf,
    Mmat,Cmat,Hmat,Gmat,nr,nz,dr,dz,dt,
    Rz_B,Srz_B,RzPE_B,
    LPML_z)
#Filling field values at j==1 (r) & k=2:nz-1
#j=1 corresponds to r=0 for variables of vz[], vphi[], tii[], tpz[]. Other values are @dr/2, not r=0
#2nd order FD (Graves, 1996)
#Supplement discussion from Mittet, 1996, Geophysics,61,1,21-33.
# m=0, monopole: vr(r)=-vr(-r),vphi=-vphi(-r),vz(r)=vz(-r), tij(r)=tij(-r)
# m=1, dipole: vr(r)=vr(-r),vphi=vphi(-r),vz(r)=-vz(-r), tij(r)=-tij(-r)
# nonzero@r=0 --> only "symmetric" components, otherwize 0.
#PML_only version


#Borehole axis (Left Edge): d/dr=0 (Symmetric BC)
j=1
#Variables@left-edge are vz[], vphi[], tii[], tpz[] Mittet
#Case m=0: vr=0, trz=0 (Randall, Mittet->fluid), trp=0 (Randall, Mittet->fluid)

 #PML(Bottom)
 for k=nz-LPML_z+1:nz #2nd order FD
    #I am accesing [k,j], but the definition of r_now changes with each component!
    r_now=(j-1)*dr #for trr,tpp,tzz,tpz Mittet

    vr_f1,vr_b1=vr[k,j],-vr[k,j] #symmtr
    vz_f1,vz_b1=vz[k,j],vz[k-1,j]

    vfr_f1,vfr_b1=vfr[k,j],-vfr[k,j] #symmtr
    vfz_f1,vfz_b1=vfz[k,j],vfz[k-1,j]

    H=Hmat[k,j] #be careful
    C=Cmat[k,j] #be careful
    G=Gmat[k,j] #be careful
    M=Mmat[k,j] #be careful

    #Special treatment @r=0
    drvr=(vr_f1-vr_b1)
    drvr=drvr/dr
    dzvz=(vz_f1-vz_b1)
    dzvz=dzvz/dz

    drvfr=(vfr_f1-vfr_b1)
    drvfr=drvfr/dr
    dzvfz=(vfz_f1-vfz_b1)
    dzvfz=dzvfz/dz

    #w/o PML then PML addition
    trr[k,j]=trr[k,j]+dt*(
                         +(H-2G)*(drvr+dzvz)
                         +H*drvr
                         +C*(2drvfr+dzvfz)
                         )
    trr[k,j]=trr[k,j]
            +dt*(H-2G)*Rz_B[k-nz+LPML_z,j]
            +dt*C*RzPE_B[k-nz+LPML_z,j]

    tpp[k,j]=tpp[k,j]+dt*(
                        +(H-2G)*(drvr+dzvz)
                        +H*drvr
                        +C*(2drvfr+dzvfz)
                        )
    tpp[k,j]=tpp[k,j]
             +dt*(H-2G)*Rz_B[k-nz+LPML_z,j]
             +dt*C*RzPE_B[k-nz+LPML_z,j]

    tzz[k,j]=tzz[k,j]+dt*(
                        +(H-2.0*G)*2.0*drvr
                        +H*dzvz
                        +C*(2.0*drvfr+dzvfz)
                        )
    tzz[k,j]=tzz[k,j]
             +dt*H*Rz_B[k-nz+LPML_z,j]
             +dt*C*RzPE_B[k-nz+LPML_z,j]

    pf[k,j]=pf[k,j]+dt*(
                        -C*(2.0*drvr+dzvz)
                        -M*(2.0*drvfr+dzvfz)
                        )
    pf[k,j]=pf[k,j]
            -dt*C*Rz_B[k-nz+LPML_z,j]
            -dt*M*RzPE_B[k-nz+LPML_z,j]


    if(k!=nz)
        #---trp and trz not on left edge!
        r_now=(j-1)*dr+dr/2.0 #for trp,trz Mittet
        trz_now=trz[k,j]
        if (Gmat[k,j]*Gmat[k,j+1]*Gmat[k+1,j]*Gmat[k+1,j+1]==0.0)
             G_av=0.0
        else
            gj1=2.0*Gmat[k,j]*Gmat[k+1,j]/(Gmat[k,j]+Gmat[k+1,j])
            gj2=2.0*Gmat[k,j+1]*Gmat[k+1,j+1]/(Gmat[k,j+1]+Gmat[k+1,j+1])
            G_av=2.0*gj1*gj2/(gj1+gj2)
        end

        vr_f1,vr_b1=vr[k+1,j],vr[k,j]
        vz_f1,vz_b1=vz[k,j+1],vz[k,j]
        trz[k,j]=update_trz(trz_now,vr_f1,vr_b1,vz_f1,vz_b1,
            dz,dr,dt,r_now,G_av,0)
        #PML addition
        trz[k,j]=trz[k,j]+dt*G_av*Srz_B[k-nz+LPML_z,j]
        #experimental: when k==LPML_z, trz is outside PML
    end # if k==nz

end #k (z)


end #function


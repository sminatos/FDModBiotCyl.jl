#Functions in testing FDTD in cylindrical coordinate
#Randall et al (1991)


#01: start. to use, include(this file)
#02: minor update
#03: explit input of "m", Tuning multithreading
#03B: checking r_now againg!
#04: explicit input/output
#04B: better convention
#Mittet01: field variable locations are from Mittet (1996)
#Mittet02: rechecking dipole periodic boundary condition (not used)
#Mittet03: minor bug fix
#Mittet04: minor bug fix in init_receiver (index)
#05tmp: test implementing NPML of Wang and Tang (2003)
#05tmp2: update_vel and stress now assumes PML layer part skip. Plus, Left B.C. m=0, Solid case.
#05tmp3: Left B.C. m=0, trz and trp are not on edge but assigned to be zero -> test fixing...
#      : minor bugfix about symetric in Left B.C. stress/vel (m=1). experimental check..
#      : Left BCs are rechecked. Contradiction between Mittet and Randall on monipole trz.
#      : If Randall is true, then trz(r)=-trz(-r). For fluid no problem.
#      : This version relies on Randall about trz (unti)symetry with m=0.
#05tmp4: Right PML
#05tmp5: naming convnention modified to avoid misfeeding (memPQRS B and T)
#checking why bottom and top is non synmetric...
#     --> fixed. because of mismach of main for loop k range (L to nz-L should be L+1 to nz-L)
#tmp6: minor fixing dz/2 difference in PML_Wz etc...
#    : bugfix in PML_IWr. Look great
#tmp7: implementing TopRight -> temporary OK.
#tmp8: testing dr/2 difference in PML_Wr etc...
#05: Released. Monopole looks OK.
#06tmp1: minor bugfix (update_vel/stress range)
#---
#06tmp2: Right PML checking try to improve
#06tmp3: init_PML_profile separately provides dr,2dr,3dr.. array and dr/2+dr,dr/2+2dr,dr/2+3dr ... array
#06tmp4: Just modified bug in init_receiver (tii receiver index)

#06tmp5(under test on LVTS): including PML at Left Edge
#      :PML_update_memRS/PQ_Top/Bottom modified to include j==1
#      :PML_update_stress_XXX, input rhomat-->lmat,mmat

#06tmp6(under test on LVTS plane wave): PML-R memPQ/RS correct filling until Left Edge
#       : also for PML-TR...

# --------------
# major change (adding new funcions, keeping old functions)
# --------------
# Based from CylFDMod_Mittet06tmp6.jl
# --PcylFD -> Biot's poroelastic media.
# --Main -> Ou 2019, GJI, doi: 10.1093/gji/ggz144. Sub->Sidler, 2014, GJI, doi: 10.1093/gji/ggt447
# Parameters definitions are better in Sidler 2014
# --radially invariant source only (m=0)
# PCylFDMod_Ou01: start
# PCylFDMod_Ou02: several minor updates for easier maintenace.
# PCylFDMod_Ou03: update_vX_2nd_Por now includes Flag_vf_zero so as to accept Acoustic/Elastic modeling
#               : --> not compatible with ealier version
# PCylFDMod_Ou04: clean
# PCylFDMod_Ou1st01: 1st order FD only. For easier mentainance
# PCylFDMod_Ou1st02tmp: rechecking PML geometry and FD grids. testing. Definition of OMEGA function (init_PML),
#                     :  and update_memPQ, update_memRS
#   : memPQRS-->Top done. Right done. TopRight done. Bottom done. BottomRight done.
#   : update vel-stress PML --> Right done. Top done. TopRight done (this was important for edge reflection).
#   :                           Bottom done. BottomRight done (this was important for edge reflection).
#   :                       --> j=1 done.
# PCylFDMod_Ou1st02: checked with Guan's Fig3
# PCylFDMod_Ou1st03: So far, Elastic to Poroelastic horizontal boundary was OK,
#                  : but Poroelastic to Elastic horizontal bournadyr was not tested (vfz should be zero, but will not be zero at the current implementation)
#                  : --> dirty check  Flag_vf_zero(k,j)==0 && D1(k+1,j)==Inf, then Flag_vf_zero activated

#several bug-fix before module



function myricker2(tvec,f0,delay,derivative_number)
    # Gaussian=exp(-sigma*(t-t0)^2), sigma=-pi^2*f0^2. Mittet (eq 25), 1996, Geophysics
    nt=length(tvec)
    y=zeros(1,nt)
    if (derivative_number==2)
        println("Myricker: 2nd derivative of Gaussian (Ricker)")
    #---2nd derivative Gaussian (Ricker)
    for it=1:nt
     y[it] = -2*pi^2*f0^2*(1.0 - 2.0*(pi^2)*(f0^2)*((tvec[it]-delay).^2)) .* exp(-(pi^2)*(f0^2)*((tvec[it]-delay).^2))
    end
    elseif (derivative_number==1)
        println("Myricker: 1st derivative of Gaussian")
    #---1st derivative Gaussian
    for it=1:nt
     y[it] = -2*pi^2*f0^2*(tvec[it]-delay) .* exp(-(pi^2)*(f0^2)*((tvec[it]-delay).^2))
    end
    elseif (derivative_number==0)
        println("Myricker: Gaussian")
        #---Gaussian
    for it=1:nt
        y[it] = exp(-(pi^2)*(f0^2)*((tvec[it]-delay).^2))
    end
    elseif (derivative_number==3)
        println("Myricker: 3rd derivative of Gaussian")
    #---3rd derivative of Gaussian
    for it=1:nt
        y[it] = 4.0*pi^4*f0^4*(-2.0*pi^2*f0^2*(tvec[it]-delay)^3+3.0*(tvec[it]-delay))*
        exp(-(pi^2)*(f0^2)*((tvec[it]-delay).^2))
    end
    else
        println("Myricker: Error. Derivative_number should be 0, 1, 2 or 3.")
    end

    return y
end


function mywavelet_TsangRader(tvec,f0,Tc)
    #Tsang and Rader wavelet: Guan2011 eq 6.1
    # 1/2*(1+cos(2pi/Tc*(t-Tc/2)))*cos(2pi*f0*(t-Tc/2))
    # Gaussian=exp(-sigma*(t-t0)^2), sigma=-pi^2*f0^2. Mittet (eq 25), 1996, Geophysics
    nt=length(tvec)
    y=zeros(1,nt)
    println("Mywavelet: Tsang and Rader")
    #---2nd derivative Gaussian (Ricker)
    nt_end=ceil.(Int,round(Tc/dt)+1)

    for it=1:minimum([nt_end nt])
     y[it] = 1/2*(1+cos(2pi/Tc*(tvec[it]-Tc/2)))*cos(2pi*f0*(tvec[it]-Tc/2))
    end

    return y
end


#function to update vr
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
#                +(trr_av-tpp_av)/(r_now) #trr/r-tpp/r
#                +m*trp_now/(r_now) #m*trp/r
                +(___trr_av-___tpp_av+___m*___trp_now)/(___r_now) #trr/r-tpp/r+m*trp/r
                +(___trz_f-___trz_b)/___dz #dztrz
                )
end
#Poroelasticity
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

    #updating solid alternative (using drpf, instead of A)
    #avrg_vfr=(vfr_now+vfr_old)/2
    #vr_old=vr_now
    #vr_now=vr_old+dt/Rhof*(-drpf-D1*dert_vfr-D2*avrg_vfr)
    #

    #==crazy test
    vr_old=vr_now
    vr_now1=vr_old+dt/Rho*(A-Rhof*dert_vfr)
    avrg_vfr=(vfr_now+vfr_old)/2
    vr_now2=vr_old+dt/Rhof*(-drpf-D1*dert_vfr-D2*avrg_vfr)
    vr_now=(vr_now1+vr_now2)/2
    if(abs((vr_now1-vr_now2)/vr_now)*100>1)
        println("vr not consistent!")
        println(vr_now1,",",vr_now2)
    end
    ==#

    return vr_now,vfr_now
end


function update_vr_1st_Por_TEST(j,k,vr_now,vfr_now,
    trr_f,trr_b,trz_f,trz_b,trr_av,tpp_av,
    pf_f,pf_b,
    dz,dr,dt,r_now,Rho,Rhof,D1,D2,Flag_vf_zero)
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
    #vr_old=vr_now
    #vr_now=vr_old+dt/Rho*(A-Rhof*dert_vfr)

    #updating solid alternative (using drpf, instead of A)
    #avrg_vfr=(vfr_now+vfr_old)/2
    #vr_old=vr_now
    #vr_now=vr_old+dt/Rhof*(-drpf-D1*dert_vfr-D2*avrg_vfr)
    #

    #crazy test
    vr_old=vr_now
    vr_now1=vr_old+dt/Rho*(A-Rhof*dert_vfr)
    avrg_vfr=(vfr_now+vfr_old)/2
    vr_now2=vr_old+dt/Rhof*(-drpf-D1*dert_vfr-D2*avrg_vfr)
    vr_now=(vr_now1+vr_now2)/2
    if(abs((vr_now1-vr_now2)/vr_now)*100>1)
        println("vr not consistent!")
        println(vr_now1,",",vr_now2,"@j=",j,"k=",k)
    end


    return vr_now,vfr_now
end


#function to update vphi
function update_vphi(___vphi_now,___trp_f,___trp_b,___tpz_f,___tpz_b,___trp_av,___tpp_now,
    ___dz,___dr,___dt,___r_now,___Rho,___m)
    #first-order difference
    #X_now: current time, current position [k,j]
    #X_f: foward position (trp[k,j+1],tpz[k+1,j])
    #X_b: backward position (trp[k,j],tpz[k,j])
    #_av: averaged values (trp_av)
    #tpp is at the same location
    #r_now: current radial position @vphi
    #m: order (0: monopole, 1:dipole)
    return ___vphi_now+___dt/___Rho*(
              +(___trp_f-___trp_b)/___dr #drtrp
#              +2.0*(trp_av)/(r_now+dr/2) #2*trp/r
#              -m*tpp_now/(r_now+dr/2) #m*tpp/r
              +(2.0*(___trp_av)-___m*___tpp_now)/(___r_now) #2*trp/r-m*tpp/r
              +(___tpz_f-___tpz_b)/___dz #dztpz
              )
end

#function to update vphi
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
#              +(trz_av)/(r_now+dr/2) #trz/r
#              +m*tpz_now/(r_now+dr/2) #m*tpz/r
              +(___trz_av+___m*___tpz_now)/(___r_now) #trz/r+m*tpz/r
              +(___tzz_f-___tzz_b)/___dz #dztzz
              )
end
#Poroelasticity
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

    #updating solid alternative (using dzpf, instead of A)
#    avrg_vfz=(vfz_now+vfz_old)/2
#    vz_old=vz_now
#    vz_now=vz_old+dt/Rhof*(-dzpf-D1*dert_vfz-D2*avrg_vfz)
    #

    return vz_now,vfz_now
end

#-
#function to update trr
function update_trr(___trr_now,___vr_f,___vr_b,___vz_f,___vz_b,
    ___vr_av,___vphi_now,___dz,___dr,___dt,___r_now,___lambda,___mu,___m)
    #first-order difference
    #X_now: current time, current position [k,j]
    #X_f: foward position (vr[k,j+1],vz[k+1,j])
    #X_b: backward position (vr[k,j],vz[k,j])
    #_av: averaged values (vr_av)
    #vphi is at the same location
    #r_now: current radial position @trr
    #m: order (0: monopole, 1:dipole)
    return ___trr_now+___dt*(
              +(___lambda+2.0*___mu)*(___vr_f-___vr_b)/___dr #(l+2mu)*drvr
#              +lambda*(vr_av)/(r_now+dr/2) #l*vr/r
#              +lambda*m*vphi_now/(r_now+dr/2) #l*m*vphi/r
              +___lambda*(___vr_av+___m*___vphi_now)/(___r_now) #l*vr/r+l*m*vphi/r
              +___lambda*(___vz_f-___vz_b)/___dz #l*dvzdz
              )
end
#Poroelasic
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


#function to update tpp
function update_tpp(___tpp_now,___vr_f,___vr_b,___vz_f,___vz_b,___vr_av,
    ___vphi_now,___dz,___dr,___dt,___r_now,___lambda,___mu,___m)
    #first-order difference
    #X_now: current time, current position [k,j]
    #X_f: foward position (vr[k,j+1],vz[k+1,j])
    #X_b: backward position (vr[k,j],vz[k,j])
    #_av: averaged values (vr_av)
    #vphi is at the same location
    #r_now: current radial position @tpp
    #m: order (0: monopole, 1:dipole)
    return ___tpp_now+___dt*(
              +(___lambda+2.0*___mu)*(___vr_av+___m*___vphi_now)/(___r_now) #(l+2mu)*(vr+m*vphi)/r
              +___lambda*(___vr_f-___vr_b)/___dr #l*dvrdr
              +___lambda*(___vz_f-___vz_b)/___dz #l*dvzdz
              )
end
#Poroelasic
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


 #function to update tzz
function update_tzz(___tzz_now,___vr_f,___vr_b,___vz_f,___vz_b,
    ___vr_av,___vphi_now,___dz,___dr,___dt,___r_now,___lambda,___mu,___m)
     #first-order difference
     #X_now: current time, current position [k,j]
     #X_f: foward position (vr[k,j+1],vz[k+1,j])
     #X_b: backward position (vr[k,j],vz[k,j])
     #_av: averaged values (vr_av)
     #vphi is at the same location
     #r_now: current radial position @tzz
     #m: order (0: monopole, 1:dipole)
     return ___tzz_now+___dt*(
               +(___lambda+2.0*___mu)*(___vz_f-___vz_b)/___dz #(l+2mu)*dzvz
               +___lambda*(___vr_f-___vr_b)/___dr #l*drvr
               +___lambda*(___vr_av+___m*___vphi_now)/(___r_now) #l*(vr+m*vphi)/r
               )
end
#Poroelasic
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


#function to update trp
function update_trp(___trp_now,___vphi_f,___vphi_b,___vphi_av,___vr_now,
    ___dz,___dr,___dt,___r_now,___mu,___m)
    #first-order difference
    #X_now: current time, current position [k,j]
    #X_f: foward position (vphi[k,j])
    #X_b: backward position (vphi[k,j-1])
    #_av: averaged values (vphi_av)
    #vr is at the same location
    #r_now: current radial position @trp
    #m: order (0: monopole, 1:dipole)
    return ___trp_now+___dt*(
            +___mu*(___vphi_f-___vphi_b)/___dr #mu*drvphi
            -___mu*(___vphi_av+m*___vr_now)/(___r_now) #-mu*(vphi+m*vr)/r
            )
end

#function to update trz
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

#function to update tpz
function update_tpz(___tpz_now,___vphi_f,___vphi_b,___vz_now,___dz,___dt,
    ___r_now,___mu,___m)
    #first-order difference
    #X_now: current time, current position [k,j]
    #X_f: foward position (vphi[k,j])
    #X_b: backward position (vphi[k-1,j])
    #_av: averaged values ()
    #vz is at the same location
    #r_now: current radial position @tpz
    #m: order (0: monopole, 1:dipole)
    return ___tpz_now+___dt*(
            +___mu*(___vphi_f-___vphi_b)/___dz #mu*dzvphi
            -___mu*___m*___vz_now/(___r_now) #mu*m*vz/r
            )
end

#Poroelasticity
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

#Poroelasticity
function update_velocity_1st_Por!(vr,vz,trr,tpp,tzz,trz,vfr,vfz,pf,
    rhomat,rhofmat,D1mat,D2mat,Flagmat_vf_zero,nr,nz,dr,dz,dt,LPML_z,LPML_r)

#PE-to-E horizontal boundary incorpolated by checking if Flagmat_vf_zero[k,j]==0 && D1[k+1,j]==Inf

#2nd order FD (Graves, 1996)
#assuming PML!
#@inbounds Threads.@threads for j=3:nr-2
#      for k=3:nz-2

@inbounds Threads.@threads for j=2:nr-LPML_r
    for k=LPML_z+1:nz-LPML_z
          #I am alywas accesing [k,j], but definition of r_now changes with components!
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

#         if(j!=8)
#             vr[k,j],vfr[k,j]=update_vr_1st_Por_TEST(j,k,vr_now,vfr_now,
#                 trr_f1,trr_b1,trz_f1,trz_b1,trr_av,tpp_av,
#                 pf_f1,pf_b1,
#                 dz,dr,dt,r_now,rho_av,rhof_av,D1_av,D2_av,Flagmat_vf_zero[k,j])
#         end

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


#Poroelasticity
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
         #I am alywas accesing [k,j], but definition of r_now changes with components!
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

    #I am alywas accesing [k,j], but definition of r_now changes with components!
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
    #I am alywas accesing [k,j], but definition of r_now changes with components!
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



#---Hicks interpolation For Source to simualte band-limited dirac function---
# Hicks 2002 Geophysics
#src_index,src_dn=src_Hicks_fz(srcgeom)
#srcgeom=[z r]
function src_Hicks_fz(srcgeom)
    an_z=zeros(size(vz))
    an_r=zeros(size(vz))
    # if src is fz (-->vz): vz[1,1]=(0,dr/2)
    for iz=1:nz
       tmpz=(iz-1)*dz #because vz[1,1]=(0,dr/2)
       tmp_an=floor((tmpz-srcgeom[1,1])/dz)+mod((tmpz-srcgeom[1,1]),dz)
       an_z[iz,:]=tmp_an*ones(1,nr)
    end
    for ir=1:nr
       tmpr=(ir-1)*dr+dr/2 #because vz[1,1]=(0,dr/2)
       tmp_an=floor((tmpr-srcgeom[1,2])/dr)+mod((tmpr-srcgeom[1,2]),dr)
       an_r[:,ir]=tmp_an*ones(nz,1)
    end

    Hicks_para_r=4
    Hicks_para_b=4.14
    src_index=[]
    src_dn=[]
    for iz=1:nz
       for ir=1:nr
          if (abs(an_z[iz,ir])<=Hicks_para_r) && (abs(an_r[iz,ir])<=Hicks_para_r)
             src_index=[src_index; iz; ir]

             dn_z=besseli(0,Hicks_para_b*sqrt(1-(an_z[iz,ir]/Hicks_para_r)^2))/besseli(0,Hicks_para_b)
             dn_z=dn_z*sinc(an_z[iz,ir])
             dn_r=besseli(0,Hicks_para_b*sqrt(1-(an_r[iz,ir]/Hicks_para_r)^2))/besseli(0,Hicks_para_b)
             dn_r=dn_r*sinc(an_r[iz,ir])
             src_dn=[src_dn; dn_z*dn_r]

          end
       end
    end
    src_index=reshape(src_index,2,Int(length(src_index)/2)) #(:,isrc_index)->(iz,ir), value is src_dn(isrc_index)
    return src_index,src_dn #Inserting these into vz grid
end

function src_Hicks_fr(srcgeom)
    an_z=zeros(size(vz))
    an_r=zeros(size(vz))
    # if src is fr (-->vr): vr[1,1]=(dz/2,0)
    for iz=1:nz
       tmpz=(iz-1)*dz+dz/2 #because vr[1,1]=(dz/2,0)
       tmp_an=floor((tmpz-srcgeom[1,1])/dz)+mod((tmpz-srcgeom[1,1]),dz)
       an_z[iz,:]=tmp_an*ones(1,nr)
    end
    for ir=1:nr
       tmpr=(ir-1)*dr #because vr[1,1]=(dz/2,0)
       tmp_an=floor((tmpr-srcgeom[1,2])/dr)+mod((tmpr-srcgeom[1,2]),dr)
       an_r[:,ir]=tmp_an*ones(nz,1)
    end

    Hicks_para_r=4
    Hicks_para_b=4.14
    src_index=[]
    src_dn=[]
    for iz=1:nz
       for ir=1:nr
          if (abs(an_z[iz,ir])<=Hicks_para_r) && (abs(an_r[iz,ir])<=Hicks_para_r)
             src_index=[src_index; iz; ir]

             dn_z=besseli(0,Hicks_para_b*sqrt(1-(an_z[iz,ir]/Hicks_para_r)^2))/besseli(0,Hicks_para_b)
             dn_z=dn_z*sinc(an_z[iz,ir])
             dn_r=besseli(0,Hicks_para_b*sqrt(1-(an_r[iz,ir]/Hicks_para_r)^2))/besseli(0,Hicks_para_b)
             dn_r=dn_r*sinc(an_r[iz,ir])
             src_dn=[src_dn; dn_z*dn_r]

          end
       end
    end
    src_index=reshape(src_index,2,Int(length(src_index)/2)) #(:,isrc_index)->(iz,ir), value is src_dn(isrc_index)
    return src_index,src_dn #Inserting these into vr grid
end

function src_Hicks_tii(srcgeom) #NOT YET
    an_z=zeros(size(vz))
    an_r=zeros(size(vz))
    # if src is tzz/tpp/trr: tii[1,1]=(dz/2,dr/2)
    for iz=1:nz
       tmpz=(iz-1)*dz+dz/2 #because tii[1,1]=(dz/2,dr/2)
       tmp_an=floor((tmpz-srcgeom[1,1])/dz)+mod((tmpz-srcgeom[1,1]),dz)
       an_z[iz,:]=tmp_an*ones(1,nr)
    end
    for ir=1:nr
       tmpr=(ir-1)*dr+dz/2 #because tii[1,1]=(dz/2,dr/2)
       tmp_an=floor((tmpr-srcgeom[1,2])/dr)+mod((tmpr-srcgeom[1,2]),dr)
       an_r[:,ir]=tmp_an*ones(nz,1)
    end

    Hicks_para_r=4
    Hicks_para_b=4.14
    src_index=[]
    src_dn=[]
    for iz=1:nz
       for ir=1:nr
          if (abs(an_z[iz,ir])<=Hicks_para_r) && (abs(an_r[iz,ir])<=Hicks_para_r)
             src_index=[src_index; iz; ir]

             dn_z=besseli(0,Hicks_para_b*sqrt(1-(an_z[iz,ir]/Hicks_para_r)^2))/besseli(0,Hicks_para_b)
             dn_z=dn_z*sinc(an_z[iz,ir])
             dn_r=besseli(0,Hicks_para_b*sqrt(1-(an_r[iz,ir]/Hicks_para_r)^2))/besseli(0,Hicks_para_b)
             dn_r=dn_r*sinc(an_r[iz,ir])
             src_dn=[src_dn; dn_z*dn_r]

          end
       end
    end
    src_index=reshape(src_index,2,Int(length(src_index)/2)) #(:,isrc_index)->(iz,ir), value is src_dn(isrc_index)
    return src_index,src_dn #Inserting these into vr grid
end

#--Apply source (multi-node application of source function <--> call after src_Hicks_fz,fr)
function srcapply!(paramarray,src_index,src_dn,srcamp)
    #--e.g. fz src: vz=vz+fz (same form as vz update, thus fz is as-is and not time-derivatives)
    #--input array can be vx,vz etc
    for isrcind=1:length(src_dn)
       isz=src_index[1,isrcind]
       isr=src_index[2,isrcind]
       paramarray[isz,isr]=paramarray[isz,isr]+srcamp*src_dn[isrcind]
    end
end

#---Receiver utilities---
function getRecInterpIndex(recgeom,nrec)
#--Receiver interpolation indexes
#wis_allrec_vr,wis_allrec_vz,wis_allrec_vphi=getRecInterpIndex(recgeom,nrec)
#wis are indexes to evaluate bilinear interpolation
#see Doc of Interpolations.jl for more details
wis_allrec_vr=Array{Any}(undef,nrec)
wis_allrec_vz=Array{Any}(undef,nrec)
wis_allrec_vphi=Array{Any}(undef,nrec)
itp=interpolate(vr,BSpline(Linear())) #change here if you change interpolation.
    for irec=1:nrec
       #--index conversion vr: [1,1] <-> (dz/2,0)
       xtmp=((recgeom[irec,1]-dz/2.0)/dz+1,recgeom[irec,2]/dr+1)
       wis = Interpolations.weightedindexes((Interpolations.value_weights,), Interpolations.itpinfo(itp)..., xtmp)
       wis_allrec_vr[irec]=wis
       #--index conversion vz: [1,1] <-> (0,dr/2)
       xtmp=(recgeom[irec,1]/dz+1,(recgeom[irec,2]-dr/2.0)/dr+1)
       wis = Interpolations.weightedindexes((Interpolations.value_weights,), Interpolations.itpinfo(itp)..., xtmp)
       wis_allrec_vz[irec]=wis
       #--index conversion vphi: [1,1] <-> (dz/2,dr/2)
       xtmp=((recgeom[irec,1]-dz/2.0)/dz+1,(recgeom[irec,2]-dr/2.0)/dr+1)
       wis = Interpolations.weightedindexes((Interpolations.value_weights,), Interpolations.itpinfo(itp)..., xtmp)
       wis_allrec_vphi[irec]=wis
    end
    return wis_allrec_vr,wis_allrec_vz,wis_allrec_vphi
end

#-Receiver field extraction/interpolation
function getRecData!(rec_vr,rec_vz,rec_vphi,wis_allrec_vr,wis_allrec_vz,wis_allrec_vphi,ii,nrec)
for irec=1:nrec
       rec_vr[ii,irec]=vr[wis_allrec_vr[irec]...]
       rec_vz[ii,irec]=vz[wis_allrec_vz[irec]...]
       rec_vphi[ii,irec]=vphi[wis_allrec_vphi[irec]...]
    end
end

function getRecInterpIndex_stress(recgeom,nrec)
#--Receiver interpolation indexes
#wis_allrec_vr,wis_allrec_vz,wis_allrec_vphi=getRecInterpIndex(recgeom,nrec)
#wis are indexes to evaluate bilinear interpolation
#see Doc of Interpolations.jl for more details
wis_allrec_trr=Array{Any}(undef,nrec)
itp=interpolate(vr,BSpline(Linear())) #change here if you change interpolation.
    for irec=1:nrec
       #--index conversion trr: [1,1] <-> (dz/2,dr/2)
       xtmp=((recgeom[irec,1]-dz/2.0)/dz+1,(recgeom[irec,2]-dr/2.0)/dr+1)
       wis = Interpolations.weightedindexes((Interpolations.value_weights,), Interpolations.itpinfo(itp)..., xtmp)
       wis_allrec_trr[irec]=wis
    end
    return wis_allrec_trr
end

#-Receiver field extraction/interpolation
function getRecData_stress!(rec_trr,wis_allrec_trr,ii,nrec)
    for irec=1:nrec
       rec_trr[ii,irec]=trr[wis_allrec_trr[irec]...]
    end
end

#---------------

#--getshots
function get_snapshots!(snapshots_vr,snapshots_vz,isnap,vrmat,vzmat,nr,nz)
@inbounds Threads.@threads for ir=1:nr
        for iz=1:nz
            snapshots_vr[iz,ir,isnap]=vrmat[iz,ir]
            snapshots_vz[iz,ir,isnap]=vzmat[iz,ir]
        end
    end
end

function get_snapshots_t!(snapshots_trr,isnap,trrmat,nr,nz)
@inbounds Threads.@threads for ir=1:nr
        for iz=1:nz
            snapshots_trr[iz,ir,isnap]=trrmat[iz,ir]
        end
    end
end

#--Snapshots (sparse version)
function get_snapshots_sp!(snapshots_vr,snapshots_vz,isnap,vrmat,vzmat,nr,nz)

ir_sp=0
@inbounds for ir=1:10:nr
        ir_sp=ir_sp+1
        iz_sp=0
        for iz=1:10:nz
            iz_sp=iz_sp+1
            snapshots_vr[iz_sp,ir_sp,isnap]=vrmat[iz,ir]
            snapshots_vz[iz_sp,ir_sp,isnap]=vzmat[iz,ir]
        end
    end
end

function get_snapshots_sp_t!(snapshots_trr,isnap,trrmat,nr,nz)
ir_sp=0
@inbounds for ir=1:10:nr
        ir_sp=ir_sp+1
        iz_sp=0
        for iz=1:10:nz
            iz_sp=iz_sp+1
            snapshots_trr[iz_sp,ir_sp,isnap]=trrmat[iz,ir]
        end
    end
end

#--Source function from Randall (1989)--
#Psi_m of equation (B-6)
function Randallsrc(volume_injection_func,dr,dz,dt,r0,m,Vf,Rhof)
    #volume_injection_func is e.g., Ricker wavelet
    #r0 is a source radius and very small
    lambda=Rhof*Vf^2
    if(m==0)
        em=1. #Neumann's factor (See Kurjkian and Chen, 1986)
    else
        em=2.
    end
    dvdt=(volume_injection_func[2:end]-volume_injection_func[1:end-1])/dt
    dvdt=[0;dvdt]

    return psi_m=-lambda/(2pi)*em*dvdt*1.0/(0.5*dr^2*dz)*(r0/dr)^m
end

function cumtrapz(X::T, Y::T) where {T <: AbstractVector}
  #--Function Y=F(X)
  # Check matching vector length
  @assert length(X) == length(Y)
  # Initialize Output
  out = similar(X)
  out[1] = 0
  # Iterate over arrays
  for i in 2:length(X)
    out[i] = out[i-1] + 0.5*(X[i] - X[i-1])*(Y[i] + Y[i-1])
  end
  # Return output
  out
end

#---Blakman-Harris wavelet to reproduce Randall
# requires DSP.jl
function myBHwave01(dt,nt,delay,fc,Wn)
    #Wn=(freq1,freq2) for Bandpass filter
    #fc: center freq (defines wavelet length)
    #--Blackman-Harris window
    #https://nl.mathworks.com/help/signal/ref/blackmanharris.html
    N=Int(ceil(1.5/fc/dt)+1)

    Blackman=zeros(N,1)
    Blackman4_coef=[0.35875 0.48829 0.14128 0.01168]
    for ii=1:N
       Blackman[ii]=Blackman4_coef[1]-Blackman4_coef[2]*cos(2pi*(ii-1)/(N-1))
          +Blackman4_coef[3]*cos(4pi*(ii-1)/(N-1))-Blackman4_coef[4]*cos(6pi*(ii-1)/(N-1))
    end

    #--2nd derivative of BH---
    BH=zeros(nt,1)
    BH[1:N]=Blackman[:]
    BH[N+1:nt]=ones(nt-N,1)*Blackman[end]

    BH2=zeros(nt,1)
    for ii=2:N-1
       BH2[ii]=BH[ii-1]-2*BH[ii]+BH[ii+1]
    end

    display(plot(BH2,title="2nd derivative of Blackman-Harris window"))

    #-Filtering--
    #responsetype = Bandpass(Wn[1],Wn[2]; fs=1/dt)
    responsetype = Lowpass(Wn[2]; fs=1/dt)
    designmethod = FIRWindow(hanning(nt))
    filter1=digitalfilter(responsetype, designmethod)
    fs=1/dt
    delay_fir=(nt-1)/(2*fs)

    BH2_filt=filt(filter1, BH2)

    display(plot(BH2_filt,title="filtered"))

    BH2_filt_f=fft(BH2_filt)
    T=(nt-1)*dt
    fvec=[0:1/T:(nt-1)/T]
    tmp_abs=map(abs,BH2_filt_f)
    tmp_abs_dB=10*map(log10,tmp_abs/maximum(tmp_abs))
    #display(plot(fvec,tmp_abs/maximum(tmp_abs[:]),yaxis=:log,xlims=(0,1E5)))
    display(plot(fvec,tmp_abs_dB,xlims=(0,5E3),title="filtered"))

    #--compensating linear delay (when FIR)
    isdelay_fir=Int(ceil(delay_fir/dt))+1
    #--including input delay
    isdelay0=Int(ceil(delay/dt))

    isdelay_total=isdelay_fir-isdelay0
    BH2_filt3=0.0*BH2_filt
    BH2_filt3[1:nt-isdelay_total+1]=BH2_filt[isdelay_total:nt]

    #--additional smoothing (to be robust to delay compensation later)
    T_window=dt*(nt-1)
    isTw=Int(ceil(T_window/dt)+1)
    ww=10
    tmpp=range(0,stop=pi/2,length=ww)
    for is=1:ww
       BH2_filt3[is]=BH2_filt3[is]*(sin(tmpp[is]))^4
       BH2_filt3[isTw-is+1]=BH2_filt3[isTw-is+1]*(sin(tmpp[is]))^4
    end
    BH2_filt3[isTw+1:nt]=BH2_filt3[isTw+1:nt]*0.0


    tvec=[0:dt:(nt-1)*dt]
    display(plot(tvec,BH2_filt3,xlims=(0,2E-3),title="filtered(delayed)"))

    println("HEY",BH2_filt3[1]/maximum(map(abs,BH2_filt3[:])))
    #--Nomralize--
    BH2_filt3[:]/maximum(map(abs,BH2_filt3[:]))
end


#---
function init_snapshots_v(nskip,nt,nz,nr)
    #---Snapshots
    #nskip=100
    itvec_snap=zeros(Float64,Int(ceil(nt/nskip))-1)
    tvec_snap=zeros(Float64,Int(ceil(nt/nskip))-1)

    cnt=1
    for itmp=1:Int(ceil(nt/nskip))-1
       cnt+=nskip
       itvec_snap[itmp]=cnt
    end

    if (maximum(itvec_snap)>nt)
       println("Error in snapshots settings!")
       return
    else
       nsnap=length(itvec_snap)
    end

    snapshots_vr=zeros(nz,nr,nsnap)
    snapshots_vz=zeros(nz,nr,nsnap)
    return snapshots_vr,snapshots_vz,nsnap,itvec_snap
end

function init_snapshots_t(nskip,nt,nz,nr)
    #---Snapshots
    #nskip=100
    itvec_snap=zeros(Float64,Int(ceil(nt/nskip))-1)
    tvec_snap=zeros(Float64,Int(ceil(nt/nskip))-1)

    cnt=1
    for itmp=1:Int(ceil(nt/nskip))-1
       cnt+=nskip
       itvec_snap[itmp]=cnt
    end

    if (maximum(itvec_snap)>nt)
       println("Error in snapshots settings!")
       return
    else
       nsnap=length(itvec_snap)
    end

    snapshots_trr=zeros(nz,nr,nsnap)
    return snapshots_trr,nsnap,itvec_snap
end

function init_fields(nz,nr)
    vr=zeros(Float64,nz,nr)
    vphi=zeros(Float64,nz,nr)
    vz=zeros(Float64,nz,nr)
    trr=zeros(Float64,nz,nr)
    tpp=zeros(Float64,nz,nr)
    tzz=zeros(Float64,nz,nr)
    trp=zeros(Float64,nz,nr)
    trz=zeros(Float64,nz,nr)
    tpz=zeros(Float64,nz,nr)
    return vr,vphi,vz,trr,tpp,tzz,trp,trz,tpz
end

#Poroelasticity (m=0 only)
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


#--Source functions--
function makesrc_multipole01(f0,delay,dz,dt,nt,T,r0,m)

    #calculating stress values for multipole source (Randall 1991)
    tvec=range(0.0,T,length=nt) # range object (no memory allocation)
    tvec=collect(tvec) # a vector

    #--Dipole/Monopole Source on the borehole axis (fluid)
    #-first define volume change function (consistent to Kurjkian), then calculate quantities for Randall-FD
    Vf=1500.0
    Rhof=1000.0
    r0=1E-6 #monopole sources on the circle of this radius to construct dipole source
    #--first define volume change function v0(t)
    tmp=myricker(tvec,f0,delay)
    #tmp=myBHwave01(dt,nt,3/f0,f0,(0.25E3,3E3))
    volume_injection_func=tmp
    #--"src_func" will be added to tpp,tzz,trr at source location (0+dr/2 in the radial direction)
    src_func=Randallsrc(volume_injection_func,dr,dz,dt,r0,m,Vf,Rhof)

    return src_func,volume_injection_func,tvec
end
function makesrc_bodyforce01(f0,delay,dt,nt,T,rho_at_src)
    #calculating source term for body force (added to velocities)
    tvec=range(0.0,T,length=nt) # range object (no memory allocation)
    tvec=collect(tvec) # a vector

#    src_func=myricker2(tvec,f0,delay,0) #when using Gaussian
    src_func=myricker2(tvec,f0,delay,1) #when using 1st derivative Gaussian
#    src_func=myricker2(tvec,f0,delay,2) #when using 2nd derivative Gaussian (Ricker)
#    src_func=myricker2(tvec,f0,delay,3) #when using 3rd derivative Gaussian (Ricker)

    src_func=src_func/rho_at_src*dt # this directly adds to velocities

    #---tmporary--
    #tmp=src_func[2:end]-src_func[1:end-1]
    #tmp=[0;tmp]
    #src_func=tmp
    #-------------
    #---tmporary--
    #tmp=cumtrapz(src_func[:],tvec[:])
    #src_func=tmp
    #-------------

    return src_func,tvec
end
function makesrc_dwi01(f0,delay,dt,nt,T,rho_at_src)
    #calculating source term for comparison with Tubman's DWI
    tvec=range(0.0,T,length=nt) # range object (no memory allocation)
    tvec=collect(tvec) # a vector

    # DWI is tested using Ricker (Gauss 2nd derivative)
    # => gives 2nd derivative Ricker (4th derivative Gauss)
    # To have the same wavelet in pressure,
    # FD may require 1st derivative Ricker (3rd derivative Gauss) as input source, added to d/dt(trr)

    src_func=myricker2(tvec,f0,delay,1) #when using 1st derivative Gaussian

#    src_func=myricker2(tvec,f0,delay,3) #when using 3rd derivative Gaussian
    src_func=-src_func

    return src_func,tvec
end
#--Source geometry--
function get_srcindex(srcgeom,dr,dz)
    #src index and amplitudes

    #--For body force, band-limited delta function from Hicks formula
    #src_index,src_dn=src_Hicks_fz(srcgeom) #--fz src---

    #--For monopole, band-limited delta function from Hicks formula
    #src_index,src_dn=src_Hicks_tii(srcgeom) #--tii src---

    #--Manual input
    izsrc_ref=Int(ceil(srcgeom[1,1]/dz)+1)
#    irsrc_ref=Int(ceil(srcgeom[1,2]/dr)+1)
#    src_index=[izsrc_ref;irsrc_ref]
#    src_dn=[1.]

    #-1 point on axis--
#    src_index=[izsrc_ref;1]
#    src_dn=[1.]
    #-2 points lateral--
    src_index=[izsrc_ref izsrc_ref;1 2]
    src_dn=[1. 1.]/2.
    #-3 points lateral--
#    src_index=[izsrc_ref izsrc_ref izsrc_ref;1 2 3]
#    src_dn=[1. 1. 1.]/3.

    #-1 point non on axis-
#    src_index=[izsrc_ref;2]
#    src_dn=[1.]

    #-4 square points
#    src_index=[izsrc_ref izsrc_ref izsrc_ref+1 izsrc_ref+1;1 2 1 2]
#    src_dn=[1. 1. 1. 1.]/4.

    #--For dipole source, try NxN grid
    #izsrc_ref=Int(ceil(srcgeom[1,1]/dz)+1)
    #Ngrid=3 # odd number!
    #
    #src_index=zeros(2,Ngrid^2)
    #for itmp=1:Ngrid
    #   src_index[1,1+Ngrid*(itmp-1):Ngrid*itmp]=ones(1,Ngrid)*(izsrc_ref-(Ngrid-1)/2+itmp-1)
    #   src_index[2,1+Ngrid*(itmp-1):Ngrid*itmp]=transpose(collect(range(1,Ngrid,length=Ngrid)))
    #end
    #src_index=ceil.(Int,src_index)
    #src_dn=ones(1,Ngrid^2)/Ngrid^2

    #src_index=[400 400 400 401 401 401 402 402 402;1 2 3 1 2 3 1 2 3]
    #src_dn=[1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0];

    return src_index,src_dn
end

function get_srcindex_monopole(srcgeom,dr,dz)
    #src index and amplitudes
    #--point input
    izsrc_ref=Int(ceil(srcgeom[1,1]/dz)+1)
    irsrc_ref=Int(ceil(srcgeom[1,2]/dr)+1)
    src_index=[izsrc_ref;irsrc_ref]
    src_dn=[1.]
    return src_index,src_dn
end
function get_srcindex_dipole(srcgeom,dr,dz)
    #src index and amplitudes
    #dipole assume on-axis source
    #--point input
    izsrc_ref=Int(ceil(srcgeom[1,1]/dz)+1)
    src_index=[izsrc_ref;1]
    src_dn=[1.]
    return src_index,src_dn
end

#--Receiver geometry--
function init_receiver_interp(recgeom,nrec,dr,dz,nr,nz,nt)
    #using Interpolations.jl to off-grids receivers

    if (minimum(recgeom[:,1])<0.) || (maximum(recgeom[:,1])>(nz-1)*dz)
       error("Check receiver geometry!")
    end
    if (minimum(recgeom[:,2])<0.) || (maximum(recgeom[:,2])>(nr-1)*dr)
       error("Check receiver geometry!")
    end

    #--Receiver interpolation indexes---
    wis_allrec_vr,wis_allrec_vz,wis_allrec_vphi=getRecInterpIndex(recgeom,nrec)
    wis_allrec_tzz=getRecInterpIndex_stress(recgeom,nrec)
    #---Initializing receiver data matrces
    rec_vr=zeros(nt,nrec)
    rec_vz=zeros(nt,nrec)
    rec_tii=zeros(nt,nrec)

    return rec_vr,rec_vz,rec_tii,wis_allrec_vr,wis_allrec_vz,wis_allrec_tii
end

function init_receiver(recgeom,nrec,dr,dz,nr,nz,nt)
    #receivers are assigned at nearest grids

    if (minimum(recgeom[:,1])<0.) || (maximum(recgeom[:,1])>(nz-1)*dz)
       error("Check receiver geometry!")
    end
    if (minimum(recgeom[:,2])<0.) || (maximum(recgeom[:,2])>(nr-1)*dr)
       error("Check receiver geometry!")
    end

    #--Receiver indexes---
    index_allrec_vr=zeros(nrec,2) #z,r
    index_allrec_vz=zeros(nrec,2) #z,r
    index_allrec_tii=zeros(nrec,2) #z,r
    for irec=1:nrec
        tmpz=recgeom[irec,1]
        tmpr=recgeom[irec,2]

        iz=round((tmpz-dz/2)/dz)+1 #vz
        ir=round(tmpr/dr)+1 #vz
        index_allrec_vz[irec,1]=iz
        index_allrec_vz[irec,2]=ir

        iz=round(tmpz/dz)+1 #vr
        ir=round((tmpr-dr/2)/dr)+1 #vr
        index_allrec_vr[irec,1]=iz
        index_allrec_vr[irec,2]=ir

#        iz=round((tmpz-dz/2)/dz)+1 #tii
#        ir=round((tmpr-dr/2)/dr)+1 #tii
        iz=round(tmpz/dz)+1 #tii
        ir=round(tmpr/dr)+1 #tii
        index_allrec_tii[irec,1]=iz
        index_allrec_tii[irec,2]=ir

    end
    #---Initializing receiver data matrces
    rec_vr=zeros(nt,nrec)
    rec_vz=zeros(nt,nrec)
    rec_tii=zeros(nt,nrec)

    return rec_vr,rec_vz,rec_tii,ceil.(Int,index_allrec_vr),ceil.(Int,index_allrec_vz),ceil.(Int,index_allrec_tii)
end

#---
function getRecData_from_index!(paramarray,rec_data,rec_index,nrec,it)
    #receivers are assigned at nearest grids

    for irec=1:nrec
        iz=rec_index[irec,1]
        ir=rec_index[irec,2]
        rec_data[it,irec]=paramarray[iz,ir]
    end
end


function check_stability01(dt,dr,Vp,m)
#From Randall 1991, cdt/dr<=1/sqrt(2)*(1+m^2dr^2/(8*rj^2))^(-1/2)
#dr=dz, rj=j*dr
    #Rough version : V*dt/dr<=1/sqrt(2), V-->maximum Vp
    Vp_max=maximum(Vp)
    if (Vp_max*dt/dr <= 1/sqrt(2))
        println("Stability OK.")
    else
        error("########Stability NG.######")
    end

end




#---PMLs---
function init_PML_profile(LPML_r,LPML_z,Vmax,dr,dz,nr)
#function init_PML_profile(LPML_r,LPML_z,Vmax)
 #creating PML profile function
 #LPML_ : PML thickness (points)
 #Vmax: maximum velocity

 #Try to look at PWL_Wz: when slope is very sharp reflection can occur.
 #                     : The highest point is controled by Ralpha,
 #                     : slope continuing from non-PML to PML is controlled by a, b
 #                     : large Ralpha -> too sharp. large a -> continuation may be discontinous

 #---Wang and Tang 2003 eqs (42-43)
# a=0.25
# b=0.75
 a=0.0
 b=1.0
 c=0.0
 Ralpha=1E-8 #required reflection magnitude order: ~log()=3pi(m+1)? m is order of OMEGA function
# Ralpha=1E-3 #required reflection magnitude order (1E-3 was better than 1E-6)
# Ralpha=0.01 #this was even better, but FSref can be large
# Ralpha=1E-4 #Wang and Tang (good with a=0.25,b=0.75, L=lambda/2)

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


 #--integral PML_IWr(r)=1/r*integral(0--r)[PML_Wr] dr
# PML_IWr=0.0*PML_Wr
 #PML_IWr[1]=PML_Wr[1]/2.0*dr
 #for ir=2:LPML_r
#     PML_IWr[ir]=PML_IWr[ir-1]+(PML_Wr[ir-1]+PML_Wr[ir])/2.0*dr
 #end
 #PML_IWr=PML_IWr./(dvec_r+((nr-LPML_r)*dr+dr/2)*ones(length(dvec_r)))

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



 return  LPML_r,LPML_z,PML_Wr,PML_Wz,PML_IWr,PML_Wr2,PML_Wz2,PML_IWr2
end

function PML_check(LPML_r,LPML_z,Vmax,dr,dz,f0)
     #Check if LPML < LAMBDA=VMAX/f0/2
     #LPML_ : PML thickness (points)
     #Vmax: maximum velocity

     LPML_r_meter=LPML_r*dr
     LPML_z_meter=LPML_z*dz

     Max_wavelength=Vmax/f0

     if (LPML_r_meter > Max_wavelength/2)
         println("PML thickness OK in R.")
     else
         tmp=round(Max_wavelength/2/dr+1)
         println("==Suggested LPML_r=",tmp)
         error("########INCREASE PML THICKNESS in R ######")
     end

     if (LPML_z_meter > Max_wavelength/2)
         println("PML thickness OK in Z.")
     else
         tmp=round(Max_wavelength/2/dz+1)
         println("==Suggested LPML_z=",tmp)
         error("########INCREASE PML THICKNESS in Z ######")
     end

end

function init_PML(LPML_r,LPML_z,nr,nz)
 #initializing PML variables (18 additional variables)
 #LPML_ : PML thickness (points)

 #Top and Bottom (z streching only)
 Prz_T=zeros(LPML_z,nr-LPML_r)
 Ppz_T=zeros(LPML_z,nr-LPML_r)
 Pzz_T=zeros(LPML_z,nr-LPML_r)
 Rz_T=zeros(LPML_z,nr-LPML_r)
 Srz_T=zeros(LPML_z,nr-LPML_r)
 Rpz_T=zeros(LPML_z,nr-LPML_r)

 Prz_B=zeros(LPML_z,nr-LPML_r)
 Ppz_B=zeros(LPML_z,nr-LPML_r)
 Pzz_B=zeros(LPML_z,nr-LPML_r)
 Rz_B=zeros(LPML_z,nr-LPML_r)
 Srz_B=zeros(LPML_z,nr-LPML_r)
 Rpz_B=zeros(LPML_z,nr-LPML_r)

 #Right (r streching only)
 Prr_R=zeros(nz-2*LPML_z,LPML_r)
 Qrp_R=zeros(nz-2*LPML_z,LPML_r)
 Ppr_R=zeros(nz-2*LPML_z,LPML_r)
 Qpp_R=zeros(nz-2*LPML_z,LPML_r)
 Pzr_R=zeros(nz-2*LPML_z,LPML_r)
 Qzp_R=zeros(nz-2*LPML_z,LPML_r)
 Rr_R=zeros(nz-2*LPML_z,LPML_r)
 Rp_R=zeros(nz-2*LPML_z,LPML_r)
 Rrp_R=zeros(nz-2*LPML_z,LPML_r)
 Srp_R=zeros(nz-2*LPML_z,LPML_r)
 Rrz_R=zeros(nz-2*LPML_z,LPML_r)
 Spz_R=zeros(nz-2*LPML_z,LPML_r)

 #TopRight and BottomRight (both z and r streching)
 Prr_TR=zeros(LPML_z,LPML_r)
 Qrp_TR=zeros(LPML_z,LPML_r)
 Prz_TR=zeros(LPML_z,LPML_r)
 Ppr_TR=zeros(LPML_z,LPML_r)
 Qpp_TR=zeros(LPML_z,LPML_r)
 Ppz_TR=zeros(LPML_z,LPML_r)
 Pzr_TR=zeros(LPML_z,LPML_r)
 Qzp_TR=zeros(LPML_z,LPML_r)
 Pzz_TR=zeros(LPML_z,LPML_r)
 Rr_TR=zeros(LPML_z,LPML_r)
 Rp_TR=zeros(LPML_z,LPML_r)
 Rz_TR=zeros(LPML_z,LPML_r)
 Rrp_TR=zeros(LPML_z,LPML_r)
 Srp_TR=zeros(LPML_z,LPML_r)
 Rrz_TR=zeros(LPML_z,LPML_r)
 Srz_TR=zeros(LPML_z,LPML_r)
 Rpz_TR=zeros(LPML_z,LPML_r)
 Spz_TR=zeros(LPML_z,LPML_r)

 Prr_BR=zeros(LPML_z,LPML_r)
 Qrp_BR=zeros(LPML_z,LPML_r)
 Prz_BR=zeros(LPML_z,LPML_r)
 Ppr_BR=zeros(LPML_z,LPML_r)
 Qpp_BR=zeros(LPML_z,LPML_r)
 Ppz_BR=zeros(LPML_z,LPML_r)
 Pzr_BR=zeros(LPML_z,LPML_r)
 Qzp_BR=zeros(LPML_z,LPML_r)
 Pzz_BR=zeros(LPML_z,LPML_r)
 Rr_BR=zeros(LPML_z,LPML_r)
 Rp_BR=zeros(LPML_z,LPML_r)
 Rz_BR=zeros(LPML_z,LPML_r)
 Rrp_BR=zeros(LPML_z,LPML_r)
 Srp_BR=zeros(LPML_z,LPML_r)
 Rrz_BR=zeros(LPML_z,LPML_r)
 Srz_BR=zeros(LPML_z,LPML_r)
 Rpz_BR=zeros(LPML_z,LPML_r)
 Spz_BR=zeros(LPML_z,LPML_r)

#----
return Prz_T,Ppz_T,Pzz_T,Rz_T,Srz_T,Rpz_T,
Prz_B,Ppz_B,Pzz_B,Rz_B,Srz_B,Rpz_B,
Prr_R,Qrp_R,Ppr_R,Qpp_R,Pzr_R,Qzp_R,Rr_R,Rp_R,Rrp_R,Srp_R,Rrz_R,Spz_R,
Prr_TR,Qrp_TR,Prz_TR,Ppr_TR,Qpp_TR,Ppz_TR,Pzr_TR,Qzp_TR,Pzz_TR,Rr_TR,Rp_TR,Rz_TR,Rrp_TR,Srp_TR,Rrz_TR,Srz_TR,Rpz_TR,Spz_TR,
Prr_BR,Qrp_BR,Prz_BR,Ppr_BR,Qpp_BR,Ppz_BR,Pzr_BR,Qzp_BR,Pzz_BR,Rr_BR,Rp_BR,Rz_BR,Rrp_BR,Srp_BR,Rrz_BR,Srz_BR,Rpz_BR,Spz_BR

end

#Poroelasticity
function init_PML_Por(LPML_r,LPML_z,nr,nz)
 #initializing PML variables (16 additional variables): m=0 only
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


function init_memory_variables(LPML_r,LPML_z,nr,nz)
    #--memory variables to calculate PML
    #for 2nd order FD, storing +2 sample layer larger than PML layer
    #Top and Bottom
    memT_vr=zeros(LPML_z+2,nr)
    memT_vphi=zeros(LPML_z+2,nr)
    memT_vz=zeros(LPML_z+2,nr)
    memT_trr=zeros(LPML_z+2,nr)
    memT_tpp=zeros(LPML_z+2,nr)
    memT_tzz=zeros(LPML_z+2,nr)
    memT_trp=zeros(LPML_z+2,nr)
    memT_trz=zeros(LPML_z+2,nr)
    memT_tpz=zeros(LPML_z+2,nr)

    memB_vr=zeros(LPML_z+2,nr)
    memB_vphi=zeros(LPML_z+2,nr)
    memB_vz=zeros(LPML_z+2,nr)
    memB_trr=zeros(LPML_z+2,nr)
    memB_tpp=zeros(LPML_z+2,nr)
    memB_tzz=zeros(LPML_z+2,nr)
    memB_trp=zeros(LPML_z+2,nr)
    memB_trz=zeros(LPML_z+2,nr)
    memB_tpz=zeros(LPML_z+2,nr)

    #Right
    memR_vr=zeros(nz,LPML_r+2)
    memR_vphi=zeros(nz,LPML_r+2)
    memR_vz=zeros(nz,LPML_r+2)
    memR_trr=zeros(nz,LPML_r+2)
    memR_tpp=zeros(nz,LPML_r+2)
    memR_tzz=zeros(nz,LPML_r+2)
    memR_trp=zeros(nz,LPML_r+2)
    memR_trz=zeros(nz,LPML_r+2)
    memR_tpz=zeros(nz,LPML_r+2)

return    memT_vr,memT_vphi,memT_vz,memT_trr,memT_tpp,memT_tzz,memT_trp,memT_trz,memT_tpz,
    memB_vr,memB_vphi,memB_vz,memB_trr,memB_tpp,memB_tzz,memB_trp,memB_trz,memB_tpz,
    memR_vr,memR_vphi,memR_vz,memR_trr,memR_tpp,memR_tzz,memR_trp,memR_trz,memR_tpz

end
#Poroelasticity
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

function PML_save_vel_Top!(memT_vr,memT_vphi,memT_vz,vr,vphi,vz,nr,nz,LPML_z,LPML_r)
    save_memory_Top!(memT_vr,vr,nr,nz,LPML_z,LPML_r)
    save_memory_Top!(memT_vphi,vphi,nr,nz,LPML_z,LPML_r)
    save_memory_Top!(memT_vz,vz,nr,nz,LPML_z,LPML_r)
end
function PML_save_vel_Bottom!(memB_vr,memB_vphi,memB_vz,vr,vphi,vz,nr,nz,LPML_z,LPML_r)
    save_memory_Bottom!(memB_vr,vr,nr,nz,LPML_z,LPML_r)
    save_memory_Bottom!(memB_vphi,vphi,nr,nz,LPML_z,LPML_r)
    save_memory_Bottom!(memB_vz,vz,nr,nz,LPML_z,LPML_r)
end
function PML_save_vel_Right!(memR_vr,memR_vphi,memR_vz,vr,vphi,vz,nr,nz,LPML_z,LPML_r)
    save_memory_Right!(memR_vr,vr,nr,nz,LPML_z,LPML_r)
    save_memory_Right!(memR_vphi,vphi,nr,nz,LPML_z,LPML_r)
    save_memory_Right!(memR_vz,vz,nr,nz,LPML_z,LPML_r)
end

#Poroelasticity
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

function PML_save_stress_Top!(memT_trr,memT_tpp,memT_tzz,memT_trp,memT_trz,memT_tpz,
    trr,tpp,tzz,trp,trz,tpz,nr,nz,LPML_z,LPML_r)
    save_memory_Top!(memT_trr,trr,nr,nz,LPML_z,LPML_r)
    save_memory_Top!(memT_tpp,tpp,nr,nz,LPML_z,LPML_r)
    save_memory_Top!(memT_tzz,tzz,nr,nz,LPML_z,LPML_r)
    save_memory_Top!(memT_trp,trp,nr,nz,LPML_z,LPML_r)
    save_memory_Top!(memT_trz,trz,nr,nz,LPML_z,LPML_r)
    save_memory_Top!(memT_tpz,tpz,nr,nz,LPML_z,LPML_r)
end

function PML_save_stress_Bottom!(memB_trr,memB_tpp,memB_tzz,memB_trp,memB_trz,memB_tpz,
    trr,tpp,tzz,trp,trz,tpz,nr,nz,LPML_z,LPML_r)
    save_memory_Bottom!(memB_trr,trr,nr,nz,LPML_z,LPML_r)
    save_memory_Bottom!(memB_tpp,tpp,nr,nz,LPML_z,LPML_r)
    save_memory_Bottom!(memB_tzz,tzz,nr,nz,LPML_z,LPML_r)
    save_memory_Bottom!(memB_trp,trp,nr,nz,LPML_z,LPML_r)
    save_memory_Bottom!(memB_trz,trz,nr,nz,LPML_z,LPML_r)
    save_memory_Bottom!(memB_tpz,tpz,nr,nz,LPML_z,LPML_r)
end

function PML_save_stress_Right!(memR_trr,memR_tpp,memR_tzz,memR_trp,memR_trz,memR_tpz,
    trr,tpp,tzz,trp,trz,tpz,nr,nz,LPML_z,LPML_r)
    save_memory_Right!(memR_trr,trr,nr,nz,LPML_z,LPML_r)
    save_memory_Right!(memR_tpp,tpp,nr,nz,LPML_z,LPML_r)
    save_memory_Right!(memR_tzz,tzz,nr,nz,LPML_z,LPML_r)
    save_memory_Right!(memR_trp,trp,nr,nz,LPML_z,LPML_r)
    save_memory_Right!(memR_trz,trz,nr,nz,LPML_z,LPML_r)
    save_memory_Right!(memR_tpz,tpz,nr,nz,LPML_z,LPML_r)
end
#Poroelasticity
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


#Poroelastic
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


#Poroelastic
function PML_update_memPQ_1st_Bottom_Por!(Prz,Pzz,PzzPE,
    memB_trr,memB_tpp,memB_tzz,memB_trz,memB_pf,
    trr,tpp,tzz,trz,pf,nr,nz,dr,dz,dt,LPML_z,LPML_r,PML_Wz,PML_Wz2)
#mod in 06tmp5
#mod in Ou1st02tmp:
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


#--updating PML variables: those for stress (Rx and Sxx)
#Poroelasticity
function PML_update_memRS_1st_Top_Por!(Rz,Srz,RzPE,
    memT_vr,memT_vz,memT_vfr,memT_vfz,
    vr,vz,vfr,vfz,nr,nz,dr,dz,dt,LPML_z,LPML_r,PML_Wz,PML_Wz2)
#modefied in 06tmp5
#modified 2: Srz (k=2) 2nd, Srz (k=1) 1st.
#mod in Ou1st02tmp:
#  LPML corresponds to the location of tii, and at tii(k=LPML), OMEGA function is zero
#  Mittet grid. origin is tii
#
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


#Poroelasic
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


#--PML: updating velocities
#Poroelasticity
function PML_update_velocity_1st_Top_Por!(vr,vz,trr,tpp,tzz,trz,vfr,vfz,pf,
    rhomat,rhofmat,D1mat,D2mat,Flagmat_vf_zero,nr,nz,dr,dz,dt,Prz_T,Pzz_T,PzzPE_T,LPML_z,LPML_r)
    #PML_variables (PQRS) index starts from first layer in PML (maybe reversed order)
    #z streching only
    #1st order
    #new in 06tmp5: k=2, vr,vphi->1st, vz->2nd
@inbounds Threads.@threads for j=2:nr-LPML_r #assumig LPML_r>2
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
              #I am alywas accesing [k,j], but definition of r_now changes with components!
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

function PML_update_velocity_Right!(vr,vphi,vz,trr,tpp,tzz,trp,trz,tpz,
    rhomat,m,nr,nz,dr,dz,dt,
    Prr_R,Qrp_R,Ppr_R,Qpp_R,Pzr_R,Qzp_R,Rr_R,Rp_R,Rrp_R,Srp_R,Rrz_R,Spz_R,
    LPML_z,LPML_r)
    #PML_variables (PQRS) index starts from first layer in PML (maybe reversed order)
    #r streching only

    #j=nr-1->vphi,vz (2nd), vr (1st)
    #j=nr->vphi,vz (1st), vr (zero)

    #2nd order
@inbounds Threads.@threads  for j=nr-LPML_r+1:nr-2
          for k=LPML_z+1:nz-LPML_z #assuming LPML_r>2
              #I am alywas accesing [k,j], but definition of r_now changes with components!
             r_now=(j-1)*dr+dr/2.0 #for vr (Mittet)

             vr_now=vr[k,j]
             trr_f1,trr_b1=trr[k,j+1],trr[k,j]
             trr_f2,trr_b2=trr[k,j+2],trr[k,j-1]
             trz_f1,trz_b1=trz[k,j],trz[k-1,j]
             trz_f2,trz_b2=trz[k+1,j],trz[k-2,j]
             trp_now=trp[k,j]
             rho_av=(rhomat[k,j]+rhomat[k,j+1])/2.0 #be carefull

             tpp_2,tpp_1=tpp[k,j+1],tpp[k,j] #see trr
             tpp_3,tpp_0=tpp[k,j+2],tpp[k,j-1]
             vr[k,j]=update_vr_2nd(vr_now,trr_f1,trr_f2,trr_b1,trr_b2,
             trz_f1,trz_f2,trz_b1,trz_b2,trp_now,dz,dr,dt,r_now,rho_av,m,tpp_0,tpp_1,tpp_2,tpp_3)
             #PML addition : caution k and j index (starts from 1)
             vr[k,j]=vr[k,j]+dt/rho_av*(
                Prr_R[k-LPML_z,j-nr+LPML_r]+1.0/r_now*Qrp_R[k-LPML_z,j-nr+LPML_r])

             r_now=(j-1)*dr #for vphi and vz (Mittet)

             vphi_now=vphi[k,j]
             trp_f1,trp_b1=trp[k,j],trp[k,j-1]
             trp_f2,trp_b2=trp[k,j+1],trp[k,j-2]
             tpz_f1,tpz_b1=tpz[k,j],tpz[k-1,j]
             tpz_f2,tpz_b2=tpz[k+1,j],tpz[k-2,j]
             tpp_now=tpp[k,j]
             rho_av=rhomat[k,j] #be carefull
             vphi[k,j]=update_vphi_2nd(vphi_now,trp_f1,trp_f2,trp_b1,trp_b2,
             tpz_f1,tpz_f2,tpz_b1,tpz_b2,tpp_now,dz,dr,dt,r_now,rho_av,m)
             #PML addition : caution k and j index (starts from 1)
             vphi[k,j]=vphi[k,j]+dt/rho_av*(
                Ppr_R[k-LPML_z,j-nr+LPML_r]+1.0/r_now*Qpp_R[k-LPML_z,j-nr+LPML_r])


             vz_now=vz[k,j]
             trz_f1,trz_b1=trz[k,j],trz[k,j-1]
             trz_f2,trz_b2=trz[k,j+1],trz[k,j-2]
             tzz_f1,tzz_b1=tzz[k+1,j],tzz[k,j]
             tzz_f2,tzz_b2=tzz[k+2,j],tzz[k-1,j]
             tpz_now=tpz[k,j]
             rho_av=(rhomat[k,j]+rhomat[k+1,j])/2.0 #be careful
             vz[k,j]=update_vz_2nd(vz_now,trz_f1,trz_f2,trz_b1,trz_b2,
             tzz_f1,tzz_f2,tzz_b1,tzz_b2,tpz_now,dz,dr,dt,r_now,rho_av,m)
             #PML addition : caution k and j index (starts from 1)
             vz[k,j]=vz[k,j]+dt/rho_av*(
                Pzr_R[k-LPML_z,j-nr+LPML_r]+1.0/r_now*Qzp_R[k-LPML_z,j-nr+LPML_r])

         end #k (z)
    end #j (r)

    #taken from CylFDMod_LVTS02.jl
    #j=nr-1->vphi,vz (2nd), vr (1st)
    #j=nr->vphi,vz (1st), vr (zero)
    #--1st order FD

 j=nr-1
@inbounds for k=LPML_z+1:nz-LPML_z
    r_now=(j-1)*dr+dr/2.0 #for vr (Mittet)

    vr_now=vr[k,j]
    trr_f,trr_b=trr[k,j+1],trr[k,j]
    trz_f,trz_b=trz[k,j],trz[k-1,j]
    trr_av=(trr[k,j+1]+trr[k,j])/2.0
    tpp_av=(tpp[k,j+1]+tpp[k,j])/2.0
    trp_now=trp[k,j]
    rho_av=(rhomat[k,j]+rhomat[k,j+1])/2.0 #be carefull
    vr[k,j]=update_vr(vr_now,trr_f,trr_b,trz_f,trz_b,trr_av,tpp_av,trp_now,dz,dr,dt,r_now,rho_av,m)
    #PML addition : caution k and j index (starts from 1)
    vr[k,j]=vr[k,j]+dt/rho_av*(
       Prr_R[k-LPML_z,j-nr+LPML_r]+1.0/r_now*Qrp_R[k-LPML_z,j-nr+LPML_r])

   r_now=(j-1)*dr #for vphi and vz (Mittet)

   vphi_now=vphi[k,j]
   trp_f1,trp_b1=trp[k,j],trp[k,j-1]
   trp_f2,trp_b2=trp[k,j+1],trp[k,j-2]
   tpz_f1,tpz_b1=tpz[k,j],tpz[k-1,j]
   tpz_f2,tpz_b2=tpz[k+1,j],tpz[k-2,j]
   tpp_now=tpp[k,j]
   rho_av=rhomat[k,j] #be carefull
   vphi[k,j]=update_vphi_2nd(vphi_now,trp_f1,trp_f2,trp_b1,trp_b2,
   tpz_f1,tpz_f2,tpz_b1,tpz_b2,tpp_now,dz,dr,dt,r_now,rho_av,m)
   #PML addition : caution k and j index (starts from 1)
   vphi[k,j]=vphi[k,j]+dt/rho_av*(
      Ppr_R[k-LPML_z,j-nr+LPML_r]+1.0/r_now*Qpp_R[k-LPML_z,j-nr+LPML_r])


   vz_now=vz[k,j]
   trz_f1,trz_b1=trz[k,j],trz[k,j-1]
   trz_f2,trz_b2=trz[k,j+1],trz[k,j-2]
   tzz_f1,tzz_b1=tzz[k+1,j],tzz[k,j]
   tzz_f2,tzz_b2=tzz[k+2,j],tzz[k-1,j]
   tpz_now=tpz[k,j]
   rho_av=(rhomat[k,j]+rhomat[k+1,j])/2.0 #be careful
   vz[k,j]=update_vz_2nd(vz_now,trz_f1,trz_f2,trz_b1,trz_b2,
   tzz_f1,tzz_f2,tzz_b1,tzz_b2,tpz_now,dz,dr,dt,r_now,rho_av,m)
   #PML addition : caution k and j index (starts from 1)
   vz[k,j]=vz[k,j]+dt/rho_av*(
      Pzr_R[k-LPML_z,j-nr+LPML_r]+1.0/r_now*Qzp_R[k-LPML_z,j-nr+LPML_r])
  end #k

  j=nr
 @inbounds for k=LPML_z+1:nz-LPML_z
     r_now=(j-1)*dr #for vphi and vz (Mittet)

     vphi_now=vphi[k,j]
     trp_f,trp_b=trp[k,j],trp[k,j-1]
     tpz_f,tpz_b=tpz[k,j],tpz[k-1,j]
     trp_av=(trp[k,j]+trp[k,j-1])/2.0
     tpp_now=tpp[k,j]
     rho_av=rhomat[k,j] #be carefull
     vphi[k,j]=update_vphi(vphi_now,trp_f,trp_b,tpz_f,tpz_b,trp_av,tpp_now,dz,dr,dt,r_now,rho_av,m)
     #PML addition : caution k and j index (starts from 1)
     vphi[k,j]=vphi[k,j]+dt/rho_av*(
        Ppr_R[k-LPML_z,j-nr+LPML_r]+1.0/r_now*Qpp_R[k-LPML_z,j-nr+LPML_r])


     vz_now=vz[k,j]
     trz_f,trz_b=trz[k,j],trz[k,j-1]
     tzz_f,tzz_b=tzz[k+1,j],tzz[k,j]
     trz_av=(trz[k,j]+trz[k,j-1])/2.0
     tpz_now=tpz[k,j]
     rho_av=(rhomat[k,j]+rhomat[k+1,j])/2.0 #be careful
     vz[k,j]=update_vz(vz_now,trz_f,trz_b,tzz_f,tzz_b,trz_av,tpz_now,dz,dr,dt,r_now,rho_av,m)
     #PML addition : caution k and j index (starts from 1)
     vz[k,j]=vz[k,j]+dt/rho_av*(
        Pzr_R[k-LPML_z,j-nr+LPML_r]+1.0/r_now*Qzp_R[k-LPML_z,j-nr+LPML_r])
   end #k



#==original
    #--1st order FD
 j=nr-1
@inbounds for k=LPML_z+1:nz-LPML_z
    r_now=(j-1)*dr+dr/2.0 #for vr (Mittet)

    vr_now=vr[k,j]
    trr_f,trr_b=trr[k,j+1],trr[k,j]
    trz_f,trz_b=trz[k,j],trz[k-1,j]
    trr_av=(trr[k,j+1]+trr[k,j])/2.0
    tpp_av=(tpp[k,j+1]+tpp[k,j])/2.0
    trp_now=trp[k,j]
    rho_av=(rhomat[k,j]+rhomat[k,j+1])/2.0 #be carefull
    vr[k,j]=update_vr(vr_now,trr_f,trr_b,trz_f,trz_b,trr_av,tpp_av,trp_now,dz,dr,dt,r_now,rho_av,m)
    #PML addition : caution k and j index (starts from 1)
    vr[k,j]=vr[k,j]+dt/rho_av*(
       Prr_R[k-LPML_z,j-nr+LPML_r]+1.0/r_now*Qrp_R[k-LPML_z,j-nr+LPML_r])

    r_now=(j-1)*dr #for vphi and vz (Mittet)

    vphi_now=vphi[k,j]
    trp_f,trp_b=trp[k,j],trp[k,j-1]
    tpz_f,tpz_b=tpz[k,j],tpz[k-1,j]
    trp_av=(trp[k,j]+trp[k,j-1])/2.0
    tpp_now=tpp[k,j]
    rho_av=rhomat[k,j] #be carefull
    vphi[k,j]=update_vphi(vphi_now,trp_f,trp_b,tpz_f,tpz_b,trp_av,tpp_now,dz,dr,dt,r_now,rho_av,m)
    #PML addition : caution k and j index (starts from 1)
    vphi[k,j]=vphi[k,j]+dt/rho_av*(
       Ppr_R[k-LPML_z,j-nr+LPML_r]+1.0/r_now*Qpp_R[k-LPML_z,j-nr+LPML_r])


    vz_now=vz[k,j]
    trz_f,trz_b=trz[k,j],trz[k,j-1]
    tzz_f,tzz_b=tzz[k+1,j],tzz[k,j]
    trz_av=(trz[k,j]+trz[k,j-1])/2.0
    tpz_now=tpz[k,j]
    rho_av=(rhomat[k,j]+rhomat[k+1,j])/2.0 #be careful
    vz[k,j]=update_vz(vz_now,trz_f,trz_b,tzz_f,tzz_b,trz_av,tpz_now,dz,dr,dt,r_now,rho_av,m)
    #PML addition : caution k and j index (starts from 1)
    vz[k,j]=vz[k,j]+dt/rho_av*(
       Pzr_R[k-LPML_z,j-nr+LPML_r]+1.0/r_now*Qzp_R[k-LPML_z,j-nr+LPML_r])
  end #k
  ==#

end

#Poroelasic
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
              #I am alywas accesing [k,j], but definition of r_now changes with components!
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


#Poroelasic
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
              #I am alywas accesing [k,j], but definition of r_now changes with components!
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


#--PML: updating stress
#Poroelastic
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
             #I am alywas accesing [k,j], but definition of r_now changes with components!
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


#Poroelastic
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
             #I am alywas accesing [k,j], but definition of r_now changes with components!
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
             #I am alywas accesing [k,j], but definition of r_now changes with components!
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
             #I am alywas accesing [k,j], but definition of r_now changes with components!
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
             #I am alywas accesing [k,j], but definition of r_now changes with components!
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

function ApplyBC_stress_AcousticMedia!(trr,tpp,tzz,trz,pf,Flag_AC,nr,nz)
    #when FlagAC==1, then set tau=-pf*eye
    @inbounds Threads.@threads for j=1:nr
        for k=1:nz
            if(Flag_AC[k,j]==1)
                trr[k,j]=-pf[k,j]
                tpp[k,j]=-pf[k,j]
                tzz[k,j]=-pf[k,j]
                trz[k,j]=0.0
            end
        end
    end
end

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

function ApplyBC_stress_ElasticMedia_Guan!(trr,tpp,tzz,pf,Flag_E,nr,nz)
#when FlagE==1, then set pf to be negative value of normal stress (Guan 2011)
#doi: 10.4208/cicp.020810.161210a
#
@inbounds Threads.@threads for j=1:nr
    for k=1:nz
        if(Flag_E[k,j]==1)
            pf[k,j]=-(trr[k,j]+tpp[k,j]+tzz[k,j])/3
        end
    end
end
end



function get_Flag_vf_zero(Flag_AC,Flag_E,nr,nz)
    #when FlagAC==1 or Flag_E==1, then set Flag_vf_zero==1
    Flag_vf_zero=zeros(nz,nr)
    @inbounds Threads.@threads for j=1:nr
        for k=1:nz
            if(Flag_AC[k,j]==1 || Flag_E[k,j]==1)
                Flag_vf_zero[k,j]=1
            end
        end
    end
    return Flag_vf_zero
end




#---
#---Activating vfr@vertical boundary separating Fluid-PoroElastic media
function update_vr_vfr_1st_vertical(vr,trr,tpp,tzz,trz,vfr,pf,
    vr_org,vfr_org,
    Rho,Rhof,D1,D2,nr,nz,dr,dz,dt,LPML_z,
    Prz_T,Prz_B,
    ir_wall,Flag_AC,Flag_E)
    #Acoustic-Poroelastic vertical boundary
    #Calculating vrw and vr at boundary

    j=ir_wall #wall
    #main region
    @inbounds Threads.@threads for k=LPML_z+1:nz-LPML_z

        #check if I need to activate vfr
        if(Flag_AC[k,j]==1 && Flag_AC[k,j+1]==0 && Flag_E[k,j+1]==0)#
              #I am alywas accesing [k,j], but definition of r_now changes with components!
             r_now=(j-1)*dr+dr/2.0 #for vr (Mittet)
             vr_now=vr_org[k,j]
             trr_f1,trr_b1=trr[k,j+1],trr[k,j]
             trz_f1,trz_b1=trz[k,j],trz[k-1,j] #should be zero

             trr_av=0.5*(trr_f1+trr_b1) #this looks ok according to BC (pf=-trr, now trr_b1 is -pf)
             tpp_av=0.5*(tpp[k,j+1]+tpp[k,j]) #see trr. I am not sure if this is ok for BC
#             tpp_av=tpp[k,j+1] #test. assuming tpp is discont

             vfr_now=vfr_org[k,j]
             pf_f1,pf_b1=pf[k,j+1],pf[k,j] #

             rho_av=0.5*(Rho[k,j]+Rho[k,j+1]) #be carefull (normal)
             rhof_av=0.5*(Rhof[k,j]+Rhof[k,j+1]) #be carefull (normal)
             D1_av=0.5*(D1[k,j]+D1[k,j+1]) #be carefull (normal)
             D2_av=0.5*(D2[k,j]+D2[k,j+1]) #be carefull (normal)

             #Normal
             vr[k,j],vfr[k,j]=update_vr_1st_Por(vr_now,vfr_now,
                 trr_f1,trr_b1,trz_f1,trz_b1,trr_av,tpp_av,
                 pf_f1,pf_b1,
                 dz,dr,dt,r_now,rho_av,rhof_av,D1_av,D2_av,0)
             #

             #== TMP
             drtrr=(trr_f1-trr_b1)/dr
             dztrz=(trz_f1-trz_b1)/dz
             A=drtrr+(trr_av-tpp_av)/(r_now)+dztrz #common term in solid and fluid phase
             drpf=(pf_f1-pf_b1)/dr

             vfr_old=vfr_now
             vfr_now=(D2_av/2+(D1_av-rhof_av^2/rho_av)/dt)^(-1) * (
                         (-D2_av/2+(D1_av-rhof_av^2/rho_av)/dt)*vfr_old-drpf-rhof_av/rho_av*A
                         )
             dert_vfr=(vfr_now-vfr_old)/dt

             vr_old=vr_now
             vr_now1=vr_old+dt/rho_av*(A-rhof_av*dert_vfr)

             avrg_vfr=(vfr_now+vfr_old)/2
             vr_now2=vr_old+dt/rhof_av*(-drpf-D1_av*dert_vfr-D2_av*avrg_vfr)
             if(vr_now1!=vr_now2)
                 println("vr is not consistent!",vr_now1,",",vr_now2)
             end
             vfr[k,j]=vfr_now
             vr[k,j]=vr_now1
             ==# #TMP


         end #if AC-PE BC
    end #k

#
    #Top PML
    @inbounds Threads.@threads for k=2:LPML_z
        #check if I need to activate vfr
        if(Flag_AC[k,j]==1 && Flag_AC[k,j+1]==0 && Flag_E[k,j+1]==0)#
              #I am alywas accesing [k,j], but definition of r_now changes with components!
             r_now=(j-1)*dr+dr/2.0 #for vr (Mittet)
             vr_now=vr_org[k,j]
             trr_f1,trr_b1=trr[k,j+1],trr[k,j]
             trz_f1,trz_b1=trz[k,j],trz[k-1,j] #should be zero

             trr_av=0.5*(trr_f1+trr_b1) #this looks ok according to BC (pf=-trr, now trr_b1 is -pf)
             tpp_av=0.5*(tpp[k,j+1]+tpp[k,j]) #see trr. I am not sure if this is ok for BC
#             tpp_av=tpp[k,j+1] #test. assuming tpp is discont

             vfr_now=vfr_org[k,j]
             pf_f1,pf_b1=pf[k,j+1],pf[k,j] #

             rho_av=0.5*(Rho[k,j]+Rho[k,j+1]) #be carefull (normal)
             rhof_av=0.5*(Rhof[k,j]+Rhof[k,j+1]) #be carefull (normal)
             D1_av=0.5*(D1[k,j]+D1[k,j+1]) #be carefull (normal)
             D2_av=0.5*(D2[k,j]+D2[k,j+1]) #be carefull (normal)

             vfr_old=vfr_now

             #Normal
             vr_dummy,vfr[k,j]=update_vr_1st_Por(vr_now,vfr_now,
                    trr_f1,trr_b1,trz_f1,trz_b1,trr_av,tpp_av,
                    pf_f1,pf_b1,
                    dz,dr,dt,r_now,rho_av,rhof_av,D1_av,D2_av,0)


              #PML addition (vfr)
              vfr[k,j]=vfr[k,j]+
                     (D2_av/2+(D1_av-rhof_av^2/rho_av)/dt)^(-1)*(
                                       -rhof_av/rho_av*Prz_T[LPML_z-k+1,j]
                                       )

               #PML (vr) calculate elastic one then vfr addition
               vr[k,j]=update_vr(vr_now,trr_f1,trr_b1,trz_f1,trz_b1,trr_av,tpp_av,
                      0,dz,dr,dt,r_now,rho_av,0)
               vr[k,j]=vr[k,j]+dt/rho_av*Prz_T[LPML_z-k+1,j] #elastic+PML
               dert_vfr=(vfr[k,j]-vfr_old)/dt
               vr[k,j]=vr[k,j]-dt*rhof_av/rho_av*dert_vfr #poroelastic PML
             #
         end #if AC-PE BC
    end #k

#
    #Bottom PML
    @inbounds Threads.@threads for k=nz-LPML_z+1:nz-1

        #check if I need to activate vfr
        if(Flag_AC[k,j]==1 && Flag_AC[k,j+1]==0 && Flag_E[k,j+1]==0)#
              #I am alywas accesing [k,j], but definition of r_now changes with components!
             r_now=(j-1)*dr+dr/2.0 #for vr (Mittet)
             vr_now=vr_org[k,j]
             trr_f1,trr_b1=trr[k,j+1],trr[k,j]
             trz_f1,trz_b1=trz[k,j],trz[k-1,j] #should be zero

             trr_av=0.5*(trr_f1+trr_b1) #this looks ok according to BC (pf=-trr, now trr_b1 is -pf)
             tpp_av=0.5*(tpp[k,j+1]+tpp[k,j]) #see trr. I am not sure if this is ok for BC
#             tpp_av=tpp[k,j+1] #test. assuming tpp is discont

             vfr_now=vfr_org[k,j]
             pf_f1,pf_b1=pf[k,j+1],pf[k,j] #

             rho_av=0.5*(Rho[k,j]+Rho[k,j+1]) #be carefull (normal)
             rhof_av=0.5*(Rhof[k,j]+Rhof[k,j+1]) #be carefull (normal)
             D1_av=0.5*(D1[k,j]+D1[k,j+1]) #be carefull (normal)
             D2_av=0.5*(D2[k,j]+D2[k,j+1]) #be carefull (normal)

             vfr_old=vfr_now

             #Normal
             vr_dummy,vfr[k,j]=update_vr_1st_Por(vr_now,vfr_now,
                    trr_f1,trr_b1,trz_f1,trz_b1,trr_av,tpp_av,
                    pf_f1,pf_b1,
                    dz,dr,dt,r_now,rho_av,rhof_av,D1_av,D2_av,0)


              #PML addition (vfr)
              vfr[k,j]=vfr[k,j]+
                     (D2_av/2+(D1_av-rhof_av^2/rho_av)/dt)^(-1)*(
                                       -rhof_av/rho_av*Prz_B[k-nz+LPML_z,j]
                                       )

               #PML (vr) calculate elastic one then vfr addition
               vr[k,j]=update_vr(vr_now,trr_f1,trr_b1,trz_f1,trz_b1,trr_av,tpp_av,
                      0,dz,dr,dt,r_now,rho_av,0)
               vr[k,j]=vr[k,j]+dt/rho_av*Prz_B[k-nz+LPML_z,j] #elastic PML
               dert_vfr=(vfr[k,j]-vfr_old)/dt
               vr[k,j]=vr[k,j]-dt*rhof_av/rho_av*dert_vfr #poroelastic PML
             #
         end #if AC-PE BC
    end #k

end



#---

function tmp_check(vr,vz,trr,tpp,tzz,trz,vfr,vfz,pf,
    rhomat,rhofmat,D1mat,D2mat,Flagmat_vf_zero,nr,nz,dr,dz,LPML_z,LPML_r)
#2nd order FD (Graves, 1996)
#assuming PML!
#@inbounds Threads.@threads for j=3:nr-2
#      for k=3:nz-2

k=LPML_z
j=5

          #I am alywas accesing [k,j], but definition of r_now changes with components!
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

         tmp_vr,tmp_vfr=update_vr_1st_Por(vr_now,vfr_now,
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

         tmp_vz,tmp_vfz=update_vz_1st_Por(vz_now,vfz_now,
             trz_f1,trz_b1,tzz_f1,tzz_b1,trz_av,
             pf_f1,pf_b1,
             dz,dr,dt,r_now,rho_av,rhof_av,D1_av,D2_av,Flagmat_vf_zero[k,j])

             println("Normal vz(LPML_z,5)=",tmp_vz)

end #function



function mycopy_mat(MATin,MATout,nz,nr)
    @inbounds Threads.@threads for j=1:nr
        for k=1:nz
            MATout[k,j]=MATin[k,j]
        end
    end
end




#---for checking purpose---
function calculate_divv_1st_Por!(divv,vr,vz,
    nr,nz,dr,dz)
#calculate divergence of velocity:dvrdr+vr/r+dvzdz

@inbounds Threads.@threads for j=2:nr
      for k=2:nz
         #I am alywas accesing [k,j], but definition of r_now changes with components!
         r_now=(j-1)*dr #for trr,tpp,tzz,tpz Mittet

         vr_f1,vr_b1=vr[k,j],vr[k,j-1]
         vz_f1,vz_b1=vz[k,j],vz[k-1,j]

         vr_av=0.5*(vr_f1+vr_b1)

         drvr=vr_f1-vr_b1
         drvr=drvr/dr
         dzvz=vz_f1-vz_b1
         dzvz=dzvz/dz

         divv[k,j]=drvr+vr_av/r_now+dzvz


   end #k (z)
end #j (r)

end #function

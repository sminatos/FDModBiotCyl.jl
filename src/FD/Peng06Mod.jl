#Peng.jl
#Solving BC, obtain coeffients , and calculate waveform
#n=0 (axisymetric) only
#01:start
#  : Scheoenebrg's A0,B0,C0 was negative sign of the current code's A0,B0,C0.
#  : This is due to the definition of V(w) -> vz=-V exp(ikz*z-wt)? Maybe OK.
#  :
#02: Important; Everything assumes V(w)=w^2, where V(w) is defined in plane P-wave displacement potential phiP
#  : as phi_P=-Vp_solid/w^2*V(w)exp(ikz*z-iwt) <--> vz=-V(w)exp(ikz*z-iwt)
#  : think carefully which quantity you want to have source function controlled!
#  : several tests pass
#03: waveform
#04: compare with Schoenberg's approximation
#Mod: jsut removing "include" and "using"
#06: modifications found in using PCylFD_planewave


#==
using CPUTime
using ProgressMeter
using Interpolations
using MAT
using SpecialFunctions
using Base64
using DSP
using FFTW
using Plots
using JLD
using Dierckx
using DelimitedFiles
using Printf

include("./CylFDMod_Mittet06tmp3.jl") #main modules
==#

function tmp_solveBC()
    rb=0.1

    Vp_fluid=1500.
    Rho_fluid=1000.
    Vp_solid=3000.
    Vs_solid=2000.
    Rho_solid=2000.

    cT=sqrt(Vp_fluid^2/(1+Rho_fluid*Vp_fluid^2/(Rho_solid*Vs_solid^2)))
    g=(Vs_solid/Vp_solid)^2

    freq=100.
    w=2pi*freq
    angle_delta=1E-2 # =0 normal incidence, but diveeges due to Hankel function.
    A0,B0,C0=SolvePengBC01(w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)

    @show A0,B0,C0

    #--check continuation
#==continuation trr is ok
@show  trr_f(rb,A0,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)
@show  trr_phi(rb,B0,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)+
      trr_xi(rb,C0,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)+
      trr_P(rb,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)
==#
#==continuation vr is ok
@show ur_f(rb,A0,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)
@show ur_phi(rb,B0,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)+
     ur_xi(rb,C0,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)+
     ur_P(rb,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)
==#
#==vanishing trz is ok
@show trz_phi(rb,B0,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)+
      trz_xi(rb,C0,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)+
      trz_P(rb,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)
==#
    #--check uz ratio
    #==uz ratio is OK
    test=uz_f(rb,A0,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)/
    (uz_phi(rb,B0,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)+
     uz_xi(rb,C0,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)+
     uz_P(rb,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid))
    #Schoenberg's ratio
@show    test_Sch=cT^2/Vs_solid^2*(1-2*g)/(1-cT^2/Vp_solid^2)
@show test
    ==#


end


function getwaveform_Peng()

    Vp_fluid=1500.
    Rho_fluid=1000.
    Vp_solid=3000.
    Vs_solid=2000.
    Rho_solid=2000.

    cT=sqrt(Vp_fluid^2/(1+Rho_fluid*Vp_fluid^2/(Rho_solid*Vs_solid^2)))
    g=(Vs_solid/Vp_solid)^2


    #---src parameters--
    ns=3001
    dt=0.5E-5
    T=(ns-1)*dt
    f0=300 #src Freq
    delay=1/f0*1.1

    tvec=range(0.0,T,length=ns) # range object (no memory allocation)
    tvec=collect(tvec) # a vector
    src_func=myricker2(tvec,f0,delay,2) #when using 2nd derivative Gaussian
    tmp_maxamp=maximum(map(abs,src_func))
    src_func=src_func/tmp_maxamp
    display(plot(tvec,src_func[:],title="src function"))

    #---freq processings--
    fvec=[0:1/T:1/T*(ns-1)]
    fvec=collect(fvec[1])

    fsrc=rfft(src_func[:])
    test_src=irfft(fsrc,ns)
    ns_half=length(fsrc)
    fvec=fvec[1:ns_half]
    plot(fvec,map(abs,fsrc),xlims=(0,1000))

    max_freq_src=3*f0
    if_max=minimum(findall(x -> x > max_freq_src,fvec)) #to be used to filter high-freq noise
    #-----

    #--Dimension parameters-
    rb=0.15
    nr=101
    dr=0.025
    ir_wall=Int(round(rb/dr))
    rb=(ir_wall-1)*dr+dr/2.0

    #---output waveforms
    #--outputs are vr[ns,nr] and vz[ns,nr]
    vr_Peng=zeros(ns,nr)
    vz_Peng=zeros(ns,nr)

    vphi_Peng=zeros(ns,nr)

    trr_Peng=zeros(ns,nr)
    tpp_Peng=zeros(ns,nr)
    tzz_Peng=zeros(ns,nr)
    trp_Peng=zeros(ns,nr)
    trz_Peng=zeros(ns,nr)
    tpz_Peng=zeros(ns,nr)


    #--repeat every r_now
    for ir=1:nr

    #---Start frequency response
    fvz=complex(zeros(1,ns_half))
    fvr=complex(zeros(1,ns_half))

    ftrr=complex(zeros(1,ns_half))
    ftpp=complex(zeros(1,ns_half))
    ftzz=complex(zeros(1,ns_half))

    ftrz=complex(zeros(1,ns_half)) #other shear component are zero (axisymmetric)

    for iw=2:if_max
        freq=fvec[iw]
        w=2pi*freq
        angle_delta=1E-4 # =0 normal incidence, but give a small number as it diverges due to Hankel function.
        A0,B0,C0=SolvePengBC01(w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)
        #    @show A0,B0,C0
        r_now=(ir-1)*dr #for vz (I keep this same in all field variables. When Staggered grid FDTD, this should be different.)
        if (ir<=ir_wall)
            fvz[iw]=(-im*w)*uz_f(r_now,A0,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)
            fvr[iw]=(-im*w)*ur_f(r_now,A0,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)

            ftrr[iw]=trr_f(r_now,A0,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)
            ftpp[iw]=trr_f(r_now,A0,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)
            ftzz[iw]=trr_f(r_now,A0,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)

            ftrz[iw]=0.0
        else
            fvz[iw]=(-im*w)*(
                uz_phi(r_now,B0,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)+
                uz_xi(r_now,C0,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)+
                uz_P(r_now,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)
                )

            fvr[iw]=(-im*w)*(
                ur_phi(r_now,B0,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)+
                ur_xi(r_now,C0,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)+
                ur_P(r_now,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)
                )

            ftrr[iw]=(
                trr_phi(r_now,B0,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)+
                trr_xi(r_now,C0,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)+
                trr_P(r_now,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)
                )

            ftpp[iw]=(
                tpp_phi(r_now,B0,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)+
                tpp_xi(r_now,C0,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)+
                tpp_P(r_now,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)
                )

            ftzz[iw]=(
                tzz_phi(r_now,B0,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)+
                tzz_xi(r_now,C0,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)+
                tzz_P(r_now,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)
                )

            ftrz[iw]=(
                trz_phi(r_now,B0,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)+
                trz_xi(r_now,C0,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)+
                trz_P(r_now,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)
                )

        end

        fvz[iw]=conj(fvz[iw]) #Aki -> Kees FT
        fvz[iw]=fvz[iw]*(fsrc[iw]/w^2) # /w^2 scaling is due to V(w)=w^2 in Peng.

        fvr[iw]=conj(fvr[iw]) #Aki -> Kees FT
        fvr[iw]=fvr[iw]*(fsrc[iw]/w^2)

        ftrr[iw]=conj(ftrr[iw]) #Aki -> Kees FT
        ftrr[iw]=ftrr[iw]*(fsrc[iw]/w^2)

        ftpp[iw]=conj(ftpp[iw]) #Aki -> Kees FT
        ftpp[iw]=ftpp[iw]*(fsrc[iw]/w^2)

        ftzz[iw]=conj(ftzz[iw]) #Aki -> Kees FT
        ftzz[iw]=ftzz[iw]*(fsrc[iw]/w^2)

        ftrz[iw]=conj(ftrz[iw]) #Aki -> Kees FT
        ftrz[iw]=ftrz[iw]*(fsrc[iw]/w^2)

    end
    vz_tmp=irfft(fvz[:],ns)
    vr_tmp=irfft(fvr[:],ns)
    trr_tmp=irfft(ftrr[:],ns)
    tpp_tmp=irfft(ftpp[:],ns)
    tzz_tmp=irfft(ftzz[:],ns)
    trz_tmp=irfft(ftrz[:],ns)

    vz_Peng[:,ir]=vz_tmp
    vr_Peng[:,ir]=vr_tmp
    trr_Peng[:,ir]=trr_tmp
    tpp_Peng[:,ir]=tpp_tmp
    tzz_Peng[:,ir]=tzz_tmp
    trz_Peng[:,ir]=trz_tmp
    end

return vr_Peng,vz_Peng,trr_Peng,tpp_Peng,tzz_Peng,trz_Peng

    #--check continuation
#==continuation trr is ok
@show  trr_f(rb,A0,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)
@show  trr_phi(rb,B0,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)+
      trr_xi(rb,C0,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)+
      trr_P(rb,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)
==#
#==continuation vr is ok
@show ur_f(rb,A0,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)
@show ur_phi(rb,B0,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)+
     ur_xi(rb,C0,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)+
     ur_P(rb,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)
==#
#==vanishing trz is ok
@show trz_phi(rb,B0,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)+
      trz_xi(rb,C0,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)+
      trz_P(rb,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)
==#
    #--check uz ratio
    #==uz ratio is OK
    test=uz_f(rb,A0,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)/
    (uz_phi(rb,B0,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)+
     uz_xi(rb,C0,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)+
     uz_P(rb,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid))
    #Schoenberg's ratio
@show    test_Sch=cT^2/Vs_solid^2*(1-2*g)/(1-cT^2/Vp_solid^2)
@show test
    ==#

#    @show    test_Sch=cT^2/Vs_solid^2*(1-2*g)/(1-cT^2/Vp_solid^2)

end

function getwaveform_Schoenberg()

    Vp_fluid=1500.
    Rho_fluid=1000.
    Vp_solid=3000.
    Vs_solid=2000.
    Rho_solid=2000.

    cT=sqrt(Vp_fluid^2/(1+Rho_fluid*Vp_fluid^2/(Rho_solid*Vs_solid^2)))
    g=(Vs_solid/Vp_solid)^2


    #---src parameters--
    ns=3001
    dt=0.5E-5
    T=(ns-1)*dt
    f0=300 #src Freq
    delay=1/f0*1.1

    tvec=range(0.0,T,length=ns) # range object (no memory allocation)
    tvec=collect(tvec) # a vector
    src_func=myricker2(tvec,f0,delay,2) #when using 2nd derivative Gaussian
    tmp_maxamp=maximum(map(abs,src_func))
    src_func=src_func/tmp_maxamp
    display(plot(tvec,src_func[:],title="src function"))

    #---freq processings--
    fvec=[0:1/T:1/T*(ns-1)]
    fvec=collect(fvec[1])

    fsrc=rfft(src_func[:])
    test_src=irfft(fsrc,ns)
    ns_half=length(fsrc)
    fvec=fvec[1:ns_half]
    plot(fvec,map(abs,fsrc),xlims=(0,1000))

    max_freq_src=3*f0
    if_max=minimum(findall(x -> x > max_freq_src,fvec)) #to be used to filter high-freq noise


    #-----

    #--Dimension parameters-
    rb=0.15
    nr=101
    dr=0.025
    ir_wall=Int(round(rb/dr))
    rb=(ir_wall-1)*dr+dr/2.0

    #---output waveforms
    #--outputs are vr[ns,nr] and vz[ns,nr]
    vr_Sch=zeros(ns,nr)
    vz_Sch=zeros(ns,nr)

#    vphi_Sch=zeros(ns,nr) #zero

    trr_Sch=zeros(ns,nr)
    tpp_Sch=zeros(ns,nr)
    tzz_Sch=zeros(ns,nr)
#    trp_Sch=zeros(ns,nr) #zero
    trz_Sch=zeros(ns,nr)
#    tpz_Sch=zeros(ns,nr) #zero

    GamL_vec=2.0*pi*fvec*rb/Vp_solid

    #--repeat every r_now
    for ir=1:nr

    #---Start frequency response
    fvz=complex(zeros(1,ns_half))
    fvr=complex(zeros(1,ns_half))

    ftrr=complex(zeros(1,ns_half))
    ftpp=complex(zeros(1,ns_half))
    ftzz=complex(zeros(1,ns_half))

    ftrz=complex(zeros(1,ns_half)) #other shear component are zero (axisymmetric)

    for iw=2:if_max
        freq=fvec[iw]
        w=2pi*freq
        #    @show A0,B0,C0
        r_now=(ir-1)*dr #for vz (I keep this same in all field variables. When Staggered grid FDTD, this should be different.)
        if (ir<=ir_wall)
            fvz[iw]=cT^2/Vs_solid^2*(1.0-2.0*g)/(1.0-cT^2/Vp_solid^2)
            fvr[iw]=im*GamL_vec[iw]*r_now/(2.0*rb)*cT^2/(g*Vp_fluid^2)*
                   (1.0-2.0*g)*(1.0-Vp_fluid^2/Vp_solid^2)/(1.0-cT^2/Vp_solid^2)

            ftrr[iw]=-Vp_solid*Rho_fluid*cT^2/Vs_solid^2*(1.0-2.0*g)/(1.0-cT^2/Vp_solid^2)
            ftpp[iw]=-Vp_solid*Rho_fluid*cT^2/Vs_solid^2*(1.0-2.0*g)/(1.0-cT^2/Vp_solid^2)
            ftzz[iw]=-Vp_solid*Rho_fluid*cT^2/Vs_solid^2*(1.0-2.0*g)/(1.0-cT^2/Vp_solid^2)

            ftrz[iw]=0.0
        else
            fvz[iw]=1.0

            fvr[iw]=im*GamL_vec[iw]*rb/(2.0*r_now)*cT^2/(g*Vp_fluid^2)*
                 (1.0-2.0*g)*(1.0-Vp_fluid^2/Vp_solid^2)/(1.0-cT^2/Vp_solid^2)

            TMP=im*GamL_vec[iw]*rb/2.0*cT^2/(g*Vp_fluid^2)*
                    (1.0-2.0*g)*(1.0-Vp_fluid^2/Vp_solid^2)/(1.0-cT^2/Vp_solid^2)

            ftrr[iw]=(-Rho_solid*Vp_solid+2.0*Rho_solid*Vs_solid^2/Vp_solid) +
                      1.0/(-im*2pi*fvec[iw]).*(-2.0*Rho_solid*Vs_solid^2*TMP/r_now^2)

            ftpp[iw]=-Rho_solid/Vp_solid*(Vp_solid^2-2.0*Vs_solid^2)-
                      1.0/(-im*2pi*fvec[iw]).*(-2.0*Rho_solid*Vs_solid^2*TMP/r_now^2)

            ftzz[iw]=-Rho_solid*Vp_solid

            ftrz[iw]=0.0

        end

        fvz[iw]=conj(fvz[iw]) #Aki -> Kees FT
        fvz[iw]=fvz[iw]*fsrc[iw] #

        fvr[iw]=conj(fvr[iw]) #Aki -> Kees FT
        fvr[iw]=fvr[iw]*fsrc[iw]

        ftrr[iw]=conj(ftrr[iw]) #Aki -> Kees FT
        ftrr[iw]=ftrr[iw]*fsrc[iw]

        ftpp[iw]=conj(ftpp[iw]) #Aki -> Kees FT
        ftpp[iw]=ftpp[iw]*fsrc[iw]

        ftzz[iw]=conj(ftzz[iw]) #Aki -> Kees FT
        ftzz[iw]=ftzz[iw]*fsrc[iw]

        ftrz[iw]=conj(ftrz[iw]) #Aki -> Kees FT
        ftrz[iw]=ftrz[iw]*fsrc[iw]

    end
    vz_tmp=irfft(fvz[:],ns)
    vr_tmp=irfft(fvr[:],ns)
    trr_tmp=irfft(ftrr[:],ns)
    tpp_tmp=irfft(ftpp[:],ns)
    tzz_tmp=irfft(ftzz[:],ns)
    trz_tmp=irfft(ftrz[:],ns)

    vz_Sch[:,ir]=vz_tmp
    vr_Sch[:,ir]=vr_tmp
    trr_Sch[:,ir]=trr_tmp
    tpp_Sch[:,ir]=tpp_tmp
    tzz_Sch[:,ir]=tzz_tmp
    trz_Sch[:,ir]=trz_tmp
    end

return vr_Sch,vz_Sch,trr_Sch,tpp_Sch,tzz_Sch,trz_Sch
end


#---vel---
#---f----
function U0f(r,Vp_fluid,kf)
    Vp_fluid*kf*(-besselj(1,kf*r))
end
function W0f(r,Vp_fluid,kf,kz)
    Vp_fluid*im*kz*besselj(0,kf*r)
end
#---phi---
function U0phi(r,Vp_solid,kp)
    Vp_solid*kp*(-hankelh1(1,kp*r))
end
function W0phi(r,Vp_solid,kp,kz)
    Vp_solid*im*kz*hankelh1(0,kp*r)
end
#---xi---
function U0xi(r,Vs_solid,kb,ks,kz)
    -Vs_solid*kz*ks/kb*(-hankelh1(1,ks*r))
end
function W0xi(r,Vs_solid,kb,ks)
    Vs_solid*im*ks^2/kb*hankelh1(0,ks*r)
end

#---stress---
#---f----
function R0f(r,Vp_fluid,Rho_fluid,w,kf)
    -Rho_fluid*Vp_fluid*w^2*besselj(0,kf*r)
end
#---phi---
#CAUTION
function R0phi(r,Vp_solid,Vs_solid,Rho_solid,w,kp,kz)
    -Rho_solid*Vp_solid*(
        (w^2-2.0*Vs_solid^2*kz^2)*hankelh1(0,kp*r)+
        2.0*Vs_solid^2*kp^2/(kp*r)*(-hankelh1(1,kp*r))
        )
end
function Z0phi(r,Vp_solid,Vs_solid,Rho_solid,w,kp,kz)
    2.0*im*Rho_solid*Vp_solid*Vs_solid^2*kz*kp*(-hankelh1(1,kp*r))
end
function L0phi(r,Vp_solid,Vs_solid,Rho_solid,w,kp,kz,ka)
    -Rho_solid*Vp_solid*(
        w^2-2.0*Vs_solid^2*(ka^2-kz^2)
        )*hankelh1(0,kp*r)
end
function M0phi(r,Vp_solid,Vs_solid,Rho_solid,w,kp,kz,ka)
    -Rho_solid*Vp_solid*(
        (Vp_solid^2-2.0*Vs_solid^2)*ka^2*hankelh1(0,kp*r)-
        2.0*Vs_solid^2*kp^2/(kp*r)*(-hankelh1(1,kp*r))
        )
end
#---xi---
function R0xi(r,Vp_solid,Vs_solid,Rho_solid,w,ks,kz,kb)
    2.0*Rho_solid*Vs_solid^3*kz*ks^2/kb*(
        hankelh1(0,ks*r)-1.0/(ks*r)*hankelh1(1,ks*r)
        )
end
function Z0xi(r,Vs_solid,Rho_solid,w,ks,kz,kb)
    -im*Rho_solid*Vs_solid^3*(kz^2-ks^2)*ks/kb*(-hankelh1(1,ks*r))
end
function L0xi(r,Vs_solid,Rho_solid,w,ks,kz,kb)
    -2.0*Rho_solid*Vs_solid^3*kz/kb*(kb^2-kz^2)*hankelh1(0,ks*r)
end
function M0xi(r,Vs_solid,Rho_solid,w,ks,kz,kb)
    -2.0*Rho_solid*Vs_solid^3*kz/kb*(
        -ks^2/(ks*r)*hankelh1(1,ks*r)
        )
end
#--plane P wave--
function U0P(r,Vp_solid,kp)
    Vp_solid*kp*(-besselj(1,kp*r))
end
function W0P(r,Vp_solid,kp,kz)
    Vp_solid*im*kz*besselj(0,kp*r)
end
function R0P(r,Vp_solid,Vs_solid,Rho_solid,w,kp,kz)
    -Rho_solid*Vp_solid*(
        (w^2-2.0*Vs_solid^2*kz^2)*besselj(0,kp*r)+
        2.0*Vs_solid^2*kp^2/(kp*r)*(-besselj(1,kp*r))
        )
end
function Z0P(r,Vp_solid,Vs_solid,Rho_solid,w,kp,kz)
    2.0*im*Rho_solid*Vp_solid*Vs_solid^2*kz*kp*(-besselj(1,kp*r))
end
function L0P(r,Vp_solid,Vs_solid,Rho_solid,w,kp,kz,ka)
    -Rho_solid*Vp_solid*(
        w^2-2.0*Vs_solid^2*(ka^2-kz^2)
        )*besselj(0,kp*r)
end
function M0P(r,Vp_solid,Vs_solid,Rho_solid,w,kp,kz,ka)
    -Rho_solid*Vp_solid*(
        (Vp_solid^2-2.0*Vs_solid^2)*ka^2*besselj(0,kp*r)-
        2.0*Vs_solid^2*kp^2/(kp*r)*(-besselj(1,kp*r))
        )
end

#---module for BC solve--
function SolvePengBC01(w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)

    kz=w/Vp_solid*cos(angle_delta*pi/180)
    ka=w/Vp_solid
    kb=w/Vs_solid
    kp=sqrt(ka^2-kz^2)
    ks=sqrt(kb^2-kz^2)
    kf=sqrt((w/Vp_fluid)^2-kz^2)


#@show R0phi(rb,Vp_solid,Vs_solid,Rho_solid,w,kp,kz)
    #BC: n=0. D*[A0;B0;C0]=Ep
    D=[U0f(rb,Vp_fluid,kf) -U0phi(rb,Vp_solid,kp) -U0xi(rb,Vs_solid,kb,ks,kz);
       R0f(rb,Vp_fluid,Rho_fluid,w,kf) -R0phi(rb,Vp_solid,Vs_solid,Rho_solid,w,kp,kz) -R0xi(rb,Vp_solid,Vs_solid,Rho_solid,w,ks,kz,kb);
       0.0  -Z0phi(rb,Vp_solid,Vs_solid,Rho_solid,w,kp,kz) -Z0xi(rb,Vs_solid,Rho_solid,w,ks,kz,kb)]
    Ep=[U0P(rb,Vp_solid,kp);
        R0P(rb,Vp_solid,Vs_solid,Rho_solid,w,kp,kz);
        Z0P(rb,Vp_solid,Vs_solid,Rho_solid,w,kp,kz)]

    Coef0=inv(D)*Ep
    A0,B0,C0=Coef0[1],Coef0[2],Coef0[3]

#==Schoenberg's ABC
    cT=sqrt(Vp_fluid^2/(1+Rho_fluid*Vp_fluid^2/(Rho_solid*Vs_solid^2)))
    g=(Vs_solid/Vp_solid)^2
    KF=kf*rb
    KL=kp*rb
    tmp=Vs_solid^2/cT^2*(1-(cT/Vp_solid)^2*(cos(angle_delta*pi/180))^2)
    Coef0_approx=[-Vp_solid/Vp_fluid*(1-2.0*g*(cos(angle_delta*pi/180))^2)/tmp;
                   im*pi/4.0*(KL^2-KF^2*(1-2.0*g*(cos(angle_delta*pi/180))^2)^2/tmp);
                  -im*pi/2.0*KF^2*cos(angle_delta*pi/180)*(1-2.0*g*(cos(angle_delta*pi/180))^2)/tmp]
==#
    return A0,B0,C0
end

#---utilities--
function get_wavenumbers(w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid)
    kz=w/Vp_solid*cos(angle_delta*pi/180)
    ka=w/Vp_solid
    kb=w/Vs_solid
    kp=sqrt(ka^2-kz^2)
    ks=sqrt(kb^2-kz^2)
    kf=sqrt((w/Vp_fluid)^2-kz^2)
    return kz,ka,kb,kp,ks,kf
end
#---Potentials Evaluation--
#---f---
function ur_f(r,A0,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)
    kz,ka,kb,kp,ks,kf=get_wavenumbers(w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid)
    return U0f(r,Vp_fluid,kf)*A0
end
function uz_f(r,A0,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)
    kz,ka,kb,kp,ks,kf=get_wavenumbers(w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid)
    return W0f(r,Vp_fluid,kf,kz)*A0
end
function trr_f(r,A0,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)
    #note: -p=trr=tpp=tzz
    kz,ka,kb,kp,ks,kf=get_wavenumbers(w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid)
    return R0f(r,Vp_fluid,Rho_fluid,w,kf)*A0
end
#---phi---
function ur_phi(r,B0,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)
    kz,ka,kb,kp,ks,kf=get_wavenumbers(w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid)
    return U0phi(r,Vp_solid,kp)*B0
end
function uz_phi(r,B0,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)
    kz,ka,kb,kp,ks,kf=get_wavenumbers(w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid)
    return W0phi(r,Vp_solid,kp,kz)*B0
end
function trr_phi(r,B0,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)
    kz,ka,kb,kp,ks,kf=get_wavenumbers(w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid)
    return R0phi(r,Vp_solid,Vs_solid,Rho_solid,w,kp,kz)*B0
end
function tzz_phi(r,B0,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)
    kz,ka,kb,kp,ks,kf=get_wavenumbers(w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid)
    return L0phi(r,Vp_solid,Vs_solid,Rho_solid,w,kp,kz,ka)*B0
end
function tpp_phi(r,B0,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)
    kz,ka,kb,kp,ks,kf=get_wavenumbers(w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid)
    return M0phi(r,Vp_solid,Vs_solid,Rho_solid,w,kp,kz,ka)*B0
end
function trz_phi(r,B0,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)
    kz,ka,kb,kp,ks,kf=get_wavenumbers(w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid)
    return Z0phi(r,Vp_solid,Vs_solid,Rho_solid,w,kp,kz)*B0
end

#---xi---
function ur_xi(r,C0,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)
    kz,ka,kb,kp,ks,kf=get_wavenumbers(w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid)
    return U0xi(r,Vs_solid,kb,ks,kz)*C0
end
function uz_xi(r,C0,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)
    kz,ka,kb,kp,ks,kf=get_wavenumbers(w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid)
    return W0xi(r,Vs_solid,kb,ks)*C0
end
function trr_xi(r,C0,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)
    kz,ka,kb,kp,ks,kf=get_wavenumbers(w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid)
    return R0xi(r,Vp_solid,Vs_solid,Rho_solid,w,ks,kz,kb)*C0
end
function tzz_xi(r,C0,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)
    kz,ka,kb,kp,ks,kf=get_wavenumbers(w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid)
    return L0xi(r,Vs_solid,Rho_solid,w,ks,kz,kb)*C0
end
function tpp_xi(r,C0,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)
    kz,ka,kb,kp,ks,kf=get_wavenumbers(w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid)
    return M0xi(r,Vs_solid,Rho_solid,w,ks,kz,kb)*C0
end
function trz_xi(r,C0,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)
    kz,ka,kb,kp,ks,kf=get_wavenumbers(w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid)
    return Z0xi(r,Vs_solid,Rho_solid,w,ks,kz,kb)*C0
end
#---plane P wave---
function ur_P(r,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)
    kz,ka,kb,kp,ks,kf=get_wavenumbers(w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid)
    return U0P(r,Vp_solid,kp)
end
function uz_P(r,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)
    kz,ka,kb,kp,ks,kf=get_wavenumbers(w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid)
    return W0P(r,Vp_solid,kp,kz)
end
function trr_P(r,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)
    kz,ka,kb,kp,ks,kf=get_wavenumbers(w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid)
    return R0P(r,Vp_solid,Vs_solid,Rho_solid,w,kp,kz)
end

function tzz_P(r,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)
    kz,ka,kb,kp,ks,kf=get_wavenumbers(w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid)
    return L0P(r,Vp_solid,Vs_solid,Rho_solid,w,kp,kz,ka)
end #my eq. --deduction done. check numerically

function tpp_P(r,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)
    kz,ka,kb,kp,ks,kf=get_wavenumbers(w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid)
    return M0P(r,Vp_solid,Vs_solid,Rho_solid,w,kp,kz,ka)
end #my eq. no deduction nor numerical

function trz_P(r,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)
    kz,ka,kb,kp,ks,kf=get_wavenumbers(w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid)
    return Z0P(r,Vp_solid,Vs_solid,Rho_solid,w,kp,kz)
end

function test_plot()
    display(heatmap(heatmap(vr_Peng,title="Peng (vr)"),heatmap(vr_Sch,title="Schoenberg (vr)"),layout=(1,2)))
    display(heatmap(heatmap(vz_Peng,title="Peng (vz)"),heatmap(vz_Sch,title="Schoenberg (vz)"),layout=(1,2)))
    display(heatmap(heatmap(trr_Peng,title="Peng (trr)"),heatmap(trr_Sch,title="Schoenberg (trr)"),layout=(1,2)))
    display(heatmap(heatmap(tpp_Peng,title="Peng (tpp)"),heatmap(tpp_Sch,title="Schoenberg (tpp)"),layout=(1,2)))
    display(heatmap(heatmap(tzz_Peng,title="Peng (tzz)"),heatmap(tzz_Sch,title="Schoenberg (tzz)"),layout=(1,2)))
    display(heatmap(heatmap(trz_Peng,title="Peng (trz)"),heatmap(trz_Sch,title="Schoenberg (trz)"),layout=(1,2)))
end

#---Entry Point----
#vr_Peng,vz_Peng,trr_Peng,tpp_Peng,tzz_Peng,trz_Peng=getwaveform_Peng()
#vr_Sch,vz_Sch,trr_Sch,tpp_Sch,tzz_Sch,trz_Sch=getwaveform_Schoenberg()

#--moved from planewave
function Peng_solution03(nr,nz,dr,dz,ir_wall,tvec,src_func,srcdepth,f0,ntaper,noffset,
    Vp,Vs,Rho)
   #important: Schoenberg/Peng -> FT Aki
   #         : rfft <-> irfft -> FT Kees
   #01: not approximate, but exact solution
   #02: option to create initial stress one step ealier for pre-update
   #03: Vp,Vs,Rho as input

   println("Peng_solution02")

   dt=tvec[2]-tvec[1]
   T=tvec[end]
   ns=length(tvec)
   fvec=[0:1/T:1/T*(ns-1)]
   fvec=collect(fvec[1])

   fsrc=rfft(src_func[:])
   test_src=irfft(fsrc,ns)
   ns_half=length(fsrc)
   fvec=fvec[1:ns_half]
   plot(fvec,map(abs,fsrc),xlims=(0,1000))
   #---calculating Schoenberg's exact solution for vz, and vr, at fluid and solid phase
   #--normal incidence Pwave only (positive z propagation)
   src_iz=Int(round(srcdepth-dz/2)/dz+1) #src defined at vz cell
   rb=(ir_wall-1)*dr+dr/2.0 # make sure Vp[1:ir_wall]=1500

   max_freq_src=3*f0
   if_max=minimum(findall(x -> x > max_freq_src,fvec)) #to be used to filter high-freq noise


   Vp_fluid=Vp[src_iz,1]
   Rho_fluid=Rho[src_iz,1]
   Vp_solid=Vp[src_iz,ir_wall+1] # make sure Vp[1:ir_wall]=1500
   Vs_solid=Vs[src_iz,ir_wall+1] # make sure Vp[1:ir_wall]=1500
   Rho_solid=Rho[src_iz,ir_wall+1] # make sure Vp[1:ir_wall]=1500

   angle_delta=1E-4 # =0 normal incidence, but give a small number as it diverges due to Hankel function.

   #---checking ir_wall value
   flag_error_irwall=0
   if (Vs_solid < 1.0)
      flag_error_irwall=1
   end
   if (Vp_solid == Vp_fluid)
      flag_error_irwall=1
   end

   if (flag_error_irwall==1)
      error("Check ir_wall in Peng_solution! Is ir_wall at the boundary??")
   end

   #---output waveforms
   #--outputs are vr[ns,nr] and vz[ns,nr]
   vr_Peng=zeros(ns,nr)
   vz_Peng=zeros(ns,nr)

   vphi_Peng=zeros(ns,nr) #zero (axisymmetric)

   trr_Peng=zeros(ns,nr)
   tpp_Peng=zeros(ns,nr)
   tzz_Peng=zeros(ns,nr)
   trp_Peng=zeros(ns,nr) #zero (axisymmetric)
   trz_Peng=zeros(ns,nr)
   tpz_Peng=zeros(ns,nr) #zero (axisymmetric)

   #setting up
   #solid: vz
   fvz_solid=fsrc #independent of r
   maxamp_fvz=maximum(map(abs,fvz_solid))
   max_freq_src=3*f0
   if_max=minimum(findall(x -> x > max_freq_src,fvec)) #to be used to filter high-freq noise


   #---creating -dz/2/Vp_solid time delayed src function for later use
   #  time rewinds because vr and tii are dz/2 upwards vz, and plane wave propagates downwards
   # plus consider stress is dt/2 time delayed in main loop
   # so, this function is used only creating stress BC
   #Option:
   #Default stress is dt/2 from velocity.
   #When enabling pre-update before main loop, stress is -dt/2 from velocity

   println("===============================")
   println("Output initial stress is +dt/2 from Velocity (default)!")
   src_func_forStress=(src_func[1:end]+[src_func[2:end]; 0])/2.0
 #  println("Output initial stress is -dt/2 from Velocity!")
 #  println("Make sure to enable pre-update!")
 #  src_func_forStress=(src_func[1:end]+[0; src_func[1:end-1]])/2.0
   println("===============================")


   fsrc_forStress=rfft(src_func_forStress[:])
   tmp_phase_delay=im*2pi*fvec*(dz/2)/Vp_solid #negative phase -> time delays -> Kees FT
   fsrc_forStress_delayed=fsrc_forStress.*map(exp,tmp_phase_delay) #I said "delay" but it actually rewinds as exlained above
   fsrc_forStress_delayed[if_max+1:end]=zeros(ns_half-if_max,1)
   src_forStress_delayed=irfft(fsrc_forStress_delayed,ns)


   #--repeat every r_now
   for ir=1:nr-noffset

      #---Start frequency response
      fvz=complex(zeros(1,ns_half))
      fvr=complex(zeros(1,ns_half))

      ftrr=complex(zeros(1,ns_half))
      ftpp=complex(zeros(1,ns_half))
      ftzz=complex(zeros(1,ns_half))

      ftrz=complex(zeros(1,ns_half)) #other shear component are zero (axisymmetric)

      for iw=2:if_max
         freq=fvec[iw]
         w=2pi*freq
         #Solving B.C. to get Coefficients of potential functions
         A0,B0,C0=SolvePengBC01(w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)

        if (ir<=ir_wall)
         #fluid:vz
         r_now=(ir-1)*dr
         fvz[iw]=(-im*w)*uz_f(r_now,A0,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)
         fvz[iw]=conj(fvz[iw]) #FT Aki->FT Kees
         fvz[iw]=fvz[iw]*fsrc[iw]/w^2 #/w^2 scaling is due to V(w)=w^2 in Peng.

         #fluid:vr (with delay=-(dz/2)/Vp_solid), time reversed because vr is dz/2 upwards vz, and plane wave downwards
         r_now=(ir-1)*dr+dr/2.0
         fvr[iw]=(-im*w)*ur_f(r_now,A0,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)
         fvr[iw]=conj(fvr[iw]) #FT Aki->FT Kees
         fvr[iw]=fvr[iw]*fsrc[iw]/w^2 #/w^2 scaling is due to V(w)=w^2 in Peng.
            #additinoal delay (positive z propagation)
            tmp_phase_delay=im*2pi*fvec[iw]*(dz/2)/Vp_solid #negative phase -> time advance -> Kees FT
         fvr[iw]=fvr[iw]*exp(tmp_phase_delay)

         #fluid: pressure (with delay=-(dz/2)/Vp_solid)
         r_now=(ir-1)*dr
         ftrr[iw]=trr_f(r_now,A0,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)
         ftrr[iw]=conj(ftrr[iw]) #FT Aki->FT Kees
         ftrr[iw]=ftrr[iw]*fsrc_forStress_delayed[iw]/w^2 #/w^2 scaling is due to V(w)=w^2 in Peng.
            #additinoal delay (positive z propagation)
            tmp_phase_delay=im*2pi*fvec[iw]*(dz/2)/Vp_solid #negative phase -> time advance -> Kees FT
         ftrr[iw]=ftrr[iw]*exp(tmp_phase_delay)

         ftpp[iw]=ftrr[iw]
         ftzz[iw]=ftrr[iw]
         ftrz[iw]=0.0
      else #(solid)
         #solid:vz
         r_now=(ir-1)*dr
         fvz[iw]=(-im*w)*(
             uz_phi(r_now,B0,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)+
             uz_xi(r_now,C0,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)+
             uz_P(r_now,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)
             )
         fvz[iw]=conj(fvz[iw]) #FT Aki->FT Kees
         fvz[iw]=fvz[iw]*fsrc[iw]/w^2 #/w^2 scaling is due to V(w)=w^2 in Peng.

         #solid:vr (with delay=-(dz/2)/Vp_solid)
         r_now=(ir-1)*dr+dr/2.0
         fvr[iw]=(-im*w)*(
             ur_phi(r_now,B0,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)+
             ur_xi(r_now,C0,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)+
             ur_P(r_now,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)
             )
         fvr[iw]=conj(fvr[iw]) #FT Aki->FT Kees
         fvr[iw]=fvr[iw]*fsrc[iw]/w^2 #/w^2 scaling is due to V(w)=w^2 in Peng.
            #additinoal delay (positive z propagation)
            tmp_phase_delay=im*2pi*fvec[iw]*(dz/2)/Vp_solid #negative phase -> time advance -> Kees FT
         fvr[iw]=fvr[iw]*exp(tmp_phase_delay)


         #solid: trr,tpp,tzz (with delay=-(dz/2)/Vp_solid)
         r_now=(ir-1)*dr
         ftrr[iw]=(
             trr_phi(r_now,B0,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)+
             trr_xi(r_now,C0,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)+
             trr_P(r_now,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)
             )
         ftpp[iw]=(
             tpp_phi(r_now,B0,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)+
             tpp_xi(r_now,C0,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)+
             tpp_P(r_now,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)
             )
         ftzz[iw]=(
             tzz_phi(r_now,B0,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)+
             tzz_xi(r_now,C0,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)+
             tzz_P(r_now,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)
             )
          ftrr[iw]=conj(ftrr[iw]) #FT Aki->FT Kees
          ftpp[iw]=conj(ftpp[iw]) #FT Aki->FT Kees
          ftzz[iw]=conj(ftzz[iw]) #FT Aki->FT Kees
          ftrr[iw]=ftrr[iw]*fsrc_forStress_delayed[iw]/w^2 #/w^2 scaling is due to V(w)=w^2 in Peng.
          ftpp[iw]=ftpp[iw]*fsrc_forStress_delayed[iw]/w^2 #/w^2 scaling is due to V(w)=w^2 in Peng.
          ftzz[iw]=ftzz[iw]*fsrc_forStress_delayed[iw]/w^2 #/w^2 scaling is due to V(w)=w^2 in Peng.
            #additinoal delay (positive z propagation)
            tmp_phase_delay=im*2pi*fvec[iw]*(dz/2)/Vp_solid #negative phase -> time advance -> Kees FT
          ftrr[iw]=ftrr[iw]*exp(tmp_phase_delay)
          ftpp[iw]=ftpp[iw]*exp(tmp_phase_delay)
          ftzz[iw]=ftzz[iw]*exp(tmp_phase_delay)

          #solid: trz
          r_now=(ir-1)*dr+dr/2.0
          ftrz[iw]=(
             trz_phi(r_now,B0,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)+
             trz_xi(r_now,C0,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)+
             trz_P(r_now,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)
             )
          ftrz[iw]=conj(ftrz[iw]) #FT Aki->FT Kees
          ftrz[iw]=ftrz[iw]*fsrc_forStress_delayed[iw]/w^2 #/w^2 scaling is due to V(w)=w^2 in Peng.
       end #if ir
    end #iw
    vz_tmp=irfft(fvz[:],ns)
    vr_tmp=irfft(fvr[:],ns)
    trr_tmp=irfft(ftrr[:],ns)
    tpp_tmp=irfft(ftpp[:],ns)
    tzz_tmp=irfft(ftzz[:],ns)
    trz_tmp=irfft(ftrz[:],ns)

    vz_Peng[:,ir]=vz_tmp
    vr_Peng[:,ir]=vr_tmp
    trr_Peng[:,ir]=trr_tmp
    tpp_Peng[:,ir]=tpp_tmp
    tzz_Peng[:,ir]=tzz_tmp
    trz_Peng[:,ir]=trz_tmp
  end #ir






    #final tapering
    if (ntaper>1)
    tmp_angle=collect(range(0,pi/2,length=ntaper))
    for ir=1:ntaper
      tmp_amp=(sin(tmp_angle[ir]))^2
      vz_Peng[:,nr-ir+1-noffset]=tmp_amp*vz_Peng[:,nr-ir+1-noffset]
      vr_Peng[:,nr-ir+1-noffset]=tmp_amp*vr_Peng[:,nr-ir+1-noffset]

      trr_Peng[:,nr-ir+1-noffset]=tmp_amp*trr_Peng[:,nr-ir+1-noffset]
      tpp_Peng[:,nr-ir+1-noffset]=tmp_amp*tpp_Peng[:,nr-ir+1-noffset]
      tzz_Peng[:,nr-ir+1-noffset]=tmp_amp*tzz_Peng[:,nr-ir+1-noffset]
      trz_Peng[:,nr-ir+1-noffset]=tmp_amp*trz_Peng[:,nr-ir+1-noffset]

    end
    end


 return vz_Peng,vr_Peng,vphi_Peng,
        trr_Peng,tpp_Peng,tzz_Peng,trp_Peng,trz_Peng,tpz_Peng,
        src_iz
end




#---Schoenberg's solution as an initial condition
function Schoenberg_solution_initialBC03(nr,nz,dr,dz,ir_wall,Vp,
         tvec_vel,tvec_stress,
         srcdepth,f0,delay,
         vz_Schoenberg,vr_Schoenberg,vphi_Schoenberg,
         trr_Schoenberg,tpp_Schoenberg,tzz_Schoenberg,trp_Schoenberg,trz_Schoenberg,tpz_Schoenberg,
         ntaper,noffset)

    #03: better interpolation scheme. Requires tvec for vel and stress

   dt=tvec_vel[2]-tvec_vel[1]
   T=tvec_vel[end]
   ns=length(tvec_vel)

   #---calculating Schoenberg's approximate solution for vz, and vr, at fluid and solid phase
   #--normal incidence Pwave only (positive z propagation)
   src_iz=Int(round(srcdepth-dz/2)/dz+1) #src defined at vz cell
   rb=(ir_wall-1)*dr+dr/2.0 # make sure Vp[1:ir_wall]=1500
   Vp_solid=Vp[src_iz,ir_wall+1] # make sure Vp[1:ir_wall]=1500

   #---Create initial field condition
   vz_init=zeros(nz,nr)
   vphi_init=zeros(nz,nr)
   vr_init=zeros(nz,nr)

   trr_init=zeros(nz,nr)
   tpp_init=zeros(nz,nr)
   tzz_init=zeros(nz,nr)
   trp_init=zeros(nz,nr)
   trz_init=zeros(nz,nr)
   tpz_init=zeros(nz,nr)


   #at iz=src_iz, input XXX_Schoenberg has been assigned
   #Now I propagate this to assign initial matrices
   #Note, XXX_Schoenberg was assumed BC src, and stress values are delayed by dt/2

   #===Folling block may not be nessesary. (not enabling following block keeps time-zero to be velocity)
   #Now, if we use initial matrices, stress matrices will be accessed first in updating velocity,
   #thus stress is time rewinded by dt/2, adding zero in the begginng
   trr_Schoenberg2=[zeros(1,nr); trr_Schoenberg[1:end-1,:]]
   tpp_Schoenberg2=[zeros(1,nr); tpp_Schoenberg[1:end-1,:]]
   tzz_Schoenberg2=[zeros(1,nr); tzz_Schoenberg[1:end-1,:]]
   trp_Schoenberg2=[zeros(1,nr); trp_Schoenberg[1:end-1,:]]
   trz_Schoenberg2=[zeros(1,nr); trz_Schoenberg[1:end-1,:]]
   tpz_Schoenberg2=[zeros(1,nr); tpz_Schoenberg[1:end-1,:]]
   ==#


   #==w/o interpotaion (fast)
   is_zero=Int(round(delay/dt)+1)
   for iz=src_iz:nz
      is_tmp=Int(round((iz-src_iz)*dz/Vp_solid/dt))+is_zero
      if (is_tmp <= ns)
      for ir=1:nr
         vz_init[iz,ir]=vz_Schoenberg[is_tmp,ir];
         vr_init[iz,ir]=vr_Schoenberg[is_tmp,ir];
         vphi_init[iz,ir]=vphi_Schoenberg[is_tmp,ir];
         trr_init[iz,ir]=trr_Schoenberg2[is_tmp,ir];
         tpp_init[iz,ir]=tpp_Schoenberg2[is_tmp,ir];
         tzz_init[iz,ir]=tzz_Schoenberg2[is_tmp,ir];
         trp_init[iz,ir]=trp_Schoenberg2[is_tmp,ir];
         trz_init[iz,ir]=trz_Schoenberg2[is_tmp,ir];
         tpz_init[iz,ir]=tpz_Schoenberg2[is_tmp,ir];
      end
      end
   end

   for iz=1:src_iz-1
      is_tmp=Int(round((src_iz-iz)*dz/Vp_solid/dt))+is_zero
      if (is_tmp > 0)
      for ir=1:nr
         vz_init[iz,ir]=vz_Schoenberg[is_tmp,ir];
         vr_init[iz,ir]=vr_Schoenberg[is_tmp,ir];
         vphi_init[iz,ir]=vphi_Schoenberg[is_tmp,ir];
         trr_init[iz,ir]=trr_Schoenberg2[is_tmp,ir];
         tpp_init[iz,ir]=tpp_Schoenberg2[is_tmp,ir];
         tzz_init[iz,ir]=tzz_Schoenberg2[is_tmp,ir];
         trp_init[iz,ir]=trp_Schoenberg2[is_tmp,ir];
         trz_init[iz,ir]=trz_Schoenberg2[is_tmp,ir];
         tpz_init[iz,ir]=tpz_Schoenberg2[is_tmp,ir];
      end
      end
   end
   ==#

   #With interpolation (slow)
   itvec_tmp=zeros(1,nz)
   for iz=src_iz:nz
      it_tmp=(iz-src_iz)*dz/Vp_solid+delay
      itvec_tmp[iz]=it_tmp
   end

   for iz=1:src_iz-1
      it_tmp=(src_iz-iz)*dz/Vp_solid+delay
      itvec_tmp[iz]=it_tmp
   end

   for iz=1:nz
#       println("itvec=",itvec_tmp[iz])
       if (itvec_tmp[iz] <= 0.0)
           itvec_tmp[iz]=0.0
       elseif (itvec_tmp[iz] >= T)
           itvec_tmp[iz]=T
       end
   end

   for ir=1:nr
        vz_spl=Spline1D(tvec_vel,vz_Schoenberg[:,ir],k=1)
        vr_spl=Spline1D(tvec_vel,vr_Schoenberg[:,ir],k=1)
        vphi_spl=Spline1D(tvec_vel,vphi_Schoenberg[:,ir],k=1)
        trr_spl=Spline1D(tvec_stress,trr_Schoenberg[:,ir],k=1)
        tpp_spl=Spline1D(tvec_stress,tpp_Schoenberg[:,ir],k=1)
        tzz_spl=Spline1D(tvec_stress,tzz_Schoenberg[:,ir],k=1)
        trp_spl=Spline1D(tvec_stress,trp_Schoenberg[:,ir],k=1)
        trz_spl=Spline1D(tvec_stress,trz_Schoenberg[:,ir],k=1)
        tpz_spl=Spline1D(tvec_stress,tpz_Schoenberg[:,ir],k=1)

        vz_init[:,ir]=vz_spl(itvec_tmp[:]);
        vr_init[:,ir]=vr_spl(itvec_tmp[:]);
        vphi_init[:,ir]=vphi_spl(itvec_tmp[:]);
        trr_init[:,ir]=trr_spl(itvec_tmp[:]);
        tpp_init[:,ir]=tpp_spl(itvec_tmp[:]);
        tzz_init[:,ir]=tzz_spl(itvec_tmp[:]);
        trp_init[:,ir]=trp_spl(itvec_tmp[:]);
        trz_init[:,ir]=trz_spl(itvec_tmp[:]);
        tpz_init[:,ir]=tpz_spl(itvec_tmp[:]);
    end


 return vz_init,vr_init,vphi_init,trr_init,tpp_init,tzz_init,
    trp_init,trz_init,tpz_init
end



#---LVTS functions
function Peng_solution_LVTS01(nz,dz,nr1,dr1,nr2,dr2,
    Vp1,Vs1,Rho1,
    ir_wall,tvec,src_func,srcdepth,f0,ntaper,noffset)
   #important: Schoenberg/Peng -> FT Aki
   #         : rfft <-> irfft -> FT Kees
   #01: not approximate, but exact solution
   #LVTS01: LVTS version (2 region)
   #Assuming: only R axis is discontinous grid, Region 1 includes borehole,
   #        : Region 2 is 1D medium,
   #        : 3 times larger grid and time in Region 2

   println("Peng_solution_LVTS01")

   dt=tvec[2]-tvec[1]
   T=tvec[end]
   ns=length(tvec)
   fvec=[0:1/T:1/T*(ns-1)]
   fvec=collect(fvec[1])

   fsrc=rfft(src_func[:])
   test_src=irfft(fsrc,ns)
   ns_half=length(fsrc)
   fvec=fvec[1:ns_half]
   plot(fvec,map(abs,fsrc),xlims=(0,1000))
   #---calculating Schoenberg's exact solution for vz, and vr, at fluid and solid phase
   #--normal incidence Pwave only (positive z propagation)
   src_iz=Int(round(srcdepth-dz/2)/dz+1) #src defined at vz cell
   rb=(ir_wall-1)*dr1+dr1/2.0 # make sure Vp[1:ir_wall]=1500

   max_freq_src=3*f0
   if_max=minimum(findall(x -> x > max_freq_src,fvec)) #to be used to filter high-freq noise


   Vp_fluid=Vp1[src_iz,1]
   Rho_fluid=Rho1[src_iz,1]
   Vp_solid=Vp1[src_iz,ir_wall+1] # make sure Vp[1:ir_wall]=1500
   Vs_solid=Vs1[src_iz,ir_wall+1] # make sure Vp[1:ir_wall]=1500
   Rho_solid=Rho1[src_iz,ir_wall+1] # make sure Vp[1:ir_wall]=1500

   angle_delta=1E-4 # =0 normal incidence, but give a small number as it diverges due to Hankel function.

   #---checking ir_wall value
   flag_error_irwall=0
   if (Vs_solid < 1.0)
      flag_error_irwall=1
   end
   if (Vp_solid == Vp_fluid)
      flag_error_irwall=1
   end

   if (flag_error_irwall==1)
      error("Check ir_wall in Peng_solution! Is ir_wall at the boundary??")
   end

   #---output waveforms (Region 1)
   #--outputs are vr[ns,nr] and vz[ns,nr]
   vr1_Peng=zeros(ns,nr1)
   vz1_Peng=zeros(ns,nr1)

   vphi1_Peng=zeros(ns,nr1) #zero (axisymmetric)

   trr1_Peng=zeros(ns,nr1)
   tpp1_Peng=zeros(ns,nr1)
   tzz1_Peng=zeros(ns,nr1)
   trp1_Peng=zeros(ns,nr1) #zero (axisymmetric)
   trz1_Peng=zeros(ns,nr1)
   tpz1_Peng=zeros(ns,nr1) #zero (axisymmetric)



   #setting up
   #solid: vz
   fvz_solid=fsrc #independent of r
   maxamp_fvz=maximum(map(abs,fvz_solid))
   max_freq_src=3*f0
   if_max=minimum(findall(x -> x > max_freq_src,fvec)) #to be used to filter high-freq noise


   #---creating -dz/2/Vp_solid time delayed src function for later use
   #  time rewinds because vr and tii are dz/2 upwards vz, and plane wave propagates downwards
   # plus consider stress is dt/2 time delayed in main loop (vel as time zero)
   # so, this function is used only creating stress BC

   src_func_forStress=(src_func[1:end]+[src_func[2:end]; 0])/2.0

   #(alternative)
   #if initial values will be created, and then pre-update region 1 is activated,
   #we need stress dt/2 prior to velocity
#   src_func_forStress=([0; src_func[1:end-1]]+src_func[1:end])/2.0
#   println("======")
#   println("stress -dt/2 is activated! Do not forget pre-update region 1!")
#   println("======")


   fsrc_forStress=rfft(src_func_forStress[:])
   tmp_phase_delay=im*2pi*fvec*(dz/2)/Vp_solid #negative phase -> time delays -> Kees FT
   fsrc_forStress_delayed=fsrc_forStress.*map(exp,tmp_phase_delay) #I said "delay" but it actually rewinds as exlained above
   fsrc_forStress_delayed[if_max+1:end]=zeros(ns_half-if_max,1)
   src_forStress_delayed=irfft(fsrc_forStress_delayed,ns)


   #--repeat every r_now
   for ir=1:nr1

      #---Start frequency response
      fvz=complex(zeros(1,ns_half))
      fvr=complex(zeros(1,ns_half))

      ftrr=complex(zeros(1,ns_half))
      ftpp=complex(zeros(1,ns_half))
      ftzz=complex(zeros(1,ns_half))

      ftrz=complex(zeros(1,ns_half)) #other shear component are zero (axisymmetric)

      for iw=2:if_max
         freq=fvec[iw]
         w=2pi*freq
         #Solving B.C. to get Coefficients of potential functions
         A0,B0,C0=SolvePengBC01(w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)

        if (ir<=ir_wall)
         #fluid:vz
         r_now=(ir-1)*dr1
         fvz[iw]=(-im*w)*uz_f(r_now,A0,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)
         fvz[iw]=conj(fvz[iw]) #FT Aki->FT Kees
         fvz[iw]=fvz[iw]*fsrc[iw]/w^2 #/w^2 scaling is due to V(w)=w^2 in Peng.

         #fluid:vr (with delay=-(dz/2)/Vp_solid), time reversed because vr is dz/2 upwards vz, and plane wave downwards
         r_now=(ir-1)*dr1+dr1/2.0
         fvr[iw]=(-im*w)*ur_f(r_now,A0,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)
         fvr[iw]=conj(fvr[iw]) #FT Aki->FT Kees
         fvr[iw]=fvr[iw]*fsrc[iw]/w^2 #/w^2 scaling is due to V(w)=w^2 in Peng.
            #additinoal delay (positive z propagation)
            tmp_phase_delay=im*2pi*fvec[iw]*(dz/2)/Vp_solid #negative phase -> time advance -> Kees FT
         fvr[iw]=fvr[iw]*exp(tmp_phase_delay)

         #fluid: pressure (with delay=-(dz/2)/Vp_solid)
         r_now=(ir-1)*dr1
         ftrr[iw]=trr_f(r_now,A0,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)
         ftrr[iw]=conj(ftrr[iw]) #FT Aki->FT Kees
         ftrr[iw]=ftrr[iw]*fsrc_forStress_delayed[iw]/w^2 #/w^2 scaling is due to V(w)=w^2 in Peng.
            #additinoal delay (positive z propagation)
            tmp_phase_delay=im*2pi*fvec[iw]*(dz/2)/Vp_solid #negative phase -> time advance -> Kees FT
         ftrr[iw]=ftrr[iw]*exp(tmp_phase_delay)

         ftpp[iw]=ftrr[iw]
         ftzz[iw]=ftrr[iw]
         ftrz[iw]=0.0
      else #(solid)
         #solid:vz
         r_now=(ir-1)*dr1
         fvz[iw]=(-im*w)*(
             uz_phi(r_now,B0,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)+
             uz_xi(r_now,C0,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)+
             uz_P(r_now,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)
             )
         fvz[iw]=conj(fvz[iw]) #FT Aki->FT Kees
         fvz[iw]=fvz[iw]*fsrc[iw]/w^2 #/w^2 scaling is due to V(w)=w^2 in Peng.

         #solid:vr (with delay=-(dz/2)/Vp_solid)
         r_now=(ir-1)*dr1+dr1/2.0
         fvr[iw]=(-im*w)*(
             ur_phi(r_now,B0,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)+
             ur_xi(r_now,C0,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)+
             ur_P(r_now,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)
             )
         fvr[iw]=conj(fvr[iw]) #FT Aki->FT Kees
         fvr[iw]=fvr[iw]*fsrc[iw]/w^2 #/w^2 scaling is due to V(w)=w^2 in Peng.
            #additinoal delay (positive z propagation)
            tmp_phase_delay=im*2pi*fvec[iw]*(dz/2)/Vp_solid #negative phase -> time advance -> Kees FT
         fvr[iw]=fvr[iw]*exp(tmp_phase_delay)


         #solid: trr,tpp,tzz (with delay=-(dz/2)/Vp_solid)
         r_now=(ir-1)*dr1
         ftrr[iw]=(
             trr_phi(r_now,B0,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)+
             trr_xi(r_now,C0,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)+
             trr_P(r_now,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)
             )
         ftpp[iw]=(
             tpp_phi(r_now,B0,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)+
             tpp_xi(r_now,C0,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)+
             tpp_P(r_now,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)
             )
         ftzz[iw]=(
             tzz_phi(r_now,B0,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)+
             tzz_xi(r_now,C0,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)+
             tzz_P(r_now,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)
             )
          ftrr[iw]=conj(ftrr[iw]) #FT Aki->FT Kees
          ftpp[iw]=conj(ftpp[iw]) #FT Aki->FT Kees
          ftzz[iw]=conj(ftzz[iw]) #FT Aki->FT Kees
          ftrr[iw]=ftrr[iw]*fsrc_forStress_delayed[iw]/w^2 #/w^2 scaling is due to V(w)=w^2 in Peng.
          ftpp[iw]=ftpp[iw]*fsrc_forStress_delayed[iw]/w^2 #/w^2 scaling is due to V(w)=w^2 in Peng.
          ftzz[iw]=ftzz[iw]*fsrc_forStress_delayed[iw]/w^2 #/w^2 scaling is due to V(w)=w^2 in Peng.
            #additinoal delay (positive z propagation)
            tmp_phase_delay=im*2pi*fvec[iw]*(dz/2)/Vp_solid #negative phase -> time advance -> Kees FT
          ftrr[iw]=ftrr[iw]*exp(tmp_phase_delay)
          ftpp[iw]=ftpp[iw]*exp(tmp_phase_delay)
          ftzz[iw]=ftzz[iw]*exp(tmp_phase_delay)

          #solid: trz
          r_now=(ir-1)*dr1+dr1/2.0
          ftrz[iw]=(
             trz_phi(r_now,B0,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)+
             trz_xi(r_now,C0,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)+
             trz_P(r_now,w,angle_delta,rb,Vp_fluid,Vp_solid,Vs_solid,Rho_fluid,Rho_solid)
             )
          ftrz[iw]=conj(ftrz[iw]) #FT Aki->FT Kees
          ftrz[iw]=ftrz[iw]*fsrc_forStress_delayed[iw]/w^2 #/w^2 scaling is due to V(w)=w^2 in Peng.
       end #if ir
    end #iw
    vz_tmp=irfft(fvz[:],ns)
    vr_tmp=irfft(fvr[:],ns)
    trr_tmp=irfft(ftrr[:],ns)
    tpp_tmp=irfft(ftpp[:],ns)
    tzz_tmp=irfft(ftzz[:],ns)
    trz_tmp=irfft(ftrz[:],ns)

    vz1_Peng[:,ir]=vz_tmp
    vr1_Peng[:,ir]=vr_tmp
    trr1_Peng[:,ir]=trr_tmp
    tpp1_Peng[:,ir]=tpp_tmp
    tzz1_Peng[:,ir]=tzz_tmp
    trz1_Peng[:,ir]=trz_tmp
  end #ir



#---output waveforms (Region 2)
#Consdiering LVTS: initial velocity of Region 2 is -2dt from that of Region 1
#                : initial stress of Region 2 is -dt from that of Region 1
#Above setting will be valid only when starting immediately main loop by accesing stress (thus vel is not used)
#2nd option: When, before entering main loop, Stress (Region 2) is updated using Vel (plus using Stress for PML memPQs)
#            , then initial stress of Region 2 is -4dt from that of Region 1

vr2_Peng=zeros(ns,nr2)
vz2_Peng=zeros(ns,nr2)

vphi2_Peng=zeros(ns,nr2) #zero (axisymmetric)

trr2_Peng=zeros(ns,nr2)
tpp2_Peng=zeros(ns,nr2)
tzz2_Peng=zeros(ns,nr2)
trp2_Peng=zeros(ns,nr2) #zero (axisymmetric)
trz2_Peng=zeros(ns,nr2)
tpz2_Peng=zeros(ns,nr2) #zero (axisymmetric)

for ir=1:nr2-noffset
    vr2_Peng[1:ns-2,ir]=vr1_Peng[3:ns,nr1]
    vz2_Peng[1:ns-2,ir]=vz1_Peng[3:ns,nr1]
end

println("===============================")
println("Output initial stress is -dt from Region1!")
for ir=1:nr2-noffset
    trr2_Peng[1:ns-1,ir]=trr1_Peng[2:ns,nr1]
    tpp2_Peng[1:ns-1,ir]=tpp1_Peng[2:ns,nr1]
    tzz2_Peng[1:ns-1,ir]=tzz1_Peng[2:ns,nr1]
    trz2_Peng[1:ns-1,ir]=trz1_Peng[2:ns,nr1]
end

#==
println("Output initial stress is -4 dt from Region1!")
println("Make sure to enable pre-update (region 2)!")
for ir=1:nr2-noffset
    trr2_Peng[1:ns-4,ir]=trr1_Peng[5:ns,nr1]
    tpp2_Peng[1:ns-4,ir]=tpp1_Peng[5:ns,nr1]
    tzz2_Peng[1:ns-4,ir]=tzz1_Peng[5:ns,nr1]
    trz2_Peng[1:ns-4,ir]=trz1_Peng[5:ns,nr1]
end
==#
println("===============================")

#final tapering
if (ntaper>1)
tmp_angle=collect(range(0,pi/2,length=ntaper))
for ir=1:ntaper
  tmp_amp=(sin(tmp_angle[ir]))^2
  vz2_Peng[:,nr2-ir+1-noffset]=tmp_amp*vz2_Peng[:,nr2-ir+1-noffset]
  vr2_Peng[:,nr2-ir+1-noffset]=tmp_amp*vr2_Peng[:,nr2-ir+1-noffset]

  trr2_Peng[:,nr2-ir+1-noffset]=tmp_amp*trr2_Peng[:,nr2-ir+1-noffset]
  tpp2_Peng[:,nr2-ir+1-noffset]=tmp_amp*tpp2_Peng[:,nr2-ir+1-noffset]
  tzz2_Peng[:,nr2-ir+1-noffset]=tmp_amp*tzz2_Peng[:,nr2-ir+1-noffset]
  trz2_Peng[:,nr2-ir+1-noffset]=tmp_amp*trz2_Peng[:,nr2-ir+1-noffset]
end
end


 return vz1_Peng,vr1_Peng,vphi1_Peng,
        vz2_Peng,vr2_Peng,vphi2_Peng,
        trr1_Peng,tpp1_Peng,tzz1_Peng,trp1_Peng,trz1_Peng,tpz1_Peng,
        trr2_Peng,tpp2_Peng,tzz2_Peng,trp2_Peng,trz2_Peng,tpz2_Peng,
        src_iz
end


#---Schoenberg's solution as an initial condition
function Schoenberg_solution_initialBC_LVTS01(nz,dz,nr1,dr1,nr2,dr2,
         ir_wall,Vp1,
         tvec_vel1,tvec_stress1,
         tvec_vel2,tvec_stress2,
         srcdepth,f0,delay,
         vz1_Schoenberg,vr1_Schoenberg,vphi1_Schoenberg,
         trr1_Schoenberg,tpp1_Schoenberg,tzz1_Schoenberg,trp1_Schoenberg,trz1_Schoenberg,tpz1_Schoenberg,
         vz2_Schoenberg,vr2_Schoenberg,vphi2_Schoenberg,
         trr2_Schoenberg,tpp2_Schoenberg,tzz2_Schoenberg,trp2_Schoenberg,trz2_Schoenberg,tpz2_Schoenberg,
         ntaper,noffset)

   #LVTS01: LVTS version (2 region)

   dt=tvec_vel1[2]-tvec_vel1[1]
   T=tvec_vel1[end]
   ns=length(tvec_vel1)
   #---calculating Schoenberg's approximate solution for vz, and vr, at fluid and solid phase
   #--normal incidence Pwave only (positive z propagation)
   src_iz=Int(round(srcdepth-dz/2)/dz+1) #src defined at vz cell
   rb=(ir_wall-1)*dr1+dr1/2.0 # make sure Vp[1:ir_wall]=1500
   Vp_solid=Vp1[src_iz,ir_wall+1] # make sure Vp[1:ir_wall]=1500

   #---Create initial field condition
   vz1_init=zeros(nz,nr1)
   vphi1_init=zeros(nz,nr1)
   vr1_init=zeros(nz,nr1)

   trr1_init=zeros(nz,nr1)
   tpp1_init=zeros(nz,nr1)
   tzz1_init=zeros(nz,nr1)
   trp1_init=zeros(nz,nr1)
   trz1_init=zeros(nz,nr1)
   tpz1_init=zeros(nz,nr1)

#-
   vz2_init=zeros(nz,nr2)
   vphi2_init=zeros(nz,nr2)
   vr2_init=zeros(nz,nr2)

   trr2_init=zeros(nz,nr2)
   tpp2_init=zeros(nz,nr2)
   tzz2_init=zeros(nz,nr2)
   trp2_init=zeros(nz,nr2)
   trz2_init=zeros(nz,nr2)
   tpz2_init=zeros(nz,nr2)


   #at iz=src_iz, input XXX_Schoenberg has been assigned
   #Now I propagate this to assign initial matrices
   #Note, XXX_Schoenberg was assumed BC src, and stress values are delayed by dt/2
   #Now, if we use initial matrices, stress matrices will be accessed first in updating velocity,
   #thus stress is time rewinded by dt/2, adding zero in the begginng
   #===Folling block may not be nessesary. (not enabling following block keeps tie-zero to be velocity)
   trr_Schoenberg2=[zeros(1,nr); trr_Schoenberg[1:end-1,:]]
   tpp_Schoenberg2=[zeros(1,nr); tpp_Schoenberg[1:end-1,:]]
   tzz_Schoenberg2=[zeros(1,nr); tzz_Schoenberg[1:end-1,:]]
   trp_Schoenberg2=[zeros(1,nr); trp_Schoenberg[1:end-1,:]]
   trz_Schoenberg2=[zeros(1,nr); trz_Schoenberg[1:end-1,:]]
   tpz_Schoenberg2=[zeros(1,nr); tpz_Schoenberg[1:end-1,:]]
   ==#

   #==w/o interpotaion (fast)
   is_zero=Int(round(delay/dt)+1)
   for iz=src_iz:nz
      is_tmp=Int(round((iz-src_iz)*dz/Vp_solid/dt))+is_zero
      if (is_tmp <= ns)
      for ir=1:nr1
         vz1_init[iz,ir]=vz1_Schoenberg[is_tmp,ir];
         vr1_init[iz,ir]=vr1_Schoenberg[is_tmp,ir];
         vphi1_init[iz,ir]=vphi1_Schoenberg[is_tmp,ir];
         trr1_init[iz,ir]=trr1_Schoenberg[is_tmp,ir];
         tpp1_init[iz,ir]=tpp1_Schoenberg[is_tmp,ir];
         tzz1_init[iz,ir]=tzz1_Schoenberg[is_tmp,ir];
         trp1_init[iz,ir]=trp1_Schoenberg[is_tmp,ir];
         trz1_init[iz,ir]=trz1_Schoenberg[is_tmp,ir];
         tpz1_init[iz,ir]=tpz1_Schoenberg[is_tmp,ir];
      end

      for ir=1:nr2
         vz2_init[iz,ir]=vz2_Schoenberg[is_tmp,ir];
         vr2_init[iz,ir]=vr2_Schoenberg[is_tmp,ir];
         vphi2_init[iz,ir]=vphi2_Schoenberg[is_tmp,ir];
         trr2_init[iz,ir]=trr2_Schoenberg[is_tmp,ir];
         tpp2_init[iz,ir]=tpp2_Schoenberg[is_tmp,ir];
         tzz2_init[iz,ir]=tzz2_Schoenberg[is_tmp,ir];
         trp2_init[iz,ir]=trp2_Schoenberg[is_tmp,ir];
         trz2_init[iz,ir]=trz2_Schoenberg[is_tmp,ir];
         tpz2_init[iz,ir]=tpz2_Schoenberg[is_tmp,ir];
      end

      end
   end

   for iz=1:src_iz-1
      is_tmp=Int(round((src_iz-iz)*dz/Vp_solid/dt))+is_zero
      if (is_tmp > 0)
      for ir=1:nr1
         vz1_init[iz,ir]=vz1_Schoenberg[is_tmp,ir];
         vr1_init[iz,ir]=vr1_Schoenberg[is_tmp,ir];
         vphi1_init[iz,ir]=vphi1_Schoenberg[is_tmp,ir];
         trr1_init[iz,ir]=trr1_Schoenberg[is_tmp,ir];
         tpp1_init[iz,ir]=tpp1_Schoenberg[is_tmp,ir];
         tzz1_init[iz,ir]=tzz1_Schoenberg[is_tmp,ir];
         trp1_init[iz,ir]=trp1_Schoenberg[is_tmp,ir];
         trz1_init[iz,ir]=trz1_Schoenberg[is_tmp,ir];
         tpz1_init[iz,ir]=tpz1_Schoenberg[is_tmp,ir];
      end

      for ir=1:nr2
         vz2_init[iz,ir]=vz2_Schoenberg[is_tmp,ir];
         vr2_init[iz,ir]=vr2_Schoenberg[is_tmp,ir];
         vphi2_init[iz,ir]=vphi2_Schoenberg[is_tmp,ir];
         trr2_init[iz,ir]=trr2_Schoenberg[is_tmp,ir];
         tpp2_init[iz,ir]=tpp2_Schoenberg[is_tmp,ir];
         tzz2_init[iz,ir]=tzz2_Schoenberg[is_tmp,ir];
         trp2_init[iz,ir]=trp2_Schoenberg[is_tmp,ir];
         trz2_init[iz,ir]=trz2_Schoenberg[is_tmp,ir];
         tpz2_init[iz,ir]=tpz2_Schoenberg[is_tmp,ir];
      end

      end
   end
   ==#

   #With interpolation (slow)
   itvec_tmp=zeros(1,nz)
   for iz=src_iz:nz
      it_tmp=(iz-src_iz)*dz/Vp_solid+delay
      itvec_tmp[iz]=it_tmp
   end

   for iz=1:src_iz-1
      it_tmp=(src_iz-iz)*dz/Vp_solid+delay
      itvec_tmp[iz]=it_tmp
   end

   for iz=1:nz
#       println("itvec=",itvec_tmp[iz])
       if (itvec_tmp[iz] <= 0.0)
           itvec_tmp[iz]=0.0
       elseif (itvec_tmp[iz] >= T)
           itvec_tmp[iz]=T
       end
   end

   for ir=1:nr1
        vz_spl=Spline1D(tvec_vel1,vz1_Schoenberg[:,ir],k=2)
        vr_spl=Spline1D(tvec_vel1,vr1_Schoenberg[:,ir],k=2)
        vphi_spl=Spline1D(tvec_vel1,vphi1_Schoenberg[:,ir],k=2)
        trr_spl=Spline1D(tvec_stress1,trr1_Schoenberg[:,ir],k=2)
        tpp_spl=Spline1D(tvec_stress1,tpp1_Schoenberg[:,ir],k=2)
        tzz_spl=Spline1D(tvec_stress1,tzz1_Schoenberg[:,ir],k=2)
        trp_spl=Spline1D(tvec_stress1,trp1_Schoenberg[:,ir],k=2)
        trz_spl=Spline1D(tvec_stress1,trz1_Schoenberg[:,ir],k=2)
        tpz_spl=Spline1D(tvec_stress1,tpz1_Schoenberg[:,ir],k=2)

        vz1_init[:,ir]=vz_spl(itvec_tmp[:]);
        vr1_init[:,ir]=vr_spl(itvec_tmp[:]);
        vphi1_init[:,ir]=vphi_spl(itvec_tmp[:]);
        trr1_init[:,ir]=trr_spl(itvec_tmp[:]);
        tpp1_init[:,ir]=tpp_spl(itvec_tmp[:]);
        tzz1_init[:,ir]=tzz_spl(itvec_tmp[:]);
        trp1_init[:,ir]=trp_spl(itvec_tmp[:]);
        trz1_init[:,ir]=trz_spl(itvec_tmp[:]);
        tpz1_init[:,ir]=tpz_spl(itvec_tmp[:]);
    end

    for ir=1:nr2
        vz_spl=Spline1D(tvec_vel2,vz2_Schoenberg[:,ir],k=2)
        vr_spl=Spline1D(tvec_vel2,vr2_Schoenberg[:,ir],k=2)
        vphi_spl=Spline1D(tvec_vel2,vphi2_Schoenberg[:,ir],k=2)
        trr_spl=Spline1D(tvec_stress2,trr2_Schoenberg[:,ir],k=2)
        tpp_spl=Spline1D(tvec_stress2,tpp2_Schoenberg[:,ir],k=2)
        tzz_spl=Spline1D(tvec_stress2,tzz2_Schoenberg[:,ir],k=2)
        trp_spl=Spline1D(tvec_stress2,trp2_Schoenberg[:,ir],k=2)
        trz_spl=Spline1D(tvec_stress2,trz2_Schoenberg[:,ir],k=2)
        tpz_spl=Spline1D(tvec_stress2,tpz2_Schoenberg[:,ir],k=2)

        vz2_init[:,ir]=vz_spl(itvec_tmp[:]);
        vr2_init[:,ir]=vr_spl(itvec_tmp[:]);
        vphi2_init[:,ir]=vphi_spl(itvec_tmp[:]);
        trr2_init[:,ir]=trr_spl(itvec_tmp[:]);
        tpp2_init[:,ir]=tpp_spl(itvec_tmp[:]);
        tzz2_init[:,ir]=tzz_spl(itvec_tmp[:]);
        trp2_init[:,ir]=trp_spl(itvec_tmp[:]);
        trz2_init[:,ir]=trz_spl(itvec_tmp[:]);
        tpz2_init[:,ir]=tpz_spl(itvec_tmp[:]);
    end



    #final vertical tapering
    ntaper=10
    tmp_angle=collect(range(0,pi/2,length=ntaper))
    wavelength_half=Vp_solid/f0/2
    wl_half_iz=Int(round(wavelength_half/dz)+1)+10
    if (src_iz-wl_half_iz-ntaper<=0)
        error("srcdepth is too shallow/src wavelength is too large!")
    elseif (src_iz+wl_half_iz+ntaper>nz)
        error("srcdepth is too deep/src wavelength is too large!")
    end

    weightvec=zeros(1,nz)
    cnt=0
    for iz=src_iz-(ntaper+wl_half_iz)+1:src_iz-wl_half_iz
        cnt=cnt+1
        weightvec[iz]=(sin(tmp_angle[cnt]))^2
    end
    for iz=src_iz-wl_half_iz+1:src_iz+wl_half_iz
        weightvec[iz]=1.0
    end
    cnt=0
    for iz=src_iz+wl_half_iz+1:src_iz+wl_half_iz+ntaper
        cnt=cnt+1
        weightvec[iz]=(sin(tmp_angle[ntaper-cnt+1]))^2
    end

    for iz=1:nz
      tmp_amp=weightvec[iz]

      vz1_init[iz,:]=tmp_amp*vz1_init[iz,:];
      vphi1_init[iz,:]=tmp_amp*vphi1_init[iz,:];
      vr1_init[iz,:]=tmp_amp*vr1_init[iz,:];
      trr1_init[iz,:]=tmp_amp*trr1_init[iz,:];
      tpp1_init[iz,:]=tmp_amp*tpp1_init[iz,:];
      tzz1_init[iz,:]=tmp_amp*tzz1_init[iz,:];
      trp1_init[iz,:]=tmp_amp*trp1_init[iz,:];
      trz1_init[iz,:]=tmp_amp*trz1_init[iz,:];
      tpz1_init[iz,:]=tmp_amp*tpz1_init[iz,:];

      vz2_init[iz,:]=tmp_amp*vz2_init[iz,:];
      vphi2_init[iz,:]=tmp_amp*vphi2_init[iz,:];
      vr2_init[iz,:]=tmp_amp*vr2_init[iz,:];
      trr2_init[iz,:]=tmp_amp*trr2_init[iz,:];
      tpp2_init[iz,:]=tmp_amp*tpp2_init[iz,:];
      tzz2_init[iz,:]=tmp_amp*tzz2_init[iz,:];
      trp2_init[iz,:]=tmp_amp*trp2_init[iz,:];
      trz2_init[iz,:]=tmp_amp*trz2_init[iz,:];
      tpz2_init[iz,:]=tmp_amp*tpz2_init[iz,:];
    end


 return vz1_init,vr1_init,vphi1_init,
    vz2_init,vr2_init,vphi2_init,
    trr1_init,tpp1_init,tzz1_init,trp1_init,trz1_init,tpz1_init,
    trr2_init,tpp2_init,tzz2_init,trp2_init,trz2_init,tpz2_init
end


#for plane wave BC on Right edge (TEST. LVTS.)
function ApplyBCRight_velocity1D01!(origin_r,vr,vphi,vz,
    trr,tpp,tzz,trp,trz,tpz,
    Rho,m,nr,nz,dr,dz,dt)
    #updating velocity on Right most edge (it was Neumann BC), assuming only 1D wave propagation
    #testing for plane wave inciden ce
    #assuming LVTS, region 2 (may be this assumption is not nessesary)

    #nonzero values at j=nr are
    #j=nr->vphi,vz (1st), vr (zero)
    #Make sure regular [update] function has been filled till the right edge!


    #Following flag is 1 when no 1D assumption, but simply zeroing blank component (i.e., vr)
    flag_zero=1

if (flag_zero==1)
    j=nr
@inbounds for k=3:nz-2
    vr[k,j]=0.0
end

else

    j=nr
@inbounds for k=3:nz-2
    r_now=(j-1)*dr+dr/2.0 #for vr (Mittet)
    r_now=r_now+origin_r #for LVTS, region 2

    vr_now=vr[k,j]
    trz_f1,trz_b1=trz[k,j],trz[k-1,j]
    trz_f2,trz_b2=trz[k+1,j],trz[k-2,j]
    trp_now=trp[k,j]

    dztrz=9.0/8.0*(trz_f1-trz_b1)-1.0/24.0*(trz_f2-trz_b2)
    dztrz=dztrz/dz

    vr[k,j]=vr_now+dt/Rho[k,j]*(
                +(trr[k,j]-tpp[k,j]+m*trp_now)/r_now
                +dztrz
                )


    #==Not used as these has been filled
    r_now=(j-1)*dr #for vphi and vz (Mittet)
    r_now=r_now+origin_r #for LVTS, region 2

    vphi_now=vphi[k,j]
    tpz_f1,tpz_b1=tpz[k,j],tpz[k-1,j]
    tpz_f2,tpz_b2=tpz[k+1,j],tpz[k-2,j]
    tpp_now=tpp[k,j]

    dztpz=9.0/8.0*(tpz_f1-tpz_b1)-1.0/24.0*(tpz_f2-tpz_b2)
    dztpz=dztpz/dz

    vphi[k,j]=vphi_now+dt/Rho[k,j]*(
              +(2.0*(trp[k,j])-m*tpp_now)/(r_now) #2*trp/r-m*tpp/r
              +dztpz
              )

    vz_now=vz[k,j]
    tzz_f1,tzz_b1=tzz[k+1,j],tzz[k,j]
    tzz_f2,tzz_b2=tzz[k+2,j],tzz[k-1,j]
    tpz_now=tpz[k,j]
    rho_av=(Rho[k,j]+Rho[k+1,j])/2.0 #be careful

    dztzz=9.0/8.0*(tzz_f1-tzz_b1)-1.0/24.0*(tzz_f2-tzz_b2)
    dztzz=dztzz/dz

    vz[k,j]=vz_now+dt/rho_av*(
              +(trz[k,j]+m*tpz_now)/(r_now) #trz/r+m*tpz/r
              +dztzz
              )
    ==# #not used until here
  end
end #if
end #function


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


function ApplyBCRight_stress1D01!(origin_r,vr,vphi,vz,
    trr,tpp,tzz,trp,trz,tpz,
    lmat,mmat,m,nr,nz,dr,dz,dt)
    #updating stress on Right most edge (it was Neumann BC), assuming only 1D wave propagation
    #testing for plane wave inciden ce
    #assuming LVTS, region 2 (may be this assumption is not nessesary)

    #nonzero values @j=nr are
    #j=nr->tii,tpz (1st), trp,trz (zero)
    #Make sure regular [update] function has been filled till the right edge!

    #Following flag is 1 when no 1D assumption, but simply zeroing blank component (i.e., trp,trz)
    flag_zero=1

if (flag_zero==1)
    j=nr
@inbounds for k=3:nz-2
    trp[k,j]=0.0
    trz[k,j]=0.0
end

else

    j=nr
@inbounds for k=3:nz-2
    #==Not used as these has been filled
    r_now=(j-1)*dr #for trr,tpp,tzz,tpz Mittet
    r_now=r_now+origin_r #for LVTS, region 2

    trr_now=trr[k,j]
    vz_f1,vz_b1=vz[k,j],vz[k-1,j]
    vz_f2,vz_b2=vz[k+1,j],vz[k-2,j]
    vphi_now=vphi[k,j]
    lambda=lmat[k,j] #be careful
    mu=mmat[k,j] #be careful
    dzvz=9.0/8.0*(vz_f1-vz_b1)-1.0/24.0*(vz_f2-vz_b2)
    dzvz=dzvz/dz
    trr[k,j]=trr_now+dt*(
              +lambda*(vr[k,j]+m*vphi_now)/(r_now) #l*vr/r+l*m*vphi/r
              +lambda*dzvz #l*dvzdz
              )

    tpp_now=tpp[k,j]
    tpp[k,j]=tpp_now+dt*(
              +(lambda+2.0*mu)*(vr[k,j]+m*vphi_now)/(r_now) #(l+2mu)*(vr+m*vphi)/r
              +lambda*dzvz #l*dvzdz
              )

    tzz_now=tzz[k,j]
    tzz[k,j]=tzz_now+dt*(
              +(lambda+2.0*mu)*dzvz #(l+2mu)*dzvz
              +lambda*(vr[k,j]+m*vphi_now)/(r_now) #l*(vr+m*vphi)/r
              )

    tpz_now=tpz[k,j]
    if (mmat[k,j]+mmat[k+1,j]==0.0)
        mu_av=0.0
    else
        mu_av=2.0*mmat[k,j]*mmat[k+1,j]/(mmat[k,j]+mmat[k+1,j]) #be careful
    end
    vphi_f1,vphi_b1=vphi[k+1,j],vphi[k,j]
    vphi_f2,vphi_b2=vphi[k+2,j],vphi[k-1,j]
    vz_now=vz[k,j]
    tpz[k,j]=update_tpz_2nd(tpz_now,vphi_f1,vphi_f2,vphi_b1,vphi_b2,
    vz_now,dz,dt,r_now,mu_av,m)
    ==# #not used until here

    r_now=(j-1)*dr+dr/2.0 #for trp,trz Mittet
    r_now=r_now+origin_r #for LVTS, region 2

    trp_now=trp[k,j]
    vr_now=vr[k,j]
    trp[k,j]=trp_now+dt*(
            -mmat[k,j]*(vphi[k,j]+m*vr_now)/(r_now) #-mu*(vphi+m*vr)/r
            )
    trz_now=trz[k,j]

    vr_f1,vr_b1=vr[k+1,j],vr[k,j]
    vr_f2,vr_b2=vr[k+2,j],vr[k-1,j]
    dzvr=9.0/8.0*(vr_f1-vr_b1)-1.0/24.0*(vr_f2-vr_b2)
    dzvr=dzvr/dz
    trz[k,j]=trz_now+dt*(
               +mmat[k,j]*dzvr #mu*dzvr
               )
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



#--
function taper_initial_field!(field,nz,nr,dz,src_iz,wavelength,ntaper)
    #tapering initial field
    #field[nz,nr]
    #flat weight 1 between src_iz-W:src_iz+W
    #ntaper after then

    #final vertical tapering
    #ntaper=10 :given
    tmp_angle=collect(range(0,pi/2,length=ntaper))
    #wavelength_half=Vp_solid/f0/2
    wavelength_half=wavelength/2
    wl_half_iz=Int(round(wavelength_half/dz)+1)+10
    if (src_iz-wl_half_iz-ntaper<=0)
        error("srcdepth is too shallow/src wavelength is too large!")
    elseif (src_iz+wl_half_iz+ntaper>nz)
        error("srcdepth is too deep/src wavelength is too large!")
    end

    weightvec=zeros(1,nz)
    cnt=0
    for iz=src_iz-(ntaper+wl_half_iz)+1:src_iz-wl_half_iz
        cnt=cnt+1
        weightvec[iz]=(sin(tmp_angle[cnt]))^2
    end
    for iz=src_iz-wl_half_iz+1:src_iz+wl_half_iz
        weightvec[iz]=1.0
    end
    cnt=0
    for iz=src_iz+wl_half_iz+1:src_iz+wl_half_iz+ntaper
        cnt=cnt+1
        weightvec[iz]=(sin(tmp_angle[ntaper-cnt+1]))^2
    end

    for iz=1:nz
      tmp_amp=weightvec[iz]
      field[iz,:]=tmp_amp*field[iz,:];
    end

    println("Initial field tapered until iz=",src_iz+wl_half_iz+ntaper)
end #function

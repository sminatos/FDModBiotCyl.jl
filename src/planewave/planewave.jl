# Setting up initial conditions for normal incidence plane wave
# Analytical plane wave solutions are based on Peng (1996)
# dependent on: rfft/irfft (FFTW)
#             : Spline1D (Dierckx)
#             : besselj, hankelh1 (SpecialFunctions)


#-------------------------------
# tapering initial field
#-------------------------------
function taper_initial_field!(field,nz,nr,dz,src_iz,wavelength,ntaper)
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

#-------------------------------
# Peng's analytical solution of borehole pressure and
# elastic field in surrounding media
#-------------------------------

function Peng_solution03(nr,nz,dr,dz,ir_wall,tvec,src_func,srcdepth,f0,ntaper,noffset,
    Vp,Vs,Rho)
   #important: Schoenberg/Peng -> FT Aki
   #         : rfft <-> irfft -> FT Kees
   #01: not approximate, but exact solution
   #02: option to create initial stress one step ealier for pre-update
   #03: Vp,Vs,Rho as input

   #println("Peng_solution02")

   dt=tvec[2]-tvec[1]
   T=tvec[end]
   ns=length(tvec)
   fvec=[0:1/T:1/T*(ns-1)]
   fvec=collect(fvec[1])

   fsrc=rfft(src_func[:])
   test_src=irfft(fsrc,ns)
   ns_half=length(fsrc)
   fvec=fvec[1:ns_half]
#   plot(fvec,map(abs,fsrc),xlims=(0,1000))
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



#-------------------------------
# Converting time-domain solution into
# Space-domain to construct initial conditions
#-------------------------------
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



#------------------------
# Solving BC problem in Peng
#------------------------
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


#------Peng's solutions----
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





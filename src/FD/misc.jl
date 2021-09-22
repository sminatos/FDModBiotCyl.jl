#---------------
# copying matrix
#---------------
function mycopy_mat(MATin,MATout,nz,nr)
    @inbounds Threads.@threads for j=1:nr
        for k=1:nz
            MATout[k,j]=MATin[k,j]
        end
    end
end


#------------------------------
# finding acoustic or elastic media,
# and assigning flags where
# vf should be zero
#------------------------------
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


#------------------------------
# Checking stability of FD modeling
#------------------------------
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





















#------------------------------------------------------------
# Following are not used, but kept for future extention
#------------------------------------------------------------
#--calculate divergence of pariticle velocity (assuming axisymetric)
function wave_div01!(divV,nr,nz,dr,dz,vr,vz)
   #div(u)=1/r*dur/dr+duz/dz: location of vphi,tii (Mittet)
#@inbounds Threads.@threads
   for j=3:nr-2
      for k=3:nz-2
         r_now=(j-1)*dr #position of tii (nonaxis to avoid singularity)
         vr_f1,vr_b1=vr[k,j],vr[k,j-1]
         vr_f2,vr_b2=vr[k,j+1],vr[k,j-2]
         vz_f1,vz_b1=vz[k,j],vz[k-1,j]
         vz_f2,vz_b2=vz[k+1,j],vz[k-2,j]
         drvr=9.0/8.0*(vr_f1-vr_b1)-1.0/24.0*(vr_f2-vr_b2)
         drvr=drvr/dr
         dzvz=9.0/8.0*(vz_f1-vz_b1)-1.0/24.0*(vz_f2-vz_b2)
         dzvz=dzvz/dz
         divV[k,j]=1.0/r_now*drvr+dzvz
      end
   end
end

#--calculate curl of pariticle velocity (assuming axisymetric)
function wave_curl01!(curV,nr,nz,dr,dz,vr,vz)
   #curl(u)=dur/dz-duz/dr: location of trz (Mittet)
#@inbounds Threads.@threads
   for j=3:nr-2
      for k=3:nz-2
         r_now=(j-1)*dr+dr/2.0 #position of trz
         vr_f1,vr_b1=vr[k+1,j],vr[k,j]
         vr_f2,vr_b2=vr[k+2,j],vr[k-1,j]
         vz_f1,vz_b1=vz[k,j+1],vz[k,j]
         vz_f2,vz_b2=vz[k,j+2],vz[k,j-1]
         dzvr=9.0/8.0*(vr_f1-vr_b1)-1.0/24.0*(vr_f2-vr_b2)
         dzvr=dzvr/dz
         drvz=9.0/8.0*(vz_f1-vz_b1)-1.0/24.0*(vz_f2-vz_b2)
         drvz=drvz/dr
         curV[k,j]=dzvr-drvz
      end
   end
end

function init_snap_sparse(nt,nz,nr)
   #---Initializing snapshots (sparse version)
   nskip=100
   nr_sp=0
   for ir=1:10:nr
      nr_sp=nr_sp+1
   end
   nz_sp=0
   for iz=1:10:nz
      nz_sp=nz_sp+1
   end

   snapshots_vr,snapshots_vz,nsnap,itvec_snap=init_snapshots_v(nskip,nt,nz_sp,nr_sp)
   snapshots_trr,nsnap,itvec_snap=init_snapshots_t(nskip,nt,nz_sp,nr_sp)
   return nskip,snapshots_trr,snapshots_vr,snapshots_vz,nsnap,itvec_snap
end

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

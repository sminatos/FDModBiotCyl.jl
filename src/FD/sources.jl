
#------------------------------
# creating source wavelets
#------------------------------
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


#------------------------------
# converting src geometry to
# index numbers of a pressure grid
#------------------------------
function get_srcindex_monopole(srcgeom,dr,dz)
    #src index and amplitudes
    #--point input
    izsrc_ref=Int(ceil(srcgeom[1,1]/dz)+1)
    irsrc_ref=Int(ceil(srcgeom[1,2]/dr)+1)
    src_index=[izsrc_ref;irsrc_ref]
    src_dn=[1.]
    return src_index,src_dn
end

#------------------------------
# applying src term into grids
#------------------------------
function srcapply!(paramarray,src_index,src_dn,srcamp)
    #--e.g. fz src: vz=vz+fz (same form as vz update, thus fz is as-is and not time-derivatives)
    #--input array can be vx,vz etc
    for isrcind=1:length(src_dn)
       isz=src_index[1,isrcind]
       isr=src_index[2,isrcind]
       paramarray[isz,isr]=paramarray[isz,isr]+srcamp*src_dn[isrcind]
    end
end


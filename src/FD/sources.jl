
#------------------------------
# creating source wavelets
#------------------------------
function myricker2(tvec,f0,delay,derivative_number)
    # Gaussian=exp(-sigma*(t-t0)^2), sigma=-pi^2*f0^2. Mittet (eq 25), 1996, Geophysics
    nt=length(tvec)
    y=zeros(1,nt)
    if (derivative_number==2)
        println("Myricker: 2nd derivative of Gaussian (Ricker with a negative sign)")
    #---2nd derivative Gaussian (Ricker with a negative sign)
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

#-----------------------------------
# converting src geometry to
# index numbers of a pressure grid
# Rubust point src by Gaussian shape
#-----------------------------------
function get_srcindex_pGauss(srcgeom,dr,dz,windowsize,sigma)
    #window size->odd
    #src index and amplitudes
    #Correct scaling factor (as a dirac delta in cylindrical coordinate) only when source is located around r=0
    #--point input
    izsrc_ref=Int(ceil(srcgeom[1,1]/dz)+1)
    irsrc_ref=Int(ceil(srcgeom[1,2]/dr)+1)

    Gauss_Wmat=Gaussian_filter(windowsize,sigma)

    halfw=Int(round((windowsize-1)/2))

    #check if gaussian window exceeds at left boundary
    if(irsrc_ref-halfw < 1)
        irvec_1row=collect(range(irsrc_ref-halfw,irsrc_ref+halfw,length=windowsize))
        ir1=findall(x->x==1,irvec_1row)

        izvec_1col=collect(range(izsrc_ref-halfw,izsrc_ref+halfw,length=windowsize))
        izvec_mcol=repeat(izvec_1col,outer=windowsize-ir1[1]+1)

        irvec_1row=collect(range(1,irsrc_ref+halfw,length=windowsize-ir1[1]+1))
        irvec_mcol=repeat(irvec_1row,inner=windowsize)

        src_index=zeros(2,windowsize*(windowsize-ir1[1]+1))
        src_index[1,:]=izvec_mcol[:]
        src_index[2,:]=irvec_mcol[:]

        #extracting symmetric components around r=0
        Wmat=Gauss_Wmat[:,ir1[1]:windowsize]
        #evaluating volume integral to get scaing factor
        V=Check_VolumeIntegral(Wmat,dr,dz)
        #saling-corrected Gaussian
        src_dn=Wmat/V

        src_dn=src_dn[:]
    else

        src_index=zeros(2,windowsize^2)

        izvec_1col=collect(range(izsrc_ref-halfw,izsrc_ref+halfw,length=windowsize))
        izvec_mcol=repeat(izvec_1col,outer=windowsize)

        irvec_1row=collect(range(irsrc_ref-halfw,irsrc_ref+halfw,length=windowsize))
        irvec_mcol=repeat(irvec_1row,inner=windowsize)

        src_index[1,:]=izvec_mcol[:]
        src_index[2,:]=irvec_mcol[:]

        #Note: scaling factor does not represent a dirac delta function
        src_dn=Gauss_Wmat[:]
    end


    src_index=map(x->Int(x),src_index)
    return src_index,src_dn
end

#Gaussian filter
function Gaussian_filter(wsize, sigma)
#wsize should be odd number
    gauss=zeros(wsize,wsize); #2D filter matrix
    for i=-round((wsize-1)/2):1:round((wsize-1)/2)
        for j=-round((wsize-1)/2):1:round((wsize-1)/2)
            i0=(wsize+1)/2;
            j0=(wsize+1)/2;
            i1=Int(round(i+i0));
            j1=Int(round(j+j0));
            gauss[j1,i1]=exp(-((i1-i0)^2+(j1-j0)^2)/(2*sigma^2));
        end
    end
gauss=gauss/sum(gauss[:]);
return gauss
end


function Check_VolumeIntegral(Wmat,dr,dz)
    #evaluate volume integral cylindrical coordinate
    #V=2pi int int r dr dz f(r,z)
    #assumes (Gaussian) matrices symmetrical around r=0
    # in case of a scalar at r=0, it results in pi*dz*dr^2*1/4
    wsize_z,wsize_r=size(Wmat)

    V=0.
    for iz=1:wsize_z
        tmpvec_r=Wmat[iz,:]
        for ir=1:wsize_r
            if (ir==1)
                V=V+pi*dz*tmpvec_r[ir]*dr^2/4
            else
                V=V+pi*dz*tmpvec_r[ir]*dr^2/4*(ir-1)*8
            end
        end
    end

    return V
end

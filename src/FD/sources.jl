
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

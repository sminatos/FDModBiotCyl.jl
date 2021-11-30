module FDModBiotCyl

    # the dependency of this module
    using FFTW, #rff,irfft (planewave.jl)
    Dierckx, #Spline1D (planewave.jl)
    SpecialFunctions #besselj,hankelh1 (planewave.jl)

    include("./FD/fdtd_export.jl")
    include("./planewave/planewave_export.jl")

end




#using Revise
#]activate .

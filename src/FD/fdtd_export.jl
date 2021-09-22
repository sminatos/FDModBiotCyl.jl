      #---FD main (fdmain.jl)----
export  init_fields_Por,
        update_velocity_1st_Por!,
        update_stress_1st_Por!,
        ApplyBCLeft_velocity_1st_Por!,
        ApplyBCLeft_stress_1st_Por!,
        ApplyBCRight_velocity1D_Por01!,
        ApplyBCRight_stress1D_Por01!, 
        ApplyBC_stress_AcousticMedia_TEST!,
        ApplyBC_stress_ElasticMedia_Ou!,
   #---Additional consideration at borehole wall (boreholewall.jl)
        update_vr_vfr_1st_vertical,
   #---PML (pml.jl)---------
           init_PML_Por,
           init_PML_profile,
           init_memory_variables_Por,
           PML_save_vel_Top_Por!,
           PML_save_vel_Bottom_Por!,
           PML_save_vel_Right_Por!,
           PML_save_stress_Top_Por!,
           PML_save_stress_Bottom_Por!,
           PML_save_stress_Right_Por!,
           PML_update_velocity_1st_Top_Por!,
           PML_update_velocity_1st_Bottom_Por!,
           PML_update_velocity_1st_Right_Por!,
           PML_update_velocity_1st_TopRight_Por!,
           PML_update_velocity_1st_BottomRight_Por!,
           PML_update_stress_1st_Top_Por!,
           PML_update_stress_1st_Bottom_Por!,
           PML_update_stress_1st_Right_Por!,
           PML_update_stress_1st_TopRight_Por!,
           PML_update_stress_1st_BottomRight_Por!,
           PML_update_memRS_1st_Top_Por!,
           PML_update_memRS_1st_Bottom_Por!,
           PML_update_memRS_1st_Right_Por!,
           PML_update_memRS_1st_TopRight_Por!,
           PML_update_memRS_1st_BottomRight_Por!,
           PML_update_memPQ_1st_Top_Por!,
           PML_update_memPQ_1st_Bottom_Por!,
           PML_update_memPQ_1st_Right_Por!,
           PML_update_memPQ_1st_TopRight_Por!,
           PML_update_memPQ_1st_BottomRight_Por!,
           ApplyBCLeft_velocity_1st_atPML_Top_Por!,
           ApplyBCLeft_velocity_1st_atPML_Bottom_Por!,
           ApplyBCLeft_stress_1st_atPML_Top_Por!,
           ApplyBCLeft_stress_1st_atPML_Bottom_Por!,
       #---snapshots (snapshots.jl)----
          init_snapshots_v,
          init_snapshots_t,
          #init_snap, #-->planewave_3L.jl
          get_snapshots!,
          get_snapshots_t!,
      #---receivers (receivers.jl)---
          init_receiver,
          getRecData_from_index!,
      #---other functions (misc.jl)---------
        mycopy_mat,
      get_Flag_vf_zero,
      check_stability01,
      #---sources (sources.jl)-----
      myricker2






include("./fdmain.jl")
include("./pml.jl")
include("./boreholewall.jl")
include("./receivers.jl")
include("./sources.jl")
include("./snapshots.jl")
include("./misc.jl")

#include("./PCylFDMod_Ou1st03.jl") #main modules. With acoustic, not good, but ok.
#include("./FD/Peng06Mod.jl") #Peng's exact solution to B.C.s of plane wave (MOD version. Use this.)

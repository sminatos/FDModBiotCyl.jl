#Additional considerations at the borehole wall


#--------------------------------------
# Activating vfr@vertical boundary
# separating Fluid-PoroElastic media
#--------------------------------------
function update_vr_vfr_1st_vertical(vr,trr,tpp,tzz,trz,vfr,pf,
    vr_org,vfr_org,
    Rho,Rhof,D1,D2,nr,nz,dr,dz,dt,LPML_z,
    Prz_T,Prz_B,
    ir_wall,Flag_AC,Flag_E)
    #Acoustic-Poroelastic vertical boundary
    #Calculating vrw and vr at the boundary

    j=ir_wall #wall
    #main region
    @inbounds Threads.@threads for k=LPML_z+1:nz-LPML_z

        #check if I need to activate vfr
        if(Flag_AC[k,j]==1 && Flag_AC[k,j+1]==0 && Flag_E[k,j+1]==0)#
             #I am accessing [k,j], but the definition of r_now changes with each component!
             r_now=(j-1)*dr+dr/2.0 #for vr (Mittet)
             vr_now=vr_org[k,j]
             trr_f1,trr_b1=trr[k,j+1],trr[k,j]
             trz_f1,trz_b1=trz[k,j],trz[k-1,j] #should be zero

             trr_av=0.5*(trr_f1+trr_b1) #this looks ok according to BC (pf=-trr, now trr_b1 is -pf)
             tpp_av=0.5*(tpp[k,j+1]+tpp[k,j]) #see trr.
#             tpp_av=tpp[k,j+1] #test. assuming tpp is discontinuous

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


         end #if AC-PE BC
    end #k

#
    #Top PML
    @inbounds Threads.@threads for k=2:LPML_z
        #check if I need to activate vfr
        if(Flag_AC[k,j]==1 && Flag_AC[k,j+1]==0 && Flag_E[k,j+1]==0)#
              #I am accessing [k,j], but the definition of r_now changes with each component!
             r_now=(j-1)*dr+dr/2.0 #for vr (Mittet)
             vr_now=vr_org[k,j]
             trr_f1,trr_b1=trr[k,j+1],trr[k,j]
             trz_f1,trz_b1=trz[k,j],trz[k-1,j] #should be zero

             trr_av=0.5*(trr_f1+trr_b1) #this looks ok according to BC (pf=-trr, now trr_b1 is -pf)
             tpp_av=0.5*(tpp[k,j+1]+tpp[k,j]) #see trr.
#             tpp_av=tpp[k,j+1] #test. assuming tpp is discontinuous

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
              #I am always accessing [k,j], but the definition of r_now changes with components!
             r_now=(j-1)*dr+dr/2.0 #for vr (Mittet)
             vr_now=vr_org[k,j]
             trr_f1,trr_b1=trr[k,j+1],trr[k,j]
             trz_f1,trz_b1=trz[k,j],trz[k-1,j] #should be zero

             trr_av=0.5*(trr_f1+trr_b1) #this looks ok according to BC (pf=-trr, now trr_b1 is -pf)
             tpp_av=0.5*(tpp[k,j+1]+tpp[k,j]) #see trr.
#             tpp_av=tpp[k,j+1] #test. assuming tpp is discontinuous

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

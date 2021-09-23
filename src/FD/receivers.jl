#-------------------------
# initializing receiver field
#-------------------------
function init_receiver(recgeom,nrec,dr,dz,nr,nz,nt)
    #vr,vz, and tii receivers
    #receivers are assigned at nearest grids

    if (minimum(recgeom[:,1])<0.) || (maximum(recgeom[:,1])>(nz-1)*dz)
       error("Check receiver geometry!")
    end
    if (minimum(recgeom[:,2])<0.) || (maximum(recgeom[:,2])>(nr-1)*dr)
       error("Check receiver geometry!")
    end

    #--Receiver indexes---
    index_allrec_vr=zeros(nrec,2) #z,r
    index_allrec_vz=zeros(nrec,2) #z,r
    index_allrec_tii=zeros(nrec,2) #z,r
    for irec=1:nrec
        tmpz=recgeom[irec,1]
        tmpr=recgeom[irec,2]

        iz=round((tmpz-dz/2)/dz)+1 #vz
        ir=round(tmpr/dr)+1 #vz
        index_allrec_vz[irec,1]=iz
        index_allrec_vz[irec,2]=ir

        iz=round(tmpz/dz)+1 #vr
        ir=round((tmpr-dr/2)/dr)+1 #vr
        index_allrec_vr[irec,1]=iz
        index_allrec_vr[irec,2]=ir

#        iz=round((tmpz-dz/2)/dz)+1 #tii
#        ir=round((tmpr-dr/2)/dr)+1 #tii
        iz=round(tmpz/dz)+1 #tii
        ir=round(tmpr/dr)+1 #tii
        index_allrec_tii[irec,1]=iz
        index_allrec_tii[irec,2]=ir

    end
    #---Initializing receiver data matrces
    rec_vr=zeros(nt,nrec)
    rec_vz=zeros(nt,nrec)
    rec_tii=zeros(nt,nrec)

    return rec_vr,rec_vz,rec_tii,ceil.(Int,index_allrec_vr),ceil.(Int,index_allrec_vz),ceil.(Int,index_allrec_tii)
end

#
function init_receiver_hydrophone(recgeom,nrec,dr,dz,nr,nz,nt)
    #tii(pressure) receivers
    #receivers are assigned at nearest grids

    if (minimum(recgeom[:,1])<0.) || (maximum(recgeom[:,1])>(nz-1)*dz)
       error("Check receiver geometry!")
    end
    if (minimum(recgeom[:,2])<0.) || (maximum(recgeom[:,2])>(nr-1)*dr)
       error("Check receiver geometry!")
    end

    #--Receiver indexes---
    index_allrec_tii=zeros(nrec,2) #z,r
    for irec=1:nrec
        tmpz=recgeom[irec,1]
        tmpr=recgeom[irec,2]

#        iz=round((tmpz-dz/2)/dz)+1 #tii
#        ir=round((tmpr-dr/2)/dr)+1 #tii
        iz=round(tmpz/dz)+1 #tii
        ir=round(tmpr/dr)+1 #tii
        index_allrec_tii[irec,1]=iz
        index_allrec_tii[irec,2]=ir

    end
    #---Initializing receiver data matrces
    rec_tii=zeros(nt,nrec)

    return rec_tii,ceil.(Int,index_allrec_tii)
end

function init_receiver_geophone(recgeom,nrec,dr,dz,nr,nz,nt)
    #vr,and vz receivers
    #receivers are assigned at nearest grids
    #borehole wall is not checked!

    if (minimum(recgeom[:,1])<0.) || (maximum(recgeom[:,1])>(nz-1)*dz)
       error("Check receiver geometry!")
    end
    if (minimum(recgeom[:,2])<0.) || (maximum(recgeom[:,2])>(nr-1)*dr)
       error("Check receiver geometry!")
    end

    #--Receiver indexes---
    index_allrec_vr=zeros(nrec,2) #z,r
    index_allrec_vz=zeros(nrec,2) #z,r
    index_allrec_tii=zeros(nrec,2) #z,r
    for irec=1:nrec
        tmpz=recgeom[irec,1]
        tmpr=recgeom[irec,2]

        iz=round((tmpz-dz/2)/dz)+1 #vz
        ir=round(tmpr/dr)+1 #vz
        index_allrec_vz[irec,1]=iz
        index_allrec_vz[irec,2]=ir

        iz=round(tmpz/dz)+1 #vr
        ir=round((tmpr-dr/2)/dr)+1 #vr
        index_allrec_vr[irec,1]=iz
        index_allrec_vr[irec,2]=ir

    end
    #---Initializing receiver data matrces
    rec_vr=zeros(nt,nrec)
    rec_vz=zeros(nt,nrec)

    return rec_vr,rec_vz,ceil.(Int,index_allrec_vr),ceil.(Int,index_allrec_vz)
end


#-------------------------
# obtaining values at receiver grid
#-------------------------
function getRecData_from_index!(paramarray,rec_data,rec_index,nrec,it)
    #receivers are assigned at nearest grids

    for irec=1:nrec
        iz=rec_index[irec,1]
        ir=rec_index[irec,2]
        rec_data[it,irec]=paramarray[iz,ir]
    end
end


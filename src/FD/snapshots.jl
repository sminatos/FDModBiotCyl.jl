
#-----------------------------
# initializing snapshots field
#-----------------------------
function init_snapshots_v(nskip,nt,nz,nr)
    #---Snapshots (velocity)
    #nskip=100
    itvec_snap=zeros(Float64,Int(ceil(nt/nskip))-1)
    tvec_snap=zeros(Float64,Int(ceil(nt/nskip))-1)

    cnt=1
    for itmp=1:Int(ceil(nt/nskip))-1
       cnt+=nskip
       itvec_snap[itmp]=cnt
    end

    if (maximum(itvec_snap)>nt)
       println("Error in snapshots settings!")
       return
    else
       nsnap=length(itvec_snap)
    end

    snapshots_vr=zeros(nz,nr,nsnap)
    snapshots_vz=zeros(nz,nr,nsnap)
    return snapshots_vr,snapshots_vz,nsnap,itvec_snap
end

function init_snapshots_t(nskip,nt,nz,nr)
    #---Snapshots (stress)
    #nskip=100
    itvec_snap=zeros(Float64,Int(ceil(nt/nskip))-1)
    tvec_snap=zeros(Float64,Int(ceil(nt/nskip))-1)

    cnt=1
    for itmp=1:Int(ceil(nt/nskip))-1
       cnt+=nskip
       itvec_snap[itmp]=cnt
    end

    if (maximum(itvec_snap)>nt)
       println("Error in snapshots settings!")
       return
    else
       nsnap=length(itvec_snap)
    end

    snapshots_trr=zeros(nz,nr,nsnap)
    return snapshots_trr,nsnap,itvec_snap
end

#----------------------------
# obtaining snapshots
#----------------------------

function get_snapshots!(snapshots_vr,snapshots_vz,isnap,vrmat,vzmat,nr,nz)
@inbounds Threads.@threads for ir=1:nr
        for iz=1:nz
            snapshots_vr[iz,ir,isnap]=vrmat[iz,ir]
            snapshots_vz[iz,ir,isnap]=vzmat[iz,ir]
        end
    end
end

function get_snapshots_t!(snapshots_trr,isnap,trrmat,nr,nz)
@inbounds Threads.@threads for ir=1:nr
        for iz=1:nz
            snapshots_trr[iz,ir,isnap]=trrmat[iz,ir]
        end
    end
end

#--Snapshots (sparse version)
function get_snapshots_sp!(snapshots_vr,snapshots_vz,isnap,vrmat,vzmat,nr,nz)

ir_sp=0
@inbounds for ir=1:10:nr
        ir_sp=ir_sp+1
        iz_sp=0
        for iz=1:10:nz
            iz_sp=iz_sp+1
            snapshots_vr[iz_sp,ir_sp,isnap]=vrmat[iz,ir]
            snapshots_vz[iz_sp,ir_sp,isnap]=vzmat[iz,ir]
        end
    end
end

function get_snapshots_sp_t!(snapshots_trr,isnap,trrmat,nr,nz)
ir_sp=0
@inbounds for ir=1:10:nr
        ir_sp=ir_sp+1
        iz_sp=0
        for iz=1:10:nz
            iz_sp=iz_sp+1
            snapshots_trr[iz_sp,ir_sp,isnap]=trrmat[iz,ir]
        end
    end
end
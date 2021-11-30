
#-----------------------------
# initializing snapshots field
#-----------------------------
function init_snap(nt,nz,nr,nskip)
   #Initializing snapshots
   snapshots_vr,snapshots_vz,nsnap,itvec_snap=init_snapshots_v(nskip,nt,nz,nr)
   snapshots_trr,nsnap,itvec_snap=init_snapshots_t(nskip,nt,nz,nr)
   return nskip,snapshots_trr,snapshots_vr,snapshots_vz,nsnap,itvec_snap
end

function init_snapshots(nt,nz,nr,nskip)
    #---Snapshots (single field variable)
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

    snapshots_field=zeros(nz,nr,nsnap)
    return snapshots_field,nsnap,itvec_snap
end

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

function get_snapshots!(snapshots_field,isnap,field,nr,nz)
@inbounds Threads.@threads for ir=1:nr
        for iz=1:nz
            snapshots_field[iz,ir,isnap]=field[iz,ir]
        end
    end
end

#--Snapshots (sparse version)
function get_snapshots_sp!(snapshots_field,isnap,field,nr,nz)

ir_sp=0
@inbounds for ir=1:10:nr
        ir_sp=ir_sp+1
        iz_sp=0
        for iz=1:10:nz
            iz_sp=iz_sp+1
            snapshots_field[iz_sp,ir_sp,isnap]=field[iz,ir]
        end
    end
end

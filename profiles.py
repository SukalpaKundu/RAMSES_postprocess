import numpy as np
import pynbody as pb


#%run -i ../plotlib/single_snapshot.ipynb



        
class simulation_ramses:
    def __init__(self,folder, name):
        
        self.folder = folder
        self.name   = name
        
        self.t_array= None
        
    def _get_time(self, snap_ids=[]):
        t_array=[]

        for snap_id in snap_ids:
            folder  = self.folder+snap_id
            snap    = pb.load(folder)
            t0      = snap.properties['time'].in_units('Myr')
            t_array = np.append(t_array,t0)
        return t_array
    
    
    def _boxlen(self,snap_id="output_00001"):
        folder  = self.folder+snap_id
        snap    = pb.load(folder)
        boxlen  = snap.properties['boxsize'].in_units('kpc')
        return boxlen
        
    
    def _load_snap(self,snap_id,verbose=False,component='star'):   
        ##------------usage----------------
        ##xx,yy,zz,vx,vy,vz,dd,PP,TT,mm,cs = sim._load_snap(snap_id,component="gas")
        ##xx,yy,zz,vx,vy,vz,mm             = sim._load_snap(snap_id,component="star/DM")
        ##----------end usage--------------

        if verbose:
            print("reading snapshot:",snap_id)

        print(snap_id)
        snap   = pb.load(self.folder+snap_id)
        
        
        if(component=="gas"):

            #load raw data
            xx     = np.array(snap.gas['x'].in_units('kpc'))
            yy     = np.array(snap.gas['y'].in_units('kpc'))
            zz     = np.array(snap.gas['z'].in_units('kpc'))
            mm     = np.array(snap.gas['mass'].in_units('Msol'))
            cs     = np.array(snap.gas['cs'].in_units('km s**-1'))

            vx     = np.array(snap.gas['vx'].in_units('km s**-1'))
            vy     = np.array(snap.gas['vy'].in_units('km s**-1'))
            vz     = np.array(snap.gas['vz'].in_units('km s**-1'))

            dd     = np.array(snap.gas['rho'].in_units('m_p cm**-3'))
            PP     = np.array(snap.gas['p'].in_units('Pa'))
            TT     = np.array(snap.gas['temp'].in_units('K'))
            #ll     = np.array(snap.gas['smooth'].in_units('kpc'))

            boxlen = snap.properties['boxsize'].in_units('kpc')

            #end load raw data



            #modify the data
            xx     = xx-boxlen/2
            yy     = yy-boxlen/2
            zz     = zz-boxlen/2


            PP     = PP * 7.24*10**16

            if verbose:
                print("snapshot loaded, gas:",snap_id)

            return xx,yy,zz,vx,vy,vz,dd,PP,TT,mm,cs
        
        elif(component=="star"):
             #load raw data
            xx     = np.array(snap.star['x'].in_units('kpc'))
            yy     = np.array(snap.star['y'].in_units('kpc'))
            zz     = np.array(snap.star['z'].in_units('kpc'))
            mm     = np.array(snap.star['mass'].in_units('Msol'))
            vx     = np.array(snap.star['vx'].in_units('km s**-1'))
            vy     = np.array(snap.star['vy'].in_units('km s**-1'))
            vz     = np.array(snap.star['vz'].in_units('km s**-1'))
            
            
            boxlen = snap.properties['boxsize'].in_units('kpc')

            #end load raw data
            #modify the data
            xx     = xx-boxlen/2
            yy     = yy-boxlen/2
            zz     = zz-boxlen/2

            if verbose:
                print("snapshot loaded, star:",snap_id)

            return xx,yy,zz,vx,vy,vz,mm
        
        elif(component=="DM"):
             #load raw data
            xx     = np.array(snap.dm['x'].in_units('kpc'))
            yy     = np.array(snap.dm['y'].in_units('kpc'))
            zz     = np.array(snap.dm['z'].in_units('kpc'))
            mm     = np.array(snap.dm['mass'].in_units('Msol'))
            vx     = np.array(snap.dm['vx'].in_units('km s**-1'))
            vy     = np.array(snap.dm['vy'].in_units('km s**-1'))
            vz     = np.array(snap.dm['vz'].in_units('km s**-1'))
            
            
            boxlen = snap.properties['boxsize'].in_units('kpc')

            #end load raw data
            #modify the data
            xx     = xx-boxlen/2
            yy     = yy-boxlen/2
            zz     = zz-boxlen/2

            if verbose:
                print("snapshot loaded, star:",snap_id)

            return xx,yy,zz,vx,vy,vz,mm
        else:
            raise("componet=star or gas or DM only")
            return    
    
    
    def _get_disk_radius(self,snap_ids= [], z_cut_for_mask=4.0,r_cut_for_mask=8.0,component="star",to_get="sigma",mass_frac=0.8):
        
        
        #****************usage******************
        #rdisk,zdisk =  sim._get_disk_radius(...)
        #***************************************
        
        
        
        
        r_disk     = []
        z_disk     = []
        folder_sim = self.folder
        for snap_id in snap_ids:

            folder                                = snap_id
            
            #print(folder)
            
            
        
            if(component=="star"):
                mm,xx,yy,zz,vx,vy,vz              = self._load_snap(folder,verbose=False,component='star')
            elif(component=="DM"):
                mm,xx,yy,zz,vx,vy,vz              = self._load_snap(folder,verbose=False,component='DM')
            elif(component=="gas"):
                xx,yy,zz,vx,vy,vz,dd,PP,TT,mm,cs  = self._load_snap(folder,verbose=False,component='gas')
            else:
                raise("component can be either 'star' or 'gas' or 'DM'")
            
            
            
            
            rr                                = np.sqrt(xx**2+yy**2)
            mask_r                            = rr<r_cut_for_mask
            mask_z                            = np.abs(zz)<z_cut_for_mask
            mask                              = mask_r*mask_z
            #print(np.shape(rr))
            if(to_get=="sigma"):
                sigma_r                       = np.sqrt(np.sum(rr**2*mm*mask)/np.sum(mm*mask))
                sigma_z                       = np.sqrt(np.sum(zz**2*mm*mask)/np.sum(mm*mask))
                

            elif(to_get=="m_to_mtot"):

                rr=np.sqrt(xx**2+yy**2) 

                mtot                          = np.sum(mm*mask)
                ii_r                          = np.argsort(rr)    
                r_sorted                      = rr[ii_r]
                m_cumsum_r                    = np.cumsum((mm[ii_r]))

                ii_z                          = np.argsort(np.abs(zz))    
                z_sorted                      = np.abs(zz[ii_z])
                m_cumsum_z                    = np.cumsum((mm[ii_z]))

                sigma_r                       = np.interp(mtot*mass_frac, m_cumsum_r, r_sorted)
                sigma_z                       = np.interp(mtot*mass_frac, m_cumsum_z, z_sorted)
                del(mtot,ii_r,ii_z,m_cumsum_r,m_cumsum_z,rr)

            else:
                raise("to get can be 'sigma' or 'm_to_mtot'")

            r_disk                        = np.append(r_disk,sigma_r)
            z_disk                        = np.append(z_disk,sigma_z)
            del(xx,yy,zz,vx,vy,vz,mm,mask,mask_r,mask_z,sigma_r,sigma_z)
            #print("loaded",snap_id,r_disk,z_disk)

        return r_disk,z_disk

    

        
   
        
        
    def _get_profiles(self,snap_id="output_00001",bins=np.linspace(0,10,1),z_cut=10,r_cut=10,z_cut_profile=10,r_cut_profile=10,component="gas",weight_by="mass",direction="r",verbose=True):
    
    
        #-----------------usage--------------------
        #dd_r,PP_r,TT_r,vp_r,cs_r,vr_r,v_sigma_r = sim._get_profiles(snap_id,bins,z_cut,r_cut,component="gas",weight_by="density/mass",direction="r")
        #dd_p,PP_p,TT_p,cs_p,vz_p                = sim._get_profiles(snap_id,bins,z_cut,r_cut,component="gas",weight_by="density/mass",direction="z")
        #dd_p,vp_p,vr_p,v_sigma_p                = sim._get_profiles(snap_id,bins,z_cut,r_cut,component="star/DM",weight_by="density/mass",direction="r")
        #dd_p,vz_p                               = sim._get_profiles(snap_id,bins,z_cut,r_cut,component="star/DM",weight_by="density/mass",direction="z")
        
        #****for star/DM, the binning must be obtained with a linear spacing
        
        #---------------end usage------------------


        if(component=="star"):
            xx,yy,zz,vx,vy,vz,mm             =  self._load_snap(snap_id,verbose=False,component='star')
        elif(component=="DM"):
            xx,yy,zz,vx,vy,vz,mm             =  self._load_snap(snap_id,verbose=False,component='DM')
        elif(component=="gas"):
            xx,yy,zz,vx,vy,vz,dd,PP,TT,mm,cs =  self._load_snap(snap_id,verbose=False,component='gas')
        else:
            raise("component can be 'gas','DM' or 'star' only")



        rr                    =  np.sqrt(xx**2+yy**2)
        vp                    =  (xx*vy-yy*vx)/rr
        vr                    =  (xx*vx+yy*vy)/rr

        if (weight_by=="mass"):
            wg = mm
        elif (weight_by=="density"):
            if(component!="gas"):
                raise("density weighting currentcly unavailable for non gas particles")
            wg = dd
        else:
            wg = np.ones_like(mm)

        


        if(direction=="r"):
            
            wg[rr>r_cut]                  =  0
            wg[np.abs(zz)>z_cut_profile]  =  0 


            id_bin                = np.digitize(rr,bins)
            vp_p                  = np.bincount(id_bin,weights=vp*wg,minlength=len(bins))/(np.bincount(id_bin,weights=wg,minlength=len(bins)))
            vr_p                  = np.bincount(id_bin,weights=vr*wg,minlength=len(bins))/(np.bincount(id_bin,weights=wg,minlength=len(bins)))
            vx_p                  = np.bincount(id_bin,weights=vx*wg,minlength=len(bins))/(np.bincount(id_bin,weights=wg,minlength=len(bins)))
            vy_p                  = np.bincount(id_bin,weights=vy*wg,minlength=len(bins))/(np.bincount(id_bin,weights=wg,minlength=len(bins)))
            vz_p                  = np.bincount(id_bin,weights=vz*wg,minlength=len(bins))/(np.bincount(id_bin,weights=wg,minlength=len(bins)))
            vx2_p                 = np.bincount(id_bin,weights=vx**2*wg,minlength=len(bins))/(np.bincount(id_bin,weights=wg,minlength=len(bins)))
            vy2_p                 = np.bincount(id_bin,weights=vy**2*wg,minlength=len(bins))/(np.bincount(id_bin,weights=wg,minlength=len(bins)))
            vz2_p                 = np.bincount(id_bin,weights=vz**2*wg,minlength=len(bins))/(np.bincount(id_bin,weights=wg,minlength=len(bins)))

            v_sigma_p             = np.sqrt((vx2_p-vx_p**2)+(vy2_p-vy_p**2)+(vz2_p-vz_p**2))

            if(component=="gas"):

                cs_p                  = np.bincount(id_bin,weights=cs*wg,minlength=len(bins))/(np.bincount(id_bin,weights=wg,minlength=len(bins)))
                dd_p                  = np.bincount(id_bin,weights=dd*wg,minlength=len(bins))/(np.bincount(id_bin,weights=wg,minlength=len(bins)))
                PP_p                  = np.bincount(id_bin,weights=PP*wg,minlength=len(bins))/(np.bincount(id_bin,weights=wg,minlength=len(bins)))
                TT_p                  = np.bincount(id_bin,weights=TT*wg,minlength=len(bins))/(np.bincount(id_bin,weights=wg,minlength=len(bins)))

                return dd_p,PP_p,TT_p,vp_p,cs_p,vr_p,v_sigma_p
            else:
                
                mm_p                  = np.bincount(id_bin,weights=mm,minlength=len(bins))
                dr_p                  = np.abs(bins[1]-bins[0])
                #dr_p                  = np.append([0],dr_p)
                VV_p                  = 2*np.pi*bins*(2*z_cut_profile)*dr_p
                
                if(component=="star"):
                    dd_p                  = mm_p/VV_p
                elif(component=="DM"):
                    dd_p                  = mm_p[1:]/VV_p
                
                
                scale_d               = 4.05e-8                   #Msol/kpc**3 to m_p/cc
                dd_p                  = dd_p*scale_d
                
                
                
                return dd_p,vp_p,vr_p,v_sigma_p

        if(direction=="z"):
            
            wg[rr>r_cut_profile]  =  0
            wg[np.abs(zz)>z_cut]  =  0
            
            id_bin                = np.digitize(zz,bins)
            vz_p                  = np.bincount(id_bin,weights=vz*wg)/(np.bincount(id_bin,weights=wg))

            if(component=="gas"):

                cs_p                  = np.bincount(id_bin,weights=cs*wg,minlength=len(bins))/(np.bincount(id_bin,weights=wg,minlength=len(bins)))
                dd_p                  = np.bincount(id_bin,weights=dd*wg,minlength=len(bins))/(np.bincount(id_bin,weights=wg,minlength=len(bins)))
                PP_p                  = np.bincount(id_bin,weights=PP*wg,minlength=len(bins))/(np.bincount(id_bin,weights=wg,minlength=len(bins)))
                TT_p                  = np.bincount(id_bin,weights=TT*wg,minlength=len(bins))/(np.bincount(id_bin,weights=wg,minlength=len(bins)))

                return dd_p,PP_p,TT_p,cs_p,vz_p
            else:
                
                mm_p                  = np.bincount(id_bin,weights=mm,minlength=len(bins))
                dz_p                  = np.abs(bins[1]-bins[0])
                #dz_p                  = np.append([0],dz_p)
                VV_p                  = np.pi*r_cut_profile**2*dz_p
                
                
                if(component=="star"):
                    dd_p                  = mm_p/VV_p
                elif(component=="DM"):
                    dd_p                  = mm_p[1:]/VV_p
                
                
                scale_d               = 4.05e-8                   #Msol/kpc**3 to m_p/cc
                dd_p                  = dd_p*scale_d
                
                return dd_p,vz_p

def amr_bins(r_points=[],z_points=[],boxlen=0,ncoarse=6,nfine=6,r_max_to_plot=0,z_max_to_plot=0):
    
    
    #rpoints and zpoints must be arrays with size 4
    #******************usage****************
    #rbin,zbin = bins(r_points=[],z_points=[],boxlen=boxlen,ncoarse=ncoarse,nfine=nfine,r_max_to_plot=..,z_max_to_plot=..)
    
    dl        = boxlen/2**nfine
    dl_coarse = boxlen/2**ncoarse

    rbin = []
    rr   = 0

    r1   = r_points[0]
    r2   = r_points[1]
    r3   = r_points[2]
    r4   = r_points[3]
    
    while rr<r_max_to_plot:

        rbin = np.append(rbin,rr)

        if(rr<r1):
            dx = boxlen/2**(nfine)
        elif(r1<rr and rr<r2):
            dx = boxlen/2**(nfine-1)
        elif(r2<rr and rr<r3):
            dx = boxlen/2**(nfine-2)
        elif(r3<rr and rr<r4):
            dx = boxlen/2**(nfine-3)
        else:
            dx = boxlen/2**ncoarse

        rr   = rr + dx





    zbin = []

    zz   = 0
    z1   = z_points[0]
    z2   = z_points[1]
    z3   = z_points[2]
    z4   = z_points[3]


    while zz<z_max_to_plot:

        zbin = np.append(zbin,zz)

        if(np.abs(zz)<z1):
            dz = boxlen/2**(nfine)
        elif(z1<np.abs(zz) and np.abs(zz)<z2):
            dz = boxlen/2**(nfine-1)
        elif(z2<np.abs(zz) and np.abs(zz)<z3):
            dz = boxlen/2**(nfine-2)
        elif(z3<np.abs(zz) and np.abs(zz)<z4):
            dz = boxlen/2**(nfine-3)
        else:
            dz = boxlen/2**ncoarse

        zz   = zz + dz


    zbin = np.append(zbin,-zbin[1:])
    zbin = sorted(zbin)
    
    return rbin,zbin
    

import numpy as np
import yt
from convert_FIRE_cubes import write_grid

field_names = [   ('PartType0','CosmicRayEnergy'),
               ('PartType0','H_nuclei_density'),
              ('PartType0', 'ElectronAbundance'),
	      ('PartType0', 'Density'), 
             ('PartType0', 'InternalEnergy'),
	     ('gas','magnetic_field_x'),
             ('gas', 'magnetic_field_y'),
	    ('gas', 'magnetic_field_z')]

write_grid('/panfs/ds09/hopkins/panopg/snapshots/m12i/cr_700/snapshot_600.0.hdf5',10,field_names,dims=[64,64,64])

#hopkins/phopkins/data1/GalaxiesOnFire/m11q_res880/output/snapdir_600/snapshot_600.0.hdf5',100,field_names,dims=[256,256,256])

#'/panfs/ds08/hopkins/sponnada/m12i/cr_700/snapdir_600/snapshot_600.0.hdf5',100,field_names,dims=[256,256,256])
#/mnt/data1/GalaxiesOnFire/multiphysics/m12i_res7000_md_mhd_cv/output/snapdir_600/snapshot_600.0.hdf5
#'/panfs/ds08/hopkins/sponnada/m12i/cr_700/snapdir_600/snapshot_600.0.hdf5',100,field_names,[512]*3)


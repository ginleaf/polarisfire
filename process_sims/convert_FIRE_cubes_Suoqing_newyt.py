import yt

import h5py

from yt.frontends.gizmo.api import GizmoDataset

from yt.utilities.physical_constants import mh

from unyt import g, cm

import numpy as np

import matplotlib.pyplot as plt

import os


@yt.derived_field(name="CosmicRayEnergy_spec", units="code_specific_energy", sampling_type = 'particle')

def _CosmicRayEnergy_spec(field,data):

    try:
        return data.ds.arr((data[("PartType0", "CosmicRayEnergy")] / data[("PartType0", "Masses")]).v, "code_specific_energy")
    except Exception as err:
        print(err)
        return data.ds.arr(np.full(data["density"].shape, np.nan), "code_specific_energy")
    
    
def write_grid(fn, radius, fields, dims=[256,256,256], overwrite=True):
    ''' Writes contents of simulation snapshot to grid using yt's 
            arbitrary_grid function within the gizmo api.
        Input:
        fn = string specifying the filename of the original snapshot you wish to grid
        radius =  float, the half side length of the resulting grid in kpc
        fields = list of tuples specifying the fields you wish to include, 
                    of the form ('PartType0','Density')
        dims = list of 3 integers, each specifying the number of pixels along a given axis
        overwrite = bool, deprecated now. By default a pre-existing grid of 
                    the same name will be overwritten
        Output: A gridded cube 
        hdf5 file named fn+'_cut%d_r%d.hdf5'%(Npix,radius)
        
    '''

    print('Requested fields:', fields)

    # Read the original file to get the names of fields it already has
    f_orig = h5py.File(fn,'r')
    print('Existing fields are',f_orig.keys())
    print(f_orig['Header'].keys())
    print(f_orig['PartType0'].keys())

    def get_field_data(field_name, grid): # Conversion for CR energy
        '''
        field_name is a tuple like ('PartType0','Density'), 
        which means the gas density (particle type 0 means gas, 4 means stars, see Gizmo user guide)
        '''
        if field_name[1] == "CosmicRayEnergy_spec":

            if not ('PartType0', 'CosmicRayEnergy') in ds.field_list: raise ValueError("No CosmicRayEnergy, skip")

            ds.add_field(('PartType0','CosmicRayEnergy_spec'),
                         function=_CosmicRayEnergy_spec,
                         units="code_specific_energy", sampling_type = 'particle')
                          
                    
            field_data = (grid[('PartType0','CosmicRayEnergy_spec')]).in_cgs()#"CosmicRayEnergy_spec"]).in_cgs()
            print('done with CR')

        else:

            field_data = grid[field_name].in_cgs()

        return field_data



    f = h5py.File('%s_cut%d_r%d.hdf5' % (fn, dims[0],radius), "w")

    ds = GizmoDataset(fn)

    v, c = ds.find_max(('gas',"density"))

    radius = ds.quan(radius, 'kpc')

    left = c - radius # yt knows length unit of kpc/h

    right = c + radius

    grid = ds.arbitrary_grid(left, right, dims)
    

    f.create_dataset("radius", [1], data=radius.v, dtype="float64")

    if hasattr(ds, "current_redshift"):
        f.create_dataset("current_redshift", [1], data=ds.current_redshift, dtype="float64")

    if hasattr(ds, "current_time"):
        f.create_dataset("current_time", [1], data=ds.current_time.in_cgs(), dtype="float64")


    for field in fields:

        if 'H_nuclei_density' in field:
            # Compute the number density of Hydrogen, will later be used for calculating number
            # density of thermal electrons (from ElectronAbundance which is per proton, not per nucleon,
            # as described in gizmo user guide)
            print('Creating H_nuclei_density field....')

            proton_mass = mh
            # we need to create this field from the mass density, asuming a helium fraction
            field_data = get_field_data(('PartType0','Density'), grid)

            # get Helium abundance
            He = get_field_data(('PartType0', 'Metallicity_01'), grid)
            Hearray = He.v.astype("float64")
            
            plt.hist(Hearray.flatten(),histtype = 'step',bins=100)
            plt.savefig('Heabundance.png')

            # Make sure you are dividing in units consistent with the proton mass (cgs)
            assert field_data.units.same_dimensions_as(g/cm**3)
            # Number density = hydrogen mass fraction * mass density / proton_mass
            n_field_data = (1-Hearray) * field_data.v.astype("float64")/proton_mass

            print('create dataset')
            f.create_dataset('H_nuclei_density', dims, 

                     data=n_field_data,

                     dtype="float64")

            f.attrs["%s_units" % 'H_nuclei_density'] = '1/cm**3'

            print("%s field interpolated in unit of %s" % ('H_nuclei_density', '1/cm**3'))

        else:

            try:


                field_name = field

                print('Getting field ', field_name)

                
                field_data = get_field_data(field_name, grid)

                f.create_dataset(field_name[1], dims, 

                                 data=field_data.v.astype("float64"),

                                 dtype="float64")

                f.attrs["%s_units" % field_name[1]] = str(field_data.units)

                print("%s field interpolated in unit of %s" % (field_name[1], str(field_data.units)))

            except Exception as err:

                print(err)
                



    f.close()
    return



def load_grid(fn, fields=None, exclude_fields=[], interval=1):

    f = h5py.File(fn, "r")

    data = {}

    if fields is None:

        fields = f.keys()
        print('Fields: ',fields)

    r = f["radius"][0]

    kpc = 3.08567758e21

    bbox = np.array([[-r, r], [-r, r], [-r, r]]) * kpc

    for field in fields:
        print('Loading field ',field)
        
        if field == 'Density':
            field_name = 'density'
        else:
            field_name = field

        if field in exclude_fields: continue

        if field not in ["radius" ,"current_redshift", "current_time"]:
            print('in if')

            if field not in f.keys():

                print("ERROR: field %s cannot be found!" % field)

                print("All fields available:", list(f.keys()))

                raise KeyError

            field_data = f[field][...]

            field_unit = f.attrs["%s_units" % field]

            
            if field_name in ["density", "temperature"]:

                if field_data.min == 0.: print("WARNING: %s field contains zero value" % field)

                field_data[field_data==0.] = field_data[field_data>0.].min()

            data[("gas", field_name)] = (field_data[::interval,::interval,::interval], field_unit)

            dims = data[("gas" ,field_name)][0].shape

            print("%s field loaded with unit of %s" % (field_name, field_unit))
            
    print('Before load')
    ds = yt.load_uniform_grid(data, dims, length_unit="cm", bbox=bbox)
    
    path, filename = os.path.split(fn)

    ds.basename = filename

    if "current_redshift" in fields: ds.current_redshift = f["current_redshift"][0]

    if "current_time" in fields: ds.current_time = f["current_time"][0]

    return ds







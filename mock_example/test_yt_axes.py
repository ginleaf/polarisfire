#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 16 23:12:07 2022

@author: gin
"""

import yt
import numpy as np
import matplotlib.pyplot as plt
import os
import h5py
import sys

from unyt import eV,erg

sys.path.append('/home/gin/projects/FIRE/')

os.chdir('/home/gin/projects/FIRE/POLARIS/polarisfire/mock_example/')

from convert_FIRE_cubes_Suoqing_newyt import write_grid
#%%

def make_mock(fn,savename,params,fields=None):
    
    '''
    modifies tha dataset to have constant B field values, 
        and two zero density planes. For testing purposes.
    '''
    
    CRdens, Bx, By, Bz = params
    
    # Read the hdf5 file you wish to modify
    f = h5py.File(fn, "r")
    f_new = h5py.File(savename,'w')
    
    if fields is None:
        fields = f.keys()
        print('Fields: ',fields)

    try:
        f_new.create_dataset("radius", [1], data=f['radius'], dtype="float64")
        f_new.create_dataset("current_redshift", [1], data=f['current_redshift'], dtype="float64")
        f_new.create_dataset("current_time", [1], data=f['current_time'], dtype="float64")
    
        dims = f['H_nuclei_density'].shape
        
        for ii,field in enumerate(fields):
            
            
            if field in ['radius','current_redshift',"current_time"]:
                continue
            
            print('Copying', field)
            
            # we will modify the following fields
            # Set Bx, By, Bz constant everywhere
            if field == 'magnetic_field_x':
                mag_field_x = np.full(dims,Bx)
                f_new.create_dataset('magnetic_field_x', dims, data = mag_field_x, dtype="float64")
                
            elif field == 'magnetic_field_y':
                mag_field_y = np.full(dims,By)
                f_new.create_dataset('magnetic_field_y', dims, data = mag_field_y, dtype="float64")
                
            elif field == 'magnetic_field_z':
                mag_field_z = np.full(dims,Bz)
                f_new.create_dataset('magnetic_field_z', dims, data = mag_field_z, dtype="float64")
                 
            else:    
                f_new.create_dataset(field, dims, data=f[field], dtype="float64")
    
            f_new.attrs["%s_units" % field] = f.attrs["%s_units" % field]
    
            print("%s field written in unit of %s" % (field, f.attrs["%s_units" % field]))
            
        # Add n_cre field:
        # Now set by hand the cr density to be constant everywhere
        n_cre = np.full(dims,CRdens)
        # Add two planes of zeros
        n_cre[:,0,:] = np.zeros_like(n_cre[:,0,:])
        n_cre[:,:,10] = np.zeros_like(n_cre[:,:,10])
            
        f_new.create_dataset('n_cre', dims, data = n_cre, dtype="float64")
        f_new.attrs["%s_units" % 'n_cre'] = '1/cm**3'
        
        # and n_e field:
        n_e = f["ElectronAbundance"][...] * f['H_nuclei_density'][...]
        f_new.create_dataset('n_e', dims, data = n_e, dtype="float64")
        f_new.attrs["%s_units" % 'n_e'] = '1/cm**3'

        f_new.close()
        f.close()
        return
                
    except Exception as err:
        print(err)
        f_new.close()
        f.close()
        return
    
    return
            
    
    

def load_grid_and_compute_nCR_ne(fn, fields=None, exclude_fields=[], interval=1, \
                          p_to_e_ratio = 50./1,compute_ncr=True):
    ''' Loads the hdf5 file aand returns yt dataset
        
        
        Also defines CRe number density and thermal electron density fields
    '''
    
    
    f = h5py.File(fn, "r")
    
    try:
        
    
        data = {}
    
        if fields is None:
            fields = f.keys()
            print('Fields: ',fields)
    
        r = f["radius"][0]
        kpc = 3.08567758e21
        bbox = np.array([[-r, r], [-r, r], [-r, r]]) * kpc
    
        for field in fields:
            
            if field in ['radius','current_redshift',"current_time"]:
                continue
            
            print('Loading field ',field)
            
            if field == 'Density':
                field_name = 'density'
            else:
                field_name = field
    
            if field in exclude_fields: continue
    
            if field not in ["radius" ,"current_redshift", "current_time"]:
    
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
                
        if compute_ncr:
            # Compute number density of CR protons
            GeV_to_erg = 1.602e-19 * 10**7 * 10**9
            
            n_cre = f["CosmicRayEnergy"][...] * f['Density'][...] / GeV_to_erg / p_to_e_ratio
            data[('gas','n_cre')] = (n_cre[::interval,::interval,::interval],  "1/cm**3")
            
            # Compute thermal electron density
            n_e = f["ElectronAbundance"][...] * f['H_nuclei_density'][...]
            data[('gas','El_number_density')] = (n_e[::interval,::interval,::interval],  "1/cm**3")
            
                

        # Load as uniform grid in yt
        ds = yt.load_uniform_grid(data, dims, length_unit="cm", bbox=bbox)
        
        path, filename = os.path.split(fn)
        ds.basename = filename
        if "current_redshift" in fields: ds.current_redshift = f["current_redshift"][0]
        if "current_time" in fields: ds.current_time = f["current_time"][0]
    
        return ds
     
        
    except Exception as err:
        print(err)
        f.close()
        
        raise
        return
        
   


#%%
    
## INPUT parameters ########
    
savename = 'mock_cube_yt.hdf5'

Npix = 32
radius = 10

By = 1e-6 #G
Bx = 0
Bz = By/10.

CRdens = 1e-10
p_to_e_ratio = 50./1

params = [CRdens,Bx,By,Bz]

fn = 'snapshot_600.0.hdf5_cut%d_r%d.hdf5'%(Npix,radius)

ds = load_grid_and_compute_nCR_ne(fn, p_to_e_ratio = p_to_e_ratio)

ds.field_list

#%%


# Modify the cube to have constant B, ncr
print('Writing to file',savename)

make_mock(fn,savename,params)

print('Load the grid --------------------')
ds_mock = load_grid_and_compute_nCR_ne(savename, p_to_e_ratio = p_to_e_ratio, compute_ncr=False)
 
ds_mock.field_list

#%%
## Check output with a slice plot

slc = yt.SlicePlot(ds, "z", ("gas", "magnetic_field_x"), center=[0., 0., 0.],origin='native')
slc.set_cmap(("gas", "magnetic_field_x"), "Blues")
slc.annotate_grids(cmap=None)
slc.save('Bx_slice.png')

slc = yt.SlicePlot(ds, "z", ("gas", "magnetic_field_y"), center=[0., 0., 0.],origin='native')
slc.set_cmap(("gas", "magnetic_field_y"), "Blues")
slc.annotate_grids(cmap=None)
slc.save('By_slice.png')

slc = yt.SlicePlot(ds, "z", ("gas", "magnetic_field_z"), center=[0., 0., 0.],origin='native')
slc.set_cmap(("gas", "magnetic_field_z"), "Blues")
slc.annotate_grids(cmap=None)
slc.save('Bz_slice.png')


slc = yt.SlicePlot(ds, "z", ("gas", "CosmicRayEnergy"), center=[0., 0., 0.],origin='native')
slc.set_cmap(("gas", "CosmicRayEnergy"), "Blues")
slc.annotate_grids(cmap=None)
slc.save('crenergy_slice.png')


slc = yt.SlicePlot(ds, "z", ("gas", "El_number_density"), center=[0., 0., 0.],origin='native')
slc.set_cmap(("gas", "El_number_density"), "Blues")
slc.annotate_grids(cmap=None)
slc.save('ne_slice.png')

#%%
# Make lots of slice plots every 1 kpc (1 pixel)
for i in np.arange(-10,10):
    
    print(i)
    
    slc = yt.SlicePlot(ds_mock, "z", ("gas", "n_cre"), center=[0., 0., i],origin='native')
    slc.set_cmap(("gas", "n_cre"), "Blues")
    slc.annotate_grids(cmap=None)
    slc.save('slices/ncre_z_slice%d.png'%i)

    slc = yt.SlicePlot(ds_mock, "y", ("gas", "n_cre"), center=[0., 0., i],origin='native')
    slc.set_cmap(("gas", "n_cre"), "Blues")
    slc.annotate_grids(cmap=None)
    slc.save('slices/ncre_y_slice%d.png'%i)

    slc = yt.SlicePlot(ds_mock, "x", ("gas", "n_cre"), center=[0., 0., i],origin='native')
    slc.set_cmap(("gas", "n_cre"), "Blues")
    slc.annotate_grids(cmap=None)
    slc.save('slices/ncre_x_slice%d.png'%i)
    
#%%
    
# Write

f_new = h5py.File(savename,'w')

radius = ds.quan(radius, 'kpc').value
print(radius)

f_new.create_dataset("radius", [1], data=radius, dtype="float64")
if hasattr(ds, "current_redshift"):
    f_new.create_dataset("current_redshift", [1], data=ds.current_redshift, dtype="float64")
if hasattr(ds, "current_time"):
    f_new.create_dataset("current_time", [1], data=ds.current_time, dtype="float64")

def get_field_data(ds,field):
    field_arr = ds.all_data()[field]
    unit = field_arr[1].units
    data = field_arr.value
    return data, unit


for ii,field in enumerate(ds.field_list):
    data, unit = get_field_data(ds,field)
    if ii == 0:
        dims = data.shape
        print(dims)
        
    f_new.create_dataset(field[1],dims,data=data.astype("float64"),

                                 dtype="float64")

    f_new.attrs["%s_units" % field[1]] = str(unit)

    print("%s field interpolated in unit of %s" % (field[1], str(unit)))
    
f_new.close()
    
    
#%%




#%%

#GeV_to_erg = 1.602e-19 * 10**7 * 10**9 * eV/erg

# Calculate cube n_CR
#@yt.derived_field(name="CosmicRayNDensity", units="1/cm**3", sampling_type = 'cell')

#def _CosmicRayNDensity(field,data):
#    #try:
#    # CR proton volume density
#    print(data[('gas','density')])
#    return data["CosmicRayEnergy"] * data[('gas','density')] / GeV_to_erg
#    #except Exception as err:
#    #    print(err)
#    #    return data.ds.arr(np.full(data["H_nuclei_density"].shape, np.nan), "1/cm**3")

# Fails because of a weird error: Could not find field ('unknown', 'Density') in snapshot_600.0.hdf5_cut32_r10.hdf5.
# I don't know why it looks for this field, since I told it to operate on 'gas','density' but anyway...
#ds.add_field('CosmicRayNDensity',
#                         function=_CosmicRayNDensity,
#                         units="1/cm**3", sampling_type = 'cell',force_override=True)

#%%
# Calculate cube thermal electron density
#@yt.derived_field(name="El_number_density", units="1/cm**3", sampling_type = 'cell')

#def _El_number_density(field,data):
#    try:
#        # CR proton volume density
#        return data.ds.arr(data["ElectronAbundance"] * data['H_nuclei_density'], "1/cm**3")
#    except Exception as err:
#        print(err)
#        return data.ds.arr(np.full(data["H_nuclei_density"].shape, np.nan), "1/cm**3")


#ds.add_field('El_number_density',
#                         function=_El_number_density,
#                         units="1/cm**3", sampling_type = 'cell',force_override=True)
# Check output


#slc = yt.SlicePlot(ds, "z", ("gas", "El_number_density"), center=[0., 0., 0.],origin='native')
#slc.set_cmap(("gas", "El_number_density"), "Blues")
#slc.annotate_grids(cmap=None)
#slc.save('ne_slice.png')


#%%
'''
By = 1e-6 #G
Bx = 0
Bz = By/10.

data = ds.all_data()

## Modify cube B field to be constant, with two zero planes


def _mag_field_x(field,data):
    return data.ds.arr(np.full(data["magnetic_field_x"].shape, Bx), "G")

def _mag_field_y(field,data):
    return data.ds.arr(np.full(data["magnetic_field_y"].shape, By), "G")

def _mag_field_z(field,data):
    return data.ds.arr(np.full(data["magnetic_field_z"].shape, Bz), "G")
    
ds.add_field(('gas','magnetic_field_x_new'),
             function = _mag_field_x,
             units='G',
             sampling_type = 'cell',
             force_override = True)
             
ds.add_field('magnetic_field_y_new',
             function = _mag_field_y,
             units='G',
             sampling_type = 'cell',
             force_override = True)

ds.add_field('magnetic_field_z_new',
             function = _mag_field_z,
             units='G',
             sampling_type = 'cell',
             force_override = True)

# Modify density to put two planes of zeros
def _dens_new(field,data):
    newarray = np.full(data[("gas","density")].shape, ds.find_max(('gas',"density"))[0])
    
    return data.ds.arr(newarray, "G")


## Check output with a slice plot

slc = yt.SlicePlot(ds, "z", ("gas", "magnetic_field_x_new"), center=[0., 0., 0.],origin='native')
slc.set_cmap(("gas", "magnetic_field_x_new"), "Blues")
slc.annotate_grids(cmap=None)
slc.save('Bx_slice.png')

slc = yt.SlicePlot(ds, "z", ("gas", "magnetic_field_y_new"), center=[0., 0., 0.],origin='native')
slc.set_cmap(("gas", "magnetic_field_y_new"), "Blues")
slc.annotate_grids(cmap=None)
slc.save('By_slice.png')

slc = yt.SlicePlot(ds, "z", ("gas", "magnetic_field_z_new"), center=[0., 0., 0.],origin='native')
slc.set_cmap(("gas", "magnetic_field_z_new"), "Blues")
slc.annotate_grids(cmap=None)
slc.save('Bz_slice.png')



write_grid(ds,'mock_cube_yt.hdf5',10,field_names,dims=[Npix,Npix,Npix])




#%%   ----- Deprecated code -----------
Npix = 32

NCR_in = 1e-6 # CR volume density in cm^-3
B_in = 1e-6 # B in Gauss

ncr = np.ones((Npix,Npix,Npix))*NCR_in
# add a line of zero ncr to be able to discern axes
ncr[:,0,:] = np.zeros_like(ncr[:,0,:])
ncr[:,:,10] = np.zeros_like(ncr[:,:,10])

Bz = np.zeros_like(ncr)+B_in/10.
Bx = np.zeros_like(ncr)
By = np.ones((Npix,Npix,Npix))*B_in
electron_abundance = np.zeros_like(ncr)+1 # fully ionized medium
dens = np.ones_like(ncr)*1.6e-24#g/cm^-3
# add a small density maximum in the center of the cube for yt to center on
dens[int(Npix/2),int(Npix/2),int(Npix/2)] = 2*dens[0,0,0]

num_dens = dens/1.6e-24

GeV_to_erg = 1.602e-19 * 10**7 * 10**9
# assume proton-to-electron ratio    
p_to_e_ratio = 50./1
n_crp = ncr*p_to_e_ratio
crenergy_density = n_crp*GeV_to_erg
cr_spec_energy = crenergy_density/dens

# Helium abundance
Hearray = np.zeros_like(ncr)+0.25
    
#data = {('gas','density'): ((dens*u.Da*zh2).to(u.g/u.cm**3), "g/cm**3"),
#            ('gas','z_velocity'): (velo[0].to(u.km/u.s).value, 'km/s'),
#            ('gas','y_velocity'): (velo[1].to(u.km/u.s).value, 'km/s'),
#            ('gas','x_velocity'): (velo[2].to(u.km/u.s).value, 'km/s'),
#           }


data = { ('gas','cr_ndensity'): (ncr,"1/cm**3") ,
            #( 'gas','H_nuclei_density') : (num_dens,"g/cm**3"),
             ( 'gas','Density') : (dens,"g/cm**3"),
            #( 'gas','El_number_density'): (nth,"1/cm**3"),
            ( 'gas','magnetic_field_x'): (Bx,"G"),
            ( 'gas','magnetic_field_y'): (By,"G"),
            ( 'gas','magnetic_field_z'): (Bz,"G"),
            ( 'gas','CosmicRayEnergy'): (cr_spec_energy,"erg/g"),
            ('gas','ElectronAbundance'):(nth,'unitless'),
            ('gas','Metallicity_01'):(Hearray,'unitless'),
            }

#data_simple_keys = { 'current_redshift': yt.YTArray(0.71, '(dimensionless)'),
#            'density': yt.YTArray(dens,"g/cm**3"),
#            'cr_ndensity': yt.YTArray(ncr,"1/cm**3") ,
#            'H_nuclei_density' : yt.YTArray(dens,"g/cm**3"),
#            'El_number_density': yt.YTArray(nth,"1/cm**3"),
#            'magnetic_field_x': yt.YTArray(Bx,"G"),
#            'magnetic_field_y': yt.YTArray(By,"G"),
#            'magnetic_field_z': yt.YTArray(Bz,"G"),
#            }



#%%
#arr = np.zeros((32,64,64)) + dens# order seems to be xyz from plotting below

#arr[0,:,:] = np.zeros_like(arr[0,:,:])

#data = dict(density = (arr, "g/cm**3"))
bbox = np.array([[-16, 16], [-16, 16], [-16, 16]])
ds = yt.load_uniform_grid(data, dens.shape, length_unit="kpc", bbox=bbox)

# Write the hdf5 file
sname = 'mock_cube_yt.hdf5'
if os.path.isfile(sname):
    os.remove(sname)
    
#yt.save_as_dataset(ds,sname,data)


if os.path.isfile(sname+'_cut32_r16.hdf5'):
    os.remove(sname+'_cut32_r16.hdf5')
    

field_names = ['H_nuclei_density', 'ElectronAbundance', 'CosmicRayEnergy',\
        'current_redshift', 'density', 'magnetic_field_x', 'magnetic_field_y', 'magnetic_field_z','radius']

write_grid_testing(ds,'mock_cube_yt.hdf5',16,field_names,dims=[int(Npix),int(Npix),int(Npix)])



#%%
# starts counting from 0, in kpc, +- radius. 
# so the zero z plane is at -6 kpc
slc = yt.SlicePlot(ds, "z", ("gas", "cr_ndensity"), center=[0., 0., -6.],origin='native')
slc.set_cmap(("gas", "cr_ndensity"), "Blues")
slc.annotate_grids(cmap=None)
slc.save('test.png')
#slc.show()


#%%


ds = load_grid('mock_cube_yt.hdf5_cut32_r16.hdf5')
'''
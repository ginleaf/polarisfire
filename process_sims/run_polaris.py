import subprocess
import os
import shutil


def write_cmd(cmd_savename, octree_filename, results_dir, r,Npix,observer_loc,external_observer_distance_kpc):
    xobs,yobs,zobs = observer_loc

    fop=open(cmd_savename,'w')
    fop.write('<task> 1\n')
    fop.write('  <cmd>	CMD_SYNCHROTRON	\n  \n <max_subpixel_lvl> 0\n    ')
    fop.write(' <path_grid>      "'+octree_filename+'"\n')
    fop.write(' <path_out>       "./results/'+results_dir+'"\n')
    fop.write(' #convesion into SI\n  <conv_dens>     1e6 #cm^-3 in m^-3	\n <mu>		      1.4\n')
    fop.write(' # 1pc = 3.0857e+16 m\n')
    fop.write(' # this is half the side length so r * 1000* 3.0857e+16\n')
    side_length_kpc = 3.0857e+16 * r * 1000
    external_obs_distance_m = 3.0857e+16* 1000*external_observer_distance_kpc
    fop.write(' <conv_len>      {:e} # the code assumes this is +- the side length of XX kpc\n'.format(side_length_kpc))
    fop.write(' <conv_mag>      1e-4 # G in T\n')
    fop.write(' <nr_threads>  -1\n  #rotation axis 100 is x axis, 001 is z axis, 010 is y axis (manual has typo)\n')
    fop.write(' <axis1> 1 0 0\n  <axis2> 0 0 1\n')
    fop.write(' #detectors at lambda m rotated arround each of the ax1 and ax2  axis\n')
    fop.write(' # detector number of pixels lamda_min lambda_max number_of_wavelenghts no_bgr_source=1 angle_axis1 angle_axis2\n')
    fop.write(' # for healpix last 3 are: r_x, r_y, r_z in meters\n')
    fop.write(' <detector_sync nr_pixel = "{:d}">	29.9	29.9	1	1	0.0	0.0	{:e}\n'.format(Npix,external_obs_distance_m))
    fop.write(' <detector_sync nr_pixel = "{:d}">	29.9	29.9	1	1	0.0	90.0	{:e}\n'.format(Npix,external_obs_distance_m))
    #fop.write(' <detector_sync nr_pixel = "{:d}">	29.9	29.9	1	1	90.0	0.0	{:e}\n'.format(Npix,side_length_kpc*10))
    #fop.write(' <detector_sync nr_pixel = "{:d}">	29.9	29.9	1	1	90.0	90.0	{:e}\n'.format(Npix,side_length_kpc*10))
    fop.write(' <detector_sync nr_pixel = "{:d}">	0.05	0.05	1	1	0.0	0.0	{:e}\n'.format(Npix,external_obs_distance_m))
    fop.write(' <detector_sync nr_pixel = "{:d}">	0.05	0.05	1	1	0.0	90.0	{:e}\n'.format(Npix,external_obs_distance_m))
    fop.write('<detector_sync_healpix nr_sides = "64"> 29.9	29.9	1	1	{:e}	{:e}	{:e}\n'.format(xobs,yobs,zobs))
    fop.write('<detector_sync_healpix nr_sides = "64"> 0.738	0.738	1	1	{:e}	{:e}	{:e}\n'.format(xobs,yobs,zobs))    
    fop.write('</task>\n')
    fop.close()
    return

Npix = 256 # the number of pixels of any axis of the cube
r = 30 # the radius (half sidelength) of the cube in kpc

x,y,z = 3.0857e15,0,0#-2.748202E+20, -1.30178E+20, 1.446422E+19#90,110,130

observer_loc_units = 'm'#'pixels'

ext_observer_distance = 3.5e3 #kpc

cmd_savename = 'cmd_file'
octree_filename = "/panfs/ds09/hopkins/panopg/m12i_cr700/snapshot_600.0.hdf5_cut{:d}_r{:d}_mock_smooth_disk3.dat".format(Npix,r)
#"/home/sponnada/octrees_2_25/octree_r30_x512_512_z512.dat"
results_dir = "./smooth_disk3/"

if observer_loc_units == 'pixels':
    # Observer location in pixels from corner of cube
    xpix_0 = x
    ypix_0 = y
    zpix_0 = z

    # Observer location in pixels from center of cube
    center_x = Npix/2
    center_y = Npix/2
    center_z = Npix/2

    xpix = xpix_0 - center_x
    ypix = ypix_0 - center_y
    zpix = zpix_0 - center_z

    # Observer location in meters
    pix2kpc = 2*r/Npix
    kpc_to_m = 3.0857e+16 * 1000
    xobs = xpix * pix2kpc * kpc_to_m
    yobs = ypix* pix2kpc * kpc_to_m
    zobs = zpix* pix2kpc * kpc_to_m
elif observer_loc_units == 'm':
    xobs = x
    yobs = y
    zobs = z
else:
    raise Exception

observer_loc = [xobs,yobs,zobs]

# Edit cmd file
write_cmd(cmd_savename, octree_filename, results_dir, r,Npix,observer_loc,ext_observer_distance)


# Call polaris
PATH_TO_POLARIS = '$POLARISPATH'+'/bin/polaris'

#subprocess.call(['sh',PATH_TO_POLARIS+' '+cmd_savename])
#output_stream = os.popen(PATH_TO_POLARIS+' '+cmd_savename)
#output_stream.read()

def run_command(command):
    p = subprocess.Popen(command,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT,shell=True)
    return iter(p.stdout.readline, b'')

command = [PATH_TO_POLARIS,cmd_savename]
for line in run_command(command):
    print(line)

# copy cmd file to results directory
shutil.copyfile(cmd_savename, './results/'+results_dir+cmd_savename)

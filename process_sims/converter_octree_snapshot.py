import numpy as np
import struct
import sys
import h5py
import os
#%%

path_input = "../../m12i_cr700/snapshot_600.0.hdf5_cut256_r30.hdf5"
path_output = "../../m12i_cr700/snapshot_600.0.hdf5_cut256_r30.dat"

# data IDs for the POLARIS header
grid_id = 20  # grid ID (20 = octree)  
# grid header see manual table 3.3
data_ids = [0, 4, 5, 6, 22, 24, 25, 26, 27]

def loadData(path):  
    f = h5py.File(path)

    print(f.keys())
    
    magx = f['magnetic_field_x'][:] # this is in G
    magy = f['magnetic_field_y'][:] # this is in G
    magz =  f['magnetic_field_z'][:] # this is in G, no need to swap z coordinate it is already in the right place for numpy
    dens = f['H_nuclei_density'][:] # this is particles/cm^3
    #temp = f['temperature'][:]
    massdens = f['Density'][:]
    #n_th = f['n_e'][:]
    # Electron abundance
    electron_abundance = f['ElectronAbundance'][:]
    # Compute thermal electron volume density
    n_th = dens*electron_abundance

    # Compute total number density of CR electrons assuming a LISM spectrum
    # where proton energy density of 1 eV/cm^3 corresponds to 0.04 eV/cm^3 electron energy density
    # and this gives 10^-5 cm^-3 number density after integrating ... (details needed here)
    # Load CR specific energy 
    CR_spec_energy = f['CosmicRayEnergy_spec'][:] # this is specific energy (erg/g)
    CR_energy_density = CR_spec_energy*massdens # erg/cm^3
    eV_to_erg = 1.602e-19 * 10**7
    # CR electron volume density
    CR_energydens_eV = CR_energydens/eV_to_erg # convert to eV/cm^3
    # CR electron number density, assuming LISM spectrum
    n_CRe = CR_energydens_eV/0.04/10**5

    # Uncomment if you are only interested in number density of 1 GeV electrons
    #GeV_to_erg = 1.602e-19 * 10**7 * 10**9 # erg - these are 1 GeV protons
    # CR proton volume density
    #n_CRp = CR_energy_density/GeV_to_erg #", units="erg/cm**3", sampling_type="cell")
    # Compute CR electrons
    # assume proton-to-electron ratio    
    #p_to_e_ratio = 50./1
    #n_CRe = n_CRp/p_to_e_ratio
    
    #n_CRe = f['n_cre'][:]
    f.close()
    
    dim=dens.shape[0]
    print("shape, dim ", dens.shape, dim)
    
    print("Density (min,max):", dens.min(),"-", dens.max(), "cm^-3")
    print("nth     (min,max):", n_th.min(),"-", n_th.max(), "cm^-3")
    print("nCR     (min,max):", n_CRe.min(),"-", n_CRe.max(), "cm^-3")
    #print("Tgas    (min,max):", temp.min(),"-", temp.max(), "K")
    print("Bx      (min,max):", magx.min(),"-", magx.max(), "Guass")
    print("By      (min,max):", magy.min(),"-", magy.max(), "Guass")
    print("Bz      (min,max):", magz.min(),"-", magz.max(), "Guass")

    return dim, dens, n_CRe, n_th, magx, magy, magz #  temp,


CLR_LINE = "                                                                  \r"
cell_counter = 0
nr_of_cells = 0

class cell_oct:
    def __init__(self, _x_min, _y_min, _z_min, _length, _level):
        self.x_min = _x_min
        self.y_min = _y_min
        self.z_min = _z_min
        
        self.length = _length
        self.level = _level
    
        self.isleaf = 0
        self.data = []
        self.branches = []      

class OcTree:
    def __init__(self, _x_min, _y_min, _z_min, _length):
        self.root = cell_oct(_x_min, _y_min, _z_min, _length, 0)   

    def initCellBoundaries(self, cell,_level):
        x_min = cell.x_min
        y_min = cell.y_min
        z_min = cell.z_min
        l = 0.5 * cell.length

        level = _level

        cell.isleaf = 0
        cell.data = []
        cell.branches = [None, None, None, None, None, None, None, None]
        cell.branches[0] = cell_oct(x_min, y_min, z_min, l, level)
        cell.branches[1] = cell_oct(x_min + l, y_min, z_min, l, level)
        cell.branches[2] = cell_oct(x_min, y_min + l, z_min, l, level)
        cell.branches[3] = cell_oct(x_min + l, y_min + l, z_min, l, level)

        cell.branches[4] = cell_oct(x_min, y_min, z_min + l, l, level)
        cell.branches[5] = cell_oct(x_min + l, y_min, z_min + l, l, level)
        cell.branches[6] = cell_oct(x_min, y_min + l, z_min + l, l, level)
        cell.branches[7] = cell_oct(x_min + l, y_min + l, z_min + l, l, level)     
        
    def insertInTree(self, cell_pos, cell, _level, _limit):    
        x_pos = cell.x_min
        y_pos = cell.y_min
        z_pos = cell.z_min
        
        if cell_pos.level == cell.level:
            cell_pos.data=cell.data  
            cell_pos.isleaf=1
        else:
            if cell_pos.level == _limit:
              
                if len(cell_pos.data)==0:
                    cell_pos.data=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            
                d_level=-float(cell_pos.level - cell.level)
                fc=8.0**(-d_level)
              
                data_len = len(data_ids)
            
                for i in range(0,data_len):
                    cell_pos.data[i]+=fc*cell.data[i]
              
                cell_pos.isleaf=1
              
            else:
                if len(cell_pos.branches)==0:
                    self.initCellBoundaries(cell_pos,_level+1)

                x_mid = cell_pos.x_min+0.5*cell_pos.length
                y_mid = cell_pos.y_min+0.5*cell_pos.length
                z_mid = cell_pos.z_min+0.5*cell_pos.length
              
                new_cell_pos = cell_pos

                if(z_pos < z_mid): #z 0 1 2 3

                    if(y_pos < y_mid): #y 0 1

                        if(x_pos < x_mid): #x 0
                            new_cell_pos = cell_pos.branches[0]
                        else: #x 1
                            new_cell_pos = cell_pos.branches[1]

                    else: #y 2 3

                        if(x_pos < x_mid): #x 2
                            new_cell_pos = cell_pos.branches[2]
                        else: #x 3
                            new_cell_pos = cell_pos.branches[3]

                else: #z 4 5 6 7

                    if(y_pos < y_mid): #y 4 5

                        if(x_pos < x_mid): #x 4
                            new_cell_pos = cell_pos.branches[4]
                        else: #x 5
                            new_cell_pos = cell_pos.branches[5]

                    else: #y 6 7

                        if(x_pos < x_mid): #x 6
                            new_cell_pos = cell_pos.branches[6]
                        else: #x 7
                            new_cell_pos = cell_pos.branches[7]

            self.insertInTree(new_cell_pos, cell, _level+1,_limit)


    def writeOcTree(self, file, cell):
        global cell_counter
        global nr_of_cells
                       
        file.write(struct.pack("H", cell.isleaf))
        file.write(struct.pack("H", cell.level))   

        if cell.isleaf == 1:    
            data_len = len(cell.data)
            
            if cell_counter % 10000 == 0:
                sys.stdout.write('-> Writing octree grid file : ' + str(100.0 * cell_counter / nr_of_cells) + ' %     \r')
                sys.stdout.flush()
                
            cell_counter += 1 
         
            for i in range(0, data_len):
                file.write(struct.pack("f", cell.data[i]))
        else:
            for i in range(8):
                self.writeOcTree(file, cell.branches[i])
                
                
    def checkOcTree(self, cell):
        global cell_counter
        global nr_of_cells
               
        if cell.isleaf == 1:    
            length = len(cell.data)
            
            if length == 0:
                return False
            
            
            if cell_counter % 10000 == 0:
                sys.stdout.write('-> Checking octree integrity : ' + str(100.0 * cell_counter / nr_of_cells) + ' %     \r')
                sys.stdout.flush()
                
            cell_counter += 1    
            
        else:
            length = len(cell.branches)
            
            if length == 0:
                return False
            
            for i in range(8):
                self.checkOcTree(cell.branches[i])                
                
        return True    

if __name__ == "__main__":
    print("Loading data from: \n", path_input, "\n")
    
    #,data_velx,data_vely,data_velz
    #dim, data_dens, data_n_th, data_n_CRe, data_temp, data_magx, data_magy, data_magz = loadData(path_input)
    dim, data_dens, data_n_CRe, data_n_th, data_magx, data_magy, data_magz = loadData(path_input)    
    
    max_level=int(np.log2(dim))
    nr_of_cells =int(8**max_level)
    
    pos_min = 0.5*dim
    pos_off = 0.5
    
    max_length = 2.0
       
    print("Pos. min.    :",pos_min)
    print("Pos. off     :",pos_off)
    print("Dimension    :",dim)
    print("Range        :",-pos_min,"-",pos_min)
    print("Max. level   :",max_level)
    print("Max. length  :",max_length)
    print("Nr. of cells :",nr_of_cells,"\n")
    
    # init. octree
    tree = OcTree(-pos_min, -pos_min, -pos_min, dim)    
    
    per_counter=0.0
    per_max=float(dim**3)
    
    # fill octree
    for ix in range(0,dim):
        for iy in range(0,dim):
            for iz in range(0,dim):
                c_x = ix-pos_min+pos_off
                c_y = iy-pos_min+pos_off
                c_z = iz-pos_min+pos_off
                
                dens = data_dens[ix,iy,iz]#[iy,ix,iz]#
                
                #Tgas=data_temp[ix,iy,iz]
                
                magx = data_magx[ix,iy,iz]
                magy = data_magy[ix,iy,iz]
                magz = data_magz[ix,iy,iz]
                
                n_th = data_n_th[ix,iy,iz]
                n_CR = data_n_CRe[ix,iy,iz]
                
                g_min = 10
                g_max = 100000
                syn_p = 3

                cell = cell_oct(c_x, c_y, c_z, 0, max_level)

                # fill single cell with the data
                #cell.data =[dens, Tgas, magx, magy, magz, n_th, n_CR, g_min, g_max, syn_p]
                cell.data =[dens, magx, magy, magz, n_th, n_CR, g_min, g_max, syn_p]
                
                per_counter+=1
                
                if per_counter %10000==0:
                    sys.stdout.write('Constructing octree: ' + str(100.0 * per_counter / per_max) + ' %    \r')
                    sys.stdout.flush()  
                    
                #insert single cell into octree
                cell_root= tree.root
                tree.insertInTree(cell_root, cell,0,100)
    
    sys.stdout.write(CLR_LINE)
    print("Constructing octree:    done   "    )
    
    #check octree integrity
    check = tree.checkOcTree(cell_root)
    sys.stdout.write(CLR_LINE)
        
    if check == False:
        print("ERROR: Octree integrity is inconsistent!   \n\n")
        exit ()
    else:
        print("Octree structure   :    OK      "    )
    
    #write octree file header
    data_len = len(data_ids)
    file = open(path_output, "wb")
        
    file.write(struct.pack("H", grid_id))
    file.write(struct.pack("H", data_len))

    for d_ids in data_ids:
        file.write(struct.pack("H", d_ids))

    file.write(struct.pack("d", max_length))
    
    #write octree
    cell_counter = 0.0
    tree.writeOcTree(file, tree.root)
    sys.stdout.write(CLR_LINE)

    print("Writing octree     :    done   \n")
    
    print("Octree successfully created")

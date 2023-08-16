import os
import glob
from matplotlib import pyplot as plt

#simulation data will be output to this directory. Graphical output will be placed here too.
output_folder = "example_dihydrogen"

#the variables below are passed into the orca.inp file
method = "PM3"
parameters = "LARGEPRINT"
atom1 = "H"
atom2 = "H"
charge = 0
multiplicity = 1

#this defines the range of interatomic distances
range_distances = list(range(2, 70, 2))
distances = [i/10 for i in range_distances]

#This object class contains all the data extracted from the .out and .xyz files
class experiment_result:
    def __init__(self, distance, energy, homo, lumo):
        self.distance = distance
        self.energy_eh = float(energy)
        self.energy_kjmol = (float(energy))*2625.5002
        self.homo = homo
        self.lumo = lumo
        self.gap = float(lumo)-float(homo)

def folder_check(subfolders):    
    x = os.getcwd()
    for i in subfolders:
        if os.path.exists(f"{x}/{i}"):
            pass
        else:
            os.mkdir(str(i))
            print(f"Folder {i} not present, folder created in {x}\n")  
            
def check_out_file(distance):
#this function confirms that a completed .out file exists    
    try:
        with open("orca.out", "r") as f:
            check = f.readlines()[-1]
        if "TOTAL" in check:
            print(f"Calculation at {distance} angstroms already exists, using...")
            return True
        else:
            print(f"Calculation at {distance} angstroms not complete, running calculation...")
            return False
    except:
        return False
                
#this function runs an optimisation experiment, automatically identifying the energy minima of the system
def run_orca_opt(atom1, atom2, charge=charge, spin=multiplicity):
    folder_check(["opt"])
    os.chdir("opt")
    if check_out_file("optimised") == False:
        print("Carrying out optimisation experiment...")
        with open("orca.inp", "w") as f:
            f.write(f"!{method} {parameters} opt\n")
            f.write(f"*xyz {charge} {spin}\n")
            f.write(f"{atom1}\t1\t0\t0\n")
            f.write(f"{atom2}\t-1\t0\t0\n*")
        os.system("orca orca.inp > orca.out")
        check = "Start"
        while "TOTAL RUN TIME" not in check:
            with open("orca.out", "r") as f:
                check = (f.readlines()[-1])
        print("Optimisation complete")
    os.chdir(f"{path}//{output_folder}")

#this runs an experiment for a single interatomic distance
def run_orca_static(atom1, atom2, distance, charge=charge, spin=multiplicity):
    folder_check([distance])
    os.chdir(str(distance))
    if check_out_file(distance) == False:
        print(f"Running calculation at a distance of {distance} angstroms...")
        with open("orca.inp", "w") as f:
            f.write(f"!{method} {parameters}\n")
            f.write(f"*xyz {charge} {spin}\n")
            f.write(f"{atom1}\t{0.5*distance}\t0\t0\n")
            f.write(f"{atom2}\t-{0.5*distance}\t0\t0\n*")
        os.system("orca orca.inp > orca.out")
        check = "Start"
        while "TOTAL RUN TIME" not in check:
            with open("orca.out", "r") as f:
                check = (f.readlines()[-1])
        print(f"Calculation at {distance} complete")
    os.chdir(f"{path}//{output_folder}")

#This function reads through the .out files and extracts what we need. Can add extra functionality here!
def orca_read_out(folder):
    os.chdir(f"{path}//{output_folder}//{str(folder)}")
    if folder == "opt":
        with open("orca.xyz", "r") as f:
            x = list(f.readlines()[-1].split(" "))
            if x[5] == "":
                bond_length = abs(float(x[4])*2)
            else:
                bond_length = abs(float(x[5])*2)
            print(f"Optimised bond length - {bond_length}")   
    else:
        with open("orca.inp", "r") as f:
            x = (f.readlines(-2)[2].split("\t")[1])
            bond_length = abs(float(x)*2)
    with open("orca.out", "r") as f:
        x = f.readlines()
        for i in x:
            if "OCC" in i:
                index = (x.index(i))
        occ = 2
        level = 1
        
        while occ == float(2):
    
            occ = float(x[index+level].split(" ")[6])
            if occ == float(2):
                level = level+1
        
        HOMO = float(x[index+level-1].split(" ")[-2])
        LUMO = float(x[index+level].split(" ")[-2])
        print(f"HOMO: {HOMO} eV, LUMO:{LUMO} eV")
        for i in x:
            if "FINAL SINGLE POINT ENERGY" in i:
                energy = (i.split()[-1])
                
        experiment_list = experiment_result(bond_length, energy, HOMO, LUMO)
    os.chdir(f"{path}//{output_folder}")
    return experiment_list

abspath = os.path.abspath("extract_data.py")
path = os.path.dirname(abspath)
    
folder_check([output_folder])
os.chdir(output_folder)

folders = glob.glob(path+"//*//")

#this list will contain the experiment_result objects
experimental_data = []

opt_run = run_orca_opt(atom1, atom2)
opt_data = orca_read_out("opt")

for i in distances:
    run_orca_static(atom1, atom2, i)
    latest_data = orca_read_out(i)
    experimental_data.append(latest_data)
    #print(latest_data.distance, latest_data.energy_eh, latest_data.gap)
    print(f"Energy of system at {i} is {latest_data.energy_kjmol} kJ/mol\n")

experimental_data.sort(key=lambda x: x.distance, reverse=True)

os.chdir(path)

figure = plt.figure(figsize=(8,6))
plt.scatter([float(i.distance) for i in experimental_data], [float(i.energy_eh) for i in experimental_data], color = "black")
plt.plot(opt_data.distance, opt_data.energy_eh, color = "red", marker = "o")
plt.xlabel("Distance between H atoms / x10$^{-10}$ m", fontsize = 14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.ylabel("Potential Energy of system / kJ mol$^{-1}$", fontsize = 14)
title = f"Potential Energy of {atom1}-{atom2} atoms as a function of distance"
plt.title(title, fontsize = 16)
plt.savefig(f"{output_folder}/Potential_energy.png")
plt.close()

figure = plt.figure(figsize=(8,6))
plt.scatter([float(i.distance) for i in experimental_data], [float(i.homo) for i in experimental_data], color = "black")
plt.scatter([float(i.distance) for i in experimental_data], [float(i.lumo) for i in experimental_data], color = "red")
plt.xlabel("Distance between H atoms / x10$^{-10}$ m", fontsize = 14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.ylabel("Energy of HOMO/LUMO / eV", fontsize = 14)
plt.title("HOMO/LUMO positions", fontsize = 16)
plt.savefig(f"{output_folder}/HOMO_LUMO.png")
plt.close()

print("Simulation completed. Graphs output.")





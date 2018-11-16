
import os
import sys
import subprocess

def truncate_output_files(init_step, natoms):
# movie.xyz energies.dat input.in state.dat snafu.out velocities.xyz PES.dat 
# init_step is read from restart file
    natoms_lines = (natoms+2)*init_step   #  header
    step_lines = init_step                
    input_files = []
    input_files.append(["movie.xyz", natoms_lines])
    input_files.append(["velocities.xyz", natoms_lines])
    input_files.append(["PES.dat", step_lines])
    input_files.append(["energies.dat", step_lines])
    input_files.append(["state.dat", step_lines])
    
    for xx in range(len(input_files)):
        nlines = input_files[xx][1]
        input_file = input_files[xx][0]
        cmd_trunc = "head -n{} {} > temp_file; mv temp_file {} ".format(nlines, input_file, input_file)    
        try:
            trim_file = subprocess.run(cmd_trunc, stdout=None, stderr=subprocess.PIPE, shell = True, check = True)	
        except subprocess.CalledProcessError as cpe: 
            print("Return code: {}\nError: {}".format(cpe.returncode, cpe.stderr))
            exit(1)
        else:
            print("File {} was truncated after {} steps.".format(input_file,init_step))

if __name__ == "__main__":
   init_step, natoms = 20, 6
   truncate_output_files(init_step, natoms)
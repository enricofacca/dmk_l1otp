#######################################################################
# Assembler for muffa inputs
#######################################################################
These python scripts are designed to build and preprocess the inputs 
required for the muffa program, in the 2d setting, i.e.
the 2d triangulation, the subgrid, the forcing term, the initial condition, etc.
( the subgrid variable are handles by the program subgrid_preprocess 
called at the end of main)

The main.py reqires 3 inputs and it has to be called as

python inputs.ctrl folder_dat folder_vtk

where
inputs.ctrl : the file with the flags and paths 
	      describing which inputs must be assembled 
folder_dat  : folder where store all inputs in files .dat
folder_vtk  : folder where store all inputs in files .vtk


the file "inputs.ctrl" has the following structure:
Most of the lines reads as

str1 extra   ! flag_name

where:
1) str1  (string/float): define what to do
2) extra_info (string) : path to fiel with more informations
3) flag       (string) : flag name. It must be uniquely defined
   	      	        the variable so that line can be searched by 
			this name

Other lines are written in the following form
str1         !  flag_name
where only str1 is read.



# README #

Preprocess programs to create/preprocess all inputs
for the solver of dynamic Monge-Kantorovich equations
with the p1-p0 scheme.

### What is this repository for? ###

* Quick summary
Python scripts and FORTRAN program to build all inputs 
required for the p0-p1 discretization of the 
Dynamic Monge-Kantorovich model.
* Version
1.0

### How do I get set up? ###
# The scripts  ../../setup.sh should have compile all the fortran code 
# required. If not use
#

for f in *_assembly/{uniform_refinement,subgrid_preprocess}/code/
do
    cd $f
    echo $f
    make dirs
    make
    cd -
done


* Configuration

* Dependencies

BLAS and LAPACK library are required for the Fortran Programs
Python script are in Python 2.0
Meshpy and pyvtk are the non-standard Pyhton library required

* Database configuration
Three folders containing the assembler for the 2d, 3d and surface examples.
Three scripts that :
1 - Create a folder
2 - Launch the scripts/program building the inputs into the created folder 

* How to run tests
# STEP 1:  Open and modify file  "inputs.ctrl"
# in folder ../../scripts_controls/
# (We use the 2d-case as example, the 3d and the surface 
# input generation are similar)
# Most of the lines of "inputs.ctrl" reads as follows

"str1 extra_info   ! flag_name"

# where:
# str1 (string/float) : define what to do
# extra_info (string) : path to file with more informations
# flag_name  (string) : flag name. It must be uniquely defined
#   	      	         the variable so that line can be searched by 
#            		 this name
#Other lines are written in the following form
#

"str1         !  flag_name"

# where only str1 is read (for example nref).
# What str1 and extra_info do in the construction of the inputs is
# described in the python scripts in 2assembly.

# STEP 2 : (run the assembler) 
# Move in fodler assembly2d and type

python main.py inputs.ctrl fname fname_vtk	

# to create the input data described by the flags 
# in inputs.ctrl in the folder named "fname"
# Input data in vtk format will be written in folder fname_vtk

* Deployment instructions

### Contribution guidelines ###

* Writing tests
* Code review
* Other guidelines

### Who do I talk to? ###

* Repo owner or admin
* Other community or team contact
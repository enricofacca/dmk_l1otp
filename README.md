# README #

## What is this repository for?

This repository contains python
programs for the design of optimal transport network via Dynamic
Dynamical Monge-Kantorovich (DMK) model [^1 , ^2 , ^3]


Version * 2.0

## How do I get set up?
### Dependencies
  * INSTALLED PROGRAMS/LIBRIERIES
    * Fortran compiler compatible with FORTRAN 2003 (gfortran v>=4.8 )
    * Cmake (version>=2.8)
  * INSTALLED PROGRAMS/LIBRIERIES
    * Blas/Lapack libraries
  * EXTERNAL REPOSITORIES:
    Clone the following repositories in the directory external
    using the command
```
git submodule init;
git submoudel update --remote
```   
    * https://gitlab.com/enrico_facca/globals.git        
    * https://gitlab.com/enrico_facca/linear_algebra.git 
    * https://gitlab.com/enrico_facca/geometry.git 
    * https://gitlab.com/enrico_facca/p1galerkin.git 
    * https://gitlab.com/enrico_facca/dmk_solver.git


  * PYTHON >=3.6 and PACKAGES 
    * numpy  (MANDATORY)  "pip3 install numpy" 
    * click  (RECOMMENDED "pip3 install click" )
    * meshpy (RECOMMENDED "pip3 install meshpy" ) 
    * f90wrap (MANDATORY for fortran-python interface "pip3 install f90wrap")
  
### Compile

In Unix-base dsystem use the classical Cmake commands to compile
libraries and programs:
```
  mkdir build/;
  cd build;
  cmake ../; 
  make
```

### How to run tests

Move into the srcs/ directory and run
```
jupyter-notebook
```

### Deployment instructions

We adpot the [GitFlow](https://nvie.com/posts/a-successful-git-branching-model/) approach 

### Contribution guidelines ###

* Writing tests
* Code review
* Other guidelines

### Who do I talk to? ###

* Repo owner or admin
  * Enrico Facca (enrico.facca@sns.it)
  * Mario Putti  (putti@math.unipd.it) 
* Other community or team contact

### References

[^1]: E.Facca, F.Cardin, and M.Putti, [*Towards a stationary
  Monge-Kantorovich dynamics: the Physarum Polycephalum experience*](
  https://doi.org/10.1137/16M1098383),
  SIAM J. Appl. Math., 78 (2018), pp.651-676
  ([ArXiv version](https://arxiv.org/abs/1610.06325)).
  
[^2]: E.Facca, S.Daneri, F.Cardin, and M.Putti, [*Numerical solution of Monge-Mantorovich equations via a dynamic formulation*](https://doi.org/10.1007/s10915-020-01170-8), J. Scient. Comput., 82 (2020), pp.1-26.

[^3]: E.Facca, F.Cardin, and M.Putti, [*Branching structures emerging
  from a continuous optimal transport model*](https://arxiv.org/abs/1811.12691)
  ArXiv preprint 1811.12691, (2018)

[^4]: E.Facca, M.Benzi [*Fast Iterative Solution of the Optimal Transport Problem on Graphs*](https://arxiv.org/abs/2009.13478) ArXiv preprint 2009.13478, (2020)
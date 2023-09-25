# qe_workflow

This repository provides a sample workflow for running quantum espresso calculations on nanoHUB.

## Installation

### Create a nanoHUB account

To use this repository, first make an account at [nanohub](nanohub.org). This is a free account and will let you access hundreds of scientific tools, all hosted in a cloud environment. 

### Use a nanoHUB terminal

To copy this workflow to your nanoHUB storage, you need to open a terminal in nanoHUB. You have two options on how to do this:

1. (recommended) open a [jupyter notebook](https://nanohub.org/tools/jupyter70) session with nanoHUB, and click on 'new terminal' from the top right

2. ssh into nanoHUB using `ssh username@nanohub.org` with *username* as your nanoHUB account username, from your choice of terminal

From here, you can copy this workflow using `git clone https://github.com/katnykiel/qe_workflow.git`. This will put the workflow and pseudopotential files in your nanoHUB storage. 

## Running this workflow

### From the command line

To run this from a command line, you first need to open a terminal in nanoHUB using one of the methods described above. Then, you can take the following steps:

1. **Create your input file**. An example input file is provided for the relaxation of a Ti3C2 MXene, `mxene_relax.in` . These files can be edited from the command line using *vi* or the text editor of your choice. Copied below is the input file for a MXene relaxation:

```
&CONTROL
  calculation = 'vc-relax',
  outdir = './',
  pseudo_dir = './pseudo/pseudo_PAW/',
  tstress = .TRUE.,
/
&SYSTEM
  ecutwfc = 50,
  ibrav = 0,
  nat = 10,
  ntyp = 2,
/
&ELECTRONS
/
&IONS
/
&CELL
/
ATOMIC_SPECIES
  C  12.0107 C.upf
  Ti  47.8670 Ti.upf
ATOMIC_POSITIONS crystal
  Ti       0.666667000   0.333333000   0.343385745
  Ti       0.333333000   0.666667000   0.843307031
  Ti       0.333333000   0.666667000   0.656614255
  Ti       0.666667000   0.333333000   0.156692969
  Ti       0.000000000   0.000000000   0.500000000
  Ti       0.000000000   0.000000000   0.000000000
  C        0.333333000   0.666667000   0.414250951
  C        0.666667000   0.333333000   0.914278916
  C        0.333333000   0.666667000   0.085721084
  C        0.666667000   0.333333000   0.585749049
K_POINTS automatic
  3 3 3 0 0 0
CELL_PARAMETERS angstrom
  1.533303 -2.655757 -0.000000
  1.533303 2.655757 0.000000
  0.000000 0.000000 15.117150
```

2. **Run the simulation**. We run the simulation using the following command, printing the output to the console as it is generated:

```
use espresso-6.2.1
pw.x -i mxene_relax.in > mxene_relax.out &
watch tail mxene_relax.out
```

This should take no more than 10 minutes.

3. **Check the output**. The outputs are stored in .out files, which you can parse using a pre-built code by running the following

```
python
get_outputs()
get_convergence_XXX()
```

### From a jupyter notebook

To run this from a jupyter notebook, open the "qe_workflow.ipynb" file and run the cells. This workflow demonstrates how to submit jobs remotely to Purdue HPC clusters. 

## Issues

Please raise any issues in the [issues](https://github.com/katnykiel/qe_workflow/issues) section of this repository. 

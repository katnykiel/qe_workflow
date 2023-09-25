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

1. **Create your input file**. An example is provided for the relaxation of a MXene, Ti$_3$C$_2$. These files can be editor from the command line using *vi* or the text editor of your choice. Quantum espresso uses .in files as its input files.

2. **Run the simulation**. We run the simulation using the following command, printing the outpout to the console as it is generated:

```
espresso XXX YYY
```

3. **Check the output**. The outputs are stored in .out files, which you can parse using a pre-built code by running the following

```
python
get_convergence_XXX()
```

### From a jupyter notebook

To run this from a jupyter notebook, open the "qe_workflow.ipynb" file and run the cells. 

## Issues

Please raise any issues in the [issues](https://github.com/katnykiel/qe_workflow/issues) section of this repository. 

# Spatial simulations for Haughey et al 2022

This repository contains the code for simulating the spatial competition of tumour subclones. The program generates an output file cells.csv containing a list of the (x,y)-coordinates of all cells in the system, and their identity (wild-type or mutated, represented by a 0 or 1 respectively). Tested using GCC v4.2.1 and Python v3.8.5.

Execute a simulation by running:

```
./subClonalMixing [-v1] [-s S] [-t T] [-q Q] [-x X]
```

where\
&nbsp; -v1 &emsp;&emsp; verbose flag (optional)\
&nbsp; -s &emsp;&emsp; selective advantage of mutant population\
&nbsp; -t &emsp;&emsp; time of mutation (fraction of final tumour size)\
&nbsp; -q &emsp;&emsp; cell pushing strength\
&nbsp; -x &emsp;&emsp; random seed

Plot the final spatial data in cells.csv by running:

```
python3 ./plot_spatial_data.py [--path PATH]
```

where [PATH] is the file path the relevant cells.csv file.

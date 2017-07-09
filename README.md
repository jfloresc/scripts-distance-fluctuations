# README #

* Scripts used for Biophys. J. 2017, 112, 2291-2300.

Required packages MDtraj, cython, numpy

Calculation of pairwise distance fluctuations (f) in the ligand binding domain dimer. f = sqrt(<d^2> - <d>^2), where d is the pairwise distance between Calpha atoms and <> is the ensemble average.

I have written some python scripts, to calculate these fluctuations and plot them in a matrix plot.

I am using mdtraj, matplotlib and cython (for trajectory analysis, plotting and fast calculations). You could install them in Maria's cluster using pip:

pip install Cython
pip install matplotlib
pip install mdtraj


After that, you can compile distance.pyx script executing the following in your working directory:

python setup.py build_ext --inplace

After that you will generate a distance.so library to be uploaded by the script distance_matrix.py. You can change the trajectory and pdb reference file for each mutant or wild type trajectory.

python distance_matrix.py

the following script generates a heatmap plot using data from two separate trajectories (arguments-> first: upper_diagonal, second: lower diagonal):

python plot_distance_2d.py $WT_DIR/distance_matrix.dat $MUTANT_DIR/distance_matrix.dat 

Written by Jose Flores-Canales.
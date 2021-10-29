# Entropic_Path_Sampling
The repository documents how to perform entropic path sampling method to evaluate the entropic path along a reaction path. The model reaction is cyclopentadiene dimerization.
![](Protocol.png)

## Entropic Path Sampling Tutorial
We will use the program [PARENT](https://pubs.acs.org/doi/abs/10.1021/acs.jctc.5b01217) (GitHub Repo: https://github.com/markusfleck/PARENT) to compute entropy from molecular dynamics trajectories. Since PARENT only takes in GROMACS input files (“.top” for topology file and “.xtc” for coordinate files), we will have to find way to convert xyz files to PARENT-readable files before doing the entropy analysis. The following instructions will help you
1. Convert xyz files to GROMACS topology files
2. Convert xyz files to GROMACS coordinate files
3. Perform entropy computations using PARENT.

### Section 1 - .top File Preparation
1. Use [openbabel](http://openbabel.org/wiki/Category:Installation) to convert optimized .xyz file to .mol2 and to .pdb

        obabel dimerization-low-hp-C2-opt.xyz -O dimerization-low-hp-C2-opt.mol2
        obabel dimerization-low-hp-C2-opt.xyz -O dimerization-low-hp-C2-opt.pdb

2. Replace openbabel's default AMBER99 labels with GAFF labels using sed

        sed -i -e "s/C.3     1  LIG1/c3      1  LIG1/g" -e "s/H       1/hc      1/g" -e "s/S.O2    1  LIG1/s2      1  LIG1/g" -e "s/O.3     1  LIG1/o       1  LIG1/g" -e "s/O.2     1  LIG1/os      1  LIG1/g" -e "s/C.2     1  LIG1/c2      1  LIG1/g" dimerization-low-hp-C2-opt.mol2 

3. Rename .frcmod (force field) file

        mv ANTECHAMBER.frcmod dimerization-low-hp-C2-opt.frcmod 

4. Generate tleap ([AMBER](http://ambermd.org/GetAmber.php) input) file

        echo "source leaprc.gaff" >>tleap 
        echo "dimerization-low-hp-C2-opt= loadmol2 dimerization-low-hp-C2-opt.mol2">>tleap 
        echo "loadamberparams dimerization-low-hp-C2-opt.frcmod">>tleap 
        echo "saveamberparm dimerization-low-hp-C2-opt dimerization-low-hp-C2-opt.prmtop dimerization-low-hp-C2-opt.inpcrd">>tleap 
        echo "quit">>tleap 


5. Run tleap - This generates prmtop and inpcrd files, which are input files for molecular simulations in AMBER

        $AMBERHOME/bin/tleap -f tleap

6. Install acpype-master from [here](https://github.com/llazzaro/acpype/blob/master/acpype/scripts/acpype.py) and run acpype.py

        python3 /INSTALLATION_PATH/acpype-master/acpype/scripts/acpype.py -p dimerization-low-hp-C2-opt.prmtop -x dimerization-low-hp-C2-opt.inpcrd 

7. Rename the resulting .gro and .top files. These will be used for analysis with PARENT. Note that the numbers preceding "\_GMX" may differ.

        mv 885_GMX.gro dimerization-low-hp-C2-opt.gro 
        mv 885_GMX.top dimerization-low-hp-C2-opt.top 


### Section 2 - .xtc File Preparation
1. Processing - run XYZtraj-Trajs.py
    - Make sure all trajectories (.xyz files) are inside Processing/ntraj/ first
    - The 6 indices should be 6 13 4 14 10 12
    - Either use interactive input or traj.conf (1st line: reset? y/n 2nd line: bifurcation mode (1 or 2) 3rd line: indices of forming bonds)
    - The function of processing trajectories including obtaining TS sampling ensembles, rearrangement, classification can be implemnted via this script. (e.g.If trajectories have been rearranged , comment line 306.)
    - It may be necessary to change the RUNPOINT_IDX on line 19 to ensure that the runpoint is found in the trajectories properly
2. Truncation - run Truncation.py 
    - Cut off the trajectory starting from the position where bond(4,14) length is 5.00Å and stopping when the bond(4,14) length is 1.58Å for the first time
3. Slicing - run Subsection.py
    - Combination all the trajectories that form product A(i.e., formation of bond(4,14))
    - The bond(4,14) length range is 2.90Å to 1.58Å, refer to cycloaddition_TS_irc_rev.log for the intrinsic reaction coordinate (IRC) calculation that bridges cycloaddtion TS1 to Cope TS2, and to cope_TS_irc_rev.log for the IRC that bridges Cope TS2 to the product with formation of bond (4,14).
    - Divide the bond length range into N(ten here) equal parts and adjust the bond lengths on lines 31-32 to select the subsection in question.
    - Count the number of points falling into the range for convergence analysis
### Section 3 - Analysis with PARENT (Code adopted from https://github.com/markusfleck/PARENT)
1. Convert the subsectioned .xyz files to .xtc

        python3 xyz2xtc.py

    - Be sure to adjust the file names inside to get all trajectories from the Subsection folder
    - The mdtraj module may need to be installed if not already present
2. Move all .xtc files and the .top file from Section 2 and 1 to test_system inside PARENT-master
3. Revise the path TRJ and TOP inside PARENT-master/parameters-xin
4. Run PARENT
    1. Make sure you have the latest versions of GCC and OpenMPI
    2. make clean in PARENT-master
    3. Export MPI environment variables

            export MPI_NUM_PROCS=1; export OMP_NUM_THREADS=4;

    5. Run run.sh with the parameters file
   
            ./run.sh parameter-xin

5. Obtain the total configuration entropy from PARENT-master_MIST.txt in each output file.(the units of numbers in output files are k<sub>b</sub>, unit conversion method from k<sub>b</sub> to kcal/mol is in the output analysis file)   



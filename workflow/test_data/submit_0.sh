#!/bin/bash
#SBATCH -J myjob
#SBATCH -o myjob.o%j
#SBATCH -e myjob.e%j
#SBATCH -p skx-normal
#SBATCH -N 4
#SBATCH --ntasks-per-node 24
#SBATCH -t 48:00:00
###SBATCH -A myproject
cd Cs2Tl1Sc1Cl6_Cs2Cu1Sc1Cl6_Tl_0.250_Cu_0.750
ibrun /home1/05018/tg843171/vasp.5.4.4_vtst/bin/vasp_std
cd ..
cd Cs2Tl1Sc1Cl6_Cs2In1Sc1Cl6_Tl_0.750_In_0.250
ibrun /home1/05018/tg843171/vasp.5.4.4_vtst/bin/vasp_std
cd ..
cd Cs2Tl1Sc1Cl6_Cs2Tl1Sc1Br6_Cl_0.750_Br_0.250
ibrun /home1/05018/tg843171/vasp.5.4.4_vtst/bin/vasp_std
cd ..

#!/bin/bash
#SBATCH -J myjob
#SBATCH -o myjob.o%j
#SBATCH -e myjob.e%j
#SBATCH -p skx-normal
#SBATCH -N 4
#SBATCH --ntasks-per-node 24
#SBATCH -t 48:00:00
###SBATCH -A myproject

cd Cs2Tl1Sc1Cl6_Cs2In1Sc1Cl6_Tl_0.250_In_0.750
cp /Users/yao/Google Drive/mmtools/workflow/test_data/custodian_job.py ./custodian_job.py
python custodian_job.py
cd ..
echo "Cs2Tl1Sc1Cl6_Cs2In1Sc1Cl6_Tl_0.250_In_0.750" >> log_1

cd Cs2Tl1Sc1Cl6_Cs2Cu1Sc1Cl6_Tl_0.750_Cu_0.250
cp /Users/yao/Google Drive/mmtools/workflow/test_data/custodian_job.py ./custodian_job.py
python custodian_job.py
cd ..
echo "Cs2Tl1Sc1Cl6_Cs2Cu1Sc1Cl6_Tl_0.750_Cu_0.250" >> log_1

#!/bin/bash
#SBATCH --job-name=create_assignment_info
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=100GB
#SBATCH --partition=regular
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zhang.ji@wehi.edu.au

module load R/4.5.2
module load ImageMagick/7.1.1
module load mariadb-connector-c/3.3.10
module load JAGS/4.3.2
module load geos/3.12.1  
module load netcdf/4.9.2   
module load hdf5/1.12.3          
module load openCV/4.10.0   
module load proj/9.4.0    
module load gdal/3.9.0     
module load gcc/14.2

Rscript --verbose --vanilla create_assignment_info.r
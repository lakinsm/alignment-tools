#!/usr/bin/env bash
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --partition=defq
#SBATCH --cpus-per-task=1
#SBATCH --job-name=nanopore_backup
#SBATCH --mem=40gb
#SBATCH --time=20:00:00
#SBATCH --output=/lustrefs/data/backup_logs/sbatch_nanopore_backup_job_%j.log

module load python37
module load lakinsm-alignment-tools

# Arguments:
# 1. Source instrument root directory
# 2. Destination backup directory

echo -e "\n$(date): Starting backup of folder $1 to $2\n"
this_datetime=$( date +"%FT%H%M" )
srun nanopore_dir_walk_backup.py "$1" "$2" --logfile "/lustrefs/data/backup_logs/${this_datetime}_python_script.log" --slurm
echo -e "\n$(date): Backup finished for folder $1 to $2\n"

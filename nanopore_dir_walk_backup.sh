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

# Arguments:
# 1. Python script executable (nanopore_dir_walk_backup.py)
# 2. Source instrument root directory
# 3. Destination backup directory

echo -e "\n$(date): Starting backup of folder $2 to $3 using script $1\n"
this_datetime=$( date +"%FT%T" )
srun python3 "$1" "$2" "$3" --logfile "/lustrefs/data/backup_logs/${this_datetime}_python_script.log"
echo -e "\n$(date): Backup finished for folder $2 to $3 using script $1\n"

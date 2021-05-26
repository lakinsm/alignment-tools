#!/usr/bin/env bash

# module load python37
# module load lakinsm-alignment-tools

# Arguments:
# 1. Source instrument root directory
# 2. Destination backup directory

echo -e "\n$(date): Starting backup of folder $1 to $2\n"
this_datetime=$( date +"%FT%H%M%S" )
/cm/local/apps/python37/bin/python3 nanopore_dir_walk_backup.py "$1" "$2" --logfile "${this_datetime}_python_script.log" -n
echo -e "\n$(date): Backup finished for folder $1 to $2\n"

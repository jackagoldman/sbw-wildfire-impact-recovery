#!/bin/zsh

# Usage: ./random_scp.sh <remote_dir> <local_dir>
# Example: ./random_scp.sh /remote/path /local/path

remote_dir=$1
local_dir=$2

# Get list of files from remote directory
files=($(ssh goldma34@jameslab "ls '$remote_dir'"))

# Check if there are files
if [ ${#files[@]} -eq 0 ]; then
    echo "No files found in $remote_dir"
    exit 1
fi

# Select a random file
selected=$files[$RANDOM % ${#files[@]} + 1]

# Transfer the file
scp goldma34@jameslab:"$remote_dir/$selected" "$local_dir/"

echo "Transferred $selected to $local_dir/"
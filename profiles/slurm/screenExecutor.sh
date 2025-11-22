#!/bin/bash
echo "Starting screen session for SLURM job..."
echo "Session Name: faeSession"
echo "Use Ctrl+A D to detach from the screen session and 'screen -r faeSession' to reattach."
screen -L -S faeSession srun --pty bash

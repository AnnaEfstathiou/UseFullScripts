#!/bin/bash

# Read the given CSV file from command
column -s, -t < "$1" | less -#2 -N -S

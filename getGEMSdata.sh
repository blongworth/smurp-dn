#!/bin/bash

# Set local and remote paths
remote_name="whoi:MachineLab/SMURP"
local_path="data/"

rclone copyto --progress --include "GEMS_2024*.txt" "$remote_name" "$local_path"

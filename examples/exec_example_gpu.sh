#!/bin/bash

python3 make_points.py sources sources_1E4.bin 10000
python3 make_points.py targets targets_1E4.bin 10000
direct-gpu sources_1E4.bin targets_1E4.bin direct_benchmark.bin direct_summary.csv 10000 10000 0.0 0 1
tree-gpu sources_1E4.bin targets_1E4.bin direct_benchmark.bin tree_summary.csv 10000 10000 0.5 5 500 5 0 0.0 1

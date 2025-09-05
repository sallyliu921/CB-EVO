#!/bin/bash

mkdir -p logs
for blif_file in ../benchmarks/EPFL/*.blif; do
    if [ -f "$blif_file" ]; then
        echo "Processing: $blif_file"
        ./cb-evo -i "$blif_file" -t 2 -s 25 | tee "logs/$(basename "$blif_file").log"
        echo "Completed: $blif_file"
        echo "----------------------------------------"
    fi
done

echo "All EPFL benchmark files processed!"
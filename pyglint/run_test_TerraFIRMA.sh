#!/bin/bash

CONFIG_DIR="configs"

for config_file in "$CONFIG_DIR"/*.yaml; do
    echo "Running with config: $config_file"
    python test_TerraFIRMA.py --config "$config_file" &
done
wait
echo "All tests completed."

echo "Making plots..."
python plot_test_TerraFIRMA.py --icesheet AIS --variable total_smb
python plot_test_TerraFIRMA.py --icesheet GrIS --variable total_smb
wait
echo "All plots completed."
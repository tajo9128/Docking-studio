#!/bin/bash
set -e

echo "BioDockify Docking Worker Starting..."
echo "Job ID: $JOB_ID"

RECEPTOR_PATH="/data/receptor.pdbqt"
LIGAND_PATH="/data/ligand.pdbqt"
OUTPUT_PATH="/data/output/output_${JOB_ID}.pdbqt"

echo "Receptor: $RECEPTOR_PATH"
echo "Ligand: $LIGAND_PATH"
echo "Output: $OUTPUT_PATH"

if [ ! -f "$RECEPTOR_PATH" ]; then
    echo "Error: Receptor file not found"
    exit 1
fi

if [ ! -f "$LIGAND_PATH" ]; then
    echo "Error: Ligand file not found"
    exit 1
fi

mkdir -p /data/output

CENTER_X=${center_x:-0}
CENTER_Y=${center_y:-0}
CENTER_Z=${center_z:-0}
SIZE_X=${size_x:-22}
SIZE_Y=${size_y:-22}
SIZE_Z=${size_z:-22}
EXHAUSTIVENESS=${exhaustiveness:-8}
NUM_MODES=${num_modes:-9}
ENERGY_RANGE=${energy_range:-3}

echo "Grid Box: Center=($CENTER_X, $CENTER_Y, $CENTER_Z) Size=($SIZE_X, $SIZE_Y, $SIZE_Z)"
echo "Exhaustiveness: $EXHAUSTIVENESS"
echo "Num Modes: $NUM_MODES"

cat > /tmp/vina_config.txt <<EOF
center_x = $CENTER_X
center_y = $CENTER_Y
center_z = $CENTER_Z
size_x = $SIZE_X
size_y = $SIZE_Y
size_z = $SIZE_Z
exhaustiveness = $EXHAUSTIVENESS
num_modes = $NUM_MODES
energy_range = $ENERGY_RANGE
EOF

if command -v vina-gpu &> /dev/null; then
    echo "Running Vina-GPU..."
    vina-gpu --config /tmp/vina_config.txt \
             --receptor $RECEPTOR_PATH \
             --ligand $LIGAND_PATH \
             --out $OUTPUT_PATH \
             --log /data/output/vina_log_${JOB_ID}.txt
elif command -v vina &> /dev/null; then
    echo "Running Vina CPU..."
    vina --config /tmp/vina_config.txt \
         --receptor $RECEPTOR_PATH \
         --ligand $LIGAND_PATH \
         --out $OUTPUT_PATH \
         --log /data/output/vina_log_${JOB_ID}.txt
elif command -v gnina &> /dev/null; then
    echo "Running GNINA..."
    gnina -r $RECEPTOR_PATH \
          -l $LIGAND_PATH \
          --center_x $CENTER_X \
          --center_y $CENTER_Y \
          --center_z $CENTER_Z \
          --size_x $SIZE_X \
          --size_y $SIZE_Y \
          --size_z $SIZE_Z \
          --num_modes $NUM_MODES \
          --cnn_scoring rescore \
          -o $OUTPUT_PATH
else
    echo "Error: No docking engine found"
    exit 1
fi

if [ -f "$OUTPUT_PATH" ]; then
    echo "Docking completed successfully"
    exit 0
else
    echo "Error: Output file not generated"
    exit 1
fi

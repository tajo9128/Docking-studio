#!/bin/bash
set -e

echo "=========================================="
echo "  BioDockify Docking Studio"
echo "=========================================="
echo "Job ID: $JOB_ID"

# Auto-detect GPU
GPU_DETECTED=false
GPU_NAME=""

if command -v nvidia-smi &> /dev/null; then
    if nvidia-smi &> /dev/null; then
        GPU_DETECTED=true
        GPU_NAME=$(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null || echo "Unknown GPU")
    fi
fi

echo ""
echo "=== Hardware Detection ==="
if [ "$GPU_DETECTED" = true ]; then
    echo "✅ GPU Detected: $GPU_NAME"
    echo "   Using: GNINA (GPU-accelerated deep learning)"
else
    echo "⚠️  No GPU detected"
    echo "   Using: Vina CPU"
fi
echo ""

# File paths
RECEPTOR_PATH="/data/receptor.pdbqt"
LIGAND_PATH="/data/ligand.pdbqt"
OUTPUT_PATH="/data/output/output_${JOB_ID}.pdbqt"

echo "Receptor: $RECEPTOR_PATH"
echo "Ligand: $LIGAND_PATH"
echo "Output: $OUTPUT_PATH"

# Validate input files
if [ ! -f "$RECEPTOR_PATH" ]; then
    echo "Error: Receptor file not found: $RECEPTOR_PATH"
    exit 1
fi

if [ ! -f "$LIGAND_PATH" ]; then
    echo "Error: Ligand file not found: $LIGAND_PATH"
    exit 1
fi

mkdir -p /data/output

# Docking parameters
CENTER_X=${center_x:-0}
CENTER_Y=${center_y:-0}
CENTER_Z=${center_z:-0}
SIZE_X=${size_x:-22}
SIZE_Y=${size_y:-22}
SIZE_Z=${size_z:-22}
EXHAUSTIVENESS=${exhaustiveness:-8}
NUM_MODES=${num_modes:-9}
ENERGY_RANGE=${energy_range:-3}

echo ""
echo "=== Docking Parameters ==="
echo "Grid Box: Center=($CENTER_X, $CENTER_Y, $CENTER_Z) Size=($SIZE_X, $SIZE_Y, $SIZE_Z)"
echo "Exhaustiveness: $EXHAUSTIVENESS"
echo "Num Modes: $NUM_MODES"
echo ""

# Run docking based on GPU availability
if [ "$GPU_DETECTED" = true ]; then
    # GPU mode: Use GNINA (best for GPU)
    if command -v gnina &> /dev/null; then
        echo "=== Running GNINA (GPU) ==="
        gnina -r "$RECEPTOR_PATH" \
              -l "$LIGAND_PATH" \
              --center_x "$CENTER_X" \
              --center_y "$CENTER_Y" \
              --center_z "$CENTER_Z" \
              --size_x "$SIZE_X" \
              --size_y "$SIZE_Y" \
              --size_z "$SIZE_Z" \
              --num_modes "$NUM_MODES" \
              --exhaustiveness "$EXHAUSTIVENESS" \
              --cnn_scoring rescore \
              -o "$OUTPUT_PATH" \
              --log /data/output/gnina_log_${JOB_ID}.txt
        
        ENGINE="GNINA"
    elif command -v vina &> /dev/null; then
        # Fallback: Vina (GPU version via pip includes GPU support)
        echo "=== Running Vina (GPU fallback) ==="
        
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
        
        vina --config /tmp/vina_config.txt \
             --receptor "$RECEPTOR_PATH" \
             --ligand "$LIGAND_PATH" \
             --out "$OUTPUT_PATH" \
             --log /data/output/vina_log_${JOB_ID}.txt
        
        ENGINE="Vina-GPU"
    else
        echo "Error: No docking engine found (tried GNINA, Vina)"
        exit 1
    fi
else
    # CPU mode: Use Vina CPU
    if command -v vina &> /dev/null; then
        echo "=== Running Vina (CPU) ==="
        
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
        
        vina --config /tmp/vina_config.txt \
             --receptor "$RECEPTOR_PATH" \
             --ligand "$LIGAND_PATH" \
             --out "$OUTPUT_PATH" \
             --log /data/output/vina_log_${JOB_ID}.txt
        
        ENGINE="Vina-CPU"
    else
        echo "Error: Vina not found"
        exit 1
    fi
fi

# Check output
if [ -f "$OUTPUT_PATH" ]; then
    echo ""
    echo "=========================================="
    echo "✅ Docking completed successfully!"
    echo "   Engine: $ENGINE"
    echo "   Output: $OUTPUT_PATH"
    echo "=========================================="
    exit 0
else
    echo "Error: Output file not generated"
    exit 1
fi

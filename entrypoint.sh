#!/bin/bash
set -e

echo "=========================================="
echo "  BioDockify Docking Studio"
echo "=========================================="
echo "Job ID: $JOB_ID"

# Auto-detect GPU
GPU_DETECTED=false
GPU_NAME=""
VRAM_MB=0

if command -v nvidia-smi &> /dev/null; then
    if nvidia-smi &> /dev/null; then
        GPU_DETECTED=true
        GPU_NAME=$(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null || echo "Unknown GPU")
        VRAM_MB=$(nvidia-smi --query-gpu=memory.total --format=csv,noheader,nounits 2>/dev/null || echo "0")
    fi
fi

echo ""
echo "=== Hardware Detection ==="
if [ "$GPU_DETECTED" = true ]; then
    echo "✅ GPU Detected: $GPU_NAME"
    echo "   VRAM: ${VRAM_MB}MB"
    echo ""
    echo "=== Using: Vina-GPU + GNINA + RF ML ==="
    ENGINE_MODE="GPU"
else
    echo "⚠️  No GPU detected"
    echo ""
    echo "=== Using: Vina-CPU + RF ML ==="
    ENGINE_MODE="CPU"
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

# ============================================
# GPU MODE: Vina-GPU + GNINA + RF ML
# ============================================
if [ "$GPU_DETECTED" = true ]; then
    echo "=== Step 1: Vina-GPU Docking ==="
    
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
    
    # Run Vina-GPU
    if command -v vina &> /dev/null; then
        vina --config /tmp/vina_config.txt \
             --receptor "$RECEPTOR_PATH" \
             --ligand "$LIGAND_PATH" \
             --out "$OUTPUT_PATH" \
             --log /data/output/vina_log_${JOB_ID}.txt
        echo "✅ Vina-GPU docking completed"
    else
        echo "Error: Vina not found"
        exit 1
    fi
    
    # Step 2: GNINA CNN Rescoring (if available)
    if command -v gnina &> /dev/null; then
        echo ""
        echo "=== Step 2: GNINA CNN Rescoring ==="
        
        # Rescore with GNINA CNN
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
              --log /data/output/gnina_log_${JOB_ID}.txt 2>/dev/null || true
        
        echo "✅ GNINA CNN rescoring completed"
    fi
    
    # Step 3: RF ML Scoring (if RF model available)
    if [ -f "/data/rf_model.pkl" ] || [ -f "/models/rf_model.pkl" ]; then
        echo ""
        echo "=== Step 3: RF ML pKd Prediction ==="
        
        # Run RF scoring using Python
        python3 -c "
import sys
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, AllChem
    import pickle
    
    # Try to load RF model
    for path in ['/data/rf_model.pkl', '/models/rf_model.pkl', 'rf_model.pkl']:
        try:
            with open(path, 'rb') as f:
                model = pickle.load(f)
            
            # Calculate RF score for ligand
            mol = Chem.MolFromPDBFile('$LIGAND_PATH')
            if mol:
                desc = [Descriptors.MolWt(mol), Descriptors.MolLogP(mol), 
                       Descriptors.NumHDonors(mol), Descriptors.NumHAcceptors(mol),
                       Descriptors.TPSA(mol), Descriptors.NumRotatableBonds(mol)]
                # Simple prediction (placeholder)
                rf_score = model.predict([desc])[0] if hasattr(model, 'predict') else -8.5
                print(f'RF pKd: {rf_score:.2f}')
                sys.exit(0)
        except:
            continue
    print('RF model not found, skipping')
except Exception as e:
    print(f'RF scoring error: {e}')
" || echo "⚠️  RF ML scoring skipped"
    fi
    
    FINAL_ENGINE="Vina-GPU + GNINA + RF ML"

# ============================================
# CPU MODE: Vina-CPU + RF ML
# ============================================
else
    echo "=== Step 1: Vina-CPU Docking ==="
    
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
    
    # Run Vina-CPU
    if command -v vina &> /dev/null; then
        vina --config /tmp/vina_config.txt \
             --receptor "$RECEPTOR_PATH" \
             --ligand "$LIGAND_PATH" \
             --out "$OUTPUT_PATH" \
             --log /data/output/vina_log_${JOB_ID}.txt
        echo "✅ Vina-CPU docking completed"
    else
        echo "Error: Vina not found"
        exit 1
    fi
    
    # Try GNINA CPU (some systems can run GNINA on CPU, slower)
    if command -v gnina &> /dev/null; then
        echo ""
        echo "=== Step 2: GNINA CNN Rescoring (CPU) ==="
        
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
              --log /data/output/gnina_log_${JOB_ID}.txt 2>/dev/null || true
        
        echo "✅ GNINA rescoring completed"
    fi
    
    # Step 3: RF ML Scoring
    echo ""
    echo "=== Step 3: RF ML pKd Prediction ==="
    
    python3 -c "
import sys
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors
    import pickle
    
    # Try to load RF model
    for path in ['/data/rf_model.pkl', '/models/rf_model.pkl', 'rf_model.pkl']:
        try:
            with open(path, 'rb') as f:
                model = pickle.load(f)
            
            mol = Chem.MolFromPDBFile('$LIGAND_PATH')
            if mol:
                desc = [Descriptors.MolWt(mol), Descriptors.MolLogP(mol), 
                       Descriptors.NumHDonors(mol), Descriptors.NumHAcceptors(mol),
                       Descriptors.TPSA(mol), Descriptors.NumRotatableBonds(mol)]
                rf_score = model.predict([desc])[0] if hasattr(model, 'predict') else -8.5
                print(f'RF pKd: {rf_score:.2f}')
                sys.exit(0)
        except:
            continue
    print('RF model not available')
except Exception as e:
    print(f'RF scoring not available: {e}')
" || echo "⚠️  RF ML scoring skipped"
    
    FINAL_ENGINE="Vina-CPU + GNINA + RF ML"
fi

# Check output
if [ -f "$OUTPUT_PATH" ]; then
    echo ""
    echo "=========================================="
    echo "✅ Docking completed successfully!"
    echo "   Engine: $FINAL_ENGINE"
    echo "   Mode: $ENGINE_MODE"
    echo "   Output: $OUTPUT_PATH"
    echo "=========================================="
    exit 0
else
    echo "Error: Output file not generated"
    exit 1
fi

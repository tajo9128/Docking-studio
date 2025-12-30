# User Guide - BioDockify Docking Studio v1.0.0

## Overview
BioDockify Docking Studio is a desktop interface for performing molecular docking simulations using AutoDock Vina. It abstracts the complexity of command-line tools into a user-friendly GUI.

## Workflow

### 1. New Project
- Click **New Job** on the dashboard.
- **Receptor**: Upload a clean `.pdb` or `.pdbqt` file.
- **Ligand**: Upload a `.sdf`, `.mol2`, or `.pdb` file.

### 2. Configuration
- **Center (X, Y, Z)**: Coordinates of the binding pocket.
- **Size (X, Y, Z)**: Dimensions of the search box in Angstroms.
- **Exhaustiveness**: Higher values (8-32) increase accuracy but take longer.

### 3. Execution
- Click **Start Docking**.
- The job is submitted to a local Docker container (`biodockify-job-{uuid}`).
- Monitor progress in the "Job Status" panel. `Agent Zero` will automatically handle minor errors.

### 4. Analysis
- Upon completion, results are displayed in the **Results** tab.
- **Binding Energy**: Predicted affinity (kcal/mol).
- **Confidence Score**: AI-driven reliability metric (0-100).
- **Interactions**: List of key residue interactions.

## Features
- **Auto-Recovery**: Agent Zero detects failures (e.g., timeouts) and retries with optimized parameters.
- **History**: View past jobs in the "Recent History" panel.
- **Logs**: Real-time logs available in the status widget.

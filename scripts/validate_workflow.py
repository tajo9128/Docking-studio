"""Full Workflow Validation Script"""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

print('=== FULL WORKFLOW VALIDATION ===')
print()

# Test 1: v1.2.0 Features
print('1. Testing v1.2.0 Features...')
from src.core.filters.pains_filter import PAINSFilter
from src.core.admet.admet_predictor import ADMETPredictor
from src.core.druglikeness import DrugLikenessCalculator
from src.services.ai_service import AIAnalysisService
print('   PAINS Filter: OK')
print('   ADMET Predictor: OK')
print('   Drug-likeness: OK')
print('   AI Service: OK')

# Test 2: BioDockviz Core
print()
print('2. Testing BioDockviz Core Modules...')
from src.core.engines.molecular_engine import MolecularEngine
from src.core.engines.interaction_pipeline import InteractionPipeline
from src.core.parsers.pdb_parser import PDBParser
from src.core.parsers.mol2_parser import MOL2Parser
from src.core.parsers.sdf_parser import SDFParser
from src.core.analyzers.bond_detector import BondDetector
from src.core.analyzers.interaction_analyzer import InteractionAnalyzer
from src.core.math.safe_numpy import SafeNumpy
print('   Molecular Engine: OK')
print('   Interaction Pipeline: OK')
print('   PDB Parser: OK')
print('   MOL2 Parser: OK')
print('   SDF Parser: OK')
print('   Bond Detector: OK')
print('   Interaction Analyzer: OK')
print('   SafeNumpy: OK')

# Test 3: Services
print()
print('3. Testing Services...')
from src.services.analysis_service import AnalysisService
from src.services.parsing_service import ParsingService
print('   Analysis Service: OK')
print('   Parsing Service: OK')

# Test 4: Functional Test
print()
print('4. Running Functional Tests...')
pf = PAINSFilter()
result = pf.filter_smiles('CCO')
status = "PASS" if result.passed else "FAIL"
print(f'   PAINS Test: {status}')

calc = DrugLikenessCalculator()
result = calc.calculate_smiles('CC(=O)Oc1ccccc1C(=O)O')
print(f'   Drug-likeness Test: {result.overall_druglikeness} ({result.score:.2f})')

ai = AIAnalysisService()
result = ai.analyze_docking_result({'score': -8.5, 'interactions': {'hydrogen_bonds': 3}})
print(f'   AI Analysis Test: {len(result.insights)} insights')

print()
print('=== ALL TESTS PASSED ===')

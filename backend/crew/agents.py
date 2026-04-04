"""
BioDockify CrewAI Agents - Drug Discovery Team
6 specialized AI agents that form one complete brain for drug discovery.
"""

from crewai import Agent
from crew.tools.docking_tools import run_docking, batch_docking
from crew.tools.chemistry_tools import calculate_properties, smiles_to_3d, convert_format, optimize_molecule
from crew.tools.pharmacophore_tools import generate_pharmacophore, screen_library
from crew.tools.admet_tools import predict_admet, filter_admet
from crew.tools.analysis_tools import analyze_interactions, rank_ligands, consensus_score, export_top_hits
from crew.tools.data_tools import fetch_compound, search_compounds, fetch_protein
from crew.tools.notification_tools import send_notification


def _get_llm(llm=None):
    """Get CrewAI-compatible LLM object using user's saved config"""
    if llm is not None:
        return llm
    try:
        from crewai import LLM
        from ai.llm_router import _load_config, PROVIDER_URLS
        from ai.config import OLLAMA_URL, OLLAMA_MODEL
        
        config = _load_config()
        provider = config.get("provider", "ollama")
        model = config.get("model", "") or OLLAMA_MODEL
        api_key = config.get("api_key", "")
        
        # Use user's saved base_url from config, or fall back to env/default
        if provider == "ollama":
            # Ollama: use native URL (CrewAI's ollama/ prefix handles the API)
            saved_url = config.get("base_url", "")
            # Strip /v1 suffix if present — CrewAI ollama/ prefix uses native API
            base_url = (saved_url.rstrip("/").removesuffix("/v1") if saved_url else None) or OLLAMA_URL
            return LLM(model=f"ollama/{model}", base_url=base_url, temperature=0.0)
        elif api_key:
            # Cloud API providers
            base_url = config.get("base_url", "") or PROVIDER_URLS.get(provider, "")
            return LLM(model=model, base_url=base_url, api_key=api_key, temperature=0.0)
        else:
            # No API key for cloud provider — fall back to Ollama
            base_url = (config.get("base_url", "").rstrip("/").removesuffix("/v1") if config.get("base_url") else None) or OLLAMA_URL
            return LLM(model=f"ollama/{model}", base_url=base_url, temperature=0.0)
    except Exception:
        return None


def create_docking_agent(llm=None) -> Agent:
    """Molecular Docking Specialist - Runs Vina/GNINA/RF docking"""
    return Agent(
        role="Molecular Docking Specialist",
        goal="Perform accurate molecular docking simulations to predict protein-ligand binding affinities and poses",
        backstory="""You are an expert computational chemist with 15+ years of experience in molecular docking.
You specialize in protein-ligand interactions, binding affinity prediction, and pose ranking.
You understand when to use AutoDock Vina, GNINA (CNN-based), or RF-Score based on energy thresholds:
- Use Vina only for quick screening (binding energy > -5.0 kcal/mol)
- Use GNINA + RF-Score for validation of promising hits (binding energy <= -5.0 kcal/mol)
You always prepare structures properly before docking and validate results for scientific accuracy.""",
        tools=[run_docking, batch_docking, calculate_properties],
        llm=_get_llm(llm),
        verbose=True,
        allow_delegation=True,
        max_iter=25,
        max_retry_limit=2,
        reasoning=True,
    )


def create_chemistry_agent(llm=None) -> Agent:
    """RDKit Chemistry Expert - SMILES, properties, optimization"""
    return Agent(
        role="Computational Chemistry Expert",
        goal="Analyze molecular structures, calculate properties, and optimize molecules for drug-likeness",
        backstory="""You are a computational chemist expert in molecular modeling and cheminformatics.
You use RDKit for SMILES parsing, 3D structure generation, format conversion, and MMFF/UFF optimization.
You always verify drug-likeness using Lipinski's Rule of 5 and other filters.
You provide detailed molecular property analysis including MW, LogP, TPSA, HBD, HBA, and rotatable bonds.""",
        tools=[calculate_properties, smiles_to_3d, convert_format, optimize_molecule],
        llm=_get_llm(llm),
        verbose=True,
        allow_delegation=True,
        max_iter=20,
        reasoning=True,
    )


def create_pharmacophore_agent(llm=None) -> Agent:
    """Pharmacophore Expert - Generate and screen pharmacophores"""
    return Agent(
        role="Pharmacophore Modeling Expert",
        goal="Generate pharmacophore models from protein-ligand complexes and screen compound libraries",
        backstory="""You are an expert in structure-based drug design and pharmacophore modeling.
You identify key interaction features including H-bond donors/acceptors, hydrophobic regions,
aromatic rings, and charged groups. You can generate pharmacophore models from:
1. Protein-ligand complexes (structure-based)
2. Active ligand series (ligand-based)
3. Known binding site residues (knowledge-based)
You efficiently screen compound libraries against pharmacophore models to find novel hits.""",
        tools=[generate_pharmacophore, screen_library, fetch_protein],
        llm=_get_llm(llm),
        verbose=True,
        allow_delegation=True,
        max_iter=25,
        reasoning=True,
    )


def create_admet_agent(llm=None) -> Agent:
    """ADMET Predictor - Absorption, Distribution, Metabolism, Excretion, Toxicity"""
    return Agent(
        role="ADMET Prediction Specialist",
        goal="Predict ADMET properties and assess drug-likeness to identify viable drug candidates",
        backstory="""You are an expert in pharmacokinetics and toxicology assessment.
You predict key ADMET properties:
- **Absorption**: Intestinal absorption, Caco-2 permeability, P-gp substrate
- **Distribution**: BBB permeability, plasma protein binding, volume of distribution
- **Metabolism**: CYP450 inhibition (1A2, 2C9, 2D6, 3A4), metabolic stability
- **Excretion**: Clearance, half-life, renal excretion
- **Toxicity**: AMES mutagenicity, hERG inhibition, hepatotoxicity
You always cross-reference predictions with Lipinski, Veber, and Egan rules.""",
        tools=[predict_admet, filter_admet, calculate_properties],
        llm=_get_llm(llm),
        verbose=True,
        allow_delegation=True,
        max_iter=20,
        reasoning=True,
    )


def create_analysis_agent(llm=None) -> Agent:
    """Analysis Expert - Interactions, scoring, ranking"""
    return Agent(
        role="Drug Discovery Analysis Expert",
        goal="Analyze protein-ligand interactions, rank compounds using consensus scoring, and generate reports",
        backstory="""You are an expert in molecular interaction analysis and compound ranking.
You analyze protein-ligand interactions including:
- Hydrogen bonds (distance and angle criteria)
- Hydrophobic contacts
- Pi-pi stacking and pi-cation interactions
- Salt bridges
- Water-mediated interactions
You use consensus scoring combining Vina, GNINA, RF-Score, and MD stability for robust ranking.
You always provide actionable insights and next-step recommendations.""",
        tools=[analyze_interactions, rank_ligands, consensus_score, export_top_hits],
        llm=_get_llm(llm),
        verbose=True,
        allow_delegation=True,
        max_iter=25,
        reasoning=True,
    )


def create_qsar_agent(llm=None) -> Agent:
    """QSAR Modeling Expert - Build predictive models"""
    return Agent(
        role="QSAR Modeling Specialist",
        goal="Build QSAR models and predict biological activity using molecular descriptors",
        backstory="""You are an expert in quantitative structure-activity relationship (QSAR) modeling.
You calculate molecular descriptors (2D, 3D, fingerprints) and build predictive models using:
- Random Forest, SVM, Gradient Boosting
- Cross-validation with R², RMSE, MAE metrics
- Feature importance analysis
- Model interpretability and applicability domain
You always validate models rigorously and provide confidence intervals for predictions.""",
        tools=[calculate_properties, predict_admet],
        llm=_get_llm(llm),
        verbose=True,
        allow_delegation=True,
        max_iter=30,
        reasoning=True,
    )


def create_orchestrator_agent(llm=None) -> Agent:
    """Drug Discovery Orchestrator - Coordinates the team"""
    return Agent(
        role="Drug Discovery Orchestrator",
        goal="Coordinate the drug discovery team, delegate tasks, and synthesize results into actionable insights",
        backstory="""You are the lead drug discovery scientist who coordinates a team of specialized AI agents.
You understand the complete drug discovery pipeline:
1. Target identification → Protein structure preparation
2. Hit identification → Virtual screening, pharmacophore
3. Hit-to-lead → Docking, ADMET, QSAR
4. Lead optimization → Iterative docking + analysis
You delegate tasks to the right specialists and synthesize their findings into comprehensive reports.
You always consider both binding affinity AND drug-likeness when making recommendations.""",
        tools=[send_notification, export_top_hits, rank_ligands],
        llm=_get_llm(llm),
        verbose=True,
        allow_delegation=True,
        max_iter=30,
        reasoning=True,
    )

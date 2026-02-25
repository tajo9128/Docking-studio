"""
Offline Assistant - Always available fallback
Deterministic responses for common docking questions
"""

import logging
from typing import Dict, List

logger = logging.getLogger(__name__)


class OfflineAssistant:
    """
    Offline fallback assistant for when Ollama is unavailable.
    Provides deterministic responses to common molecular docking questions.
    """
    
    KNOWLEDGE_BASE = {
        "vina": {
            "keywords": ["vina", "autodock", "binding affinity", "kcal/mol"],
            "response": "AutoDock Vina calculates binding affinity in kcal/mol. More negative values indicate stronger predicted binding. Typical drug-like molecules bind with -5 to -12 kcal/mol."
        },
        "gnina": {
            "keywords": ["gnina", "cnn", "convolutional", "deep learning", "neural network"],
            "response": "GNINA uses deep learning (CNN) to evaluate pose quality. CNN scores range from 0 to 1, where higher values indicate better predicted binding. It combines traditional scoring with AI-powered pose assessment."
        },
        "rf": {
            "keywords": ["rf", "random forest", "score", "scoring"],
            "response": "Random Forest (RF) scoring uses an ensemble of decision trees to predict binding affinity. It provides robust scoring by averaging predictions from multiple tree models."
        },
        "consensus": {
            "keywords": ["consensus", "combine", "average", "ensemble"],
            "response": "Consensus scoring combines Vina, GNINA CNN, and RF scores to provide more reliable predictions. The final score is typically an average or weighted combination of individual scores."
        },
        "hbond": {
            "keywords": ["hydrogen bond", "hbond", "h-bond", "donor", "acceptor"],
            "response": "Hydrogen bonds require: 1) Distance ≤ 3.5 Å between donor and acceptor, 2) Angle ≥ 120° at the hydrogen atom. They are crucial for ligand-receptor specificity."
        },
        "hydrophobic": {
            "keywords": ["hydrophobic", "hydrophobicity", "lipophilicity", "pi stacking"],
            "response": "Hydrophobic interactions occur between non-polar regions. They drive the burial of hydrophobic ligand groups into hydrophobic receptor pockets, contributing significantly to binding free energy."
        },
        "rmsd": {
            "keywords": ["rmsd", "root mean square", "deviation", "similarity"],
            "response": "RMSD (Root Mean Square Deviation) measures structural similarity. RMSD < 2 Å indicates similar poses. It's calculated by averaging the squared distance between corresponding atoms."
        },
        "pocket": {
            "keywords": ["pocket", "binding site", "active site", "cavity", " cleft"],
            "response": "Binding pockets are regions on the protein surface where ligands bind. They often contain key residues for ligand recognition. Tools like AutoSite or Fpocket can identify potential binding sites."
        },
        "grid": {
            "keywords": ["grid", "autogrid", "maps", "energy grid"],
            "response": "Grid maps pre-compute interaction energies between the receptor and probe atoms. They speed up docking by allowing rapid energy lookups during ligand pose evaluation."
        },
        "exhaustiveness": {
            "keywords": ["exhaustiveness", "exhaust", "search", "monte carlo"],
            "response": "Exhaustiveness controls docking search depth. Higher values (8-32) explore more pose space but take longer. For virtual screening, 8-16 is usually sufficient."
        },
        "pose": {
            "keywords": ["pose", "conformation", "orientation", "geometry"],
            "response": "A pose is a specific 3D conformation and orientation of a ligand in the binding site. Docking generates multiple poses which are ranked by score."
        },
        "docking": {
            "keywords": ["dock", "docking", "virtual screening"],
            "response": "Molecular docking predicts how a small molecule (ligand) binds to a protein target. It's used in drug discovery for virtual screening and understanding ligand-protein interactions."
        },
        "help": {
            "keywords": ["help", "what can you", "commands"],
            "response": "I can help with questions about: Vina scoring, GNINA CNN, Random Forest, consensus scoring, hydrogen bonds, hydrophobic interactions, RMSD, binding pockets, grid parameters, exhaustiveness, and docking basics."
        }
    }
    
    def __init__(self):
        logger.info("OfflineAssistant initialized")
    
    def respond(self, message: str) -> str:
        """
        Generate a deterministic response based on message keywords.
        
        Args:
            message: User input message
            
        Returns:
            Response string
        """
        message_lower = message.lower()
        
        best_match = None
        best_score = 0
        
        for topic, data in self.KNOWLEDGE_BASE.items():
            keywords = data["keywords"]
            score = sum(1 for kw in keywords if kw in message_lower)
            
            if score > best_score:
                best_score = score
                best_match = data["response"]
        
        if best_match:
            return best_match
        
        return (
            "I'm running in offline mode. I can help with:\n"
            "• Vina scoring and binding affinity\n"
            "• GNINA CNN deep learning scores\n"
            "• Random Forest scoring\n"
            "• Consensus scoring methods\n"
            "• Hydrogen bonds and hydrophobic interactions\n"
            "• RMSD and pose comparison\n"
            "• Binding pocket analysis\n"
            "• Grid and exhaustiveness parameters\n\n"
            "Ask me about any of these topics!"
        )
    
    def get_available_topics(self) -> List[str]:
        """Return list of available topics"""
        return list(self.KNOWLEDGE_BASE.keys())


if __name__ == "__main__":
    assistant = OfflineAssistant()
    
    test_messages = [
        "What is Vina?",
        "How does GNINA work?",
        "What is consensus scoring?",
        "Tell me about hydrogen bonds"
    ]
    
    for msg in test_messages:
        print(f"Q: {msg}")
        print(f"A: {assistant.respond(msg)}\n")

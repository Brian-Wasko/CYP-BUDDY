
# BUDDY 🧬
### Bioinformatic Utility for Diagnostic Discovery in Yeast

**BUDDY** is a React-based bioinformatics pipeline designed to bridge the gap between human medical genetics and yeast (*Saccharomyces cerevisiae*) model organism research. 
It streamlines the process of assessing the functional impact of human Variants of Unknown Significance (VUS) by mapping them to yeast orthologs, visualizing structural conservation, designing CRISPR editing strategies, and generating AI-assisted experimental plans.
![React](https://img.shields.io/badge/built%20with-React-61DAFB.svg)
![Gemini](https://img.shields.io/badge/AI-Google%20Gemini-8E75B2.svg)

## Key Features

*   **Ortholog Mapping:** Instantly maps Human genes to Yeast orthologs (and vice-versa) using the **DIOPT** integrative score to ensure evolutionary relevance.
*   **Sequence Alignment:** Performs Pairwise Sequence Alignment (Global/Local) to map human amino acid residues to their yeast equivalents.
*   **Variant Analysis:** 
    *   Fetches ClinVar data and **AlphaMissense** pathogenicity scores.
    *   Calculates "Local Homology Scores" to determine if a variant occurs in a conserved region.
    *   Allows manual entry of specific mutations.
*   **3D Structure Visualization:** 
    *   Interactive 3D viewer using **3Dmol.js**.
    *   Fetches structures from **AlphaFold** and **RCSB PDB**.
    *   **Structural Superposition:** Overlays Human and Yeast protein structures to visualize spatial conservation of specific residues.
*   **CRISPR/Cas9 Design:** 
    *   Automated design of repair templates and sgRNAs for yeast gene editing.
    *   Optimized for the **pML104** plasmid system.
    *   Includes synonymous mutation generation for PAM disruption.
*   **AI Experimental Assistant:** 
    *   Uses **Google Gemini 3** models.
    *   Scrapes literature for specific yeast phenotypes.
    *   Generates a tailored, safety-conscious AI-generated experimental assay protocol.

## Getting Started

### Prerequisites

*   Node.js (v16 or higher)
*   npm or yarn
*   A Google Gemini API Key (Get one [here](https://aistudio.google.com/))

### Installation

1.  **Clone the repository**
    ```bash
    git clone https://github.com/yourusername/buddy-pipeline.git
    cd buddy-pipeline
    ```

2.  **Install dependencies**
    ```bash
    npm install
    ```

3.  **Configure API Key**
    Create a `.env` file in the root directory and add your Google Gemini API key:
    ```env
    API_KEY=your_actual_api_key_here
    ```
    *Note: The application uses `process.env.GEMINI_API_KEY` to authenticate with Google's GenAI SDK.*

4.  **Run the application**
    ```bash
    npm start
    ```
    Open [http://localhost:1234](http://localhost:1234) (or the port specified by your bundler) to view it in the browser.

##  Usage Guide

### 1. Select Input Mode
*   **Manual Input:** Enter a specific Human or Yeast gene symbol (e.g., `ATP6V1B1` or `VMA2`).
*   **Rare Diseases:** Browse a curated list of rare disease genes with high-confidence yeast orthologs.
*   **Search by Topic:** Use AI to find genes related to specific biological concepts (e.g., "Autophagy", "Mitochondria").

### 2. Analyze Conservation
*   Review the **DIOPT Score** to assess orthology confidence.
*   Examine the **Sequence Alignment** view. Variants are marked on the sequence.
*   Use the **3D Structure Viewer** to overlay human and yeast proteins and see where the mutation lies in 3D space.

### 3. Design Experiment
*   **Select a Variant:** Click a variant in the table or alignment view.
*   **Generate CRISPR Oligos:** The app fetches genomic DNA from Ensembl and generates oligonucleotide sequences for editing the yeast genome (w/pML104 or similar plasmid) to mimic the human variant.
*   **Generate AI Plan:** Click "Generate AI-assisted Experimental Plan". Gemini will propose a specific functional assay (e.g., growth on specific media) based on known null phenotypes.

##  Tech Stack & APIs

**Frontend:**
*   React 19
*   TypeScript
*   Tailwind CSS (Styling & Dark Mode)
*   Lucide React (Icons)

**Bioinformatics APIs:**
*   **MyGene.info:** Gene metadata and ID conversion.
*   **DIOPT (DRSC):** Ortholog prediction (a local file with genes fit to a best ortholog (via DIOPT score) with human/yeast gene pairs with DIOPT >2 is now used).
*   **UniProt:** Protein sequences.
*   **Ensembl REST:** Yeast genomic DNA sequences.
*   **MyVariant.info:** ClinVar and dbNSFP (AlphaMissense) data.
*   **Alliance of Genome Resources:** Yeast phenotypes.
*   **AlphaFold / RCSB:** 3D Protein structures.

**AI:**
*   **Google GenAI SDK:** Powered by current Gemini 3 Pro/Flash models.

## ⚠️ Disclaimers

**AI Use:**
Generative AI was used to help make BUDDY, to help make this readme, is used within BUDDY.

**For Research Use Only.** 
BUDDY is a research utility. The experimental plans, primer designs, and variant interpretations are generated using bioinformatics algorithms and AI models. They should **not** be used for clinical diagnosis or treatment decisions without independent verification. Practice appropriate laboratory safety and have a qualified scientist supervise experiements.

##  License

TBD

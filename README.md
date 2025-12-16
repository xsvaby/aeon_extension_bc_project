# BC_project -- Gene ontlogy enrichment anaylsing on partially specified boolean network

A computational pipeline for analyzing biological signaling pathways
(in form partially specified boolean networks) using **AEON models**, 
with automated **enrichment analysis**, **Excel-based reporting**,
and with **Postprocessing visualization**.

------------------------------------------------------------------------
## Installation

This project was developed and tested on Google Colab. Follow the instructions below to run it either on Colab (recommended) or locally.

**Option 1: Run on Google Colab (Recommended)**
Google Colab already provides most system dependencies and works out of the box.
Step 1: Open a New Colab Notebook - https://colab.research.google.com
Step 2: Install Dependencies
- !pip install biodivine_aeon==1.3.0a3 -qq
- !pip install SPARQLWrapper networkx matplotlib -qq
- !apt-get install graphviz graphviz-dev -y -qq
- !pip install pygraphviz -qq

Step 3: Import Modules
from biodivine_aeon import *

from EnrichmentClasses import EnrichmentPSBN

from Visualization import (
    get_quickgo_terms_batch,
    set_nodes_for_graph,
    make_graph,
    visualize_subgraphs,
    get_roots_and_leafs,
    sort_roots_and_leafs,
    visualize_subgraphs_on_whole_net,
    print_roots_and_leafs_on_whole_net,
    visualize_subgraphs_on_each_instance,
    print_roots_and_leafs_per_instance,
    visualize_unmapped_nodes_frequencies,
    visualize_mapped_nodes_frequencies,
)

from Pipeline import pipeline
from CheckFixedPoints import print_fixed_points_on_new_stg


You are now ready to run the pipeline.

**Option 2: Local Installation (Linux / macOS / WSL)**
Note: Native Windows is not supported due to pygraphviz.
Use WSL2, Docker, or Google Colab instead.

1. System Requirements
- Python 3.9 – 3.11
- pip
Graphviz system libraries

2. Install System Dependencies
Ubuntu / Debian / WSL:
- sudo apt update
- sudo apt install -y graphviz graphviz-dev python3-dev

macOS (Homebrew):
- brew install graphviz

3. (Optional) Create a Virtual Environment
- python3 -m venv venv
- source venv/bin/activate
- pip install --upgrade pip

4. Install Python Dependencies
- pip install biodivine_aeon==1.3.0a3
- pip install SPARQLWrapper networkx matplotlib
- pip install pygraphviz

5. Verify Installation
- import biodivine_aeon
- import pygraphviz
- import networkx

If no errors occur, the installation was successful.

**Troubleshooting**
pygraphviz Installation Fails

This usually means Graphviz headers are missing.

Fix (Ubuntu / WSL):

- sudo apt install graphviz graphviz-dev
- pip install --no-cache-dir pygraphviz

------------------------------------------------------------------------
## Project Structure

    BC_project/
    │
    ├── case_studies/
    │   ├── deathReceptorSignaling_unspecified_pathway.ipynb
    │   ├── interferon_pathway.ipynb
    │   ├── mapk_unspecified_pathway.ipynb
    │   └── test_case.ipynb
    │
    ├── data/
    │   ├── deathReceptorSignaling-unspecified.aeon
    │   ├── interferon(virus_true_drugs_false).aeon
    │   ├── mapk-unspecified.aeon
    │   └── test_case.aeon
    │
    ├── output_excel_files/
    │   ├── dr_OnAllInstances.xlsx
    │   ├── dr_OnAttractors.xlsx
    │   ├── dr_OnInstance.xlsx
    │   ├── interferon_OnAllInstances.xlsx
    │   ├── interferon_OnAttractors.xlsx
    │   ├── interferon_OnInstance.xlsx
    │   ├── mapk_OnAllInstances.xlsx
    │   ├── mapk_OnAttractors.xlsx
    │   ├── mapk_OnInstance.xlsx
    │   ├── test_case_OnAllInstances.xlsx
    │   ├── test_case_OnAttractors.xlsx
    │   └── test_case_OnInstance.xlsx
    │
    ├── src/
    │   ├── CheckFixedPoints.py
    │   ├── Enrichment.py
    │   ├── EnrichmentClasses.py
    │   ├── Pipeline.py
    │   └── Visualization.py
    |
    ├── LICENSE
    |
    └── README.md

------------------------------------------------------------------------
## Folder Overview

### `case_studies/`

Interactive **Jupyter notebooks** used to run and validate analyses: -
Death receptor pathway - Interferon signaling - MAPK pathway - Test pathway

------------------------------------------------------------------------

### `data/`

All **AEON biological network models**: - Used as input to the
computational pipeline - Each model represents a regulatory signaling
pathway

------------------------------------------------------------------------

### `output_excel_files/`

Pipeline generated **Excel reports**:
- `*_OnAllInstances.xlsx` → Whole network intersection of GO terms
- `*_OnInstance.xlsx` → Single network instance intersections of GO terms; 1 column = 1 instance
- `*_OnAttractors.xlsx` → Per attractor per instance GO terms; 1 column = 1 attractor, in the column name -> number of color (instance) + number of attractor

"*" stands for shortcut name of input pathways in case studies

------------------------------------------------------------------------

### `src/`

Core Python implementation:

  ------------------------ ----------------------------
  - `Pipeline.py`            Main flow
  - `CheckFixedPoints.py`    Fixed-point computation
  - `Enrichment.py`          Enrichment analysis
  - `EnrichmentClasses.py`   Supporting data structures
  - `Visualization.py`       Graphing & visualizing results

------------------------------------------------------------------------
## Workflow

1. Load AEON model from `data/`
2. Run a notebook from `case_studies/`
3. Pipeline performs:
    -   Instance evaluation
    -   Attractor identification
    -   Enrichment analysis
4. Results exported to `output_excel_files/`
5. Results postprocessed and printed in the python notebook

------------------------------------------------------------------------
## License

Distributed under the terms defined in the `LICENSE` file.

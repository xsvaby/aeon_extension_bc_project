# BC_project -- Gene ontlogy enrichment anaylsing on partially specified boolean network

A computational pipeline for analyzing biological signaling pathways
(in form partially specified boolean networks) using **AEON models**, 
with automated **enrichment analysis**, **Excel-based reporting**,
and with **Postprocessing visualization**.

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

Automatically generated **Excel reports**: - `*_OnAllInstances.xlsx` →
Whole network intersetion - `*_OnAttractors.xlsx` → Only attractor
states per instance - `*_OnInstance.xlsx` → Single network instances intersections

------------------------------------------------------------------------

### `src/`

Core Python implementation:

  File                     Purpose
  ------------------------ ----------------------------
  `Pipeline.py`            Main logic
  `CheckFixedPoints.py`    Fixed-point computation
  `Enrichment.py`          Enrichment analysis
  `EnrichmentClasses.py`   Supporting data structures
  `Visualization.py`       Graphing & visualizing results

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

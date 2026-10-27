\# Execution environment



This folder contains the Spatial-mIFlow source code and the local Sopa workflow used during development. Container images, execution caches, and generated analysis outputs are excluded.



The Sopa modules reference versioned software containers, including:



\- `quentinblampey/sopa:2.1.11`

\- `quentinblampey/sopa:2.1.11-cellpose`



Running the workflow requires Nextflow, Singularity, input data, and suitable computing resources. With the appropriate configuration and registry access, Nextflow can obtain and cache the referenced containers. On restricted clusters, obtain the containers beforehand and configure their location.



Before running, review the launcher and configuration files. Adapt input paths, output directories, scheduler settings, and container cache locations to your environment.



This workflow was developed on an institutional cluster. Execution from a fresh clone in another environment has not yet been validated.


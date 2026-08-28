Spatial-mIFlow workflow

A modular Nextflow workflow for large-scale multiplex immunofluorescence image analysis from OME-TIFF images to SOPA segmentation, single-cell AnnData objects, RAPIDS clustering, SCIMAP-style gating, contingency analysis, segmentation export, and spatial analysis.

This project was designed for COMET / Lunaphore-like multiplex IF images, but the architecture is intended to be portable across HPC systems when the input images are OME-TIFF and the channel configuration is correctly described in params.yml.

1. Project purpose

The workflow automates the main technical and computational steps required to process large multiplex IF images:

Prepare raw OME-TIFF images for SOPA.

Add a composite_MIP channel for segmentation support.

Run SOPA / Cellpose segmentation.

Collect SOPA single-cell measurements as AnnData.

Perform DAPI-based quality control and RAPIDS/Scanpy clustering.

Export SOPA segmentation polygons to .npz for downstream visualization.

Run SCIMAP-style automatic marker gating.

Explore clusters, marker signatures, and UMAPs.

Compare clusters with gated populations.

Run spatial analysis on clustered or gated cell populations.

The current project is organized around a reusable pipeline code folder, Spatial-mIFlow/, and a user-facing project root that contains the sample sheet, parameters, output folder, launcher, and deployment-level configuration.

2. Current repository layout

Development_COMET_imaging_pipeline/
├── master_samplesheet.csv              # User sample sheet
├── params.yml                          # User-facing scientific/runtime parameters
├── nextflow.config                     # Deployment-level engine/cache configuration
├── run_pipeline_slurm.sh               # SLURM launcher
├── data/                               # Raw input images
├── outputs/                            # Generated results
├── logs/                               # SLURM stdout/stderr logs
├── Spatial-mIFlow/                          # Reusable pipeline code
│   ├── main.nf                         # Main DSL2 workflow
│   ├── nextflow.config                 # Pipeline-level process/resource config
│   ├── modules/                        # Independent workflow modules
│   └── programs_hierarchy/             # Embedded SOPA pipeline/configs
└── archive/                            # Old code / backups, not part of runtime

The important design decision is:

project root      = user deployment layer
Spatial-mIFlow/        = reusable pipeline code
params.yml        = user-editable parameters
outputs/          = generated analysis outputs
scratch/work dir  = Nextflow work, SOPA work, caches, temporary runtime logs

3. Main execution command

Submit the pipeline with SLURM:

sbatch run_pipeline_slurm.sh

Resume a previous run:

sbatch run_pipeline_slurm.sh -resume

Use a custom Nextflow work directory:

sbatch run_pipeline_slurm.sh -w /scratch/$USER/comet_pipeline/work

The launcher calls:

nextflow run Spatial-mIFlow/main.nf \
  -profile slurm \
  -params-file params.yml \
  -w /scratch/$USER/comet_pipeline/work

The launcher also sets environment variables to reduce local Python/Conda conflicts:

export PYTHONNOUSERSITE=1
export PIP_NO_CACHE_DIR=1

4. User-facing configuration: params.yml

The file params.yml is the main file a user should edit. It contains sample paths, output directory, image-preparation options, SOPA parameters, and notebook parameters.

4.1 Required paths

samplesheet: "master_samplesheet.csv"
outdir: "./outputs"

The sample sheet should contain:

sample,data_path
LIVER_region_B,data/Region-B_abII-subtracted.ome.tiff
LIVER_region_D,data/Region-D_abII-subtracted.ome.tiff

Relative data_path values are resolved against the directory where the workflow is launched. In the current launcher, this is the project root.

4.2 Image preparation parameters

image_preparation:
  mip_end: 21
  drop_channel_names: "TRITC,Cy5"
  write_clean_channel_names_to_ome: false
  prepare_skip_if_complete: true

The image-preparation module also supports default values for omitted parameters. For the liver images, DAPI is the first channel and the marker channels are C2–C21, so mip_end: 21 includes all marker channels in Python slicing convention where the end index is excluded.

4.3 SOPA adapter parameters

sopa_adapter:
  required_channels: "composite_MIP,DAPI"
  skip_if_complete: true
  resume: true
  profile: "singularity"

This controls the outer SOPA adapter, not the internal SOPA segmentation parameters themselves.

4.4 SOPA segmentation parameters

sopa:
  technology: ome_tif
  patch_width_pixel: 3000
  patch_overlap_pixel: 80
  use_cellpose: true
  cellpose_model_type: "sam"
  cellpose_channels: "composite_MIP,DAPI"
  cellpose_diameter: 30
  cellprob_threshold: 0
  cellpose_use_gpu: true
  cellpose_max_forks: 8
  aggregate_channels: true
  min_intensity_ratio: 0.2
  min_area_pixels2: 200
  expand_radius_ratio: "0.1"
  pixel_size: 0.28

The user should not put SOPA input: or outdir: here. These are generated automatically by the pipeline into:

outputs/00_runtime/sopa_inputs/sopa_params.generated.yml

This avoids hard-coded cluster paths and makes the workflow portable.

4.5 Notebook parameters

notebooks:
  rapids:
    dapi_qc_percentile: 15
    render_quarto_html: true
    render_quarto_pdf: false

  scimap:
    gmm_components: 3
    log: true
    method: "all"
    random_state: 0
    positivity_threshold: 0.5
    top_n_populations: 30
    top_n_spatial_populations: 15

  cluster_exploration:
    selected_normalizations: "arcsinh,log1p"
    selected_resolutions: "0.5,0.1"
    selected_umap_min_dists: "0.3"
    top_n_markers: 30
    marker_colormap: "RdYlBu_r"

  contingency:
    top_n_populations: 30

  spatial:
    cluster_normalization: "arcsinh"
    cluster_resolution: "0.5"

Notebook parameters are injected at runtime using Papermill. For example, the cluster-exploration process passes selected_resolutions into the notebook as:

-p SELECTED_RESOLUTIONS "${params.notebooks.cluster_exploration.selected_resolutions}"

Changing a notebook parameter and running with -resume should rerun the corresponding notebook process if the injected command changes.

5. Deployment-level nextflow.config

The root nextflow.config is intentionally small. It contains deployment and engine settings, such as Conda and Singularity cache locations:

conda.cacheDir = "/scratch/$USER/comet_pipeline/nextflow_conda"
singularity.cacheDir = "/scratch/$USER/comet_pipeline/singularity"

This file should not contain biological parameters or notebook parameters. Those belong in params.yml.

6. Pipeline code: Spatial-mIFlow/

The reusable workflow code lives in:

Spatial-mIFlow/
├── main.nf
├── nextflow.config
├── modules/
└── programs_hierarchy/

6.1 Main workflow

Spatial-mIFlow/main.nf imports DSL2 modules:

include { PREPARE_COMPOSITE_MIP } from './modules/prepare_composite_mip/main.nf'
include { PUBLISH_PREPARED_IMAGE } from './modules/publish_prepared_image/main.nf'
include { BUILD_SOPA_SAMPLESHEET } from './modules/build_sopa_samplesheet/main.nf'
include { RUN_SOPA_AND_COLLECT_ADATA } from './modules/run_sopa_and_collect_adata/main.nf'
include { RAPIDS_CLUSTERING_NOTEBOOK } from './modules/rapids_clustering_notebook/main.nf'
include { SCIMAP_AUTOGATING_NOTEBOOK } from './modules/scimap_autogating_notebook/main.nf'
include { CLUSTER_EXPLORATION_NOTEBOOK } from './modules/cluster_exploration_notebook/main.nf'
include { CONTINGENCY_CLUSTER_POPULATIONS_NOTEBOOK } from './modules/contingency_cluster_populations_notebook/main.nf'
include { EXPORT_SOPA_CELL_BOUNDARIES_NPZ } from './modules/export_sopa_cell_boundaries_npz/main.nf'
include { SPATIAL_ANALYSIS_CLUSTERS_NOTEBOOK } from './modules/spatial_analysis_clusters_notebook/main.nf'

The workflow dataflow is:

flowchart TD
    A[master_samplesheet.csv] --> B[PREPARE_COMPOSITE_MIP]
    B --> C[PUBLISH_PREPARED_IMAGE]
    C --> D[BUILD_SOPA_SAMPLESHEET]
    D --> E[RUN_SOPA_AND_COLLECT_ADATA]
    E --> F[RAPIDS_CLUSTERING_NOTEBOOK]
    E --> G[EXPORT_SOPA_CELL_BOUNDARIES_NPZ]
    F --> H[CLUSTER_EXPLORATION_NOTEBOOK]
    F --> I[SCIMAP_AUTOGATING_NOTEBOOK]
    I --> J[CONTINGENCY_CLUSTER_POPULATIONS_NOTEBOOK]
    I --> K[SPATIAL_ANALYSIS_CLUSTERS_NOTEBOOK]

7. Module documentation

7.1 PREPARE_COMPOSITE_MIP

Location:

Spatial-mIFlow/modules/prepare_composite_mip/
├── main.nf
├── environment.yml
└── bin/prepare_composite_mip.py

Environment:

name: image_preparation
channels:
  - conda-forge
dependencies:
  - python=3.10
  - numpy
  - pandas
  - tifffile>=2023.7.10
  - imagecodecs

Purpose:

Read large OME-TIFF images in streaming mode.

Avoid loading the full image into RAM.

Optionally drop technical channels using name tokens.

Compute a composite_MIP channel using a one-pass streaming generator.

Write a SOPA-compatible OME-TIFF.

Validate the written OME-TIFF.

Save a channel map CSV and summary JSON.

Important output files:

outputs/01_prepared_images/<sample>/
├── <sample>_plus_composite_MIP.ome.tif
├── <sample>_channel_name_map.csv
└── <sample>_image_preparation_summary.json

The script writes tiled, compressed OME-TIFF output with zstd compression and horizontal predictor. The module validates axes, OME metadata, channel count, and duplicate channel names before considering the file valid.

The image-preparation module can reuse an existing prepared image if prepare_skip_if_complete is true and the existing files pass validation.

7.2 PUBLISH_PREPARED_IMAGE

Purpose:

Publish prepared OME-TIFFs and metadata to outputs/01_prepared_images/.

Provide a manifest for SOPA runtime input generation.

This module separates the preparation work directory from the permanent output directory.

7.3 BUILD_SOPA_SAMPLESHEET

Purpose:

Build a SOPA-compatible sample sheet from prepared image manifests.

Generate a SOPA params file using the sopa: section of params.yml.

Dynamically add SOPA input and outdir paths.

Generated runtime files:

outputs/00_runtime/sopa_inputs/
├── sopa_samplesheet.csv
├── resolved_sopa_channels.json
└── sopa_params.generated.yml

This prevents generated runtime files from being written inside pipeline code folders.

7.4 RUN_SOPA_AND_COLLECT_ADATA

Purpose:

Launch the embedded nf-core/SOPA pipeline.

Use the generated SOPA sample sheet and generated SOPA params file.

Run SOPA with Singularity and SLURM.

Collect the resulting adata.h5ad paths into a manifest.

Output manifest:

sopa_adata_manifest.csv

The manifest contains:

sample_index
sample_id
sample_dir_name
adata_path
sopa_explorer_dir

The downstream pipeline uses this manifest to pass one sample at a time into downstream modules.

7.5 Inner SOPA configuration

Location:

Spatial-mIFlow/programs_hierarchy/sopa_setup/configs/sopa.config

Purpose:

Configure inner SOPA resource allocation.

Separate CPU and GPU queues.

Enable Singularity GPU support with --nv.

Pass cache variables into containers.

Control Cellpose parallelism and retries.

Important features:

withName: '.*CELLPOSE.*' {
    queue = GPU_QUEUE
    cpus = 4
    memory = 8.GB
    time = '24h'
    maxForks = (params.cellpose_max_forks ?: 4)
    errorStrategy = 'retry'
    maxRetries = 3
}

This makes SOPA/Cellpose execution more predictable on HPC.

7.6 RAPIDS_CLUSTERING_NOTEBOOK

Purpose:

Load SOPA AnnData.

Apply DAPI-based quality control.

Keep the DAPI-QC-filtered cell population.

Apply marker transformations such as arcsinh/log1p.

Compute PCA/neighbors/UMAP/Leiden clusters.

Save clustered AnnData.

Main output:

outputs/03_RAPIDS_clustering/<sample_dir>/
├── clustering/rapids_clustered.h5ad
├── dapi_qc/
└── report/

Conceptually, rapids_clustered.h5ad is the clean downstream analysis object after DAPI QC. Later notebooks should use this object when they need the quality-controlled cell population.

7.7 EXPORT_SOPA_CELL_BOUNDARIES_NPZ

Purpose:

Read SOPA SpatialData .zarr output.

Extract cell-boundary polygons.

Save a compact .npz file for later visualization.

Main output:

outputs/03_SEGMENTATION_export/<sample_dir>/segmentation/cell_boundaries_for_napari.npz

The expected .npz keys are:

polygons
cell_ids
shape_key
bounds

This file supports downstream visualization of real SOPA cell polygons without reopening the full SpatialData object.

7.8 SCIMAP_AUTOGATING_NOTEBOOK

Purpose:

Use the RAPIDS clustered AnnData as input.

Estimate marker gates with SCIMAP-style automatic gating.

Add positive/negative marker information and population combinations.

Save SCIMAP-gated AnnData.

Generate marker summaries and population plots.

Main output:

outputs/04_SCIMAP_autogating/<sample_dir>/
├── clustering/scimap_gated.h5ad
├── plots/
├── tables/
└── report/

This object becomes the input for contingency and spatial analysis modules.

7.9 CLUSTER_EXPLORATION_NOTEBOOK

Purpose:

Explore RAPIDS/Scanpy clusters.

Plot available clustering choices.

Visualize selected normalizations, resolutions, and UMAP min-distance values.

Summarize marker expression by cluster.

Main inputs:

rapids_clustered.h5ad
params.notebooks.cluster_exploration.*

Output:

outputs/05_CLUSTER_exploration/<sample_dir>/

The module uses Papermill to inject sample-specific paths and user-selected plotting parameters. It serializes Jupyter kernel installation with flock to avoid race conditions when multiple samples run in parallel.

7.10 CONTINGENCY_CLUSTER_POPULATIONS_NOTEBOOK

Purpose:

Compare RAPIDS Leiden clusters with SCIMAP-gated populations.

Generate raw count heatmaps for top gated populations across selected clusterings.

Output:

outputs/06_CONTINGENCY_clusters_populations/<sample_dir>/plots/contingency/

This module helps interpret how unsupervised clusters relate to automatic marker-gated cell populations.

7.11 SPATIAL_ANALYSIS_CLUSTERS_NOTEBOOK

Purpose:

Use the SCIMAP-gated AnnData or clustered AnnData with spatial coordinates.

Run spatial analysis of selected cluster/gate labels.

Generate PDF plots for spatial organization.

Typical analyses include:

spatial maps,

neighborhood relationships,

interaction matrices,

co-occurrence patterns,

centrality or graph-based summaries,

spatial autocorrelation where configured.

Output:

outputs/07_SPATIAL_cluster_analysis/<sample_dir>/

8. AnnData contracts between modules

The pipeline relies on a clear AnnData contract.

8.1 SOPA AnnData

Produced by SOPA:

outputs/02_sopa/<sample>.explorer/adata.h5ad

Expected content:

one row per segmented cell/object,

marker intensity features in .X or layers,

spatial coordinates in .obsm['spatial'] or equivalent,

cell identifiers in .obs or .obs_names.

8.2 RAPIDS clustered AnnData

Produced by RAPIDS_CLUSTERING_NOTEBOOK:

outputs/03_RAPIDS_clustering/<sample_dir>/clustering/rapids_clustered.h5ad

Expected content:

DAPI-QC-filtered cells only,

marker transformations,

UMAP embeddings,

Leiden cluster columns such as:

leiden_arcsinh_res_0_5
leiden_log1p_res_0_5

The cells removed by DAPI QC are not considered part of the clean clustered analysis population.

8.3 SCIMAP-gated AnnData

Produced by SCIMAP_AUTOGATING_NOTEBOOK:

outputs/04_SCIMAP_autogating/<sample_dir>/clustering/scimap_gated.h5ad

Expected content:

DAPI-QC-filtered cells,

RAPIDS cluster columns preserved,

SCIMAP-style marker positivity and population labels,

spatial coordinates preserved.

9. Output directory structure

Typical output tree:

outputs/
├── 00_runtime/
│   └── sopa_inputs/
│       ├── sopa_samplesheet.csv
│       ├── resolved_sopa_channels.json
│       └── sopa_params.generated.yml
├── 01_prepared_images/
├── 02_sopa/
├── 03_RAPIDS_clustering/
├── 03_SEGMENTATION_export/
├── 04_SCIMAP_autogating/
├── 05_CLUSTER_exploration/
├── 06_CONTINGENCY_clusters_populations/
└── 07_SPATIAL_cluster_analysis/

10. Runtime and cache files

Do not commit or treat these as pipeline code:

.nextflow/
.nextflow.log*
logs/
outputs/
archive/
**/__pycache__/
*.pyc
Spatial-mIFlow/programs_hierarchy/sopa_setup/null/
Spatial-mIFlow/programs_hierarchy/sopa_setup/singularity/

Large Singularity images and old SOPA reports should be stored in scratch/cache directories, not inside reusable pipeline code.

Recommended .gitignore entries:

.nextflow/
.nextflow.log*
logs/
outputs/
archive/
**/__pycache__/
*.pyc
Spatial-mIFlow/programs_hierarchy/sopa_setup/null/
Spatial-mIFlow/programs_hierarchy/sopa_setup/singularity/

11. How -resume behaves

Nextflow resumes based on task hashes. If a user changes a parameter that appears in a process command, the affected process should rerun.

Examples:

Changing notebooks.cluster_exploration.selected_resolutions reruns CLUSTER_EXPLORATION_NOTEBOOK.

Changing SCIMAP parameters reruns SCIMAP and downstream contingency/spatial modules.

Changing image-preparation parameters reruns image preparation and everything downstream.

Changing SOPA segmentation parameters reruns SOPA and everything downstream, unless skip_if_complete causes existing SOPA outputs to be reused.

For example:

notebooks:
  cluster_exploration:
    selected_resolutions: "0.5,1.0"
    selected_umap_min_dists: "0.3"

This asks the cluster-exploration notebook to use resolutions 0.5 and 1.0 with UMAP min-dist 0.3 only. These cluster columns and embeddings must already exist in the RAPIDS AnnData.

12. Manual thresholding annotation side analysis

A separate manual-thresholding notebook package can be used beside this pipeline. It is not part of the automated pipeline because manual gating is user-driven, but it is designed to consume pipeline outputs.

Expected inputs:

outputs/03_RAPIDS_clustering/<sample_dir>/clustering/rapids_clustered.h5ad
outputs/03_SEGMENTATION_export/<sample_dir>/segmentation/cell_boundaries_for_napari.npz

Recommended standalone notebooks:

manual_thresholding_annotation/
├── 01_manual_gate_definition.ipynb
├── 02_manual_gate_statistics_spatial_visualization.ipynb
└── 03_manual_gate_spatial_analysis.ipynb

Manual thresholding should use only the DAPI-QC-filtered RAPIDS AnnData. In segmentation visualizations, SOPA polygons that do not match RAPIDS/DAPI-QC cells can be shown as white or light grey context.

13. Article-ready software architecture summary

The workflow was refactored into a modular Nextflow architecture named Spatial-mIFlow, where each computational step is encapsulated as an independent module with its own script, notebook, and software environment. User-editable parameters are centralized in a YAML configuration file, while deployment-specific Nextflow engine settings remain in a small root nextflow.config. The pipeline code is isolated from user data and generated outputs, improving portability across HPC systems.

A key technical improvement is the streaming OME-TIFF preparation module. This module avoids loading full multiplex images into memory, computes a composite MIP channel, writes a tiled and compressed SOPA-compatible OME-TIFF, validates the written image, and records channel-level provenance. SOPA runtime inputs, including the sample sheet and generated SOPA parameters, are created dynamically under outputs/00_runtime/, avoiding hard-coded paths. Downstream single-cell analyses are implemented as Papermill-executed notebooks, allowing sample-specific paths and analysis parameters to be injected while preserving executed notebooks for provenance.

This design improves reproducibility, maintainability, and scalability by separating user configuration from workflow code, isolating each computational step as a module, validating intermediate image products, and organizing outputs by analysis stage.

14. Minimal checklist before running

Before launching a fresh run:

# 1. Check sample sheet
cat master_samplesheet.csv

# 2. Check input images exist
ls -lh data/

# 3. Check user parameters
cat params.yml

# 4. Check Nextflow is available
nextflow -version

# 5. Submit
sbatch run_pipeline_slurm.sh

After the run starts:

# Check runtime SOPA inputs
ls -lh outputs/00_runtime/sopa_inputs/
cat outputs/00_runtime/sopa_inputs/sopa_params.generated.yml

# Check prepared images
ls -lh outputs/01_prepared_images/

# Check stage outputs
ls -lh outputs/

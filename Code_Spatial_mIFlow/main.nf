nextflow.enable.dsl = 2


// =============================================================================
// Scratch defaults
// =============================================================================
//
// workDir isn't readable from nextflow.config (config parses before workDir
// is resolved), so these default to sitting alongside it here instead. Only
// fills in values not already set via params.yml/-params-file.
//
// Must run before the `include` statements below: DSL2 modules capture
// `params` at include-time, so a module can only see values assigned to
// params.* earlier in this file, not later.

params.sopa_work = params.sopa_work ?: "${workDir}/sopa_work"
params.tmp_dir = params.tmp_dir ?: "${workDir}/tmp"
params.matplotlib_cache_dir = params.matplotlib_cache_dir ?: "${workDir}/matplotlib"
params.cellpose_model_cache = params.cellpose_model_cache ?: "${workDir}/cellpose_models"
params.sopa_logs_dir = params.sopa_logs_dir ?: "${workDir}/logs/sopa"


// =============================================================================
// MODULES
// =============================================================================

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


// =============================================================================
// WORKFLOW
// =============================================================================

workflow {

    Channel
        .fromPath(params.samplesheet)
        .ifEmpty {
            error "Samplesheet not found: ${params.samplesheet}"
        }
        .splitCsv(header: true)
        .map { row ->
            def sample_id = row.sample.toString()
            def data_path = row.data_path.toString().trim()

            // Relative data_path entries are resolved against launchDir (where
            // you invoke from), not projectDir (Spatial-mIFlow/, the pipeline code
            // itself) -- data/ lives at the project root alongside the
            // samplesheet, not inside Spatial-mIFlow/.
            def resolved_input = data_path.startsWith("/")
                ? file(data_path, checkIfExists: true)
                : file("${launchDir}/${data_path}", checkIfExists: true)

            tuple(
                sample_id,
                resolved_input
            )
        }
        .set { raw_samples_ch }


    PREPARE_COMPOSITE_MIP(raw_samples_ch)


    PUBLISH_PREPARED_IMAGE(
        PREPARE_COMPOSITE_MIP.out.prepared
    )


    BUILD_SOPA_SAMPLESHEET(
        PUBLISH_PREPARED_IMAGE.out.manifest.collect()
    )


    RUN_SOPA_AND_COLLECT_ADATA(
        BUILD_SOPA_SAMPLESHEET.out.sopa_samplesheet,
        BUILD_SOPA_SAMPLESHEET.out.resolved_channels,
        BUILD_SOPA_SAMPLESHEET.out.sopa_params
    )


    RUN_SOPA_AND_COLLECT_ADATA.out.adata_manifest
        .splitCsv(header: true)
        .map { row ->
            tuple(
                row.sample_index.toInteger(),
                row.sample_id.toString(),
                row.sample_dir_name.toString(),
                file(row.adata_path.toString())
            )
        }
        .set { sopa_adata_ch }

    RAPIDS_CLUSTERING_NOTEBOOK(sopa_adata_ch)


    EXPORT_SOPA_CELL_BOUNDARIES_NPZ(sopa_adata_ch)


    CLUSTER_EXPLORATION_NOTEBOOK(
        RAPIDS_CLUSTERING_NOTEBOOK.out.rapids_adata
    )


    SCIMAP_AUTOGATING_NOTEBOOK(
        RAPIDS_CLUSTERING_NOTEBOOK.out.rapids_adata
    )


    CONTINGENCY_CLUSTER_POPULATIONS_NOTEBOOK(
        SCIMAP_AUTOGATING_NOTEBOOK.out.scimap_adata
    )

    SPATIAL_ANALYSIS_CLUSTERS_NOTEBOOK(
        SCIMAP_AUTOGATING_NOTEBOOK.out.scimap_adata
    )

}

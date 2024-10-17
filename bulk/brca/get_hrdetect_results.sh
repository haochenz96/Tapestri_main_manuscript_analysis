# Get all the tsv files under the specified directory and create symlinks in the current directory
for tsv_file in /juno/work/iacobuzc/haochen/bulk_analysis/tempo_pipeline_cmo/results/somatic/*/hrdetect*/*.tsv; do
    ln -s "$tsv_file" "$(basename "$tsv_file")"
done

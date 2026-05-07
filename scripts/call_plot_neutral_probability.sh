

rm result_summaries/forward_validation_*_summary.Rds results/forward_validation_*.Rds

#rm results/group_fit_df.Rds core_data/simulation_intervals.Rds core_data/karyotypes.Rds  core_data/grouped_intervals.Rds

#Rscript scripts/run_forward_validation_all.R     --grouped_intervals_path=core_data/grouped_intervals.Rds     --group_fit_path=results/group_fit_df.Rds     --output_dir=results     --summary_output_dir=result_summaries    --n_reps=30 --n_null_pairs=100  --n_cores=4  --bottleneck_size=1000  --expansion_factor=32 --distance_metric=chrom_weighted_wasserstein

Rscript scripts/run_forward_validation_all.R --grouped_intervals_path=core_data/grouped_intervals.Rds --group_fit_path=results/group_fit_df.Rds --output_dir=results --summary_output_dir=result_summaries --n_reps=200 --n_null_pairs=1000 --n_cores=7 --bottleneck_size=3000 --expansion_factor=32 --distance_metric=wasserstein

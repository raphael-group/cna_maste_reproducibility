# Preprocessing: CalicoST vs cnamaste

**TL;DR:** spot removal (steps 1-4) runs before gene and SNP removal in both tools. Steps 1-10 match except the
spot floor (3-4); 12 of 21 steps differ, 5 of them (17-21) after the BAF-only clones: normal-spot selection and the
filters it drives.

Compared: CalicoST `c1abcae` (`arg_parse.py` defaults with the shipped `configuration_cna` /
`configuration_cna_multi`; "fixed" marks a value hard-coded, not configurable) and cnamaste `126ff79`. A cnamaste value
reads `zenodo_sim_config.yaml` / `config.yaml` where the two differ. **Same?** is true where the mechanism and the value
agree under both cnamaste configurations.

| Step | Parameter | CalicoST | cnamaste | Same? |
|---|---|---|---|---|
| 1. Keep in-tissue spots in the count matrix | `in_tissue` | `in_tissue == True`, intersected with count barcodes | same | true |
| 2. Keep spots with SNP counts | `barcodes.txt` | count barcodes intersected with SNP barcodes | same | true |
| 3. Drop spots by total UMI | `min_snpumi_perspot` / `quality.spot_min_snp_umis` | `>= 50` multi-slice, `> 50` single-slice | `>= 50` / `>= 10` | false |
| 4. Drop spots by SNP-covering UMI | as step 3 | `>= 50` | `>= 50` / `>= 10` | false |
| 5. Drop genes expressed in few spots | `min_percent_expressed_spots` | count > 0 in >= 0.005 x kept spots | same | true |
| 6. Drop listed genes | `filtergenelist_file` | `ig_gene_list.txt` | same | true |
| 7. Drop SNPs in listed regions | `filterregion_file` | `HLA_regions.bed`, Start <= pos < End | same | true |
| 8. Zero outlier genes by total UMI | `quality.local_outlier_filter` | always, `LocalOutlierFactor(n_neighbors=200)` | `True`, same estimator | true |
| 9. Thin genes above a UMI cap | `quality.max_umi_per_gene` | n/a | n/a (`None`) | true |
| 10. Drop genes off chr1-22 or absent from the gene table | `hgtable_file` | chr1-22 genes present in the counts | same | true |
| 11. Drop SNPs with no overlapping gene | rows searched back (fixed) | 50 | 100 | false |
| 12. Merge genes into blocks to a SNP-UMI floor | `initial_min_umi` / `quality.phasing_min_snp_umis` | 15 (fixed) | 50 / 150 | false |
| 13. Merge blocks into bins: SNP-UMI floor | `secondary_min_umi` / `quality.secondary_min_snp_umi` | 300 | 200 / 300 | false |
| 14. Merge blocks into bins: total-UMI floor | `quality.secondary_min_umi` | n/a | n/a (0) | true |
| 15. Merge blocks into bins: length cap | `max_binlength` | 5e6 (fixed); ends a bin below its floor | 5e6; applied only once every floor is met | false |
| 16. Merge a segment's last bin into the previous | rule (fixed) | SNP UMI < 0.5 x floor and length < 0.5 x cap | any floor unmet | false |
| 17. Select normal spots | `EPS_BAF` 0.05, `PERCENT_NORMAL` 40 (fixed) | near-normal BAF clone, lowest-spread 40%, +10% until normal UMI > 200 x bins; spread on 7-spot pooled RDR | same start; stop rule never fires, ends at 100%; spread on unpooled RDR | false |
| 18. Drop bins with unbalanced normal-spot BAF | `quality.normal_allele_specific_confidence`; tau floor 30 (fixed) | (0.05, 0.95) (fixed) | (0.01, 0.99) | false |
| 19. Zero DE genes, normal vs tumor/unsure | `quality.filter_normal_diffexp`; logFC 2 / 4, UMI 80th percentile (fixed) | always | n/a (`False`) | false |
| 20. Re-merge bins to a normal-spot UMI floor | `quality.secondary_min_normal_umi` | n/a | 100, with the floors of 13-15, broken at contigs only | false |
| 21. Zero read depth of bins with low normal count | `MIN_NORMAL_COUNT_PERBIN` / `quality.min_normal_count_perbin` | 20 (fixed) | 100 | false |

> [!NOTE]
> Order: in both tools the floor of step 5 is a share of the spots kept by 1-4, and the per-spot SNP UMI of step 4
> counts every SNP, before steps 7 and 11 drop any. Phasing (between 12 and 13) and the BAF-only clones (between 16
> and 17) drop no spot, gene or SNP and are not listed.

> [!WARNING]
> Step 17, cnamaste: the stop test `stdthreshold > 1.5 * prior_stdthreshold` compares against `prior_stdthreshold =
> np.inf`, never updated, and the 200 x bins test sits under `if False`. The normal spots are the near-normal clone
> less its highest-spread spot(s).

## Sources

CalicoST paths are under
[`src/calicost` @ `c1abcae`](https://github.com/raphael-group/CalicoST/tree/c1abcae/src/calicost); `utils_IO`
lines are `load_joint_data`'s, `load_data`'s (single-slice) in brackets.

| Step | CalicoST | cnamaste |
|---|---|---|
| 1-2 | [utils_IO.py#L173](https://github.com/raphael-group/CalicoST/blob/c1abcae/src/calicost/utils_IO.py#L173), [#L194](https://github.com/raphael-group/CalicoST/blob/c1abcae/src/calicost/utils_IO.py#L194) ([#L46](https://github.com/raphael-group/CalicoST/blob/c1abcae/src/calicost/utils_IO.py#L46), [#L57](https://github.com/raphael-group/CalicoST/blob/c1abcae/src/calicost/utils_IO.py#L57)) | [io.py#L277](https://github.com/michaelJwilson/cna_maste/blob/42b96ff/python/cnamaste/io.py#L277), [#L504](https://github.com/michaelJwilson/cna_maste/blob/42b96ff/python/cnamaste/io.py#L504), [#L543](https://github.com/michaelJwilson/cna_maste/blob/42b96ff/python/cnamaste/io.py#L543) |
| 3-4 | [utils_IO.py#L231-L243](https://github.com/raphael-group/CalicoST/blob/c1abcae/src/calicost/utils_IO.py#L231-L243) ([#L66-L75](https://github.com/raphael-group/CalicoST/blob/c1abcae/src/calicost/utils_IO.py#L66-L75)) | [io.py#L580-L592](https://github.com/michaelJwilson/cna_maste/blob/42b96ff/python/cnamaste/io.py#L580-L592) |
| 5-8 | [utils_IO.py#L247-L283](https://github.com/raphael-group/CalicoST/blob/c1abcae/src/calicost/utils_IO.py#L247-L283) ([#L78-L115](https://github.com/raphael-group/CalicoST/blob/c1abcae/src/calicost/utils_IO.py#L78-L115)) | [io.py#L619-L752](https://github.com/michaelJwilson/cna_maste/blob/42b96ff/python/cnamaste/io.py#L619-L752) |
| 9 | n/a | [io.py#L767](https://github.com/michaelJwilson/cna_maste/blob/42b96ff/python/cnamaste/io.py#L767) |
| 10-11 | [utils_IO.py#L504-L537](https://github.com/raphael-group/CalicoST/blob/c1abcae/src/calicost/utils_IO.py#L504-L537) | [reference.py#L45](https://github.com/michaelJwilson/cna_maste/blob/42b96ff/python/cnamaste/reference.py#L45), [omics.py#L118](https://github.com/michaelJwilson/cna_maste/blob/42b96ff/python/cnamaste/omics.py#L118), [#L143](https://github.com/michaelJwilson/cna_maste/blob/42b96ff/python/cnamaste/omics.py#L143) |
| 12 | [utils_IO.py#L540](https://github.com/raphael-group/CalicoST/blob/c1abcae/src/calicost/utils_IO.py#L540), [parse_input.py#L94](https://github.com/raphael-group/CalicoST/blob/c1abcae/src/calicost/parse_input.py#L94) | [omics.py#L476](https://github.com/michaelJwilson/cna_maste/blob/42b96ff/python/cnamaste/omics.py#L476), [run_cnamaste.py#L164](https://github.com/michaelJwilson/cna_maste/blob/42b96ff/python/cnamaste/scripts/run_cnamaste.py#L164) |
| 13-16 | [utils_IO.py#L827-L866](https://github.com/raphael-group/CalicoST/blob/c1abcae/src/calicost/utils_IO.py#L827-L866), [parse_input.py#L110](https://github.com/raphael-group/CalicoST/blob/c1abcae/src/calicost/parse_input.py#L110) | [omics.py#L16-L95](https://github.com/michaelJwilson/cna_maste/blob/42b96ff/python/cnamaste/omics.py#L16-L95), [run_cnamaste.py#L410](https://github.com/michaelJwilson/cna_maste/blob/42b96ff/python/cnamaste/scripts/run_cnamaste.py#L410) |
| 17 | [calicost_main.py#L113-L125](https://github.com/raphael-group/CalicoST/blob/c1abcae/src/calicost/calicost_main.py#L113-L125), [arg_parse.py#L42](https://github.com/raphael-group/CalicoST/blob/c1abcae/src/calicost/arg_parse.py#L42) | [normal_spot.py#L22-L91](https://github.com/michaelJwilson/cna_maste/blob/42b96ff/python/cnamaste/normal_spot.py#L22-L91), [run_cnamaste.py#L789](https://github.com/michaelJwilson/cna_maste/blob/42b96ff/python/cnamaste/scripts/run_cnamaste.py#L789) |
| 18 | [utils_IO.py#L970](https://github.com/raphael-group/CalicoST/blob/c1abcae/src/calicost/utils_IO.py#L970), [calicost_main.py#L139](https://github.com/raphael-group/CalicoST/blob/c1abcae/src/calicost/calicost_main.py#L139) | [normal_spot.py#L372](https://github.com/michaelJwilson/cna_maste/blob/42b96ff/python/cnamaste/normal_spot.py#L372), [run_cnamaste.py#L827](https://github.com/michaelJwilson/cna_maste/blob/42b96ff/python/cnamaste/scripts/run_cnamaste.py#L827) |
| 19 | [utils_IO.py#L1078](https://github.com/raphael-group/CalicoST/blob/c1abcae/src/calicost/utils_IO.py#L1078), [calicost_main.py#L145](https://github.com/raphael-group/CalicoST/blob/c1abcae/src/calicost/calicost_main.py#L145) | [normal_spot.py#L231](https://github.com/michaelJwilson/cna_maste/blob/42b96ff/python/cnamaste/normal_spot.py#L231), [run_cnamaste.py#L856](https://github.com/michaelJwilson/cna_maste/blob/42b96ff/python/cnamaste/scripts/run_cnamaste.py#L856) |
| 20 | n/a | [run_cnamaste.py#L877](https://github.com/michaelJwilson/cna_maste/blob/42b96ff/python/cnamaste/scripts/run_cnamaste.py#L877) |
| 21 | [calicost_main.py#L146-L151](https://github.com/raphael-group/CalicoST/blob/c1abcae/src/calicost/calicost_main.py#L146-L151) | [normal_spot.py#L171-L209](https://github.com/michaelJwilson/cna_maste/blob/42b96ff/python/cnamaste/normal_spot.py#L171-L209), [run_cnamaste.py#L930](https://github.com/michaelJwilson/cna_maste/blob/42b96ff/python/cnamaste/scripts/run_cnamaste.py#L930) |

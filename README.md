# mand v29a source patch

This is the first post-v28 development patch. Numeric image-processing conditions remain frozen, while tissue multicoat QC intentionally uses the supplied `src.r` lineage through an isolated adapter.

## Important

- Do not modify the frozen v28 release.
- Keep the existing dependency files in `mand_patch/R`; this patch is layered on top of them.
- The generated multicoat PNG filenames now use the actual `case_id`; this is an intentional v29a path change.
- Numeric NIfTI outputs should match v28. QC PNG hashes are not expected to match.
- Copy `R/*.R` into `mand_patch/R`, rebuild/load the package, then run `scripts/run_v29a_pilot6.R`.

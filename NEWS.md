# GRETTA 3.4.1
- Fix for situations when dep scores are identical across all samples. ie. there is no variation. This issue occurs in `GI_screen()`. 

# GRETTA 3.4.0
- Fixes issues in `GI_screen_perms()` that caused control and mutant groups to be mismatched due to in efficient fitering.

# GRETTA 3.3.0
- Fixes issue #44 compatibility with R-4.4.

# GRETTA 3.2.0
- In `GI_screen_perms()` NAs in dep probs have been filtered early on improve speed.

# GRETTA 3.1.0
- Changing random sampling method in `GI_screen_perms()` to ensure sample size of groups are the same. Otherwise, errors are common.

# GRETTA 3.0.0
GRETTA now supports RNAi screns! 
- Formatted RNAi data can be downloaded from https://www.bcgsc.ca/downloads/ytakemon/GRETTA/RNAi/
- A script used to format is availale in [DepMap_data_versions/RNAi/](https://github.com/ytakemon/GRETTA/tree/main/DepMap_data_versions/RNAi).
- Turn on the option to perform RNAi screens using the following optioons `select_cell_lines(rnai_screen = TRUE)`, `GI_screen(rnai_screen = TRUE)`, and `GI_screen_perms(rnai_screen = TRUE)`. Details can be found in `help(GI_screen)`. 
- *Note:* The DepMap Consortium no longer performs RNAi screens therefore the data available are fixed. 

# GRETTA 2.7.0
Bugs in `get_GeneNameID ` were fixed.

# GRETTA 2.6.1
Bugs in `select_cell_lines()` were fixed.

# GRETTA 2.6.0
Bugs in `select_cell_lines()` were fixed.

# GRETTA 2.5.1

# GRETTA 2.5.0
Bugs in `list_mutations()` were fixed.

# GRETTA 2.4.0
- Latest DepMap dataset 23Q4 is now available. See https://github.com/ytakemon/GRETTA/tree/main/DepMap_data_versions/ for more details.
- `list_mutations()`, `list_cancer_types()`, `extract_rna()`, `select_cell_lines()`, and `list_mutations()` have been updated to accommodate different column names used in 23Q4. 
	
# GRETTA 2.3.0

# GRETTA 2.2.0

# GRETTA 2.1.0
`plot_screen()` now has new option to plot specific genes using the `gene_list =` argument. Please read `?plot_screen` for more info. 

# GRETTA 2.0.0

Major bump. New function `GI_screen_perms()` has been added to perform permutation on GI screen results for p-value correction/

# GRETTA 1.0.6

Minor patch to for Connibear lab requests

# GRETTA 1.0.5

Minor patch to remove duplicates from common screens

# GRETTA 1.0.4

Minor patch to fix typo

# GRETTA 1.0.3

Minor patch to fix protein expression extraction.

# GRETTA 1.0.2

Minor patch to fix typo.

# GRETTA 1.0.1

Minor patch to ensure top/bottoms are selected when annotating.

# GRETTA 1.0.0

Due to slow turn around by Bioconductor, it has been decided that GRETTA will not be uploaded to their repository. 

## Introducing new functions!

- `protien_coexpress()` to perform co-expression analysis for protein data.
- `common_coefs_prot()` to map the Pearson's coefficient between input proteins.
- `common_coefs_rna()` to map the Pearson's coefficient between input RNA.
- `common_coefs_coess()` to map the Pearson's coefficient between input co/anti-essential genes.

# GRETTA 0.99.5

- `anntate_coess()` was renamed to `anntate_df()`
- `coessential_map()`, `get_inflection_point()`, and `anntate_coess()` can now handle multiple input genes.

## Introducing new functions!:

- `common_coefs()` to map the Pearson's coefficient between all input genes.
- `rna_coexpress()` to perform co-expression analysis for mRNA.

# GRETTA 0.99.4

- Fixing path for figures in README to pass bioconductor warnings

# GRETTA 0.99.3

- Untracked files unrelated to R package devel that caused errors in bioconductor checks.

# GRETTA 0.99.2

- Removed all instances of lab path to pass bioconductor pre-check. 

# GRETTA 0.99.1

- Added smaller example data to pass bioconductor pre-check. 
- Example data is now loaded using `download_example_data()`

# GRETTA 0.99.0

- Pre-release for Bioconductor submission

# GRETTA 0.6.0

- Version 0.5.0 was submitted to Bioinformatics and version 0.6.0 addresses several reviewer comments (detailed below).
- GRETTA is now GRETTA (with two T's) due conflict with existing package on CRAN.
- Upon loading GRETTA into R, a welcome message is displayed to indicate package version and the latest DepMap data set that is available. 
- Function arguments are now all lower case.
- `GI_screen()` has a new argument `gene_list = ` to allow small-scale GI screens and reduce computational time. 
- Default pan-cancer `coessential_map()` now uses pre-computed cor matrix to reduce computational time. This matrix is provided along with the DepMap data.
- A vignette is now available and in a Bioconductor format to prepare for submission (same as readme). 
- Examples are now mostly run-able to comply with Bioconductor.

# GRETTA 0.5.0

- Default output file name from `GI_screen()` has been renamed.
- Singularity definition file and tutorial file has been included as a supplement.

# GRETTA 0.4.1

-   Users can now custom name their output files in `GI_screen()` and `coessential_map()` with the `output_filename =` option. 
-   Technically I think this should be a minor version update, but it seemed too insignificant to go to v0.5.0. 

# GRETTA 0.4.0

-   All GRETTA functions are now compatible with multiple DepMap data versions (20Q1, 21Q4, and 22Q2)! 
-   An FAQ section is now available [here](https://github.com/ytakemon/GRETTA/wiki/Frequently-Asked-Questions).

# GRETTA 0.3.2

-   GINIR package name has now changed to GRETTA (Genetic inteRaction and EssenTially mApper).
-   `GINI_screen()` is now `GI_screen()`
-   Added a `NEWS.md` file to track changes to the package.
-   Back filled log for previous releases

# GINIR 0.3.1

-   The patch fixes a bug that occurred when NAs were found in dependency probability columns (#33)

# GINIR 0.3.0

-   This minor feature improvement now shows cell line IDs while listing all available mutations using `list_available_mutations()`.

# GINIR 0.2.0

-   First stable release of GINIR
-   Contains citable DOI generated by Zenodo
-   Users can now refer to the README page for a full tutorial of the tool and most up to date news.

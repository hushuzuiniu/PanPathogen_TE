# Input data

This directory separates small versioned reference files from large analysis
inputs. Do not commit patient-identifiable or controlled-access data.

## Versioned reference files

### `reference/sample_info_7_datasets_with_timecourse.csv`

Required columns:

| Column | Description |
| --- | --- |
| `Dataset` | Dataset identifier used in expression filenames |
| `Sample_Name` | Original count-matrix column name |
| `Infection_State` | Pre/post infection label |
| `Time_Point` | Sampling time, such as `0d` or `6h` |

### `reference/ruTE_other_intergenic_n595.txt`

One ruTE locus identifier per line. The included list contains 595 loci.

## Large files not tracked by Git

Place these files in `data/reference/`:

| File | Expected format | Size in the source analysis |
| --- | --- | ---: |
| `Genes_0based.bed` | Five tab-separated columns: chromosome, start, end, GTF-style attributes, strand | 10,500,914 bytes |
| `TEs_all_other_intergenic_0based.bed` | At least six tab-separated columns; TE attributes in column 4 and strand in column 6 for the source file | 302,740,397 bytes |

Place raw count matrices in `data/raw_expression/`. Each of the following
dataset prefixes requires one `_readscounts_matrix_Gene.txt` file and one
`_readscounts_matrix_TE_Loci.txt` file:

- `Virus_EBV_2`
- `Virus_EBV_5`
- `Fungi_Af_4`
- `Fungi_Af_6`
- `Virus_IAV_11`
- `Virus_RSV_1`
- `Virus_HCV_1`

Gene matrices contain an ID column followed by sample count columns. TE
matrices contain six annotation columns (`Geneid`, `Chr`, `Start`, `End`,
`Strand`, `Length`) followed by sample count columns.

## Integrity checks

SHA-256 checksums for the files used in the source analysis are recorded in
[`input_manifest.tsv`](input_manifest.tsv). On macOS or Linux, verify a file
with `shasum -a 256 <file>`.

Step 1 writes standardized metadata to `data/processed/` and normalized
expression matrices to `data/normalized_expression/`. Both directories are
generated and excluded from Git.

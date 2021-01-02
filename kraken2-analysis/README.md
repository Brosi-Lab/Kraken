
#### Directory Layout

```
kraken2-analysis
├── commands.sh
├── public-samples
│   ├── ${SPECIES}-${TAXID}
│   │   ├── ${BIOPROJECT}-downloaded.txt
│   │   ├── ${BIOPROJECT}-random.txt
│   │   └── ${BIOPROJECT}-runs.txt
│   └── ${SPECIES}-${TAXID}
│       └── ${ASSEMBLY_ACCESSION}.txt
└── results
    ├── fastq-stats
    │   └── ${SAMPLE}
    │       ├── ${SAMPLE}_R(1|2)-(final|original)_fastqc.html
    │       └── ${SAMPLE}_R(1|2)-(final|original).json
    └── reports
        └── (paired|single-end)
            └── ${SAMPLE}-report.txt

```

#### Base Directory

| Filename    | Description                                           |
|-------------|-------------------------------------------------------|
| commands.sh | The commands used for the Kraken2 setup and analysis. |
| kraken2-dbinfo.txt | A summary of the Kraken2 database base used in this study (`kraken2-inspect`). |

#### `public-samples` Directory
The `public-samples` directory includes information of additional sequences that were added to the default `plant` database for Kraken2.


| Filename                     | Description                                              |
|------------------------------|----------------------------------------------------------|
| ${ASSEMBLY_ACCESSION}.txt    | The NCBI Assembly accession that was downloaded.         |
| ${BIOPROJECT}-downloaded.txt | The SRA Runs which FASTQs were actually downloaded from. |
| ${BIOPROJECT}-random.txt     | A random subset of SRA Runs to download.                 |
| ${BIOPROJECT}-runs.txt       | The full list of SRA Runs associated with a BioProject.  |

#### `results` Directory
The `results` directory contains the summary statistics of the input FASTQs and the Kraken2 report files.

| Extension    | Description                                                                            |
|--------------|----------------------------------------------------------------------------------------|
| _fastqc.html | FASTQC HTML report of the raw FASTQs (original) or the quality filtered FASTQs (final) |
| .json        | Summary statistics of the raw FASTQs (original) or the quality filtered FASTQs (final) |
| -report.txt  | A Kraken2 report for a given sample                                                    |

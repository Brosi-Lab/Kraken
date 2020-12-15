# Create Conda environment
conda create -n plant-analysis -c conda-forge -c bioconda bactopia bbmap illumina-cleanup kraken2 ncbi-genome-download sra-tools
conda activate plant-analysis

# Download assemblies, add taxid to header
# GCA_008642245.1 (Bassia scoparia, taxid: 83154)
ncbi-genome-download  plant -s genbank  -F fasta -A GCA_008642245.1 
mv genbank/plant/GCA_008642245.1/GCA_008642245.1_CSU_KoSco_1.0_genomic.fna.gz  ./
gunzip GCA_008642245.1_CSU_KoSco_1.0_genomic.fna.gz 
mv GCA_008642245.1_CSU_KoSco_1.0_genomic.fna GCA_008642245.1.fna
cat GCA_008642245.1.fna | sed "s/^>/>kraken:taxid|83154 /" > public-samples/bassia-scoparia-83154/GCA_008642245.fasta

# GCA_011037805.1 (Carya illinoinensis taxid: 32201)
ncbi-genome-download  plant -s genbank  -F fasta --verbose -A GCA_011037805.1
mv genbank/plant/GCA_011037805.1/GCA_011037805.1_ASM1103780v1_genomic.fna.gz ./
gunzip GCA_011037805.1.fna.gz 
cat GCA_011037805.1.fna |  sed "s/^>/>kraken:taxid|32201 /" > public-samples/carya-illinoinensis-32201/GCA_011037805.fasta
rm -rf genbank/

# Prepare FASTQS
# Randomly select 10 Runs per project and download fastq
cd public samples
mkdir fastqs
find | grep runs | sed 's/-runs.txt//' | xargs -I {} sh -c 'shuf -n 10 {}-runs.txt > {}-random.txt'
find | grep random | xargs -I {} cat {} | grep SRR | xargs -I {} fasterq-dump {} --split-files -p -O fastqs/

# Merge FASTQs
cat ./ambrosia-artemisiifolia-4212/PRJNA449949-random.txt | \
    xargs -I {} cat fastqs/{}_1.fastq > ambrosia-artemisiifolia-4212/PRJNA449949_R1.fastq.gz
cat ./ambrosia-artemisiifolia-4212/PRJNA449949-random.txt | \
    xargs -I {} cat fastqs/{}_2.fastq > ambrosia-artemisiifolia-4212/PRJNA449949_R2.fastq.gz

cat ./broussonetia-papyrifera-172644/PRJNA437223-random.txt | \
    xargs -I {} cat fastqs/{}_1.fastq > ./broussonetia-papyrifera-172644/PRJNA437223_R1.fastq.gz
cat ./broussonetia-papyrifera-172644/PRJNA437223-random.txt | \
    xargs -I {} cat fastqs/{}_2.fastq > ./broussonetia-papyrifera-172644/PRJNA437223_R2.fastq.gz

cat ./populus-tremuloides-3693/PRJNA299390-random.txt | \
    xargs -I {} cat fastqs/{}_1.fastq > ./populus-tremuloides-3693/PRJNA299390_R1.fastq.gz
cat ./populus-tremuloides-3693/PRJNA299390-random.txt | \
    xargs -I {} cat fastqs/{}_2.fastq > ./populus-tremuloides-3693/PRJNA299390_R2.fastq.gz

cat ./artemisia-tridentata-55611/PRJNA258169-random.txt | \
    xargs -I {} cat fastqs/{}_1.fastq > ./artemisia-tridentata-55611/PRJNA258169_R1.fastq.gz
cat ./artemisia-tridentata-55611/PRJNA258169-random.txt | \
    xargs -I {} cat fastqs/{}_2.fastq > ./artemisia-tridentata-55611/PRJNA258169_R2.fastq.gz

cat ./poa-pratensis-1753/PRJNA517968-random.txt | \
    xargs -I {} cat fastqs/{}_1.fastq > ./poa-pratensis-1753/PRJNA517968_R1.fastq.gz
cat ./poa-pratensis-1753/PRJNA517968-random.txt | \
    xargs -I {} cat fastqs/{}_2.fastq > ./poa-pratensis-1753/PRJNA517968_R2.fastq.gz

cat ./populus-deltoides-3696/SAMN08381122-random.txt | \
    xargs -I {} cat fastqs/{}_1.fastq > ./populus-deltoides-3696/SAMN08381122_R1.fastq.gz
cat ./populus-deltoides-3696/SAMN08381122-random.txt | \
    xargs -I {} cat fastqs/{}_2.fastq > ./populus-deltoides-3696/SAMN08381122_R2.fastq.gz

# Convert FASTQ to FASTQ, add taxid to header
reformat.sh in1=./ambrosia-artemisiifolia-4212/PRJNA449949_R1.fastq.gz \
            in2=./ambrosia-artemisiifolia-4212/PRJNA449949_R2.fastq.gz \
            out=stdout.fasta | sed "s/^>/>kraken:taxid|4212 /" > ./ambrosia-artemisiifolia-4212/PRJNA449949.fasta

reformat.sh in1=./artemisia-tridentata-55611/PRJNA258169_R1.fastq.gz \
            in2=./artemisia-tridentata-55611/PRJNA258169_R2.fastq.gz \
            out=stdout.fasta | sed "s/^>/>kraken:taxid|55611 /" > ./artemisia-tridentata-55611/PRJNA258169.fasta

reformat.sh in1=./broussonetia-papyrifera-172644/PRJNA437223_R1.fastq.gz \
            in2=./broussonetia-papyrifera-172644/PRJNA437223_R2.fastq.gz \
            out=stdout.fasta | sed "s/^>/>kraken:taxid|172644 /" > ./broussonetia-papyrifera-172644/PRJNA437223.fasta

reformat.sh in1=./populus-deltoides-3696/SAMN08381122_R1.fastq.gz \
            in2=./populus-deltoides-3696/SAMN08381122_R2.fastq.gz \
            out=stdout.fasta | sed "s/^>/>kraken:taxid|3696 /" > ./populus-deltoides-3696/SAMN08381122.fasta

reformat.sh in1=./poa-pratensis-1753/PRJNA517968_R1.fastq.gz \
            in2=./poa-pratensis-1753/PRJNA517968_R2.fastq.gz \
            out=stdout.fasta | sed "s/^>/>kraken:taxid|1753 /" > ./poa-pratensis-1753/PRJNA517968.fasta

reformat.sh in1=./populus-tremuloides-3693/PRJNA299390_R1.fastq.gz \
            in2=./populus-tremuloides-3693/PRJNA299390_R2.fastq.gz \
            out=stdout.fasta | sed "s/^>/>kraken:taxid|3693 /" > ./populus-tremuloides-3693/PRJNA299390.fasta


# Merge all fastas
find -name "*.fasta" | xargs -I {} cat {} > all-taxons.fasta

# Create Kraken2 Plant database, add our samples to database
kraken2-build --download-taxonomy --db kraken-db/k2-plants
kraken2-build --download-library plant --db kraken-db/k2-plants
kraken2-build --add-to-library all-taxons.fasta --db kraken-db/k2-plants/
kraken2-build --build --db kraken-db/k2-plants

# Prepare and QC Sample FASTQS
bactopia prepare data/fastqs/ | cut -f1-3 > wgs-fastqs.txt
mkdir qc-fastqs
cd qc-fastqs
illumina-cleanup --fastqs ../wgs-fastqs.txt --max_cpus 70 --cpus 8

# Run kraken2
mkdir kraken2
ls | grep -v kraken | \
    xargs -I {} -P 10 -n 1 kraken2 --db ../kraken-db/k2-plants/ \
                                   --threads 7 \
                                   --unclassified-out kraken2/{}/{}-unclassified_R#.fastq \
                                   --classified-out kraken2/{}/{}-classified_R#.fastq \
                                   --output kraken2/{}/{}-kraken2.txt \
                                   --report kraken2/{}/{}-report.txt \
                                   --use-names \
                                   --gzip-compressed \
                                   --report-zero-counts \
                                   --paired {}/{}_R1.fastq.gz {}/{}_R2.fastq.gz

mkdir ../reports
find -n "*-report.txt" | xargs -I {} cp {} ../reports

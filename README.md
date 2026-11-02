# domains2genome

domains2genome places protein domain calls on the genome, so a genome browser
can draw them over the exons that encode them.

- **In:** gene models (GFF3 or GTF) and domain calls on their proteins
  (InterProScan TSV or GFF3, or HMMER `--domtblout`)
- **Out:** BED12, one row per domain, a block per exon it spans; or with
  `--gff3`, each coding transcript with its domains as children
- JBrowse's
  [protein domains tutorial](https://jbrowse.org/jb2/docs/tutorials/protein_domains/)
  runs it with `npx`

## Install

```bash
npm install -g domains2genome
```

or run it once with `npx domains2genome`.

## From a genome to its domains

Translate the gene models, call domains on the proteins with InterProScan or
HMMER, then place the calls back on the genome.

### 1. Translate the gene models

```bash
# one protein per coding transcript, named by its transcript ID; gffread
# reads an uncompressed FASTA
gffread -y proteins.faa -g genome.fa genes.gff3

# write each stop codon as X: InterProScan refuses `*` and `.`, and X keeps
# every residue at its position
sed -i '/^>/!s/[.*]/X/g' proteins.faa
```

### 2a. Call domains with InterProScan

```bash
nextflow run ebi-pf-team/interproscan6 -r 6.0.3.0 -profile docker \
  --datadir data --input proteins.faa \
  --applications Pfam --formats gff3
```

- InterProScan writes `proteins.faa.gff3`. Its GFF3 names each domain by Pfam's
  short name (`Gp_dh_N`), which labels a gene track better than the description
  its TSV gives (`--formats tsv` reads the same otherwise)
- Leave out `--applications` to run every analysis, then pick the ones to place
  with `--source` in step 3
- Without Docker, `-profile apptainer` runs the same containers, and Apptainer
  installs without root (`tools/install-unprivileged.sh` in its repository)
- On a cluster, `-profile singularity,slurm` runs each batch as a job, as
  InterProScan's scaling docs describe
- InterProScan 5 takes
  `interproscan.sh -i proteins.faa -f tsv -appl Pfam -o proteins.faa.tsv`

### 2b. Or call Pfam domains with HMMER

```bash
# once: Pfam's models, pressed for hmmscan (about 1.5 GB)
curl -fLO https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
gunzip Pfam-A.hmm.gz
hmmpress Pfam-A.hmm

hmmscan --cut_ga --cpu 4 --noali --domtblout pfam.domtblout \
  Pfam-A.hmm proteins.faa > /dev/null
```

`--cut_ga` applies Pfam's curated threshold for each family, as InterPro does.
hmmscan keeps about two cores busy whatever `--cpu` asks, so for a whole
proteome split the proteins and scan the parts at once:

```bash
awk -v n=8 '/^>/ { i = (i + 1) % n } { print > ("chunk." i ".faa") }' proteins.faa
for f in chunk.*.faa; do
  hmmscan --cut_ga --cpu 2 --noali --domtblout "$f.domtblout" \
    Pfam-A.hmm "$f" > /dev/null &
done
wait
cat chunk.*.faa.domtblout > pfam.domtblout
```

### 3. Place the calls on the genome

```bash
domains2genome genes.gff3 proteins.faa.gff3 > domains.bed  # InterProScan
domains2genome genes.gff3 pfam.domtblout > domains.bed     # HMMER

# bgzipped and indexed, for a browser to read by region
bgzip domains.bed
tabix -p bed domains.bed.gz
```

- `--source Pfam,SMART` keeps only the analyses named, as InterProScan names
  them in its fourth column
- `--out file.bed` writes the BED to a file instead of to stdout
- `--gff3` writes GFF3 instead of BED12 (see [Output](#output))

The calls must come from proteins translated from the same gene models. Each
call names its protein by the transcript ID, the CDS `protein_id` or the CDS ID,
which covers `gffread -y`, NCBI and Ensembl protein FASTA. Where HMMER hits on
one protein overlap by more than half, the one with the lowest E-value is kept,
whatever their Pfam clans; `pfam_scan` and InterProScan resolve only hits in one
clan.

## Output

BED12 sorted by contig and start, with a `#` header line naming the columns.
After the twelve BED columns come `accession`, `description`, `source`,
`transcript`, `proteinStart`, `proteinEnd` and `evalue`. The score column is ten
times −log10(E-value), capped at 1000.

- `transcript` is the label JBrowse gives the transcript by default: its GFF3
  `Name`, else its `ID`, read with the same parser
  ([@gmod/gff-nostream](https://github.com/GMOD/gff-nostream)). Where two coding
  transcripts would write one label, each writes its ID. From a GTF it is the
  `transcript_id`, as JBrowse labels GTF transcripts.

IGV, the UCSC browser and JBrowse all draw BED12.

### GFF3

With `--gff3`, each coding transcript is an `mRNA` with its `exon` and `CDS`
lines, and each domain placed on it a `protein_match` child, so the file carries
which transcript a domain sits on. Every coding transcript is written, with or
without a domain on it.

```bash
domains2genome genes.gff3 pfam.domtblout --gff3 | bgzip > domains.gff3.gz
tabix -p gff domains.gff3.gz
```

- The `mRNA` keeps the annotation's `ID` and its `Name`, as above, with the
  CDS's `protein_id` and `product`.
- A `protein_match`'s `Name` is the domain's name, its source column the
  analysis and its score the E-value; `accession`, `description`, `proteinStart`
  and `proteinEnd` are attributes. A domain crossing an intron is written on a
  line per exon, the lines sharing one `ID`.
- Lines are sorted by contig and start, a parent before the children starting
  where it does, as tabix needs.

JBrowse's multi-row display draws each transcript of this file on a row of its
own, its exons with the domains over them, as the
[protein domains tutorial](https://jbrowse.org/jb2/docs/tutorials/protein_domains/)
shows.

## As a library

```ts
import { placeProteinFeatures } from 'domains2genome'

const { lines, unmatched } = await placeProteinFeatures(
  'genes.gff3',
  'interpro.tsv',
  { sources: ['Pfam'] },
)
```

## From interproscan2genome

domains2genome was called `interproscan2genome` before 2.0. Version 1.x read
InterProScan's GFF3 only and printed GFF3 with one span per domain, introns
included. The `interproscan2genome` command and its `--genomeGFF` and
`--interproGFF` flags still run, and now print BED12 a block per exon. Positions
differ as well: 1.x divided residue positions by three, which put each domain
about a third of the way along its protein.

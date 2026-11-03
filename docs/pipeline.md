# From a genome to its domains

Translate the gene models, call domains on the proteins with InterProScan or
HMMER, then place the calls. JBrowse's
[protein domains tutorial](https://jbrowse.org/jb2/docs/tutorials/protein_domains/)
runs the same steps on the honey bee genome.

## 1. Translate the gene models

```bash
# one protein per coding transcript, named by its transcript ID; gffread
# reads an uncompressed FASTA
gffread -y proteins.faa -g genome.fa genes.gff3

# write each stop codon as X: InterProScan refuses `*` and `.`, and X keeps
# every residue at its position
sed -i '/^>/!s/[.*]/X/g' proteins.faa
```

The calls must come from proteins translated from the same gene models, named by
the transcript ID, the CDS `protein_id` or the CDS ID, as `gffread -y`, NCBI and
Ensembl protein FASTA are.

## 2a. Call domains with InterProScan

```bash
nextflow run ebi-pf-team/interproscan6 -r 6.0.3.0 -profile docker \
  --datadir data --input proteins.faa \
  --applications Pfam --formats gff3
```

- This writes `proteins.faa.gff3`. Prefer it to `--formats tsv`: the GFF3 names
  each domain by its short name (`Gp_dh_N`), the TSV by its description, and the
  name labels a browser track better
- Leave out `--applications` to run every analysis, then choose with `--source`
- Without Docker, use `-profile apptainer`; on a cluster,
  `-profile singularity,slurm` runs each batch as a job
- InterProScan 5 takes
  `interproscan.sh -i proteins.faa -f tsv -appl Pfam -o proteins.faa.tsv`

## 2b. Or call Pfam domains with HMMER

```bash
# once: Pfam's models, pressed for hmmscan (about 1.5 GB)
curl -fLO https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
gunzip Pfam-A.hmm.gz
hmmpress Pfam-A.hmm

hmmscan --cut_ga --cpu 4 --noali --domtblout pfam.domtblout \
  Pfam-A.hmm proteins.faa > /dev/null
```

`--cut_ga` applies Pfam's curated threshold for each family, as InterPro does.
HMMER reports every family that clears it, so related families stack on one
stretch of protein; where two hits overlap by more than half, domains2genome
keeps the lower E-value. (`pfam_scan` and InterProScan only resolve hits within
one Pfam clan.) `hmmsearch --domtblout` output reads too.

hmmscan keeps about two cores busy whatever `--cpu` asks, so for a whole
proteome, split the proteins and scan the parts at once:

```bash
awk -v n=8 '/^>/ { i = (i + 1) % n } { print > ("chunk." i ".faa") }' proteins.faa
for f in chunk.*.faa; do
  hmmscan --cut_ga --cpu 2 --noali --domtblout "$f.domtblout" \
    Pfam-A.hmm "$f" > /dev/null &
done
wait
cat chunk.*.faa.domtblout > pfam.domtblout
```

## 3. Place the calls on the genome

```bash
domains2genome genes.gff3 proteins.faa.gff3 > domains.bed  # InterProScan
domains2genome genes.gff3 pfam.domtblout > domains.bed     # HMMER

bgzip domains.bed
tabix -p bed domains.bed.gz
```

See [output formats](output.md) for the columns, and for `--gff3`.

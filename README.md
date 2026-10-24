# domains2genome

domains2genome places protein domain calls on the genome, so a genome browser
can draw them over the exons that encode them.

- **In:** gene models (GFF3 or GTF) and domain calls on their proteins
  (InterProScan TSV or GFF3, or HMMER `--domtblout`)
- **Out:** BED12, one row per domain, a block per exon it spans
- `domains2genome` needs no JBrowse install. JBrowse's
  `jbrowse make-protein-features` calls this package to place the same calls

## Install

```bash
npm install -g domains2genome
```

or run it once with `npx domains2genome`.

## Usage

```bash
# Pfam rows of an InterProScan run
domains2genome genes.gff3 interpro.tsv --source Pfam > pfam.bed

# hmmscan against Pfam
domains2genome genes.gff3 pfam.domtblout > pfam.bed

# bgzipped and indexed, for a browser to read by region
domains2genome genes.gff3 interpro.tsv | bgzip > domains.bed.gz
tabix -p bed domains.bed.gz
```

- `--source Pfam,SMART` keeps only the analyses named, as InterProScan names
  them in its fourth column
- `--out file.bed` writes the BED to a file instead of to stdout

The calls must come from proteins translated from the same gene models. Each
call names its protein by the transcript ID, the CDS `protein_id` or the CDS ID,
which covers `gffread -y`, NCBI and Ensembl protein FASTA. Where HMMER hits on
one protein overlap by more than half, the one with the lowest E-value is kept,
as `pfam_scan` and InterProScan do.

## Output

BED12 sorted by contig and start, with a `#` header line naming the columns.
After the twelve BED columns come `accession`, `description`, `source`,
`transcript`, `proteinStart`, `proteinEnd` and `evalue`. The score column is ten
times −log10(E-value), capped at 1000.

IGV, the UCSC browser and JBrowse all draw BED12. In JBrowse 5, a gene track
whose `proteinFeatures` names the BED's track draws each domain on the isoform
it was called on, as the
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

This package was `interproscan2genome`. Its 1.x read InterProScan's GFF3 only
and printed GFF3 with one span per domain, introns included. The
`interproscan2genome` command and its `--genomeGFF` and `--interproGFF` flags
still run, and now print BED12 a block per exon. Positions differ as well: 1.x
divided residue positions by three, which put each domain about a third of the
way along its protein.

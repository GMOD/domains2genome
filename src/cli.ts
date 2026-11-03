#!/usr/bin/env node
import { writeFile } from 'node:fs/promises'
import { parseArgs } from 'node:util'

import { placeProteinFeatures, placementWarnings } from './place.ts'

const usage = `Usage: domains2genome <genes.gff3|gtf> <calls> [options]

Places protein features (InterProScan TSV or GFF3, or HMMER --domtblout) on
the genome through each transcript's CDS, as BED12 with a block per exon.

Options:
  --out <file>       write the BED here instead of to stdout
  --source <a,b>     keep only these analyses, as InterProScan names them
                     (e.g. Pfam,SMART)
  --gff3             write GFF3: each coding transcript as an mRNA with its
                     exons and CDS, and each feature on it as a protein_match
                     child, so a browser draws them together
  -h, --help         show this help

Examples:
  domains2genome genes.gff3 interpro.tsv --source Pfam > pfam.bed
  domains2genome genes.gff3 pfam.domtblout | bgzip > pfam.bed.gz
  tabix -p bed pfam.bed.gz

The calls name each protein by its transcript ID, its CDS protein_id, or its
CDS ID, as gffread -y, NCBI and Ensembl protein FASTA all do. Where HMMER hits
on one protein overlap by more than half, the lowest E-value is kept.`

async function main() {
  const { values, positionals } = parseArgs({
    allowPositionals: true,
    options: {
      out: { type: 'string' },
      source: { type: 'string' },
      gff3: { type: 'boolean' },
      help: { type: 'boolean', short: 'h' },
      // the 1.x spelling, kept so existing commands still run
      genomeGFF: { type: 'string' },
      interproGFF: { type: 'string' },
    },
  })
  const annotation = positionals[0] ?? values.genomeGFF
  const calls = positionals[1] ?? values.interproGFF
  if (values.help) {
    process.stdout.write(`${usage}\n`)
    process.exit(0)
  }
  if (!annotation || !calls) {
    console.error(usage)
    process.exit(1)
  }

  // an empty --source, as from an unset shell variable, keeps every analysis
  const sources = values.source
    ?.split(',')
    .map(s => s.trim())
    .filter(s => s !== '')
  const placed = await placeProteinFeatures(annotation, calls, {
    sources: sources?.length ? sources : undefined,
    gff3: values.gff3,
  })
  const text = [placed.header, ...placed.lines, ''].join('\n')
  if (values.out) {
    await writeFile(values.out, text)
  } else {
    process.stdout.write(text)
  }
  for (const warning of placementWarnings(placed)) {
    console.warn(`Warning: ${warning}`)
  }
  const skipped =
    placed.skippedSource > 0
      ? ` (${placed.skippedSource} from other analyses left out)`
      : ''
  console.error(
    `Placed ${placed.features} protein features on ${placed.transcripts} coding transcripts${skipped}`,
  )
}

main().catch((e: unknown) => {
  console.error(`Error: ${e instanceof Error ? e.message : String(e)}`)
  process.exit(1)
})

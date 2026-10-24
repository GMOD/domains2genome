#!/usr/bin/env node
import { createWriteStream } from 'node:fs'
import { parseArgs } from 'node:util'

import { placeProteinFeatures, placementWarnings } from './place.ts'
import { BED_HEADER } from './project.ts'

const usage = `Usage: domains2genome <genes.gff3|gtf> <calls> [options]

Places protein features (InterProScan TSV or GFF3, or HMMER --domtblout) on
the genome through each transcript's CDS, as BED12 with a block per exon.

Options:
  --out <file>       write the BED here instead of to stdout
  --source <a,b>     keep only these analyses, as InterProScan names them
                     (e.g. Pfam,SMART)
  -h, --help         show this help

Examples:
  domains2genome genes.gff3 interpro.tsv --source Pfam > pfam.bed
  domains2genome genes.gff3 pfam.domtblout | bgzip > pfam.bed.gz
  tabix -p bed pfam.bed.gz

The calls name each protein by its transcript ID, its CDS protein_id, or its
CDS ID, as gffread -y, NCBI and Ensembl protein FASTA all do. Where HMMER hits
on one protein overlap by more than half, the lowest E-value is kept.`

const { values, positionals } = parseArgs({
  allowPositionals: true,
  options: {
    out: { type: 'string' },
    source: { type: 'string' },
    help: { type: 'boolean', short: 'h' },
    // the 1.x spelling, kept so existing commands still run
    genomeGFF: { type: 'string' },
    interproGFF: { type: 'string' },
  },
})

const annotation = positionals[0] ?? values.genomeGFF
const calls = positionals[1] ?? values.interproGFF
if (values.help || !annotation || !calls) {
  console.error(usage)
  process.exit(values.help ? 0 : 1)
}

async function main(annotation: string, calls: string) {
  const placed = await placeProteinFeatures(annotation, calls, {
    sources: values.source?.split(',').map(s => s.trim()),
  })
  const out = values.out ? createWriteStream(values.out) : process.stdout
  out.write(`#${BED_HEADER.join('\t')}\n`)
  for (const line of placed.lines) {
    out.write(`${line}\n`)
  }
  if (values.out) {
    out.end()
  }
  for (const warning of placementWarnings(placed)) {
    console.warn(`Warning: ${warning}`)
  }
  const skipped =
    placed.skippedSource > 0
      ? ` (${placed.skippedSource} from other analyses left out)`
      : ''
  console.error(
    `Placed ${placed.lines.length} protein features on ${placed.transcripts} coding transcripts${skipped}`,
  )
}

main(annotation, calls).catch((e: unknown) => {
  console.error(`Error: ${e instanceof Error ? e.message : String(e)}`)
  process.exit(1)
})

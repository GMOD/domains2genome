import {
  GFF3_HEADER,
  proteinMatchGff3Rows,
  transcriptGff3Rows,
} from './gff3.ts'
import {
  BED_HEADER,
  bedLine,
  bestOfOverlappingHmmerHits,
  findTranscript,
  parseCall,
  projectResidues,
  readCodingTranscripts,
  readLines,
} from './project.ts'

import type { Gff3Row, PlacedCall } from './gff3.ts'
import type { ProteinCall } from './project.ts'

export interface PlacedFeatures {
  /** the line the file opens with: BED's column names, or GFF3's version */
  header: string
  /** BED12 or GFF3 rows, sorted by contig then start, without the header */
  lines: string[]
  /** protein features placed */
  features: number
  /** coding transcripts the annotation holds */
  transcripts: number
  /** protein names that match no transcript ID, protein_id or CDS ID */
  unmatched: Set<string>
  /** calls starting past the end of their transcript's CDS */
  pastCds: number
  /** calls `sources` left out */
  skippedSource: number
}

/**
 * Places the protein features `calls` holds (InterProScan TSV or GFF3, or
 * HMMER `--domtblout`) on the coding transcripts of `annotation` (GFF3 or
 * GTF, plain or gzipped), one BED12 row per feature with a block per exon.
 * `sources` keeps only the analyses it names, as InterProScan names them
 * (`Pfam`, `SMART`), ignoring case. `gff3` writes GFF3 instead: each coding
 * transcript as an mRNA with its exons and CDS, and each feature placed on
 * it as a `protein_match` child.
 */
export async function placeProteinFeatures(
  annotation: string,
  calls: string,
  { sources, gff3 = false }: { sources?: string[]; gff3?: boolean } = {},
): Promise<PlacedFeatures> {
  const { lookup, transcripts, codingTranscripts } =
    await readCodingTranscripts(annotation)
  if (transcripts === 0) {
    throw new Error(
      `${annotation} has no CDS rows to place protein features on`,
    )
  }
  const kept = sources?.map(s => s.toLowerCase())
  const isProtein = (id: string) => findTranscript(lookup, id) !== undefined
  const parsed: ProteinCall[] = []
  const skipped = new Set<string>()
  let skippedSource = 0
  for await (const line of readLines(calls)) {
    const call = parseCall(line, isProtein)
    if (call) {
      if (kept && !kept.includes(call.source.toLowerCase())) {
        skippedSource++
        skipped.add(call.source)
      } else {
        parsed.push(call)
      }
    }
  }
  if (parsed.length === 0 && skippedSource > 0) {
    throw new Error(
      `No call in ${calls} comes from ${sources?.join(', ')}; its calls come from ${[...skipped].join(', ')}`,
    )
  }

  const placedCalls: PlacedCall[] = []
  const unmatched = new Set<string>()
  let pastCds = 0
  for (const call of bestOfOverlappingHmmerHits(parsed)) {
    const transcript = findTranscript(lookup, call.protein)
    if (!transcript) {
      unmatched.add(call.protein)
      continue
    }
    const blocks = projectResidues(transcript, call.start, call.end)
    if (blocks.length === 0) {
      pastCds++
      continue
    }
    placedCalls.push({ transcript, call, blocks })
  }
  if (placedCalls.length === 0) {
    throw new Error(
      `No protein features in ${calls} matched a transcript in ${annotation}${
        unmatched.size > 0
          ? `; protein names such as ${[...unmatched].slice(0, 3).join(', ')} match no transcript ID, protein_id or CDS ID`
          : ''
      }`,
    )
  }
  const rows: Gff3Row[] = gff3
    ? [
        ...codingTranscripts.flatMap(t => transcriptGff3Rows(t)),
        ...placedCalls.flatMap(p => proteinMatchGff3Rows(p)),
      ]
    : placedCalls.map(({ transcript, call, blocks }) => ({
        refName: transcript.refName,
        start: blocks[0]!.start,
        depth: 0,
        line: bedLine(transcript, call, blocks),
      }))
  // the whole line breaks a tie, so the same calls in any order (hmmscan
  // chunks concatenated, say) write the same file
  const compare = (a: string, b: string) => (a < b ? -1 : a > b ? 1 : 0)
  rows.sort(
    (a, b) =>
      compare(a.refName, b.refName) ||
      a.start - b.start ||
      a.depth - b.depth ||
      compare(a.line, b.line),
  )
  return {
    header: gff3 ? GFF3_HEADER : `#${BED_HEADER.join('\t')}`,
    lines: rows.map(r => r.line),
    features: placedCalls.length,
    transcripts,
    unmatched,
    pastCds,
    skippedSource,
  }
}

/** The warnings a run of `placeProteinFeatures` earns, one string each */
export function placementWarnings({ unmatched, pastCds }: PlacedFeatures) {
  return [
    ...(unmatched.size > 0
      ? [
          `${unmatched.size} protein(s) match no transcript, e.g. ${[...unmatched].slice(0, 3).join(', ')}`,
        ]
      : []),
    ...(pastCds > 0
      ? [
          `${pastCds} feature(s) start past the end of their transcript's CDS; is the annotation the one the proteins came from?`,
        ]
      : []),
  ]
}

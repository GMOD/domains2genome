import {
  bedLine,
  bestOfOverlappingHmmerHits,
  findTranscript,
  parseCall,
  projectResidues,
  readCodingTranscripts,
  readLines,
} from './project.ts'

import type { ProteinCall } from './project.ts'

export interface PlacedFeatures {
  /** BED12 rows, sorted by contig then start, without the header */
  lines: string[]
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
 * (`Pfam`, `SMART`), ignoring case.
 */
export async function placeProteinFeatures(
  annotation: string,
  calls: string,
  { sources }: { sources?: string[] } = {},
): Promise<PlacedFeatures> {
  const { lookup, transcripts } = await readCodingTranscripts(annotation)
  if (transcripts === 0) {
    throw new Error(
      `${annotation} has no CDS rows to place protein features on`,
    )
  }
  const kept = sources?.map(s => s.toLowerCase())
  const isProtein = (id: string) => findTranscript(lookup, id) !== undefined
  const parsed: ProteinCall[] = []
  let skippedSource = 0
  for await (const line of readLines(calls)) {
    const call = parseCall(line, isProtein)
    if (call) {
      if (kept && !kept.includes(call.source.toLowerCase())) {
        skippedSource++
      } else {
        parsed.push(call)
      }
    }
  }

  const rows: { refName: string; start: number; line: string }[] = []
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
    rows.push({
      refName: transcript.refName,
      start: blocks[0]!.start,
      line: bedLine(transcript, call, blocks),
    })
  }
  if (rows.length === 0) {
    throw new Error(
      `No protein features in ${calls} matched a transcript in ${annotation}${
        unmatched.size > 0
          ? `; protein names such as ${[...unmatched].slice(0, 3).join(', ')} match no transcript ID, protein_id or CDS ID`
          : ''
      }`,
    )
  }
  rows.sort((a, b) =>
    a.refName === b.refName
      ? a.start - b.start
      : a.refName < b.refName
        ? -1
        : 1,
  )
  return {
    lines: rows.map(r => r.line),
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

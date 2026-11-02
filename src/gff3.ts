import type { CodingTranscript, ProteinCall } from './project.ts'

export const GFF3_HEADER = '##gff-version 3'

/** The `protein_match` a placed protein feature is written as */
export const PROTEIN_MATCH = 'protein_match'

export interface PlacedCall {
  transcript: CodingTranscript
  call: ProteinCall
  blocks: { start: number; end: number }[]
}

export interface Gff3Row {
  refName: string
  start: number
  /** a parent sorts before the children that start where it does */
  depth: 0 | 1
  line: string
}

const percentEncode = (c: string) =>
  `%${c.charCodeAt(0).toString(16).toUpperCase().padStart(2, '0')}`

// GFF3 reserves these characters in a column-9 value; the rest stand as
// written, so a description keeps its spaces
const attributeValue = (v: string) =>
  v.replaceAll(/[%;=&,\t\n\r]/g, percentEncode)

const columnValue = (v: string) => v.replaceAll(/[%\t\n\r]/g, percentEncode)

const isStated = (v: string) => v !== '' && v !== '-'

function attributes(pairs: [string, string | number][]) {
  return pairs
    .filter(([, v]) => isStated(String(v)))
    .map(([k, v]) => `${k}=${attributeValue(String(v))}`)
    .join(';')
}

function gff3Row(
  transcript: CodingTranscript,
  {
    source,
    type,
    start,
    end,
    score,
    phase,
    depth,
    pairs,
  }: {
    source: string
    type: string
    start: number
    end: number
    score: string
    phase: string
    depth: 0 | 1
    pairs: [string, string | number][]
  },
): Gff3Row {
  return {
    refName: transcript.refName,
    start,
    depth,
    line: [
      columnValue(transcript.refName),
      columnValue(isStated(source) ? source : '.'),
      type,
      start + 1,
      end,
      isStated(score) ? score : '.',
      transcript.strand === 1 ? '+' : '-',
      phase,
      attributes(pairs),
    ].join('\t'),
  }
}

/**
 * A coding transcript as GFF3: an mRNA spanning its exons, with its exons
 * and its CDS, so a browser draws the transcript its protein features are
 * placed on.
 */
export function transcriptGff3Rows(transcript: CodingTranscript): Gff3Row[] {
  const parent: [string, string] = ['Parent', transcript.id]
  const { exons } = transcript
  return [
    gff3Row(transcript, {
      source: '.',
      type: 'mRNA',
      start: Math.min(transcript.start, ...exons.map(e => e.start)),
      end: Math.max(transcript.end, ...exons.map(e => e.end)),
      score: '.',
      phase: '.',
      depth: 0,
      pairs: [
        ['ID', transcript.id],
        ['Name', transcript.name],
        ['protein_id', transcript.protein],
        ['product', transcript.product],
      ],
    }),
    ...exons.map(exon =>
      gff3Row(transcript, {
        source: '.',
        type: 'exon',
        start: exon.start,
        end: exon.end,
        score: '.',
        phase: '.',
        depth: 1,
        pairs: [parent],
      }),
    ),
    ...transcript.subfeatures.map(cds =>
      gff3Row(transcript, {
        source: '.',
        type: 'CDS',
        start: cds.start,
        end: cds.end,
        score: '.',
        phase: String(cds.phase ?? 0),
        depth: 1,
        pairs: [parent],
      }),
    ),
  ]
}

/**
 * A protein feature placed on its transcript as GFF3: a `protein_match`
 * child of the mRNA, one line per block sharing one ID, as GFF3 writes a
 * feature split across introns.
 */
export function proteinMatchGff3Rows({
  transcript,
  call,
  blocks,
}: PlacedCall): Gff3Row[] {
  const id = `${transcript.id}:${call.accession}:${call.start}-${call.end}`
  return blocks.map(block =>
    gff3Row(transcript, {
      source: call.source,
      type: PROTEIN_MATCH,
      start: block.start,
      end: block.end,
      score: call.evalue,
      phase: '.',
      depth: 1,
      pairs: [
        ['ID', id],
        ['Parent', transcript.id],
        ['Name', call.name],
        ['accession', call.accession],
        ['description', call.description],
        ['proteinStart', call.start],
        ['proteinEnd', call.end],
      ],
    }),
  )
}

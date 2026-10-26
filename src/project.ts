import { createReadStream } from 'node:fs'
import { createInterface } from 'node:readline'
import { createGunzip } from 'node:zlib'

import { getAttribute, parseFeatureLazy } from '@gmod/gff-nostream'
import { proteinRangeToGenome } from 'g2p_mapper'

import type { Feat } from 'g2p_mapper'

export interface CodingTranscript extends Feat {
  id: string
  /** JBrowse's default label for the transcript, its `Name` else its ID, or its ID where another transcript writes that Name */
  name: string
  strand: 1 | -1
  subfeatures: Feat[]
}

export interface ProteinCall {
  protein: string
  start: number
  end: number
  accession: string
  name: string
  description: string
  source: string
  evalue: string
}

export async function* readLines(file: string) {
  const raw = createReadStream(file)
  const input = /\.b?gz$/i.test(file) ? raw.pipe(createGunzip()) : raw
  yield* createInterface({ input, crlfDelay: Infinity })
}

interface AnnotationRow {
  type: string | null
  refName: string
  start: number
  end: number
  strand: 1 | -1
  phase: number
  attribute: (key: string) => string[]
}

const asList = (v: string | string[] | undefined) =>
  v === undefined ? [] : Array.isArray(v) ? v : [v]

// gff-nostream is the parser JBrowse reads the gene track with, so a tag's
// case, escapes and multiple values come out as the gene track's labels do
function gff3Row(line: string): AnnotationRow {
  const feature = parseFeatureLazy(line)
  return {
    type: feature.type,
    refName: feature.refName,
    start: feature.start,
    end: feature.end,
    strand: feature.strand === -1 ? -1 : 1,
    phase: feature.phase ?? 0,
    attribute: key => asList(getAttribute(feature, key)),
  }
}

function gtfRow(cols: string[]): AnnotationRow {
  const attrs = new Map<string, string>()
  for (const m of cols[8]!.matchAll(/(\S+)\s+"([^"]*)"/g)) {
    attrs.set(m[1]!.toLowerCase(), m[2]!)
  }
  return {
    type: cols[2]!,
    refName: cols[0]!,
    start: +cols[3]! - 1,
    end: +cols[4]!,
    strand: cols[6] === '-' ? -1 : 1,
    phase: cols[7] === '.' ? 0 : +cols[7]!,
    attribute: key => asList(attrs.get(key)),
  }
}

const withoutPrefix = (id: string) =>
  id.replace(/^(cds|rna|transcript|CDS)[-:]/, '')
const withoutVersion = (id: string) => id.replace(/\.\d+$/, '')

/**
 * Every spelling a protein FASTA might use for a transcript's product: the
 * transcript ID (gffread, funannotate, BRAKER), the CDS's `protein_id` (NCBI,
 * Ensembl) or its ID, each also without an Ensembl/NCBI `type:`/`type-`
 * prefix and without a `.N` version.
 */
function proteinKeys(ids: string[]) {
  return new Set(
    ids.flatMap(id => {
      const bare = withoutPrefix(id)
      return [id, bare, withoutVersion(id), withoutVersion(bare)]
    }),
  )
}

/** Reads the CDS rows of a GFF3 or GTF into coding transcripts, keyed by every name their protein might go by. */
export async function readCodingTranscripts(file: string) {
  const byTranscript = new Map<
    string,
    { transcript: CodingTranscript; aliases: string[] }
  >()
  const gff3Names = new Map<string, string>()
  for await (const line of readLines(file)) {
    const cols = line.split('\t')
    if (line.startsWith('#') || cols.length < 9 || cols[2] === 'exon') {
      continue
    }
    const isGtf = /^\s*\w+\s+"/.test(cols[8]!)
    const row = isGtf ? gtfRow(cols) : gff3Row(line)
    if (row.type !== 'CDS') {
      const [id] = row.attribute('id')
      const [name] = row.attribute('name')
      if (!isGtf && id && name) {
        gff3Names.set(id, name)
      }
      continue
    }
    const parents = row.attribute(isGtf ? 'transcript_id' : 'parent')
    const aliases = [
      ...row.attribute('protein_id'),
      ...(isGtf ? [] : row.attribute('id')),
    ]
    const cds = {
      refName: row.refName,
      start: row.start,
      end: row.end,
      type: 'CDS',
      phase: row.phase,
    }
    for (const parent of parents) {
      const entry = byTranscript.get(parent) ?? {
        transcript: {
          id: parent,
          name: parent,
          refName: cds.refName,
          start: cds.start,
          end: cds.end,
          strand: row.strand,
          subfeatures: [],
        },
        aliases: [],
      }
      byTranscript.set(parent, entry)
      const { transcript } = entry
      transcript.start = Math.min(transcript.start, cds.start)
      transcript.end = Math.max(transcript.end, cds.end)
      transcript.subfeatures.push(cds)
      entry.aliases.push(...aliases)
    }
  }
  // a Name another transcript writes, as its Name or its ID, is not the
  // transcript's alone
  const nameCounts = new Map<string, number>()
  for (const id of byTranscript.keys()) {
    const name = gff3Names.get(id)
    for (const written of name && name !== id ? [id, name] : [id]) {
      nameCounts.set(written, (nameCounts.get(written) ?? 0) + 1)
    }
  }
  const lookup = new Map<string, CodingTranscript>()
  for (const { transcript, aliases } of byTranscript.values()) {
    const name = gff3Names.get(transcript.id)
    if (name && nameCounts.get(name) === 1) {
      transcript.name = name
    }
    for (const key of proteinKeys([transcript.id, ...aliases])) {
      lookup.set(key, transcript)
    }
  }
  return { lookup, transcripts: byTranscript.size }
}

export function findTranscript(
  lookup: Map<string, CodingTranscript>,
  protein: string,
) {
  for (const key of proteinKeys([protein])) {
    const hit = lookup.get(key)
    if (hit) {
      return hit
    }
  }
  return undefined
}

/**
 * A feature row of InterProScan's GFF3. Version 6 writes the signature's short
 * name in `Name` and its accession in `Alias`; version 5 writes the accession
 * in `Name` and the description in `signature_desc`. A `polypeptide` row is
 * the protein itself.
 */
function interProGffCall(
  line: string,
  cols: string[],
): ProteinCall | undefined {
  const [protein, source, type, start, end, score] = cols
  if (type === 'polypeptide') {
    return undefined
  }
  const { attribute } = gff3Row(line)
  const [name = ''] = attribute('name')
  const [alias] = attribute('alias')
  const [description = ''] = attribute('signature_desc')
  return {
    protein: protein!,
    start: +start!,
    end: +end!,
    accession: alias ?? name,
    name: alias ? name : description || name,
    description,
    source: source!,
    evalue: score === '.' ? '-' : score!,
  }
}

/**
 * One line of InterProScan's TSV or GFF3, or of HMMER's `--domtblout` from
 * either hmmscan (query is the protein) or hmmsearch (target is the protein).
 * `isProtein` decides which side of a HMMER row names the protein.
 */
export function parseCall(
  line: string,
  isProtein: (id: string) => boolean,
): ProteinCall | undefined {
  if (line.startsWith('#') || line.trim() === '') {
    return undefined
  }
  const tabs = line.split('\t')
  if (tabs.length === 9) {
    return interProGffCall(line, tabs)
  }
  if (tabs.length >= 11) {
    const [protein, , , source, accession, description, start, end, evalue] =
      tabs
    return {
      protein: protein!,
      start: +start!,
      end: +end!,
      accession: accession!,
      name: description && description !== '-' ? description : accession!,
      description: tabs[12] && tabs[12] !== '-' ? tabs[12] : description!,
      source: source!,
      evalue: evalue!,
    }
  }
  const cols = line.trim().split(/\s+/)
  if (cols.length < 23) {
    return undefined
  }
  const hmmIsTarget = isProtein(cols[3]!) || !isProtein(cols[0]!)
  const [name, accession, protein] = hmmIsTarget
    ? [cols[0]!, cols[1]!, cols[3]!]
    : [cols[3]!, cols[4]!, cols[0]!]
  return {
    protein,
    start: +cols[19]!,
    end: +cols[20]!,
    accession,
    name,
    description: cols.slice(22).join(' '),
    source: 'HMMER',
    evalue: cols[12]!,
  }
}

/**
 * HMMER reports every family that clears its threshold, so related families
 * (SH3_1, SH3_2, SH3_9) stack on one stretch of protein. Where two hits on a
 * protein share more than half the shorter one, the lower E-value wins,
 * whatever their Pfam clans: pfam_scan and InterProScan resolve only hits in
 * one clan, which needs Pfam's clan file. Calls from InterProScan pass
 * through, already resolved.
 */
export function bestOfOverlappingHmmerHits(calls: ProteinCall[]) {
  const kept: ProteinCall[] = []
  const byProtein = new Map<string, ProteinCall[]>()
  for (const call of calls) {
    if (call.source === 'HMMER') {
      const list = byProtein.get(call.protein) ?? []
      list.push(call)
      byProtein.set(call.protein, list)
    } else {
      kept.push(call)
    }
  }
  for (const hits of byProtein.values()) {
    const accepted: ProteinCall[] = []
    for (const hit of hits.sort((a, b) => +a.evalue - +b.evalue)) {
      const clashes = accepted.some(a => {
        const shared =
          Math.min(a.end, hit.end) - Math.max(a.start, hit.start) + 1
        const shorter = Math.min(a.end - a.start, hit.end - hit.start) + 1
        return shared * 2 > shorter
      })
      if (!clashes) {
        accepted.push(hit)
      }
    }
    kept.push(...accepted)
  }
  return kept
}

/**
 * The genome blocks under residues `startAa..endAa` (1-based, inclusive, as
 * domain callers report them), in genome order. Callers run on gffread, NCBI
 * or Ensembl proteins, which start at the first complete codon.
 */
export function projectResidues(
  transcript: CodingTranscript,
  startAa: number,
  endAa: number,
) {
  return proteinRangeToGenome(transcript, startAa - 1, endAa, {
    partialCodon: 'dropped',
  }).map(([start, end]) => ({ start, end }))
}

export const BED_HEADER = [
  'chrom',
  'chromStart',
  'chromEnd',
  'name',
  'score',
  'strand',
  'thickStart',
  'thickEnd',
  'itemRgb',
  'blockCount',
  'blockSizes',
  'chromStarts',
  'accession',
  'description',
  'source',
  'transcript',
  'proteinStart',
  'proteinEnd',
  'evalue',
]

/** BED's 0..1000 score as ten times -log10(E), so 1e-50 scores 500 */
export function evalueScore(evalue: string) {
  const score = Math.round(-10 * Math.log10(Number(evalue)))
  return Number.isNaN(score) ? 0 : Math.min(1000, Math.max(0, score))
}

const tsvField = (s: string) => s.replaceAll(/[\t\n\r]/g, ' ')

export function bedLine(
  transcript: CodingTranscript,
  call: ProteinCall,
  blocks: { start: number; end: number }[],
) {
  const start = blocks[0]!.start
  const end = blocks.at(-1)!.end
  return [
    transcript.refName,
    start,
    end,
    tsvField(call.name),
    evalueScore(call.evalue),
    transcript.strand === 1 ? '+' : '-',
    start,
    end,
    0,
    blocks.length,
    blocks.map(b => b.end - b.start).join(','),
    blocks.map(b => b.start - start).join(','),
    tsvField(call.accession),
    tsvField(call.description),
    tsvField(call.source),
    tsvField(transcript.name),
    call.start,
    call.end,
    call.evalue,
  ].join('\t')
}

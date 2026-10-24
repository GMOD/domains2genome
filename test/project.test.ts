import { mkdtempSync, writeFileSync } from 'node:fs'
import { tmpdir } from 'node:os'
import path from 'node:path'

import { expect, test } from 'vitest'

import {
  bestOfOverlappingHmmerHits,
  evalueScore,
  findTranscript,
  parseCall,
  projectResidues,
  readCodingTranscripts,
  readLines,
} from '../src/index.ts'
import { gff3 } from './fixtures.ts'

function tmpFile(name: string, text: string) {
  const file = path.join(mkdtempSync(path.join(tmpdir(), 'ips2g-')), name)
  writeFileSync(file, text)
  return file
}

test('reads CDS parts and matches a protein by any of its names', async () => {
  const { lookup, transcripts } = await readCodingTranscripts(
    tmpFile('genes.gff3', gff3),
  )
  expect(transcripts).toBe(3)
  const ncbi = findTranscript(lookup, 'XP_1.1')!
  expect(ncbi.id).toBe('rna-XM_1.1')
  // residues 8..14 are CDS bases 21..42: 9 in the first part, 12 in the second
  expect(projectResidues(ncbi, 8, 14)).toEqual([
    { start: 121, end: 130 },
    { start: 200, end: 212 },
  ])
  expect(findTranscript(lookup, 'XM_1')).toBe(ncbi)
  // a minus-strand CDS translates from its rightmost part
  expect(projectResidues(findTranscript(lookup, 'minus')!, 1, 12)).toEqual([
    { start: 324, end: 330 },
    { start: 400, end: 430 },
  ])
  // the leading phase base belongs to no residue
  expect(projectResidues(findTranscript(lookup, 'partial')!, 1, 1)).toEqual([
    { start: 11, end: 14 },
  ])
  expect(findTranscript(lookup, 'XP_2')).toBeUndefined()
})

test('reads a GTF', async () => {
  const { lookup } = await readCodingTranscripts(
    tmpFile(
      'genes.gtf',
      'chr1\t.\tCDS\t11\t40\t.\t+\t0\tgene_id "g"; transcript_id "T1"; protein_id "P1";\n',
    ),
  )
  expect(findTranscript(lookup, 'P1')!.id).toBe('T1')
})

const hmmscanRow =
  'Pkinase PF00069.29 264 XP_1.1 - 600 1e-50 170.2 0.1 1 1 1e-53 2e-50 169.9 0.1 1 264 4 9 3 10 0.95 Protein kinase domain'
const hmmsearchRow =
  'XP_1.1 - 600 Pkinase PF00069.29 264 1e-50 170.2 0.1 1 1 1e-53 2e-50 169.9 0.1 1 264 4 9 3 10 0.95 Protein kinase domain'

test('a HMMER row names the protein on whichever side it sits', () => {
  const isProtein = (id: string) => id === 'XP_1.1'
  const expected = {
    protein: 'XP_1.1',
    start: 3,
    end: 10,
    accession: 'PF00069.29',
    name: 'Pkinase',
    description: 'Protein kinase domain',
    source: 'HMMER',
    evalue: '2e-50',
  }
  expect(parseCall(hmmscanRow, isProtein)).toEqual(expected)
  expect(parseCall(hmmsearchRow, isProtein)).toEqual(expected)
  expect(parseCall('# comment', isProtein)).toBeUndefined()
})

// rows from interproscan6's own tests/data/output (Apache-2.0)
test('reads an InterProScan TSV row', () => {
  const pfam =
    'sp|Q9QXX0|JAG1_MOUSE\t0001C5D32DCCB8E03267C370872E1694\t1218\tPfam\tPF07657\tN terminus of Notch ligand C2-like domain\t32\t185\t9.0E-55\tPfam\t12-04-2026\tIPR011651\tNotch ligand, N-terminal domain\t-\t-'
  expect(parseCall(pfam, () => false)).toEqual({
    protein: 'sp|Q9QXX0|JAG1_MOUSE',
    start: 32,
    end: 185,
    accession: 'PF07657',
    name: 'N terminus of Notch ligand C2-like domain',
    description: 'Notch ligand, N-terminal domain',
    source: 'Pfam',
    evalue: '9.0E-55',
  })
})

test('reads an InterProScan 6 GFF3 row, its short name the name', () => {
  const row =
    'sp|Q9QXX0|JAG1_MOUSE\tPfam\tpolypeptide_domain\t32\t185\t9.0E-55\t.\t.\tName=MNNL;Alias=PF07657;Dbxref=InterPro:IPR011651;type=Domain;representative=false'
  expect(parseCall(row, () => false)).toEqual({
    protein: 'sp|Q9QXX0|JAG1_MOUSE',
    start: 32,
    end: 185,
    accession: 'PF07657',
    name: 'MNNL',
    description: '',
    source: 'Pfam',
    evalue: '9.0E-55',
  })
})

test('reads an InterProScan 5 GFF3 row and skips the protein itself', () => {
  const row =
    'P51587\tPfam\tprotein_match\t1002\t1036\t1.1E-6\t+\t.\tTarget=P51587 1002 1036;ID=match$1_1002_1036;signature_desc=BRCA2 repeat;Name=PF00634;status=T'
  expect(parseCall(row, () => false)).toMatchObject({
    protein: 'P51587',
    accession: 'PF00634',
    name: 'BRCA2 repeat',
    evalue: '1.1E-6',
  })
  expect(
    parseCall(
      'P51587\t.\tpolypeptide\t1\t3418\t.\t+\t.\tID=P51587',
      () => false,
    ),
  ).toBeUndefined()
})

// InterProScan scores PANTHER differently in its two formats, so E-values are
// compared on the Pfam rows alone
test("InterProScan's TSV and GFF3 of one run give the same calls", async () => {
  const keys = async (file: string) => {
    const found: string[] = []
    for await (const line of readLines(file)) {
      const c = parseCall(line, () => false)
      if (c) {
        const evalue = c.source === 'Pfam' ? c.evalue : ''
        found.push(
          [c.protein, c.source, c.accession, c.start, c.end, evalue].join(' '),
        )
      }
    }
    return found.sort()
  }
  const data = path.join(import.meta.dirname, 'data')
  const tsv = await keys(path.join(data, 'jag1.tsv'))
  expect(tsv).toHaveLength(169)
  expect(await keys(path.join(data, 'jag1.gff3'))).toEqual(tsv)
})

test('the best of overlapping HMMER hits wins; InterProScan calls pass through', () => {
  const hit = (name: string, start: number, end: number, evalue: string) => ({
    protein: 'p',
    start,
    end,
    accession: name,
    name,
    description: '',
    source: 'HMMER',
    evalue,
  })
  const kept = bestOfOverlappingHmmerHits([
    hit('SH3_1', 10, 60, '1e-10'),
    hit('SH3_2', 12, 62, '1e-12'),
    hit('SH2', 58, 150, '1e-5'),
    { ...hit('PF1', 10, 60, '1'), source: 'Pfam' },
  ])
  expect(kept.map(k => k.name).sort()).toEqual(['PF1', 'SH2', 'SH3_2'])
})

test('the BED score is ten times -log10(E), within 0..1000', () => {
  expect(evalueScore('1e-50')).toBe(500)
  expect(evalueScore('0')).toBe(1000)
  expect(evalueScore('2')).toBe(0)
  expect(evalueScore('-')).toBe(0)
})

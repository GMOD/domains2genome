import { spawnSync } from 'node:child_process'
import { mkdtempSync, readFileSync, writeFileSync } from 'node:fs'
import { tmpdir } from 'node:os'
import path from 'node:path'

import { expect, test } from 'vitest'

import { gff3 } from './fixtures.ts'

const cli = path.join(import.meta.dirname, '..', 'src', 'cli.ts')

function run(args: string[], files: Record<string, string>) {
  const dir = mkdtempSync(path.join(tmpdir(), 'ips2g-cli-'))
  for (const [name, text] of Object.entries(files)) {
    writeFileSync(path.join(dir, name), text)
  }
  const result = spawnSync(process.execPath, [cli, ...args], {
    cwd: dir,
    encoding: 'utf8',
  })
  return { ...result, dir }
}

const calls = [
  'XP_1.1\tm\t20\tPfam\tPF00069\tProtein kinase domain\t8\t14\t1E-5\tT\td\t-\t-',
  'XP_1.1\tm\t20\tSMART\tSM00220\tS_TKc\t8\t14\t1E-5\tT\td\t-\t-',
  'nobody\tm\t20\tPfam\tPF1\tX\t1\t2\t1\tT\td\t-\t-',
].join('\n')

const kinase = [
  'ctgA',
  121,
  212,
  'Protein kinase domain',
  50,
  '+',
  121,
  212,
  0,
  2,
  '9,12',
  '0,79',
  'PF00069',
  'Protein kinase domain',
  'Pfam',
  'XM_1.1',
  8,
  14,
  '1E-5',
].join('\t')

test('writes BED12 a block per exon to stdout, warnings to stderr', () => {
  const r = run(['genes.gff3', 'calls.tsv', '--source', 'pfam'], {
    'genes.gff3': gff3,
    'calls.tsv': calls,
  })
  expect(r.status).toBe(0)
  const lines = r.stdout.trim().split('\n')
  expect(lines[0]!.split('\t').slice(0, 4)).toEqual([
    '#chrom',
    'chromStart',
    'chromEnd',
    'name',
  ])
  expect(lines.slice(1)).toEqual([kinase])
  expect(r.stderr).toContain('nobody')
})

test('--gff3 writes each coding transcript with its exons and CDS, and its features as protein_match children', () => {
  const genes = [
    '##gff-version 3',
    'ctgA\t.\tmRNA\t101\t230\t.\t+\t.\tID=rna-XM_1.1;Name=XM_1.1',
    'ctgA\t.\texon\t101\t140\t.\t+\t.\tParent=rna-XM_1.1',
    'ctgA\t.\texon\t191\t260\t.\t+\t.\tParent=rna-XM_1.1',
    'ctgA\t.\tCDS\t201\t230\t.\t+\t0\tID=cds-XP_1.1;Parent=rna-XM_1.1;protein_id=XP_1.1;product=kinase%3B isoform X1',
    'ctgA\t.\tCDS\t101\t130\t.\t+\t0\tID=cds-XP_1.1;Parent=rna-XM_1.1;protein_id=XP_1.1;product=kinase%3B isoform X1',
  ].join('\n')
  const r = run(['genes.gff3', 'calls.tsv', '--source', 'pfam', '--gff3'], {
    'genes.gff3': genes,
    'calls.tsv': calls,
  })
  expect(r.status).toBe(0)
  const kinaseMatch =
    'ID=rna-XM_1.1:PF00069:8-14;Parent=rna-XM_1.1;Name=Protein kinase domain;accession=PF00069;description=Protein kinase domain;proteinStart=8;proteinEnd=14'
  expect(r.stdout.trim().split('\n')).toEqual([
    '##gff-version 3',
    'ctgA\t.\tmRNA\t101\t260\t.\t+\t.\tID=rna-XM_1.1;Name=XM_1.1;protein_id=XP_1.1;product=kinase%3B isoform X1',
    'ctgA\t.\tCDS\t101\t130\t.\t+\t0\tParent=rna-XM_1.1',
    'ctgA\t.\texon\t101\t140\t.\t+\t.\tParent=rna-XM_1.1',
    `ctgA\tPfam\tprotein_match\t122\t130\t1E-5\t+\t.\t${kinaseMatch}`,
    'ctgA\t.\texon\t191\t260\t.\t+\t.\tParent=rna-XM_1.1',
    'ctgA\t.\tCDS\t201\t230\t.\t+\t0\tParent=rna-XM_1.1',
    `ctgA\tPfam\tprotein_match\t201\t212\t1E-5\t+\t.\t${kinaseMatch}`,
  ])
  expect(r.stderr).toContain(
    'Placed 1 protein features on 1 coding transcripts',
  )
})

test('a --source matching no call names the sources the calls come from', () => {
  const r = run(['genes.gff3', 'calls.tsv', '--source', 'gene3d'], {
    'genes.gff3': gff3,
    'calls.tsv': calls,
  })
  expect(r.status).not.toBe(0)
  expect(r.stderr).toContain(
    'No call in calls.tsv comes from gene3d; its calls come from Pfam, SMART',
  )
})

test('the 1.x flags still run, and --out writes a file', () => {
  const r = run(
    [
      '--genomeGFF',
      'genes.gff3',
      '--interproGFF',
      'calls.tsv',
      '--out',
      'o.bed',
    ],
    { 'genes.gff3': gff3, 'calls.tsv': calls },
  )
  expect(r.status).toBe(0)
  expect(r.stdout).toBe('')
  expect(readFileSync(path.join(r.dir, 'o.bed'), 'utf8')).toContain(kinase)
})

test('a calls file holding no call says so, naming the formats it reads', () => {
  const r = run(['genes.gff3', 'proteins.faa'], {
    'genes.gff3': gff3,
    'proteins.faa': '>XP_1.1\nMKVLA\n',
  })
  expect(r.status).toBe(1)
  expect(r.stderr).toContain(
    'No protein features could be read from proteins.faa; it should be InterProScan TSV or GFF3, or HMMER --domtblout',
  )
})

test('an empty --source keeps every analysis', () => {
  const r = run(['genes.gff3', 'calls.tsv', '--source', ''], {
    'genes.gff3': gff3,
    'calls.tsv': calls,
  })
  expect(r.status).toBe(0)
  expect(r.stderr).toContain('Placed 2 protein features')
})

test('--help prints the usage to stdout', () => {
  const r = run(['--help'], {})
  expect(r.status).toBe(0)
  expect(r.stdout).toContain('Usage: domains2genome')
  expect(r.stderr).toBe('')
})

test('refuses calls that match no transcript', () => {
  const r = run(['genes.gff3', 'calls.tsv'], {
    'genes.gff3': gff3,
    'calls.tsv': 'nobody\tm\t20\tPfam\tPF1\tX\t1\t2\t1\tT\td\t-\t-\n',
  })
  expect(r.status).toBe(1)
  expect(r.stderr).toContain('nobody')
})

test('the same calls in any order write the same BED', () => {
  const rows = calls.split('\n')
  const forward = run(['genes.gff3', 'calls.tsv'], {
    'genes.gff3': gff3,
    'calls.tsv': rows.join('\n'),
  })
  const reversed = run(['genes.gff3', 'calls.tsv'], {
    'genes.gff3': gff3,
    'calls.tsv': rows.toReversed().join('\n'),
  })
  expect(forward.stdout.split('\n').length).toBeGreaterThan(3)
  expect(reversed.stdout).toBe(forward.stdout)
})

test('a protein named two ways is one protein: its calls are placed once, and its HMMER hits filtered together', () => {
  const pfam = calls.split('\n')[0]!
  const tsv = run(['genes.gff3', 'calls.tsv', '--gff3'], {
    'genes.gff3': gff3,
    'calls.tsv': [pfam, pfam.replace('XP_1.1', 'rna-XM_1.1')].join('\n'),
  })
  expect(
    tsv.stdout.split('\n').filter(l => l.includes('protein_match')),
  ).toHaveLength(2)
  const domtblout = [
    'Pkinase PF00069.29 264 XP_1.1 - 20 1 1 1 1 1 1 1e-10 1 1 1 264 4 9 8 14 0.9 Protein kinase domain',
    'Pkinase_Tyr PF07714.21 264 rna-XM_1.1 - 20 1 1 1 1 1 1 1e-20 1 1 1 264 4 9 8 14 0.9 Protein tyrosine kinase',
  ].join('\n')
  const hmmer = run(['genes.gff3', 'calls.domtblout'], {
    'genes.gff3': gff3,
    'calls.domtblout': domtblout,
  })
  const names = hmmer.stdout
    .trim()
    .split('\n')
    .slice(1)
    .map(l => l.split('\t')[3])
  expect(names).toEqual(['Pkinase_Tyr'])
})

test.each([
  [
    ['genes.gff3', 'calls.tsv', '--soruce', 'Pfam'],
    "Unknown option '--soruce'",
  ],
  [['genes.gff3', 'calls.tsv', '--out', 'no/such/dir/o.bed'], 'ENOENT'],
  [['genes.gff3', 'missing.tsv.gz'], 'ENOENT'],
])('%j fails with one line of error, not a stack trace', (args, message) => {
  const r = run(args, { 'genes.gff3': gff3, 'calls.tsv': calls })
  expect(r.status).toBe(1)
  expect(r.stderr).toMatch(new RegExp(`^Error: .*${message}`))
  expect(r.stderr).not.toContain('    at ')
})

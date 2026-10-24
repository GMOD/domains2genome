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
  'rna-XM_1.1',
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

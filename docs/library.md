# Using domains2genome as a library

The package exports what the command line runs, for both `import` and `require`,
with TypeScript types.

## placeProteinFeatures

```ts
import { writeFile } from 'node:fs/promises'
import { placeProteinFeatures, placementWarnings } from 'domains2genome'

const placed = await placeProteinFeatures('genes.gff3', 'interpro.tsv', {
  sources: ['Pfam'],
})
await writeFile('domains.bed', [placed.header, ...placed.lines, ''].join('\n'))
for (const warning of placementWarnings(placed)) {
  console.warn(warning)
}
```

`placeProteinFeatures(annotation, calls, options)` reads both files (plain or
gzipped) and takes two options:

| Option    | Effect                                                      |
| --------- | ----------------------------------------------------------- |
| `sources` | keep only these analyses, ignoring case, as `--source` does |
| `gff3`    | return GFF3 rows instead of BED12, as `--gff3` does         |

It resolves to a `PlacedFeatures`:

| Field           | Holds                                                                |
| --------------- | -------------------------------------------------------------------- |
| `header`        | the first line of the file: BED's column names, or `##gff-version 3` |
| `lines`         | the rows, sorted by contig then start, without the header            |
| `features`      | how many protein features were placed                                |
| `transcripts`   | how many coding transcripts the annotation holds                     |
| `unmatched`     | the protein names that matched no transcript                         |
| `pastCds`       | how many calls start past the end of their transcript's CDS          |
| `skippedSource` | how many calls `sources` left out                                    |

`placementWarnings(placed)` turns `unmatched` and `pastCds` into the warnings
the command line prints.

It throws, as the command line exits with an error, when the annotation has no
CDS, when no call can be read from `calls`, when `sources` leaves out every
call, or when no call matches a transcript.

## Building blocks

To place calls one at a time, or from a source other than a file, use the steps
`placeProteinFeatures` is made of:

```ts
import {
  bedLine,
  findTranscript,
  parseCall,
  projectResidues,
  readCodingTranscripts,
} from 'domains2genome'

const { lookup } = await readCodingTranscripts('genes.gff3')
const isProtein = (id: string) => findTranscript(lookup, id) !== undefined

const call = parseCall(domtbloutLine, isProtein)
const transcript = call && findTranscript(lookup, call.protein)
if (call && transcript) {
  const blocks = projectResidues(transcript, call.start, call.end)
  if (blocks.length > 0) {
    console.log(bedLine(transcript, call, blocks))
  }
}
```

`projectResidues` returns no blocks for a call that starts past the end of the
CDS, which usually means the proteins came from another annotation.

| Export                       | Does                                                                                            |
| ---------------------------- | ----------------------------------------------------------------------------------------------- |
| `readCodingTranscripts`      | reads a GFF3 or GTF into coding transcripts, and a `lookup` by every name their protein goes by |
| `findTranscript`             | finds a protein's transcript in that `lookup`, also without a `.N` version or a `cds-` prefix   |
| `parseCall`                  | reads one line of InterProScan TSV or GFF3, or of HMMER `--domtblout`, into a `ProteinCall`     |
| `bestOfOverlappingHmmerHits` | keeps the lowest E-value where HMMER hits on one protein overlap by more than half              |
| `projectResidues`            | the genome blocks under a 1-based, inclusive residue range, in genome order                     |
| `bedLine`, `BED_HEADER`      | a placed call as a BED12 row, and the column names                                              |
| `transcriptGff3Rows`         | a coding transcript as GFF3 `mRNA`, `exon` and `CDS` rows                                       |
| `proteinMatchGff3Rows`       | a placed call as GFF3 `protein_match` rows                                                      |
| `evalueScore`                | BED's 0..1000 score for an E-value                                                              |
| `readLines`                  | the lines of a plain or gzipped file                                                            |

`parseCall` returns `undefined` for a comment, a blank line, or a row that is
not a call, such as InterProScan's `polypeptide` rows. Its `isProtein` argument
tells it which side of a HMMER row names the protein, so it reads both hmmscan
and hmmsearch output.

# emm_typing
emm types downloaded from CDC - https://www2a.cdc.gov/ncidod/biotech/strepblast.asp

```
usage: emm_typing.py [-h] -f FASTA [FASTA ...] [--db DB] [-o OUTDIR] [-v]

Group A streptococci emm-typer, version 0.3

optional arguments:
  -h, --help            show this help message and exit
  -f FASTA [FASTA ...], --fasta FASTA [FASTA ...]
                        FASTA file to type.
  --db DB               Database for trimmed emm types. (If using non-
                        default). It must be blastn database. Only provide the
                        file that do not end with ".n*" something (do not use
                        for example /blast_db.sequences.fasta.nhr)
  -o OUTDIR, --outdir OUTDIR
                        Output directory where to write all results.
  -v, --version         Show version and exit.
```

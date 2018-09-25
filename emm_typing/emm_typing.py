#!/usr/bin/env python

# emm-typing

import os
import sys
import re
import argparse
import csv
from pkg_resources import resource_filename
from subprocess import call
from .__init__ import __version__
from Bio.Blast.Applications import NcbiblastnCommandline

EMM_VERSION = __version__

def EmmArgumentParser():
    parser = argparse.ArgumentParser(
        description = 'Group A streptococci emm-typer, version %s' % EMM_VERSION)
    parser.add_argument('-f', '--fasta',
        help='FASTA file to type.',nargs='+', required=True, type=argparse.FileType('r'))
    parser.add_argument('--db',
        help='Database for trimmed emm types. (If using non-default)',
        default=os.path.join(resource_filename(__name__,'data/trimmed_emm_types.tfa')))
    parser.add_argument('-o', '--outdir',
        help='Output directory where to write all results.',
        default='.')
    parser.add_argument('-v', '--version',
        help='Show version and exit.',
        action='version',
        version=EMM_VERSION)
    args = parser.parse_args()

    args.fasta = [os.path.abspath(fasta.name) for fasta in args.fasta]
    args.outdir = os.path.abspath(args.outdir)

    return args

def ChooseBestMatch(lines):
    # Verify EMM <= 124
    # Length = 180
    # Pident = 100.0
    # For multiple bests, report in list
    matches = []
    unvalmatches = []
    for row in lines:
        contig = row[0]
        allele = row[1]
        pident = row[2]
        length = row[3]
        if allele.startswith("EMM"):
            alleleclean = re.match("^EMMG?(\d+)\.\d+$", allele).group(1)
            if not int(alleleclean) <= 124:
                # NOT a verified type
                if pident == "100.000" and length == "180":
                    unvalmatches.append([contig, allele, pident, length])
            else:
                if pident == "100.000" and length == "180":
                    newbest = [contig, allele, pident, length]
                    #matches.insert(0, newbest)
                    if len(matches) == 0:
                        # If no matches so far
                        matches.append(newbest)
                    else:
                        unvalmatches.append(newbest)
        else:
            if pident == "100.000" and length == "180":
                matches.append([contig, allele, pident, length])

    if len(matches) > 0 or len(unvalmatches) > 0:
        return matches, unvalmatches
    else:
        return None, None

def main():
    args = EmmArgumentParser()

    # Test if write permission
    if not os.access('.', os.W_OK):
        sys.exit("Need write permission to current directory")

    if not os.path.isdir("emm_typing"):
        os.mkdir("emm_typing")

    # Sequentially BLAST all fastas against database
    for fasta in args.fasta:
        try:
            isolatename = re.match("^([\w_\-]+)\.(fasta|fa)$", os.path.basename(fasta)).group(1)
        except AttributeError as e:
            print(e)
            sys.exit("Could not understand isolatename: %s. Only a-z, A-Z, numbers, dash and underscore is allowed")
        else:
            assert isolatename

        blastn_cline = NcbiblastnCommandline(query=fasta, db=args.db, perc_identity=100, outfmt=6, max_target_seqs=10,
                                             out=os.path.join(args.outdir,
                                                              'emm_typing',
                                                              '{}_emmresults'.format(isolatename)))
        print(blastn_cline)
        call(blastn_cline(), shell=True)

    # Write all results to communal file (or alternatively, to stdout)
    with open(args.outfile,'w') as communalfile:
        header = ["Isolate", "contig", "emm-type", "pident", "length", "unvalidatedmatches"]
        communal = csv.writer(communalfile, delimiter="\t")
        communal.writerow(header)
        #communal.write("\t".join(header) + "\n")
        for fasta in args.fasta:
            isolatename = re.match("^([\w_\-]+)\.(fasta|fa)$", os.path.basename(fasta)).group(1)
            resfile = os.path.join(args.outdir, 'emm_typing', '{}_emmresults'.format(isolatename))
            with open(resfile, 'rU') as individual:
                ilines = csv.reader(individual, delimiter="\t")
                matches, unvalmatches = ChooseBestMatch(ilines)
                if unvalmatches is not None:
                    unvalidatedmatches = [",".join([u[1] for u in unvalmatches])]
                else:
                    unvalidatedmatches = []
                if matches is not None:
                    #isolate = re.match("^emm_typing/(.*)_emmresults",isolatename).group(1)
                    communal.writerow([isolatename] + matches[0] + unvalidatedmatches)


if __name__ == '__main__':
    main()

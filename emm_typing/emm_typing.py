#!/usr/bin/env python

# emm-typing

import os
import sys
import re
import argparse
import csv
from pkg_resources import resource_filename
from subprocess import call
from __init__ import __version__
from Bio.Blast.Applications import NcbiblastnCommandline

EMM_VERSION = __version__

def EmmArgumentParser():
    parser = argparse.ArgumentParser(
        description = 'Group A streptococci emm-typer, version %s' % EMM_VERSION)
    parser.add_argument('fastafile',
        help='FASTA file to type.',nargs='+', default=None)
    parser.add_argument('--db',
        help='Database for trimmed emm types. (If using non-default)',
        default=os.path.join(resource_filename(__name__,'data/trimmed_emm_types.tfa')))
    parser.add_argument('-o', '--outfile',
        help='Output file to write all results to.',
        default='./emm_typing/emm_results.txt')
    parser.add_argument('-v', '--version',
        help='Show version and exit.',
        action='version',
        version=EMM_VERSION)
    args = parser.parse_args()
    return args

def ChooseBestMatch(lines):
    # Verify EMM <= 124
    # Length = 180
    # Pident = 100.0
    # For multiple bests, report in list
    matches = []
    for row in lines:
        contig = row[0]
        allele = row[1]
        pident = row[2]
        length = row[3]
        if allele.startswith("EMM"):
            alleleclean = re.match("^EMM(\d+)\.\d+$", allele).group(1)
            if not int(alleleclean) <= 124:
                # NOT a verified type
                if pident == "100.000" and length == "180":
                    matches.append([contig, allele, pident, length])
            else:
                if pident == "100.000" and length == "180":
                    newbest = [contig, allele, pident, length]
                    matches.insert(0, newbest)
        else:
            if pident == "100.000" and length == "180":
                matches.append([contig, allele, pident, length])

    if len(matches) > 0:
        return matches[0]
    else:
        return None

def main():
    args = EmmArgumentParser()

    # Test if write permission
    if not os.access('.', os.W_OK):
        sys.exit("Need write permission to current directory")

    if not os.path.isdir("emm_typing"):
        os.mkdir("emm_typing")

    # Sequentially BLAST all fastas against database
    for fasta in args.fastafile:
        isolatename = re.match("^([\w_\-]+)\.(fasta|fa)$", fasta).group(1)
        try:
            assert isolatename
        except:
            sys.exit("Could not understand isolatename: %s. Only a-z, A-Z, numbers, dash and underscore is allowed")
        blastn_cline = NcbiblastnCommandline(query=fasta, db=args.db, perc_identity=100, outfmt=6, max_target_seqs=10, out='./emm_typing/%s_emmresults' % isolatename)
        print(blastn_cline)
        call(blastn_cline(), shell=True)

    # Write all results to communal file (or alternatively, to stdout)
    with open(args.outfile,'w') as communalfile:
        header = ["Isolate", "contig", "emm-type", "pident", "length"]
        communal = csv.writer(communalfile, delimiter="\t")
        communal.writerow(header)
        #communal.write("\t".join(header) + "\n")
        for fasta in args.fastafile:
            isolatename = re.match("^([\w_\-]+)\.(fasta|fa)$", fasta).group(1)
            resfile = "emm_typing/%s_emmresults" % isolatename
            with open(resfile, 'rU') as individual:
                ilines = csv.reader(individual, delimiter="\t")
                bestmatch = ChooseBestMatch(ilines)
                if bestmatch is not None:
                    #isolate = re.match("^emm_typing/(.*)_emmresults",isolatename).group(1)
                    communal.writerow([isolatename] + bestmatch)


if __name__ == '__main__':
    main()
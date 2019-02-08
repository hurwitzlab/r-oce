#!/usr/bin/env python3
"""
Author : Ken Youens-Clark <kyclark@email.arizona.edu>
Date   : 2019-02-01
Purpose: BLAST hits to CMAP visualization
"""

import argparse
import csv
import os
import sys
import pandas as pd
from opedia import db
from subprocess import getstatusoutput


# --------------------------------------------------
def get_args():
    """get command-line arguments"""
    parser = argparse.ArgumentParser(
        description='BLAST hits to CMAP visualization',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument(
        '-q',
        '--query',
        help='Query file for BLAST',
        metavar='str',
        type=str,
        required=True)

    parser.add_argument(
        '-b',
        '--blast_db',
        help='BLAST db',
        metavar='str',
        type=str,
        default='blast')

    parser.add_argument(
        '-p',
        '--blast_program',
        help='BLAST program',
        metavar='str',
        type=str,
        default='blastn')

    parser.add_argument(
        '-i',
        '--perc_identity',
        help='BLAST percent identity',
        metavar='float',
        type=float,
        default=0.)

    parser.add_argument(
        '-c',
        '--qcov_hsp_perc',
        help='BLAST percent query coverage per hsp',
        metavar='float',
        type=float,
        default=0.)

    parser.add_argument(
        '-o',
        '--outdir',
        help='Output directory',
        metavar='str',
        type=str,
        default='out')

    return parser.parse_args()


# --------------------------------------------------
def warn(msg):
    """Print a message to STDERR"""
    print(msg, file=sys.stderr)


# --------------------------------------------------
def die(msg='Something bad happened'):
    """warn() and exit with error"""
    warn(msg)
    sys.exit(1)


# --------------------------------------------------
def main():
    """Make a jazz noise here"""
    args = get_args()
    query = args.query
    blast_db = args.blast_db
    blast_prg = args.blast_program
    out_dir = args.outdir

    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)

    print('Running BLAST')
    hits = run_blast(
        blast_db=blast_db,
        blast_prg=blast_prg,
        perc_identity=args.perc_identity,
        qcov_hsp_perc=args.qcov_hsp_perc,
        query=query,
        out_dir=out_dir)

    print('CMAP query')
    centroids_file = cmap_query(blast_hits=hits, out_dir=out_dir)

    print('Plotting')
    plot(centroids_file=centroids_file, out_dir=out_dir)

    print('Done')


# --------------------------------------------------
def plot(centroids_file, out_dir):
    """Given CMAP location data, plot distribution"""

    cwd = os.path.abspath(os.path.dirname(sys.argv[0]))
    plot = os.path.join(cwd, 'plot.r')

    if os.path.isfile(plot):
        rv, out = getstatusoutput('{} {}'.format(plot, centroids_file))
        print(out)

        if rv == 0:
            return 1
        else:
            die('Error plotting ({}):\n{}\n'.format(rv, out))
    else:
        die('Cannot find "{}"'.format(plot))


# --------------------------------------------------
def cmap_query(blast_hits, out_dir):
    """Given BLAST hits, query CMAP for location"""

    out_file = os.path.join(out_dir, 'oce-input.csv')
    out_fh = open(out_file, 'wt')
    out_fh.write(','.join([
        'latitude', 'longitude', 'depth', 'Relative_Abundance', 'temperature',
        'salinity'
    ]) + '\n')

    seen = set()
    with open(blast_hits) as csvfile:
        flds = [
            'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
            'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'
        ]
        reader = csv.DictReader(csvfile, fieldnames=flds, delimiter='\t')

        query = "SELECT * FROM tblesv WHERE centroid='{}'"

        for hit_num, row in enumerate(reader):
            seq_id = row['sseqid']
            print('{:5}: {}'.format(hit_num + 1, seq_id))

            if seq_id in seen:
                continue
            seen.add(seq_id)

            df = db.dbFetch(query.format(seq_id))
            for i, row in df.iterrows():
                out_fh.write(','.join(
                    map(str, [
                        row['lat'], row['lon'], row['depth'],
                        row['relative_abundance'], row['relative_abundance'],
                        row['esv_salinity']
                    ])) + '\n')

    out_fh.close()
    return out_file


# --------------------------------------------------
def run_blast(blast_db, blast_prg, query, perc_identity, qcov_hsp_perc,
              out_dir):
    """Given user query and params, run BLAST"""

    hits_file = os.path.join(out_dir, 'hits.tab')

    cmd = '{} -query {} -db {} -out {} -outfmt 6'.format(
        blast_prg, query, blast_db, hits_file)

    if perc_identity > 0.:
        cmd += ' -perc_identity {}'.format(perc_identity)

    if qcov_hsp_perc > 0.:
        cmd += ' -qcov_hsp_perc {}'.format(qcov_hsp_perc)

    (rv, blast_out) = getstatusoutput(cmd)

    if rv == 0:
        with open(hits_file) as fh:
            if len(fh.readlines()) > 0:
                return hits_file
            else:
                die('No hits from BLAST')
    else:
        die('Failed to run BLAST ({}):\n{}'.format(rv, blast_out))


# --------------------------------------------------
if __name__ == '__main__':
    main()

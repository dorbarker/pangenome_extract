import argparse
from pathlib import Path
import pandas as pd
from Bio import SeqIO
from concurrent.futures import ProcessPoolExecutor
import functools


def arguments():

    parser = argparse.ArgumentParser()

    parser.add_argument('-o', '--output',
                        required=True,
                        type=Path)

    parser.add_argument('-i', '--input',
                        required=True,
                        type=Path)

    parser.add_argument('--gffs',
                        required=True,
                        type=Path)

    parser.add_argument('--genomes',
                        required=True,
                        type=Path)

    parser.add_argument('--cpus',
                        type=int,
                        default=1)

    return parser.parse_args()


def main():

    args = arguments()

    gene_presence_absence = pd.read_csv(args.input)

    process = functools.partial(process_locus,
                                gffs=args.gffs,
                                genomes=args.genomes,
                                output=args.output)

    with ProcessPoolExecutor(max_workers=args.cpus) as ppe:

        rows = (row for i, row in gene_presence_absence.iterrows())

        ppe.map(process, rows)



def process_locus(row, gffs, genomes, output):

    genome, locus = locus_tag(row)

    gff_path = generate_path(gffs, genome, '.gff')
    genome_path = generate_path(genomes, genome, '.fasta')

    locus_params = locus_parameters(gff_path, locus)

    sequence = extract_sequence(genome_path, locus_params)

    write_record(genome, locus, sequence, output)


def locus_tag(row):

    loci = row[14:]

    genomes = loci.index.values

    not_na = ~loci.isna()

    genome, *rest = genomes[not_na]
    locus, *rest = loci[not_na]

    return genome, locus


def generate_path(parent, genome, suffix):

    return (parent / genome).with_suffix(suffix)


def locus_parameters(gff, locus):

    with gff.open('r') as f:

        for line in f:

            if locus not in line:
                continue

            contig, _, _, start, stop, _, strand, _, _ = line.split('\t')

    values = {'contig': contig,
              'start': int(start),
              'stop': int(stop),
              'strand': strand}

    return values

def extract_sequence(fasta_path, values):

    with fasta_path.open('r') as f:

        for record in SeqIO.parse(f, 'fasta'):

            if record.name != values['contig']:
                continue

            sequence = record.seq[values['start'] - 1 : values['stop']]

            if values['strand'] != '+':

                out = sequence.reverse_complement()

            else:

                out = sequence

            return out


def write_record(genome, locus, sequence, outdir):

    name = f'{genome}_{locus}'

    out = generate_path(outdir, name, '.fasta')

    record = f'>{name}\n{sequence}\n'

    out.write_text(record)

if __name__ == '__main__':
    main()


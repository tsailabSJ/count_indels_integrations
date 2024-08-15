from __future__ import print_function
import re
import gzip
import argparse
import pyfaidx
import pysam
from collections import Counter, defaultdict
from skbio.alignment import StripedSmithWaterman

def get_sequence(reference_genome, chromosome, start, end, strand="+"):
    if strand == "+":
        seq = reference_genome[chromosome][start:end]
    elif strand == "-":
        seq = reference_genome[chromosome][start:end].reverse.complement
    return str(seq)

def alignSequences(reference, read):
    query = StripedSmithWaterman(reference, gap_open_penalty=10)
    alignment = query(read)
    return alignment

def readBed(bedfile, reference_genome):
    # Input is 0-based, half-open BED file (includes start but not end)
    # Sequence is always returned where 200 is the breakpoint and in the forward orientation
    BED_dict = {}
    with open(bedfile, 'r') as f:
        for line in f:
            if line[0] != '#':
                chr, start, end, name, CHANGEseq_reads, strand = line.rstrip('\n').split('\t')
                start = int(start)
                end = int(end)
                if strand == "-" :
                    sequence = get_sequence(reference_genome, chr, start + 6 - 200, start + 6 + 200)
                else:
                    sequence = get_sequence(reference_genome, chr, end - 6 - 200 , end - 6 + 200)
                BED_dict[name] =  { 'chr' : chr,
                                    'start' : start,
                                    'end' : end,
                                    'strand' : strand,
                                    'seq' : sequence
                                    }
    return BED_dict

def readBAM(bam, BED_dict, site, log, filterlog):
        integrations = 0
        indels = 0
        total = 0
        breakpoint_pos = 200
        window = 5

        # alignment_score_dict = Counter()
        # cigar_dict = Counter()
        fwd_query= StripedSmithWaterman('GTTTAATTGAGTTGTCATATGTTAATAACGGTAT', gap_open_penalty=10)
        rev_query = StripedSmithWaterman('ATACCGTTATTAACATATGACAACTCAATTAAAC', gap_open_penalty=10)
        read_query = StripedSmithWaterman(BED_dict[site]['seq'], gap_open_penalty=10, gap_extend_penalty=1, zero_index=True)

        # with open(logfile, 'a') as o2:

        for read in bam.fetch(BED_dict[site]['chr'], BED_dict[site]['start'], BED_dict[site]['end']):

            # Align read
            read_aln = read_query(read.query_sequence)
            m = re.findall(r'(\d+)([A-Z]{1})', read_aln.cigar)

            current_pos = read_aln.query_begin
            match_min_pos = 1000 # set to number above possible current positions
            match_max_pos = -1 # set to number below current max)
            match_window = 30
            indel_detected = False

            for length, type in m:
                if type == "M":
                    match_min_pos = min(current_pos, match_min_pos) # Get minimum
                    current_pos += int(length) # Move position counter
                    match_max_pos = max(current_pos, match_max_pos) # Get maximum
                elif type == "I":
                    if breakpoint_pos + window >= current_pos and breakpoint_pos - window <= current_pos + int(length):
                        indel_detected = True
                    current_pos += int(length)
                elif type == "D":
                    if breakpoint_pos + window >= current_pos and breakpoint_pos - window <= current_pos:
                        indel_detected = True

            # Align to forward and reverse oligo
            fwd_aln = fwd_query(read.query_sequence)
            rev_aln = rev_query(read.query_sequence)

            # Tabulate integrations
            if (fwd_aln.optimal_alignment_score > 30 or rev_aln.optimal_alignment_score > 30):
                if breakpoint_pos - match_window >= match_min_pos and breakpoint_pos + match_window <= match_max_pos:
                    integrations += 1
                    print(site, 'fwd', 'query', fwd_aln.aligned_query_sequence, sep="\t", file=log)
                    print(site, 'fwd', 'target', fwd_aln.aligned_target_sequence, sep="\t", file=log)
                    print(site, 'rev', 'query', rev_aln.aligned_query_sequence, sep="\t", file=log)
                    print(site, 'rev', 'target', rev_aln.aligned_target_sequence, sep="\t", file=log)
                    print(site, 'read', 'query', read_aln.aligned_query_sequence, sep="\t", file=log)
                    print(site, 'read', 'target', read_aln.aligned_target_sequence, sep="\t", file=log)
                    print(site, 'read', 'cigar', read_aln.cigar, sep="\t", file=log)
                    print(site, 'read', 'bam_sequence', read.query_sequence, sep="\t", file=log)
                else:
                    print(site, 'fwd', 'query', fwd_aln.aligned_query_sequence, file=filterlog)
                    print(site, 'fwd', 'target', fwd_aln.aligned_target_sequence, file=filterlog)
                    print(site, 'rev', 'query', rev_aln.aligned_query_sequence, file=filterlog)
                    print(site, 'rev', 'target', rev_aln.aligned_target_sequence, file=filterlog)
                    print(site, 'read', 'query', read_aln.aligned_query_sequence, file=filterlog)
                    print(site, 'read', 'target', read_aln.aligned_target_sequence, file=filterlog)
                    print(site, 'read', 'cigar', read_aln.cigar, file=filterlog)
                    print(site, 'read', 'query_name', read.query_name, file=filterlog)
                    print(site, 'read', 'bam_sequence', read.query_sequence, file=filterlog)
            if indel_detected:
                indels += 1
            total += 1

        return (indels, integrations, total)


def main():
    parser = argparse.ArgumentParser(description='Count sequence positions')
    parser.add_argument('--bed', help='BED filename', required=True)
    parser.add_argument('--ref', help='Reference filename', required=True)
    parser.add_argument('--bam', help='BAM filename', required=True)
    parser.add_argument('--out', help='Output filename', required=True)
    parser.add_argument('--site', help='Site (optional)')
    args = parser.parse_args()

    # Load reference sequence
    print("Loading reference sequence...")
    reference_genome = pyfaidx.Fasta(args.ref)

    # Load BED file
    print("Loading BED file and get sequences...")
    BED_dict = readBed(args.bed, reference_genome)

    # Allow analysis of subset if specified
    if(args.site in BED_dict.keys()):
        BED_dict = {args.site: BED_dict[args.site]}
    else:
        print("Site not found. Processing entire BED file")

    # Load BAM file
    print("Loading BAM file...")
    bam = pysam.AlignmentFile(args.bam, "rb")
    bam_file_name = args.bam.split("/")[-1]
    with open(args.out, 'w') as o1, open(args.out + ".log", 'w') as o2, open(args.out + ".filtered.log", 'w') as o3:
        for site in BED_dict:
            indels, integrations, total = readBAM(bam, BED_dict, site, o2, o3)
            print(site, BED_dict[site]['chr'], BED_dict[site]['start'], BED_dict[site]['end'], BED_dict[site]['seq'],
                  indels, integrations, total, sep="\t")
            print(bam_file_name,site, BED_dict[site]['chr'], BED_dict[site]['start'], BED_dict[site]['end'], BED_dict[site]['seq'],
                  indels, integrations, total, sep="\t", file=o1)


if __name__ == "__main__":
    main()
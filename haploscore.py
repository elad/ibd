"""
haploscore.py
Copyright (C) 2013 23andMe, Inc.

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

If you have questions, please contact 23andMe, Inc. at
1390 Shorebird Way, Mountain View, CA 94043.

Authors: 
Eric Y. Durand  <edurand@23andme.com>
Cory Y. McLean  <cmclean@23andme.com>
"""

import numpy as np
import sys
from operator import attrgetter


class PlinkReader(object): 
    """
    Class for fast mapping from rsid to SNP and iid to Individual
    """
    def __init__(self, pedfile, mapfile): 
        self._iid_to_individual = {}
        self._rsid_to_snp = {}
        self._read_map(mapfile)
        self._read_ped(pedfile)
        
    def _read_map(self, filename):
        chrom_rsids = {}
        # Create dictionary from rsid to SNP object
        with open(filename, "r") as ifs:
            for mapix, line in enumerate(ifs): 
                chrom, rsid, dist, pos = line.strip().split()
                pos = int(pos)
                assert rsid not in self._rsid_to_snp, "Duplicate SNP found: %s\n" % rsid
                if chrom not in chrom_rsids:
                    chrom_rsids[chrom] = []
                snp = SNP(rsid, pos, mapix)
                chrom_rsids[chrom].append(snp)
                self._rsid_to_snp[rsid] = snp
        self._num_snps = len(self._rsid_to_snp)
        # Create mask from index in genome order to index in map file
        # and add genome index information to each SNP
        genome_order = []
        genomeix = 0
        for chrom in sorted(chrom_rsids.keys()):
            for snp in sorted(chrom_rsids[chrom], key=attrgetter('pos')):
                genome_order.append(snp.map_index)
                snp.genome_index = genomeix
                genomeix += 1
        self._genome_order_mask = np.array(genome_order)

    def _read_ped(self, filename):
        with open(filename, "r") as ifs: 
            for line in ifs: 
                tokens = line.strip().split()
                fid, iid, pid, mid, sex, phen = tokens[:6]
                assert iid not in self._iid_to_individual, "Duplicate individual found: %s\n" % iid
                haplotypes = np.array(tokens[6:], dtype=int)
                assert len(haplotypes) == 2 * self._num_snps, "Invalid haplotype entry: %s\n" % iid
                map_ordered_hap1 = haplotypes[::2]
                map_ordered_hap2 = haplotypes[1::2]
                genome_ordered_hap1 = map_ordered_hap1[self._genome_order_mask]
                genome_ordered_hap2 = map_ordered_hap2[self._genome_order_mask]
                ind = Individual(iid, genome_ordered_hap1, genome_ordered_hap2)
                self._iid_to_individual[iid] = ind

    def get_individual(self, iid): 
        return self._iid_to_individual[iid]

    def get_snp(self, rsid): 
        return self._rsid_to_snp[rsid]


class IBDSegment(object): 
    """
    Class representing a GERMLINE IBD segment
    """
    def __init__(self, iid1, iid2, rstart, rend): 
        self.iid1 = iid1
        self.iid2 = iid2 
        self.rstart = rstart  # rsid of first SNP in IBD segment
        self.rend = rend      # rsid of last SNP in IBD segment

    @classmethod
    def from_line(cls, line):
        tokens = line.strip().split()
        assert len(tokens) == 15, "Misformatted IBD segment: %s" % line
        # Segment created from iid1, iid2, segment start snp, segment end snp
        # See http://www1.cs.columbia.edu/~gusev/germline/
        segment = cls(tokens[1], tokens[3], tokens[7], tokens[8])
        return segment


class Individual(object): 
    """
    A diploid individual is represented by (iid, hap1, hap2).
    Haplotypes are the same length and are sorted in genome order.
    """
    def __init__(self, iid, hap1, hap2): 
        self.iid = iid
        self.haplotypes = np.array([hap1, hap2])


class SNP(object): 
    """
    A SNP is represented by (rsid, pos, map_index)
    """
    def __init__(self, rsid, pos, map_index):
        self.rsid = rsid
        self.pos = pos
        self.map_index = map_index


def compute_haploscore(segment, plinkdata, geno_penalty, switch_penalty):
    start, end = [plinkdata.get_snp(rsid).genome_index
                  for rsid in [segment.rstart, segment.rend]]
    individual1 = plinkdata.get_individual(segment.iid1).haplotypes[:, start:(end+1)]
    individual2 = plinkdata.get_individual(segment.iid2).haplotypes[:, start:(end+1)]

    def _get_current_genotype_penalty(i):
        _penalty = []
        for h1 in xrange(2):
            m1 = individual1[h1, i]
            for h2 in xrange(2): 
                m2 = individual2[h2, i]
                _penalty.append((m1 != m2)*geno_penalty)
        return np.array(_penalty, dtype=np.float)

    prevscore = _get_current_genotype_penalty(0)
    # -- Calculating true number of snps as the field reported by GERMLINE is occasionally incorrect.
    nsnp = end - start + 1
    for i in xrange(1, nsnp): 
        nextswitch = prevscore + switch_penalty
        nextscore  = _get_current_genotype_penalty(i) + nextswitch.min(1)
        prevscore  = nextscore

    return min(prevscore)/float(nsnp)


def main(args):
    if args.verbose:
        sys.stderr.write("Loading PLINK data ...\n")
    plinkdata = PlinkReader(args.ped, args.map)
    genotype_penalty = 1./args.genotype_error
    switch_penalty   = 1./args.switch_error
    switch_pen_mat   = np.array([[0,1,1,2],
                                 [1,0,2,1],
                                 [1,2,0,1],
                                 [2,1,1,0]], dtype=np.float) * switch_penalty

    if args.verbose:
        sys.stderr.write("Scoring segments ...\n")
    with open(args.ibd, "r") as ifs, open(args.out, "w") as ofs: 
        for i, line in enumerate(ifs): 
            if args.verbose and (i+1) % 1000 == 0:
                sys.stderr.write("Processed [%d] segments\r" % (i+1))
            segment = IBDSegment.from_line(line)
            score = compute_haploscore(segment, plinkdata, genotype_penalty, switch_pen_mat)
            ofs.write("%s\t%s\n" % (line.strip(), "%.3f" % score))


if __name__ == '__main__':
    import argparse
    usage = "%(prog)s [options] match ped_file map_file out_file"
    description = "Compute haploscore on provided ibd segments.\
                   The score is added as an additional column to the segment file."
    
    parser = argparse.ArgumentParser(usage=usage, description=description, 
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("ibd")
    parser.add_argument("ped")
    parser.add_argument("map")
    parser.add_argument("out")
    parser.add_argument("--genotype_error", type=float, default=0.0075, 
                        help="Genotyping error rate (#errors per marker)")
    parser.add_argument("--switch_error", type=float, default=0.003, 
                        help="Switch error rate (#errors per marker)")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()
    main(args)


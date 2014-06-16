"""
haploscore.py
Copyright (C) 2014 23andMe, Inc.

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
import os.path
import sys
from haplotypes import PlinkMap, PlinkPED
from ibd import GermlineIBDSegment
from threshold import ThresholdMatrix, DEFAULT_THRESHOLD_MATRIX


def get_score_thresholder(thresholdfn, filterval):
    """
    Return a Thresholds object that will filter IBD segments based on their
    genetic length and HaploScore to the desired filtering threshold.
    """
    if thresholdfn is None:
        threshmat = DEFAULT_THRESHOLD_MATRIX
    else:
        threshmat = ThresholdMatrix.from_file(thresholdfn)
    return threshmat.threshold(filterval)


def compute_haploscore(segment, snpdata, inddata, geno_penalty, switch_penalty):
    """
    Return the HaploScore of the given IBD segment based on the haplotypes of 
    the individuals and the expected genotype and switch error penalties.
    """
    start, end = [snpdata.get_snp(rsid).index
                  for rsid in [segment.rstart, segment.rend]]
    individual1 = inddata.get_individual(segment.iid1)
    individual2 = inddata.get_individual(segment.iid2)
    haps1 = individual1.haplotypes[:, start:(end+1)]
    haps2 = individual2.haplotypes[:, start:(end+1)]

    def _get_current_genotype_penalty(i):
        _penalty = []
        for h1 in xrange(len(haps1)):
            m1 = haps1[h1, i]
            for h2 in xrange(len(haps2)):
                m2 = haps2[h2, i]
                _penalty.append((m1 != m2) * geno_penalty)
        return np.array(_penalty, dtype=np.float)

    prevscore = _get_current_genotype_penalty(0)
    # Calculate number of snps as the field 
    # reported by GERMLINE is occasionally incorrect.
    nsnp = end - start + 1
    for i in xrange(1, nsnp): 
        nextswitch = prevscore + switch_penalty
        nextscore  = _get_current_genotype_penalty(i) + nextswitch.min(1)
        prevscore  = nextscore

    return min(prevscore) / float(nsnp)


def main(args):
    if args.verbose:
        sys.stderr.write("Loading PLINK data ...\n")
    unordered_snpdata = PlinkMap.from_file(args.map)
    unordered_inddata = PlinkPED.from_file(args.ped)
    genome_order = unordered_snpdata.genome_order_mask()
    snpdata = unordered_snpdata.ordered(genome_order)
    inddata = unordered_inddata.ordered(genome_order)
    genotype_penalty = 1./args.genotype_error
    switch_penalty   = 1./args.switch_error
    switch_pen_mat   = np.array([[0,1,1,2],
                                 [1,0,2,1],
                                 [1,2,0,1],
                                 [2,1,1,0]], dtype=np.float) * switch_penalty
    minscore = float("infinity")
    nofilter = args.filter is None
    if not nofilter:
        thresholder = get_score_thresholder(args.threshold_file, args.filter)
    
    if args.verbose:
        sys.stderr.write("Scoring segments ...\n")
    kept = 0
    with open(args.ibd, "r") as ifs, open(args.out, "w") as ofs: 
        for i, line in enumerate(ifs): 
            segment = GermlineIBDSegment.from_line(line)
            score = compute_haploscore(segment, 
                                       snpdata, 
                                       inddata,
                                       genotype_penalty,
                                       switch_pen_mat)
            minscore = min(minscore, score)
            if (nofilter or 
                thresholder.passes_threshold(segment.genetic_length, score)):
                ofs.write("%s\t%s\n" % (line.strip(), "%.3f" % score))
                kept += 1
            if args.verbose and (i+1) % 1000 == 0:
                sys.stderr.write("Processed [%d] segments\r" % (i+1))
    if args.verbose:
        sys.stderr.write("Kept %d of %d total input segments.\n" % (kept, i+1))

    # Warn users if input produces worse results than expected
    if minscore > 5:
        msg = ("### Warning: Minimum HaploScore seen in %d segments=%f. ###\n"
               "### Ensure your PED file contains phased genotype data. ###\n")
        sys.stderr.write(msg % (i+1, minscore))


if __name__ == '__main__':
    import argparse
    usage = "%(prog)s [options] match_file ped_file map_file out_file"
    description = ("Compute haploscore on provided ibd segments.\n"
                   "The score is added as an additional column to the segment "
                   "file.")
    
    parser = argparse.ArgumentParser(usage=usage, description=description) 
    parser.add_argument("ibd")
    parser.add_argument("ped")
    parser.add_argument("map")
    parser.add_argument("out")
    calcgrp = parser.add_argument_group("HaploScore calculation options")
    calcgrp.add_argument("--genotype_error", type=float, default=0.0075, 
                         help="Genotyping error rate (#errors per marker) "
                              "(default: %(default)s)")
    calcgrp.add_argument("--switch_error", type=float, default=0.003,
                         help="Switch error rate (#errors per marker) "
                              "(default: %(default)s)")
    filtgrp = parser.add_argument_group("Segment filtering options",
                                        "(Only use if you want to "
                                        "simultaneously calculate HaploScore "
                                        "and filter out poorly-scoring "
                                        "reported IBD segments)")
    filtgrp.add_argument("--filter", type=float, default=None,
                         help="Filter out segments not passing this mean "
                              "overlap threshold (default: all segments kept)")
    filtgrp.add_argument("--threshold_file", default=None,
                         help="Use this file of HaploScore thresholds "
                              "(default: use thresholds from original "
                              "HaploScore publication)")
    parser.add_argument("-v", "--verbose", default=False, action="store_true")
    args = parser.parse_args()
    
    # Validate arguments
    for argname, val in [('Genotype error', args.genotype_error),
                         ('Switch error', args.switch_error),
                         ('Filter threshold', args.filter)]:
        if val is not None and not (0 < val <= 1):
            raise ValueError("%s must be in range (0, 1]: %f" % (argname, val))
    for fn in [args.ibd, args.ped, args.map]:
        if not os.path.isfile(fn):
            raise ValueError("%s does not exist." % fn)
    if (args.threshold_file is not None and
        not os.path.isfile(args.threshold_file)):
        raise ValueError("%s does not exist." % args.threshold_file)
    
    main(args)


"""
haplotypes.py
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
from operator import attrgetter


class SNP(object): 
    """
    A SNP is represented by (rsid, chrom, pos, index).

    rsid is the unique name of the SNP.
    chrom and pos indicate its genomic location.
    index is the SNP's location within the haplotypes array of an Individual --
      ind.haplotypes[:,index] give the alleles for individual ind for this SNP.
    """
    def __init__(self, rsid, chrom, pos, index):
        self.rsid = rsid
        self.chrom = chrom
        self.pos = pos
        self.index = index


class Individual(object): 
    """
    A diploid individual is represented by (iid, hap1, hap2).
    """
    def __init__(self, iid, hap1, hap2):
        assert len(hap1) == len(hap2)
        self.iid = iid
        self.haplotypes = np.array([hap1, hap2])
        self.haplen = len(hap1)


class SNPContainer(object):
    
    """
    Class to hold fast mapping from rsid to SNP.
    Should not be instantiated directly, but have subclasses use the from_file
    classmethod instead.
    """

    def __init__(self, snps):
        self._snps = snps
        self._num_snps = len(self._snps)
        self._rsid_to_snp = dict((snp.rsid, snp) for snp in snps)
    
    @classmethod
    def from_file(cls, fn):
        raise NotImplementedError()

    def get_snp(self, rsid):
        return self._rsid_to_snp[rsid]

    def genome_order_mask(self):
        """
        Return a numpy mask to order all SNPs in genome order.
        """
        chrom_rsids = {}
        for fileix, snp in enumerate(self._snps):
            if snp.chrom not in chrom_rsids:
                chrom_rsids[snp.chrom] = []
            chrom_rsids[snp.chrom].append(snp)
            assert snp.index == fileix
        genome_order = []
        genomeix = 0
        for chrom in sorted(chrom_rsids.keys()):
            for snp in sorted(chrom_rsids[chrom], key=attrgetter('pos')):
                genome_order.append(snp.index)
                genomeix += 1
        return np.array(genome_order)

    def ordered(self, mask):
        """
        Return a new object with the same SNPs but sorted in the order
        specified by the mask.
        """
        assert len(mask) == self._num_snps
        assert sorted(mask) == range(self._num_snps)
        ordered_snps = []
        for ix, snp in enumerate(self._snps[mask]):
            newsnp = SNP(snp.rsid, snp.chrom, snp.pos, ix)
            ordered_snps.append(newsnp)
        return self.__class__(np.array(ordered_snps))


class IndividualContainer(object):

    """
    Class to hold fast mapping from iid to Individual.
    Should not be instantiated directly, but have subclasses use the from_file
    classmethod instead.
    """

    def __init__(self, inds):
        self._iid_to_individual = dict((ind.iid, ind) for ind in inds)
        self._num_individuals = len(self._iid_to_individual)
        haplens = list(set([ind.haplen for ind in inds]))
        assert len(haplens) == 1
        self.haplen = haplens[0]

    @classmethod
    def from_file(cls, fn):
        raise NotImplementedError()

    def get_individual(self, iid):
        return self._iid_to_individual[iid]

    def individuals(self):
        return self._iid_to_individual.values()

    def ordered(self, mask):
        """
        Return a new container with all Individuals' haplotypes reordered as 
        specified by the mask.
        """
        assert len(mask) == self.haplen
        assert sorted(mask) == range(self.haplen)
        ordered_inds = []
        for ind in self.individuals():
            h1, h2 = ind.haplotypes
            ordered_ind = Individual(ind.iid, h1[mask], h2[mask])
            ordered_inds.append(ordered_ind)
        return self.__class__(ordered_inds)        


class PlinkMap(SNPContainer):

    """
    Concrete class to load SNPs from a PLINK .map file.
    """
    @classmethod
    def from_file(cls, fn):
        snps = []
        seen = set()
        with open(fn, "r") as ifs:
            for fileix, line in enumerate(ifs):
                chrom, rsid, gendist, pos = line.strip().split()
                pos = int(pos)
                assert rsid not in seen, "Duplicate SNP found: %s\n" % rsid
                snp = SNP(rsid, chrom, pos, fileix)
                snps.append(snp)
                seen.add(rsid)
        return cls(np.array(snps))


class PlinkPED(IndividualContainer):

    """
    Concrete class to load Individuals from a PLINK .ped file.
    """
    @classmethod
    def from_file(cls, fn):
        inds = []
        seen = set()
        with open(fn, "r") as ifs: 
            for line in ifs: 
                tokens = line.strip().split()
                fid, iid, pid, mid, sex, phen = tokens[:6]
                assert iid not in seen, "Duplicate individual found: %s\n" % iid
                haplotypes = np.array(tokens[6:], dtype=int)
                hap1 = haplotypes[::2]
                hap2 = haplotypes[1::2]
                ind = Individual(iid, hap1, hap2)
                inds.append(ind)
                seen.add(iid)
        return cls(inds)


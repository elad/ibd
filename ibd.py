"""
ibd.py
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

class IBDSegment(object): 
    """
    Class representing an IBD segment.
    The fields needed to be able to calculate HaploScore and filter segments
    are the individuals between whom the segment is specified,
    the start and end SNPs of the segment,
    and the genetic length of the segment.
    """
    def __init__(self, iid1, iid2, rstart, rend, genlen): 
        self.iid1 = iid1
        self.iid2 = iid2 
        self.rstart = rstart  # rsid of first SNP in IBD segment
        self.rend = rend      # rsid of last SNP in IBD segment
        self.genetic_length = genlen

    @classmethod
    def from_line(cls, line):
        raise NotImplementedError()


class GermlineIBDSegment(IBDSegment):

    @classmethod
    def from_line(cls, line):
        # GERMLINE reports each IBD segment match on a single 15-token line.
        # See http://www1.cs.columbia.edu/~gusev/germline/ for details.
        tokens = line.strip().split()
        assert len(tokens) == 15, "Misformatted GERMLINE IBD segment: %s" % line
        iid1, iid2 = tokens[1], tokens[3]
        startrsid, endrsid = tokens[7], tokens[8]
        genetic_length = float(tokens[10])
        segment = cls(iid1, iid2, startrsid, endrsid, genetic_length)
        return segment


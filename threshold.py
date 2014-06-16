"""
threshold.py
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


class Thresholds(object):
    """
    This class holds the HaploScore threshold values for segments of various
    genetic lengths at a single particular mean overlap.
    """
    def __init__(self, genlens, thresholds):
        assert len(genlens) == len(thresholds)
        assert len(genlens) >= 2
        self._lengths = genlens
        self._thresholds = thresholds

    def passes_threshold(self, genlen, score):
        """
        Return True iff a segment with genetic length genlen
        and HaploScore score passes the threshold.

        Note: Any segments shorter than the minimum genetic length
              are assumed to be false and any segments longer than
              the maximum genetic length are assumed to be true.
        """
        if genlen < self._lengths[0]:
            return False
        stepsize = self._lengths[-1] - self._lengths[-2]
        if genlen > (self._lengths[-1] + stepsize):
            return True
        ix = np.searchsorted(self._lengths, genlen, side='right') - 1
        return score <= self._thresholds[ix]


class ThresholdMatrix(object):
    """
    Class that holds threshold vectors for all genetic lengths and desired mean
    overlap threshold values.
    """
    def __init__(self, genlens, threshold_keys, threshold_values):
        # Ensure lengths are in sorted order, step sizes are equal,
        # and inputs are shaped correctly
        stepsizes = genlens[1:] - genlens[:-1]
        assert np.all(stepsizes > 0)
        assert np.all(np.abs(stepsizes - stepsizes[0]) < 1e-10)
        assert len(threshold_keys) > 0 and threshold_keys[0] == 1
        assert np.all((threshold_keys[1:] - threshold_keys[:-1]) < 0)
        assert len(threshold_values.shape) == 2
        tvx, tvy = threshold_values.shape
        assert len(genlens) == tvx
        assert len(threshold_keys) == tvy
        self._lengths = genlens
        self._threshold_keys = threshold_keys
        self._threshold_values = threshold_values

    @classmethod
    def from_file(cls, fn):
        """
        Create a threshold matrix from a file.
        Assumes all mean overlap bins are equally sized.
        """
        matrix = np.loadtxt(fn)
        genlens = matrix[:,0]
        threshold_values = matrix[:,3:]
        step = 1 / float(threshold_values.shape[1])
        threshold_keys = np.arange(1, 0, -step)
        return cls(genlens, threshold_keys, threshold_values)

    def threshold(self, mean_overlap):
        """
        Return a Thresholds object corresponding to the requested mean overlap.
        """
        assert 0 < mean_overlap <= 1
        if mean_overlap == 1:
            threshix = 0
        else:
            threshix = np.flatnonzero(self._threshold_keys > mean_overlap)[-1]
        return Thresholds(self._lengths, self._threshold_values[:,threshix])


# Threshold matrix calculated on the 23andMe cohort on chromosome 21
# as part of the HaploScore publication:
# Durand EY, Eriksson N, and McLean CY. Reducing pervasive false positive
# identical-by-descent segments detected by large-scale pedigree analysis.
# Molecular Biology and Evolution, 2014.
# http://mbe.oxfordjournals.org/lookup/pmid?view=long&pmid=24784137
#
# Thresholds were built in the following way:
#   Bins of length 0.1 cM from [2.0, 2.1), [2.1, 2.2), ..., [10.0, 10.1) were
#   created. All reported child-other IBD segments on chromosome 21 in the
#   23andMe cohort were placed into their appropriate 0.1 cM bin. For each slice
#   of 1/100th of the segments in increasing order by HaploScore, the mean 
#   parent-other overlap in that slice of segments was calculated. That overlap
#   threshold value was then updated with that maximum segment score. Finally,
#   monotonicity of scores was enforced. See the below example code for a more
#   precise specification of the threshold generation.
#
#### Inefficient code for a single length bin ####
#
# def thresh_ix(ovlp, num_thresholds):
#     return max(0, num_thresholds - 1 - int(floor(ovlp*num_thresholds)))
# num_thresholds = 100
# threshold_values = np.repeat(-1, num_thresholds)
# ordered_segs = sorted(bin_segs, key=attrgetter('score'))
# num_segs = len(ordered_segs)
# unique_scores = sorted(set([seg.score for seg in ordered_segs]))
# minscore = -1
# for score in unique_scores:
#     maxscore = score
#     segs = [seg for seg in ordered_segs if minscore < seg.score <= maxscore]
#     if ((len(segs) < num_segs / float(num_thresholds)) and
#         (score < unique_scores[-1])):
#         continue
#     overlap_threshold = np.mean([seg.parent_other_overlap for seg in segs])
#     ix = thresh_ix(overlap_threshold, num_thresholds)
#     threshold_values[ix] = maxscore
#     minscore = maxscore
# # Ensure monotonicity
# for i in range(1, len(threshold_values)):
#     threshold_values[i] = max(threshold_values[i], threshold_values[i-1])

_fn = os.path.join(os.path.dirname(__file__), 'data', 'chr21.scorethresh.txt')
DEFAULT_THRESHOLD_MATRIX = ThresholdMatrix.from_file(_fn)


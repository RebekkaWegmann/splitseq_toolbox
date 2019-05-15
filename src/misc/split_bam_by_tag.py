#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This code splits a bam file by a user-specified tag. Tag values can either be exact macthes, or matches with a certain number of allowed mismatches. In the latter case, all possible tag values must have the same length.
It produces two output files: 
- file one: all of the sequences that have the specified tag and its value matches the specified value.
- file two: all of the sequences where the tag is missing or not matching the specified value. To only keep sequences where the tag is present, set the --discard_missing_tags flag.

Copyright: Rebekka Wegmann, Snijderlab, ETH Zurich, 2019

# MIT License
#
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
"""

#%% libraries
from __future__ import division 
import sys, re, os
import pysam
import numpy as np
from argparse import ArgumentParser

#%% Functions
def hamming(s1, s2):
    """Calculate the Hamming distance between two strings"""
    assert len(s1) == len(s2)
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))

#%%Setup input parser
parser = ArgumentParser()
parser.add_argument("tag_id", type=str, help="The tag ID. Needs to be a valid BAM tag, such as one of the default tags (e.g. RG, MQ, ...) or a user-defined tag (e.g. XA,XB,...).")
parser.add_argument("tag_value", type=str, help="The tag value to filter by. E.g. if tag_id = RG and tag_value = A, all reads where the value of the RG tag equals A will be put in one file, and all reads where RG is missing or has a value other than A will be put in asecond file.")
parser.add_argument("-d" "--max_edit_dist", action="store", dest="max_edit_dist", default=0, help="Maximum allowed edit distance between desired and actual tag value", type=int)
parser.add_argument("-i" "--input_bam", action="store", dest="input_bam", default="-", help="Specify the input bam file. Defaults to stdin.")
parser.add_argument("-o1" "--output_bam_1", action="store", dest="output_bam_1", default="with_specified_tag.bam", help="Specify the name of the first output bam file. Defaults to with_specified_tag.bam.")
parser.add_argument("-o2" "--output_bam_2", action="store", dest="output_bam_2", default="without_specified_tag.bam", help="Specify the name of the second output bam file. Defaults to without_specified_tag.bam.")
parser.add_argument("--discard_missing_tags",action="store_true",help="Store names of discarded reads?")

args = parser.parse_args()

#%% Read in BAM file
    
infile = pysam.AlignmentFile(args.input_bam, 'rb', check_sq=False, check_header=False)
outfile_1 = pysam.AlignmentFile(args.output_bam_1, 'wb', template=infile)
outfile_2 = pysam.AlignmentFile(args.output_bam_2, 'wb', template=infile)

n = [0, 0, 0]

for entry in infile.fetch(until_eof=True):
    if not entry.has_tag(args.tag_id):
        if args.discard_missing_tags:
            n[0]+=1
            continue
        outfile_2.write(entry)
        n[1]+=1
        continue
    
    value_of_tag = entry.get_tag(args.tag_id)

    if args.max_edit_dist > 0:
        edit_dist = hamming(args.tag_value, value_of_tag)
        if edit_dist <= args.max_edit_dist:
            outfile_1.write(entry)
            n[2]+=1
        else:
            outfile_2.write(entry)
            n[1]+=1
    else:
        if value_of_tag == args.tag_value:
            outfile_1.write(entry)
            n[2]+=1
        else:
            outfile_2.write(entry)
            n[1]+=1

infile.close()

print('Processed %d lines' % sum(n))
print('Found %d lines without the tag %s' % (n[0],args.tag_id))
print('Found %d lines where %s != %s' % (n[1],args.tag_id,args.tag_value))
print('Found %d lines where %s == %s with an edit distance of maximum %d' % (n[2],args.tag_id,args.tag_value, args.max_edit_dist))


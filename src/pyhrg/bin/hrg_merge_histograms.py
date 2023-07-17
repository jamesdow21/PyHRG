"""
MIT License

Copyright (c) 2013 Nicholas Dronen

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

Modified 2023 James Dow <jamesdow21@gmail.com>
"""

import networkx as nx
import pickle
import collections

import os

from pyhrg.hrg import ConsensusDendrogramBuilder

from optparse import OptionParser


def main():
    parser = OptionParser(
        description="Merges the split histograms of two or more HRG consensus dendrograms.  Creates a new consensus dendrogram from the merged histograms.  Saves the new consensus dendrogram to a graph markup language (GML) file.",
        prog='hrg-merge-histograms.py',
        usage='%prog [options] GRAPH_EDGELIST_FILE PICKLED_HISTOGRAM ... PICKLED_HISTOGRAM OUTPUT_FILE')

    parser.add_option('-f', '--force', action='store', type=int, default=10000,
        help='Allow overwriting of existing GML dendrogram files')

    (options, args) = parser.parse_args()

    if len(args) < 4:
        parser.print_help()
        return 1

    graph_edgelist = args[0]
    G = nx.read_edgelist(graph_edgelist, nodetype=int)
    filename = os.path.basename(graph_edgelist)
    G.name = os.path.splitext(filename)[0]
    args.remove(graph_edgelist)

    outfile=args.pop()

    if os.path.exists(outfile) and not options.force:
        raise Exception("Output file " + outfile +
            " exists.  Won't overwrite without --force option.")

    n = 0

    histograms = []
    for histfile in args:
        f = open(histfile, 'rb')
        histogram = pickle.load(f)
        if not isinstance(histogram, collections.Mapping):
            raise Exception('Object in ' + histfile +
                ' is not a dictionary: ' + str(type(histogram)))
        if n == 0:
            n = histogram['num_samples']

        if histogram['num_samples'] != n:
            raise Exception('inconsistent number of samples, '
                'expected ' + str(n) + ', actual ' +
                histogram['num_samples'])

        del histogram['num_samples']

        histograms.append(histogram)

    nodes = G.nodes()
    nodes.sort()

    builder = ConsensusDendrogramBuilder()
    C = builder.build(nodes, histograms, n)

    # Save the consensus dendrogram to a GML file.
    nx.write_gml(C, outfile)
    print("Saved merged consensus dendrogram to " + outfile + ".")

    return 0

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


from __future__ import print_function

import networkx as nx
import random
import pickle

import os
import sys

from pyhrg.hrg import Dendrogram, ConsensusDendrogramBuilder

from optparse import OptionParser


def print_status(*l):
    print("\t".join([str(x) for x in l]))

def main():
    parser = OptionParser(
        description="Finds a consensus dendrogram from an HRG model of a network.  Saves the consensus dendrogram to a graph markup language (GML) file.  Saves the histogram of splits in the consensus dendrogram to a file in Python's pickle format.",
        prog='hrg-consensus.py',
        usage='%prog [options] GRAPH_EDGELIST_FILE DENDROGRAM_GML_FILE')

    parser.add_option('-s', '--num-samples', action='store', type=int,
        default=10000, help='The number of times to sample the dendrogram\'s splits (default=10000).')

    parser.add_option('-t', '--temperature', action='store', type=float,
        default=2.0, help='The temperature at which to run (default=2.0).')

    (options, args) = parser.parse_args()

    if len(args) != 2:
        parser.print_help()
        return 1

    graph_edgelist=args[0]
    G=nx.read_edgelist(graph_edgelist, nodetype=int)
    filename=os.path.basename(graph_edgelist)
    G.name=os.path.splitext(filename)[0]

    gml_file=args[1]
    D=Dendrogram.from_gml_file(gml_file, G)

    bestL=initL=D.graph['L']
    bestI=0

    print_status("step", "L", "best L", "% complete", "consensus size")

    threshold = 1/(50.0*G.number_of_nodes())
    burnin = 200*G.number_of_nodes()
    i=1

    out = os.path.splitext(graph_edgelist)[0]
    out += '-consensus-temp-%0.2f' % options.temperature
    dendro_file = out + '-dendrogram.gml'
    hist_file = out + '-histogram.dat'
    print("HRG consensus dendrogram will be saved as " + dendro_file)
    print("Split histogram will be saved as " + hist_file)

    while D.num_samples < options.num_samples:
        taken=D.monte_carlo_move(T=options.temperature, debug=False)

        if i > burnin and random.random() < threshold:
            D.sample_splits()

        t = ''
        if taken:
            t = '*'
        if D.graph['L'] > bestL:
            bestL=D.graph['L']

        if i % 4096 == 0:
            nsplits = D.num_samples
            pct_complete = 100 * D.num_samples / float(options.num_samples)
            print_status(
                "[" + str(i) + "]",
                "%.3f" % D.graph['L'],
                "%.3f" % bestL,
                "%8.2f" % pct_complete,
                "%10d" % nsplits)

        if i % 10 == 0:
            sys.stdout.flush()

        i+=1

    # Save the histogram to a file.
    D.split_histogram['num_samples'] = D.num_samples
    pickle.dump(D.split_histogram, open(hist_file, mode='wb'))
    del D.split_histogram['num_samples']
    print("Saved split histogram to " + hist_file)

    # Build the consensus dendrogram, save it to a file.
    builder = ConsensusDendrogramBuilder()
    C = builder.build(D.graph_nodes_list, D.split_histogram, D.num_samples)
    nx.write_gml(C, out + '-dendrogram.gml')
    print("Saved consensus dendrogram to " + dendro_file)

    return 0

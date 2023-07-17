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

import optparse
import os
import sys

from pyhrg.hrg import Dendrogram

def print_status(*l):
    print("\t".join([str(x) for x in l]))

def main():
    parser = optparse.OptionParser(
        description='Fits a hierarchical random graph (HRG) model to a network.  Saves the model to a file in graph markup language (GML) format.',
        prog='hrg-fit.py',
        usage='%prog [options] GRAPH_EDGELIST_FILE')

    parser.add_option('-s', '--num-steps', action='store', type=int,
        default=100000,
        help='The number of MCMC steps to take (default=100000).')

    parser.add_option('-t', '--nodetype', action='store', type='choice',
        choices=[int,str],
        default=int,
        help='The type of the nodes in the edgelist file; "int" or "str" (default="int")')

    (options, args) = parser.parse_args()

    if len(args) != 1:
        parser.print_help()
        return 1

    G=nx.read_edgelist(args[0], nodetype=options.nodetype)
    name=os.path.splitext(args[0])[0]
    hrg_file = name + '-hrg.gml'
    print("HRG model will be saved as " + hrg_file + ".")

    D=Dendrogram.from_graph(G)

    bestL=initL=D.graph['L']
    prevL=bestL
    bestI=0

    print_status("step", "L", "best L", "MC step", "deltaL")

    for i in range(1, options.num_steps):
        taken=D.monte_carlo_move()
        t = ''
        if taken:
            t = '*'
        if D.graph['L'] > bestL:
            bestL=D.graph['L']
            bestI=i
            nx.write_gml(D, hrg_file)
            print_status("["+str(i)+"]", "%.3f" % bestL, "%.3f" % bestL, t, "%.3f"%D.deltaL)
        elif i % 4096 == 0:
            print_status("["+str(i)+"]", "%.3f" % D.graph['L'], "%.3f" % bestL, t, "%.3f"%D.deltaL)

        prevL=D.graph['L']

        if i % 10 == 0:
            sys.stdout.flush()

    print("Step number of last best fit "+str(bestI) + ".")
    print("HRG model was saved as " + hrg_file + ".")

    return 0

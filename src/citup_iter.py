import csv
import sys
import logging
import os
import ConfigParser
import itertools
import argparse
import string
import gzip
import ConfigParser
import shutil
from collections import *
import numpy as np
import pandas as pd

import pypeliner
import pypeliner.managed as mgd

import citup

if __name__ == '__main__':

    argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    pypeliner.app.add_arguments(argparser)

    argparser.add_argument('input_freqs',
                           help='Input Mutation Frequencies')

    argparser.add_argument('output_results',
                           help='Output Results HDF5 Store')

    argparser.add_argument('--min_nodes', type=int, default=1,
                           help='Output For All Trees')

    argparser.add_argument('--max_nodes', type=int, default=8,
                           help='Output For All Trees')

    argparser.add_argument('--max_children_per_node', type=int, default=100,
                           help='Output For All Trees')

    args = vars(argparser.parse_args())

    pyp = pypeliner.app.Pypeline(modules=[citup], config=args)

    citup_bin_directory = os.path.join(os.path.dirname(citup.__file__))
    citup_iter_tool = os.path.join(citup_bin_directory, 'citupiter')

    lowmem = {'mem': 1, 'ncpu': 1}

    workflow = pypeliner.workflow.Workflow()

    workflow.transform(
        name='create_trees',
        ctx=lowmem,
        func=citup.create_trees,
        ret=mgd.TempOutputObj('trees', 'tree'),
        args=(
            int(args['min_nodes']),
            int(args['max_nodes']),
            int(args['max_children_per_node']),
        ),
    )
    
    workflow.commandline(
        name='run_citup',
        axes=('tree',),
        ctx=lowmem,
        args=(
            citup_iter_tool,
            mgd.TempInputObj('trees', 'tree').prop('unlabeled_tree_string'),
            mgd.InputFile(args['input_freqs']),
            mgd.TempOutputFile('results', 'tree'),
        ),
    )

    workflow.transform(
        name='select_optimal_tree',
        ctx=lowmem,
        func=citup.select_optimal_tree,
        args=(
            mgd.InputFile(args['input_freqs']),
            mgd.TempInputObj('trees', 'tree'),
            mgd.TempInputFile('results', 'tree'),
            mgd.OutputFile(args['output_results']),
        ),
    )

    pyp.run(workflow)


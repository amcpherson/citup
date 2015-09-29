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
from sklearn import mixture

import pypeliner
import pypeliner.managed as mgd

import citup
import BeyerHedetmieni
import treenode

if __name__ == '__main__':

    argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    pypeliner.app.add_arguments(argparser)

    argparser.add_argument('input_freqs',
                           help='Input Mutation Frequencies')

    argparser.add_argument('input_clusters',
                           help='Input Mutation Clusters')

    argparser.add_argument('output_results',
                           help='Output Results HDF5 Store')

    argparser.add_argument('--min_nodes', type=int, default=1,
                           help='Output For All Trees')

    argparser.add_argument('--max_nodes', type=int, default=8,
                           help='Output For All Trees')

    argparser.add_argument('--max_children_per_node', type=int, default=100,
                           help='Output For All Trees')

    args = vars(argparser.parse_args())

    pyp = pypeliner.app.Pypeline([citup], args)

    citup_bin_directory = os.path.abspath(os.path.dirname(citup.__file__))
    citupqip_tool = os.path.join(citup_bin_directory, 'citupqip')

    lowmem = {'mem':1}
    himem = {'mem':8}

    pyp.sch.transform('create_trees', (), lowmem,
        citup.create_trees,
        mgd.TempOutputObj('trees', 'tree'),
        int(args['min_nodes']),
        int(args['max_nodes']),
        int(args['max_children_per_node']))
    
    pyp.sch.commandline('run_citup', ('tree',), himem,
        citupqip_tool,
        mgd.TempInputObj('trees', 'tree').prop('unlabeled_tree_string'),
        mgd.InputFile(args['input_freqs']),
        mgd.InputFile(args['input_clusters']),
        mgd.TempOutputFile('results', 'tree'))

    pyp.sch.transform('select_optimal_tree', (), lowmem, 
        citup.select_optimal_tree,
        None,
        mgd.InputFile(args['input_freqs']),
        mgd.TempInputObj('trees', 'tree'),
        mgd.TempInputFile('results', 'tree'),
        mgd.OutputFile(args['output_results']))

    pyp.run()

else:

    class TreeInfo(object):
        def __init__(self, num_nodes, tree_index, tree):
            self.num_nodes = num_nodes
            self.tree_index = tree_index
            self.tree = tree
        @property
        def unlabeled_tree_string(self):
            return self.tree.create_unlabeled_tree_string()
        @property
        def labeled_tree_string(self):
            return self.tree.create_labeled_tree_string()
        def __eq__(self, other):
            return self.num_nodes == other.num_nodes and self.tree_index == other.tree_index

    def generate_trees(min_nodes, max_nodes, max_children_per_node):
        for num_nodes in xrange(min_nodes, max_nodes + 1):
            for tree_index, parent_array in enumerate(BeyerHedetmieni.getParentArrays(num_nodes, max_children_per_node)):
                yield TreeInfo(num_nodes, tree_index, treenode.create_from_parent_array(parent_array))

    def create_trees(min_nodes, max_nodes, max_children_per_node):
        return dict(enumerate(generate_trees(min_nodes, max_nodes, max_children_per_node)))

    def read_dataset(filename):
        dataset = dict()
        with open(filename, 'r') as f:
            for chunk in f.read().split('#'):
                if chunk == '':
                    continue
                metadata, data = chunk.split('\n', 1)
                metadata = metadata.split()
                name = metadata[0]
                dtype = eval(metadata[1])
                shape = tuple([int(a) for a in metadata[2:]])
                dataset[name] = np.array(data.split(), dtype=dtype).reshape(shape)
        return dataset

    def create_results_entry(tree_id, tree, freq, results_data):
        results_info = dict()
        results_info['tree_id'] = tree_id
        results_info['num_nodes'] = tree.num_nodes
        results_info['tree_index'] = tree.tree_index
        results_info['tree_string'] = tree.labeled_tree_string
        results_info['num_mutations'] = freq.shape[0]
        results_info['num_samples'] = freq.shape[1]
        results_info['objective_value'] = results_data['objective_value']
        try:
            results_info['cplex_status'] = results_data['cplex_status']
            results_info['cplex_time'] = results_data['cplex_hours']
        except KeyError:
            pass
        return results_info
    
    def estimate_error_rate(freq):
        gmm_bics = list()
        for n in range(1, min(21, freq.shape[0]+1)):
            gmm = mixture.GMM(n_components=n, covariance_type='spherical')
            gmm.fit(freq)
            gmm_bics.append((gmm.bic(freq), n))
        n = sorted(gmm_bics, key=lambda a: a[0])[0][1]
        gmm = mixture.GMM(n_components=n, covariance_type='spherical')
        gmm.fit(freq)
        return np.sqrt(np.mean(gmm.covars_))
    
    def select_optimal_tree(freq_filename, trees, results_filenames, output_filename):
        with pd.HDFStore(output_filename, 'w') as store:
            freq = np.loadtxt(freq_filename)
            error_rate = estimate_error_rate(freq)
            results_table = list()
            for tree_id in trees.keys():
                results_data = read_dataset(results_filenames[tree_id])
                for key, value in results_data.iteritems():
                    table_key = 'trees/{0}/{1}'.format(tree_id, key)
                    if len(value.shape) == 2:
                        store[table_key] = pd.DataFrame(value)
                    elif len(value.shape) == 1:
                        store[table_key] = pd.Series(value)
                    elif len(value.shape) == 0:
                        store[table_key] = pd.Series([value])
                results_table.append(pd.DataFrame(create_results_entry(tree_id, trees[tree_id], freq, results_data), index=[0]))
            results_table = pd.concat(results_table, ignore_index=True)
            results_table['error_rate'] = error_rate
            results_table['likelihood'] = results_table['objective_value'] / (2.0 * error_rate * error_rate)
            results_table['bic'] = 2.0 * results_table['likelihood'] + results_table['num_samples'] * (results_table['num_nodes'] - 1.0) * np.log(results_table['num_mutations'])
            results_table.sort('bic', inplace=True)
            results_table['optimal'] = False
            results_table['optimal'].iloc[0] = True
            for name, column in results_table.iteritems():
                store['results/{0}'.format(name)] = column


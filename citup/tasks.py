import numpy as np
import pandas as pd
import sklearn.mixture


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
        gmm = sklearn.mixture.GMM(n_components=n, covariance_type='spherical')
        gmm.fit(freq)
        gmm_bics.append((gmm.bic(freq), n))
    n = sorted(gmm_bics, key=lambda a: a[0])[0][1]
    gmm = sklearn.mixture.GMM(n_components=n, covariance_type='spherical')
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


import argparse
import random
import numpy as np

import BeyerHedetmieni
import treenode

if __name__ == '__main__':

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('frequencies', help='Output mutation frequency matrix')
    parser.add_argument('samples', help='Sample by node matrix')
    parser.add_argument('siminfo', help='Simulation info')
    parser.add_argument('-e', '--error_rate', type=float, default=0.01, help='Frequency measurement error')
    parser.add_argument('-m', '--num_mutations', type=int, default=10, help='Number of mutations to simulate')
    parser.add_argument('-s', '--num_samples', type=int, default=3, help='Numner of tumour samples')
    parser.add_argument('-n', '--num_nodes', type=int, default=4, help='Number of nodes in the tree')
    parser.add_argument('-c', '--max_child', type=int, default=2, help='Maximum number of children per node')
    parser.add_argument('-d', '--dir_alpha', type=float, default=1.0, help='Node distribution dirichlet alpha')
    parser.add_argument('-r', '--rng_seed', type=int, default=1234, help='Random number generator seed')

    args = parser.parse_args()

    random.seed(args.rng_seed)
    np.random.seed(args.rng_seed)

    parent_array = random.choice(list(BeyerHedetmieni.getParentArrays(args.num_nodes, args.max_child)))

    tree = treenode.create_from_parent_array(parent_array)

    # Sample by node
    node_sampling = np.random.dirichlet([args.dir_alpha] * args.num_nodes, size=args.num_samples)

    # Subtree by node
    gamma_matrix = np.zeros((args.num_nodes, args.num_nodes))
    tree.fill_gamma_matrix(gamma_matrix)

    # Mutation by subtree
    mutation_assignment = np.random.multinomial(1, [1/float(args.num_nodes)]*args.num_nodes, size=args.num_mutations)

    # Mutation by sample
    frequencies = np.mat(mutation_assignment) * np.mat(gamma_matrix) * np.mat(node_sampling.T)

    # Noise
    frequencies += args.error_rate * np.random.randn(*frequencies.shape)

    np.savetxt(args.frequencies, frequencies)

    np.savetxt(args.samples, node_sampling)

    with open(args.siminfo, 'w') as f:
        f.write('tree: {0}\n'.format(tree.create_labeled_tree_string()))


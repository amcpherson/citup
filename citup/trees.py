import citup.BeyerHedetmieni
import citup.treenode


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
        for tree_index, parent_array in enumerate(citup.BeyerHedetmieni.getParentArrays(num_nodes, max_children_per_node)):
            yield TreeInfo(num_nodes, tree_index, citup.treenode.create_from_parent_array(parent_array))


def create_trees(min_nodes, max_nodes, max_children_per_node):
    return dict(enumerate(generate_trees(min_nodes, max_nodes, max_children_per_node)))


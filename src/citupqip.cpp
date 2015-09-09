#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <time.h>
#include <vector>
#include <assert.h>

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

#include "trees.h"

#include <ilcplex/ilocplex.h>

using namespace std;
using namespace boost;


template <typename TType>
void ReadMatrixData(istream& file, vector<vector<TType> >& matrix)
{
	string line;
	while (getline(file, line))
	{
		vector<string> fields;
		split(fields, line, is_any_of(" \t"));

		matrix.push_back(vector<TType>());
		for (vector<string>::const_iterator iter = fields.begin(); iter != fields.end(); iter++)
		{
			matrix.back().push_back(lexical_cast<TType>(*iter));
		}
	}
}

template <typename TType>
void ReadVectorData(istream& file, vector<TType>& vec)
{
	string line;
	while (getline(file, line))
	{
		vec.push_back(lexical_cast<TType>(line));
	}
}

template <typename TType>
void WriteMatrixData(ostream& out, const vector<vector<TType> >& mat)
{
	for(int i = 0; i < mat.size(); i++)
	{
		for(int j = 0; j < mat[i].size(); j++)
		{
			out << mat[i][j] << " ";
		}
		out << endl;
	}
}

template <typename TType>
void WriteVectorData(ostream& out, const vector<TType>& vec)
{
	for(int i = 0; i < vec.size(); i++)
	{
		out << vec[i] << " ";
	}
	out << endl;
}

template <typename TType>
const char* ConvertTypeToString()
{
	throw;
}

template <>
const char* ConvertTypeToString<double>()
{
	return "float";
};

template <>
const char* ConvertTypeToString<float>()
{
	return "float";
};

template <>
const char* ConvertTypeToString<int>()
{
	return "int";
};

template <typename TType>
void WriteMatrix(ostream& out, const string& name, const vector<vector<TType> >& mat)
{
	out << "#" << name << " " << ConvertTypeToString<TType>() << " " << mat.size() << " " << mat[0].size() << endl;
	WriteMatrixData(out, mat);
}

template <typename TType>
void WriteVector(ostream& out, const string& name, const vector<TType>& vec)
{
	out << "#" << name << " " << ConvertTypeToString<TType>() << " " << vec.size() << endl;
	WriteVectorData(out, vec);
}

template <typename TType>
void WriteValue(ostream& out, const string& name, const TType& val)
{
	out << "#" << name << " " << ConvertTypeToString<TType>() << endl;
	out << val << endl;
}

vector<int> RandomAssignment(int size, int choices)
{
	vector<int> assignments;

	for (int idx = 0; idx < size; idx++)
	{
		assignments.push_back(rand() % choices);
	}

	return assignments;
}

struct IloEnvSafe
{
	~IloEnvSafe()
	{
		env.end();
	}

	IloEnv env;
};

double solve_tree(const vector<vector<int> >& gamma, const vector<vector<double> >& variant_freq,
	const vector<int>& variant_cluster, vector<vector<double> >& clone_freq,
	vector<vector<double> >& clade_freq, vector<int>& cluster_assignment,
	bool verbose=false)
{
	try
	{
		// Cplex environment
		IloEnvSafe safe_env;
		IloEnv& env = safe_env.env;

		// Each column of the variant frequency matrix specifies frequences of variants in a sample
		int num_variants = variant_freq.size();

		// Each row of the variant frequency matrix are the frequencies in each sample
		int num_samples = variant_freq[0].size();

		// Number of nodes from gamma matrix (one row/column per node)
		int num_nodes = gamma.size();

		// Number of mutation clusters
		int num_clusters = *(max_element(variant_cluster.begin(), variant_cluster.end())) + 1;

		// Node frequency vars, matrix of node by sample
		IloArray<IloFloatVarArray> clone_freq_vars(env, num_nodes);
		for (int node = 0; node < num_nodes; node++)
		{
			clone_freq_vars[node] = IloFloatVarArray(env, num_samples, 0, 1);
		}

		// Node frequency vars, matrix of clade by sample
		IloArray<IloFloatVarArray> clade_freq_vars(env, num_nodes);
		for (int node = 0; node < num_nodes; node++)
		{
			clade_freq_vars[node] = IloFloatVarArray(env, num_samples, 0, 1);
		}

		// Cluster assignment matrix
		IloArray<IloIntVarArray> delta(env, num_clusters);
		for(int cluster = 0; cluster < num_clusters; cluster++)
		{
			delta[cluster] = IloIntVarArray(env, num_nodes, 0, 1);
		}

		// Optimization variables
		IloArray<IloArray<IloFloatVarArray> > x(env, num_variants);
		IloArray<IloArray<IloFloatVarArray> > y(env, num_variants);
 
		for (int variant = 0; variant < num_variants; variant++)
		{
			x[variant] = IloArray<IloFloatVarArray>(env,num_nodes);
			y[variant] = IloArray<IloFloatVarArray>(env,num_nodes);
		}
		
		for (int variant = 0; variant < num_variants; variant++)
		{
			for (int node = 0; node < num_nodes; node++)
			{
				x[variant][node] = IloFloatVarArray(env, num_samples, 0, 1);
				y[variant][node] = IloFloatVarArray(env, num_samples, 0, 1);
			}
		}
	
		// Constraint for node frequencies summing to one for each sample
		IloConstraintArray node_constraints(env);
		for (int sample = 0; sample < num_samples; sample++)
		{
			IloExpr node_expr(env);
			for (int node = 0; node < num_nodes; node++)
			{
				node_expr += clone_freq_vars[node][sample];
			}

			node_constraints.add(node_expr == 1.0);
		}

		// Constraint for relationship between node frequency and clade frequency
		IloConstraintArray clade_constraints(env);
		for (int sample = 0; sample < num_samples; sample++)
		{
			for (int clade = 0; clade < num_nodes; clade++)
			{
				assert(clade < gamma.size());

				IloExpr clade_expr(env);
				for (int node = 0; node < num_nodes; node++)
				{
					assert(node < gamma[clade].size());

					if (gamma[clade][node])
					{
						clade_expr += clone_freq_vars[node][sample];
					}
				}

				clade_constraints.add(clade_freq_vars[clade][sample] == clade_expr);
			}
		}

		// Constraint for single assignment of a cluster to a node
		IloConstraintArray assignment_constraints(env);
		for(int cluster = 0; cluster < num_clusters; cluster++)
		{
			IloIntExpr cluster_assignment_expr(env);

			for(int node = 0; node < num_nodes; node++)
			{ 
				cluster_assignment_expr += delta[cluster][node];
			}
			
			assignment_constraints.add(cluster_assignment_expr == 1);
		}

		// Squared error between clade and variant freq as objective
		IloConstraintArray objective_constraints(env);
		IloExpr objective_expr(env);
		for (int sample = 0; sample < num_samples; sample++)
		{
			for (int variant = 0; variant < num_variants; variant++)
			{
				for (int node = 0; node < num_nodes; node++)
				{
					double var_freq = variant_freq[variant][sample];

					objective_constraints.add(x[variant][node][sample] >= (var_freq - clade_freq_vars[node][sample]));
					objective_constraints.add(x[variant][node][sample] >= (clade_freq_vars[node][sample] - var_freq));

					objective_constraints.add(y[variant][node][sample] >= (1)*(delta[variant_cluster[variant]][node] - 1) + x[variant][node][sample]);

					objective_expr += y[variant][node][sample] * y[variant][node][sample];
				}
			}
		}

		// Create model
		IloModel model(env);
		model.add(IloMinimize(env, objective_expr));
		model.add(clade_constraints);
		model.add(node_constraints);
		model.add(assignment_constraints);
		model.add(objective_constraints);
		
		// Create optimizer
		IloCplex cplex(model);
		
		// Set cplex parameters
		int cplex_time_limit = 3600;
		int cplex_threads = 1;
		int cplex_trelimit = 4;
		cplex.setParam(IloCplex::TiLim, cplex_time_limit);
		cplex.setParam(IloCplex::EpGap, 0.002);
		cplex.setParam(IloCplex::Threads, cplex_threads);
		cplex.setParam(IloCplex::TreLim, cplex_trelimit*1024);
		cplex.setParam(IloCplex::ParallelMode, 1);
		cplex.setParam(IloCplex::DataCheck, true);

		// Optionally print output to stderr
		if (verbose)
		{
			cplex.setOut(env.error());
			cplex.setWarning(env.error());
		}
		else
		{
			cplex.setOut(env.getNullStream());
			cplex.setWarning(env.getNullStream());
		}
		
		// Solve
		if (!cplex.solve())
		{
			ostringstream status;
			status << cplex.getStatus();
			cerr << "Warning: failed to optimize LP: " << status.str() << endl;
			return -1.0;
		}
		
		// DEBUG
		//cplex.exportModel("model.lp");
		//cplex.writeSolution("solution.lp");

		// Objective value
		double objective_value = cplex.getObjValue();
		
		// Solution values
		clone_freq.clear();
		clade_freq.clear();
		for (int sample = 0; sample < num_samples; sample++)
		{
			clone_freq.push_back(vector<double>());
			clade_freq.push_back(vector<double>());
			for (int clade = 0; clade < num_nodes; clade++)
			{
				try
				{
					clone_freq.back().push_back(cplex.getValue(clone_freq_vars[clade][sample]));
				}
				catch (IloException& e)
				{
					cerr << "Warning: unable to extract node frequency " << clade << " " << sample << endl;
				}

				try
				{
					clade_freq.back().push_back(cplex.getValue(clade_freq_vars[clade][sample]));
				}
				catch (IloException& e)
				{
					cerr << "Warning: unable to extract clade frequency " << clade << " " << sample << endl;
				}
			}
		}

		cluster_assignment.clear();
		cluster_assignment.resize(num_clusters);
		for(int cluster = 0; cluster < num_clusters; cluster++)
		{
			for (int node = 0; node < num_nodes; node++)
			{
				try
				{
					int assigned = cplex.getValue(delta[cluster][node]);

					if (assigned)
					{
						cluster_assignment[cluster] = node;
					}
				}
				catch (IloException& e)
				{
					cerr << "Warning: unable to extract cluster assignment " << cluster << " " << node << endl;
				}
			}
		}

		return objective_value;
	}
	catch (IloException& e)
	{
		string message = e.getMessage();
		e.end();
		throw runtime_error("Cplex Exception: " + message);
	}
}

void ReadFrequencies(istream& file, vector<vector<double> >& variant_freq)
{
	variant_freq.clear();

	ReadMatrixData(file, variant_freq);

	int num_variants = variant_freq.size();

	if (num_variants == 0)
	{
		cerr << "Error: no mutations" << endl;
		exit(1);
	}

	int num_samples = variant_freq[0].size();

	for (int variant = 0; variant < num_variants; variant++)
	{
		if (variant_freq[variant].size() != num_samples)
		{
			cerr << "Error: inconsistent variant matrix" << endl;
			exit(1);
		}
	}
}

void ReadClusters(istream& file, vector<int>& variant_cluster)
{
	variant_cluster.clear();
	ReadVectorData(file, variant_cluster);
}

struct Solution
{
	Solution() : objective_value(-1.0) {}

	vector<int> cluster_assignment;
	vector<vector<double> > clone_freq;
	vector<vector<double> > clade_freq;
	double objective_value;
};

void WriteSolution(ostream& file, const Node* tree, const Solution& solution)
{
	// Write objective
	WriteValue(file, "objective_value", solution.objective_value);

	// Write adjacency matrix
	WriteMatrix(file, "adjacency_list", tree->get_adjacency_list());

	// Write gamma matrix
	WriteMatrix(file, "gamma_matrix", tree->get_gamma_matrix());

	// Write clone frequencies
	WriteMatrix(file, "clone_freq", solution.clone_freq);

	// Write clade frequencies
	WriteMatrix(file, "clade_freq", solution.clade_freq);

	// Write assignments
	WriteVector(file, "cluster_assignment", solution.cluster_assignment);
}

int main(int argc, char** argv)
{
	if (argc < 5)
	{
		cerr << "Usage: tree_string variant_freq_filename variant_cluster_filename solution_filename" << endl;

		return 1;
	}

	string tree_string = argv[1];
	string variant_freq_filename = argv[2];
	string variant_cluster_filename = argv[3];
	string solution_filename = argv[4];

	Node* tree = interpret_tree_string(tree_string);

	int num_nodes = tree->count_nodes();

	vector<vector<int> > gamma_matrix = tree->get_gamma_matrix();

	ifstream variant_freq_file(variant_freq_filename.c_str());
	assert(variant_freq_file.good());
	vector<vector<double> > variant_freq;
	ReadFrequencies(variant_freq_file, variant_freq);

	ifstream variant_cluster_file(variant_cluster_filename.c_str());
	assert(variant_cluster_file.good());
	vector<int> variant_cluster;
	ReadClusters(variant_cluster_file, variant_cluster);

	Solution solution;

	solution.objective_value = solve_tree(gamma_matrix, variant_freq, variant_cluster, solution.clone_freq, solution.clade_freq, solution.cluster_assignment);

	ofstream solution_file(solution_filename.c_str());
	WriteSolution(solution_file, tree, solution);
}




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

double solve_tree_clone_freq(const vector<vector<int> >& gamma, const vector<vector<double> >& variant_freq,
	const vector<int>& variant_assignment, vector<vector<double> >& node_freq,
	vector<vector<double> >& clone_freq, bool verbose=false)
{
	try
	{
		// Cplex environment
		IloEnvSafe safe_env;
		IloEnv& env = safe_env.env;

		// Each row of the variant frequency matrix are the frequencies in each sample
		int num_samples = variant_freq[0].size();

		// Number of nodes from gamma matrix (one row/column per node)
		int num_nodes = gamma.size();

		// Node frequency vars, matrix of node by sample
		IloNumVarArray node_freq_vars(env, num_nodes * num_samples, 0.0, 1.0);

		// Node frequency vars, matrix of clone by sample
		IloNumVarArray clone_freq_vars(env, num_nodes * num_samples, 0.0, 1.0);

		// Constraint for node frequencies summing to one for each sample
		IloConstraintArray node_constraints(env);
		for (int sample = 0; sample < num_samples; sample++)
		{
			IloExpr node_expr(env);
			for (int node = 0; node < num_nodes; node++)
			{
				node_expr += node_freq_vars[node + sample * num_nodes];
			}

			node_constraints.add(node_expr == 1.0);
		}

		// Constraint for relationship between node frequency and clone frequency
		IloConstraintArray clone_constraints(env);
		for (int sample = 0; sample < num_samples; sample++)
		{
			for (int clone = 0; clone < num_nodes; clone++)
			{
				assert(clone < gamma.size());

				IloExpr clone_expr(env);
				for (int node = 0; node < num_nodes; node++)
				{
					assert(node < gamma[clone].size());

					if (gamma[clone][node])
					{
						clone_expr += node_freq_vars[node + sample * num_nodes];
					}
				}

				clone_constraints.add(clone_freq_vars[clone + sample * num_nodes] == clone_expr);
			}
		}

		// Squared error between clone and variant freq as objective
		IloExpr objective_expr(env);
		for (int sample = 0; sample < num_samples; sample++)
		{
			for (int variant = 0; variant < variant_freq.size(); variant++)
			{
				assert(variant < variant_assignment.size());

				objective_expr += (clone_freq_vars[variant_assignment[variant] + sample * num_nodes] - variant_freq[variant][sample]) *
				                  (clone_freq_vars[variant_assignment[variant] + sample * num_nodes] - variant_freq[variant][sample]);
			}
		}

		// Create model
		IloModel model(env);
		model.add(IloMinimize(env, objective_expr));
		model.add(clone_constraints);
		model.add(node_constraints);
		
		// Create optimizer
		IloCplex cplex(model);
		
		// Optionally print output to stderr
		if (verbose)
		{
			cplex.setOut(env.error());
		}
		else
		{
			cplex.setOut(env.getNullStream());
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
		node_freq.clear();
		clone_freq.clear();
		for (int sample = 0; sample < num_samples; sample++)
		{
			node_freq.push_back(vector<double>());
			clone_freq.push_back(vector<double>());
			for (int clone = 0; clone < num_nodes; clone++)
			{
				try
				{
					node_freq.back().push_back(cplex.getValue(node_freq_vars[clone + sample * num_nodes]));
				}
				catch (IloException& e)
				{
					cout << "Warning: unable to extract node frequency " << clone << " " << sample << endl;
				}

				try
				{
					clone_freq.back().push_back(cplex.getValue(clone_freq_vars[clone + sample * num_nodes]));
				}
				catch (IloException& e)
				{
					cout << "Warning: unable to extract clone frequency " << clone << " " << sample << endl;
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

double reassign_variants(const vector<vector<double> >& variant_freq, const vector<vector<double> >& clade_freq,
	vector<int>& variant_assignment)
{
	double squared_error = 0.0;

	// Each row of the variant frequency matrix are the frequencies in each sample
	int num_samples = variant_freq[0].size();

	// Each row of the clade frequency matrix are the frequencies for each clade in that sample
	int num_nodes = clade_freq[0].size();

	for (int variant = 0; variant < variant_freq.size(); variant++)
	{
		vector<double> clone_squared_error;
		for (int clone = 0; clone < num_nodes; clone++)
		{
			double squared_error = 0.0;
			for (int sample = 0; sample < num_samples; sample++)
			{
				double abs_error = variant_freq[variant][sample] - clade_freq[sample][clone];
				squared_error += abs_error * abs_error;
			}

			clone_squared_error.push_back(squared_error);
		}

		vector<double>::const_iterator min_clone_iter = min_element(clone_squared_error.begin(), clone_squared_error.end());

		squared_error += *min_clone_iter;
		variant_assignment[variant] = min_clone_iter - clone_squared_error.begin();
	}

	return squared_error;
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

struct Solution
{
	Solution() : objective_value(-1.0) {}

	vector<int> variant_assignment;
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
	WriteVector(file, "variant_assignment", solution.variant_assignment);
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
	string solution_filename = argv[3];

	int max_iterations = 100;
	int num_restarts = 1000;
	int seed_start = 1;
	double objective_tol = 1e-9;

	srand(seed_start);

	Node* tree = interpret_tree_string(tree_string);

	int num_nodes = tree->count_nodes();

	vector<vector<int> > gamma_matrix = tree->get_gamma_matrix();

	ifstream variant_freq_file(variant_freq_filename.c_str());
	assert(variant_freq_file.good());
	vector<vector<double> > variant_freq;
	ReadFrequencies(variant_freq_file, variant_freq);

	// DEBUG
	//string variant_assignment_filename = "assign.txt";
	//ifstream variant_assignment_file(variant_assignment_filename.c_str());
	//assert(variant_assignment_file.good());
	//vector<int> variant_assignment;
	//ReadVector(variant_assignment_file, variant_assignment);

	Solution best_solution;
	
	int cplex_fail_count = 0;

	for (int restart = 0; restart < num_restarts; restart++)
	{
		Solution solution;

		solution.variant_assignment = RandomAssignment(variant_freq.size(), num_nodes);

		double previous_objective = -1.0;
		for (int i = 0; i < max_iterations; i++)
		{
			double objective_value = solve_tree_clone_freq(gamma_matrix, variant_freq, solution.variant_assignment, solution.clone_freq, solution.clade_freq);
			
			if (objective_value < 0.0)
			{
				cplex_fail_count++;
				break;
			}

			solution.objective_value = reassign_variants(variant_freq, solution.clade_freq, solution.variant_assignment);

			if (abs(solution.objective_value - previous_objective) < objective_tol)
			{
				break;
			}

			previous_objective = solution.objective_value;
		}

		if (best_solution.objective_value < 0.0)
		{
			best_solution = solution;
		}
		else if (solution.objective_value < best_solution.objective_value)
		{
			best_solution = solution;
		}
	}
	
	assert(cplex_fail_count < num_restarts / 100 + 2);
	
	ofstream solution_file(solution_filename.c_str());

	WriteSolution(solution_file, tree, best_solution);
}

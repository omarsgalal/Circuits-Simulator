#include<iostream>
#include<string>
#include<windows.h>

using namespace std;

// Structs declaration:
struct Element;
struct Node;
struct Current;


// Functions prototypes:
bool Read_elements_node(Node *& node, int index);
bool Read_line(Node * node, char check);
void solve_equations(double **coefficients, int num_nodes, int num_essential_nodes, Node** nodes);
void Make_Equations(Node **node, int num_essential_nodes, int num_nodes, bool is_one_loop);
void get_currents(Current **I_Branches, int I_Branches_counter);
void GetVoltage_NE(Node *nodes[], int num_nodes, bool is_one_loop);
void GetVoltage_NE_(Node *essential, Node *pre, Element *E, double I);
double get_In(Node**nodes, Element* R, int num_essential, int num_nodes, bool is_one_loop);
int get_nearest_essential(Element* element, Node* node, bool is_one_loop);
void OhmLaw(bool is_one_loop);
void Calculate_For_R(int i);
void Calculate_For_E(int i, bool is_one_loop);
void Calculate_For_J(int i);
void Balancecheck(double &P_dissipated, double &P_supplied);
void disable_sources(Element*onlysource);
void enable_sources();
double get_Rth(Node **node, int num_essential_nodes, int num_nodes, Element* need_Rth, bool is_one_loop);
double P_max(double In, double Rth);
void solve_one_loop(double sum_E, double sum_R);
double super_position(char required, Element* element, Element* source, Node** nodes, int num_essential_nodes, int num_nodes, bool is_one_loop, int V1, int V2);
void check_all_elements();
void print_outputs(Node ** nodes, int num_essential_nodes, int num_nodes, bool is_one_loop);
void frist_print();
void second_print();


// Global variavles:
Element *R_indexes[100];
Element  *E_indexes[100];
Element *J_indexes[100];
bool Exit_prog = false;



struct Current
{
	double I, sum_R, sum_E;
	Node *previous, *next;
};


struct Node
{
	Element *Branches[50];
	double voltage;
	int Branches_counter;
	int index;
	bool Is_essential;
	int index_essential;
	bool is_calculated;
};


struct Element
{
	Node *P_to_previous, *P_to_next;
	char type;
	int index;
	double I, R, V, original;
	Current *I_in_Element;
	int Pve_node, Nve_node;

	double P;
};


int main()
{
	Node *Nodes[50];
	bool is_one_loop = false;

	for (int i = 0; i < 50; i++)
	{
		Nodes[i] = NULL;
	}

	frist_print();

	int num_node, num_essintial_node = 0;
	for (num_node = 1; !Exit_prog && Read_elements_node(Nodes[num_node], num_node); num_node++);

	if (Exit_prog)
		return 0;

	check_all_elements();

	if (Exit_prog)
		return 0;

	num_node--;

	for (int j = 1; j <= num_node; j++)
	{
		Nodes[j]->is_calculated = false;
		if (Nodes[j]->Is_essential)
		{
			num_essintial_node++;
			Nodes[j]->index_essential = num_essintial_node;
			Nodes[j]->is_calculated = true;
		}

	}

	if (num_essintial_node == 0)
	{
		is_one_loop = true;
		Nodes[1]->Is_essential = true;
		num_essintial_node = 1;
		Nodes[1]->index_essential = 1;
	}


	Make_Equations(Nodes, num_essintial_node, num_node, is_one_loop);
	OhmLaw(is_one_loop);

	double P_dissipated, P_supplied;
	Balancecheck(P_dissipated, P_supplied);
	cout << "\nP_diss = " << P_dissipated << "\t" << "P_supp = " <<  P_supplied << endl << endl;
	
	second_print();
	print_outputs(Nodes, num_essintial_node, num_node, is_one_loop);
	
	
	return 0;
}


bool Read_line(Node * node, char check)
{
	short index;
	double value;
	cin >> index >> value;
	switch (check)
	{
	case 'E':
	{
		if (E_indexes[index] == NULL)
		{
			Element *V = new Element;
			V->P_to_next = NULL;
			V->P_to_previous = NULL;
			V->I_in_Element = NULL;
			V->V = value;
			V->type = check;
			if (value > 0)	V->Pve_node = node->index;
			else
				V->Nve_node = node->index;


			V->index = index;
			V->P_to_previous = node;
			node->Branches[node->Branches_counter] = V;
			(node->Branches_counter)++;
			E_indexes[index] = V;
		}
		else
		{
			Element *V = E_indexes[index];
			if (-1 * value != V->V)
			{
				cout << "Error: You entered the value of that element with the same sign or with another value\n";
				Exit_prog = true;
				return false;
			}
			else if (V->P_to_previous != NULL && V->P_to_next != NULL)
			{
				cout << "Error: You entered that element connected to two other nodes\n";
				Exit_prog = true;
				return false;
			}
			else
			{
				if (value > 0)
				{
					V->Pve_node = node->index;
					V->V *= -1;
				}
				else
					V->Nve_node = node->index;
				V->P_to_next = node;
				node->Branches[node->Branches_counter] = V;
				(node->Branches_counter)++;
			}
		}

		break;
	}
	case 'J':
	{
		if (J_indexes[index] == NULL)
		{
			Element  *J = new Element;
			J->P_to_next = NULL;
			J->P_to_previous = NULL;
			J->I_in_Element = NULL;
			J->I = value;
			J->type = check;
			if (value > 0)
				J->Pve_node = node->index;
			else
				J->Nve_node = node->index;
			J->index = index;
			J->P_to_previous = node;
			node->Branches[node->Branches_counter] = J;
			(node->Branches_counter)++;
			J_indexes[index] = J;
		}
		else
		{
			Element *J = J_indexes[index];
			if (-1 * value != J->I)
			{
				cout << "Error: You entered the value of that element with the same sign or with another value\n";
				Exit_prog = true;
				return false;
			}
			else if (J->P_to_previous != NULL && J->P_to_next != NULL)
			{
				cout << "Error: You entered that element connected to two other nodes\n";
				Exit_prog = true;
				return false;
			}
			else
			{
				if (value > 0)
				{
					J->Pve_node = node->index;
					J->I *= -1;
				}
				else
					J->Nve_node = node->index;
				J->P_to_next = node;
				node->Branches[node->Branches_counter] = J;
				(node->Branches_counter)++;
			}
		}

		break;
	}

	case 'R':
	{
		if (R_indexes[index] == NULL)
		{
			Element *R = new Element;
			R->P_to_next = NULL;
			R->P_to_previous = NULL;
			R->I_in_Element = NULL;
			R->R = value;
			R->type = check;
			R->index = index;
			R->P_to_previous = node;
			node->Branches[node->Branches_counter] = R;
			(node->Branches_counter)++;
			R_indexes[index] = R;
		}
		else
		{
			Element *R = R_indexes[index];
			if (value != R->R)
			{
				cout << "Error: You entered that element before with another value\n";
				Exit_prog = true;
				return false;
			}
			else if (R->P_to_previous != NULL && R->P_to_next != NULL)
			{
				cout << "Error: You entered that element connected to two other nodes\n";
				Exit_prog = true;
				return false;
			}
			else
			{
				R->P_to_next = node;
				node->Branches[node->Branches_counter] = R;
				(node->Branches_counter)++;
			}
		}

		break;
	}
	}
}


bool Read_elements_node(Node *& node, int index)
{
	char check;
	cin >> check;
	if (check == 'k')	return false;
	node = new Node;
	node->Branches_counter = 1;
	node->Is_essential = false;
	node->index = index;
	do
	{
		if (!Read_line(node, check))
			return false;
		cin >> check;
	} while (check != 'k');
	node->Branches_counter--;
	if (node->Branches_counter > 2)	node->Is_essential = true;
	return true;
}


void Make_Equations(Node **node, int num_essential_nodes, int num_nodes, bool is_one_loop)
{
	double **coefficients = new double*[num_essential_nodes - 1];
	for (int i = 0; i <= num_essential_nodes - 1; i++)
		coefficients[i] = new double[num_essential_nodes];

	Current **I_Branches = new Current*[100];
	int I_Branches_counter = 0;

	for (int i = 0; i<100; i++)
		I_Branches[i] = new Current;

	
	for (int i = 1; R_indexes[i] != NULL; i++)
	{
		R_indexes[i]->I_in_Element = NULL;
	}
	int i = 1;
	while (E_indexes[i] != NULL)
	{
		E_indexes[i]->I_in_Element = NULL;
		i++;
	}
	i = 1;
	while (J_indexes[i] != NULL)
	{
		J_indexes[i]->I_in_Element = NULL;
		i++;
	}

	for (int i = 1; i <= num_nodes; i++)
	{
		node[i]->is_calculated = false;
	}

	for (int i = 0; i < num_essential_nodes - 1; i++)
	{
		for (int j = 0; j < num_essential_nodes; j++)
		{
			coefficients[i][j] = 0;
		}
	}


	int index_refrence_node;
	for (int i = 1; i <= num_nodes; i++)
	{
		if (node[i]->Is_essential == true)
		{
			node[i]->voltage = 0;
			index_refrence_node = i;
			break;
		}

	}

	if (is_one_loop)
		index_refrence_node = 0;

	bool I_NonRepeated;
	for (int i = index_refrence_node + 1; i <= num_nodes; i++)
	{

		if (node[i]->Is_essential)
		{
			for (int b = 1; b <= node[i]->Branches_counter; b++)
			{
				switch (node[i]->Branches[b]->type)
				{
				case 'R':
				{

					Element *R = node[i]->Branches[b];
					double sum_R = R->R;
					double sum_E = 0;
					if (R->I_in_Element == NULL)
						I_NonRepeated = true;
					else
						I_NonRepeated = false;

					if (I_NonRepeated)
					{
						R->I_in_Element = I_Branches[I_Branches_counter];
						I_Branches[I_Branches_counter]->previous = node[i];
					}

					Node *second_node;
					if (R->P_to_next != node[i])
						second_node = R->P_to_next;
					else
						second_node = R->P_to_previous;

					if (second_node->Is_essential)
					{
						coefficients[node[i]->index_essential - 2][node[i]->index_essential - 2] += 1 / sum_R;
						if (second_node->index != index_refrence_node)
							coefficients[node[i]->index_essential - 2][second_node->index_essential - 2] -= 1 / sum_R;
						if (I_NonRepeated)
						{
							I_Branches[I_Branches_counter]->sum_R = sum_R;
							I_Branches[I_Branches_counter]->sum_E = 0;
							I_Branches[I_Branches_counter]->next = second_node;
							I_Branches_counter++;

						}
					}
					else
					{
						Element * next_element = R;
						Node *to_next_node = second_node;
						do
						{
							Node* node = to_next_node;


							if (node->Branches[1] != next_element)
								next_element = node->Branches[1];
							else
								next_element = node->Branches[2];

							switch (next_element->type)
							{
							case 'R':
							{
								Element *R = next_element;
								sum_R += R->R;
								if (R->I_in_Element == NULL)
									I_NonRepeated = true;
								else
									I_NonRepeated = false;

								if (I_NonRepeated) R->I_in_Element = I_Branches[I_Branches_counter];
								if (R->P_to_previous != node)
									to_next_node = R->P_to_previous;
								else
									to_next_node = R->P_to_next;
								break;
							}
							case 'E':
							{
								Element *E = next_element;
								if (E->I_in_Element == NULL)
									I_NonRepeated = true;
								else
									I_NonRepeated = false;

								if (I_NonRepeated) E->I_in_Element = I_Branches[I_Branches_counter];
								if (to_next_node->index == E->Pve_node)
									sum_E += E->V;
								else
									sum_E -= E->V;


								if (E->P_to_previous != node)
									to_next_node = E->P_to_previous;
								else
									to_next_node = E->P_to_next;
								break;
							}
							}
						} while (!to_next_node->Is_essential);


						if (I_NonRepeated)
						{
							I_Branches[I_Branches_counter]->sum_R = sum_R;
							I_Branches[I_Branches_counter]->sum_E = sum_E;
							I_Branches[I_Branches_counter]->next = to_next_node;
							I_Branches_counter++;
						}

						if (is_one_loop)
						{
							solve_one_loop(-1 * sum_E, sum_R);
							GetVoltage_NE(node, num_nodes, is_one_loop);
							return;
						}

						coefficients[node[i]->index_essential - 2][node[i]->index_essential - 2] += 1 / sum_R;
						if (to_next_node->index != index_refrence_node)
							coefficients[node[i]->index_essential - 2][to_next_node->index_essential - 2] -= 1 / sum_R;
						coefficients[node[i]->index_essential - 2][num_essential_nodes - 1] += sum_E / sum_R;
					}

					break;
				}

				case 'E':
				{
					Element *E = node[i]->Branches[b];
					double sum_R = 0;
					double sum_E;
					if (E->I_in_Element == NULL)
						I_NonRepeated = true;
					else
						I_NonRepeated = false;

					if (I_NonRepeated)
					{
						E->I_in_Element = I_Branches[I_Branches_counter];
						I_Branches[I_Branches_counter]->previous = node[i];
					}
					if (node[i]->index == E->Pve_node)
						sum_E = E->V;
					else
						sum_E = -1 * E->V;

					Node *second_node;
					if (E->P_to_next != node[i])
						second_node = E->P_to_next;
					else
						second_node = E->P_to_previous;

					Node *to_next_node = second_node;
					Element * next_element = E;
					do
					{
						Node* node = to_next_node;


						if (node->Branches[1] != next_element)
							next_element = node->Branches[1];
						else
							next_element = node->Branches[2];

						switch (next_element->type)
						{
						case 'R':
						{
							Element *R = next_element;
							sum_R += R->R;
							if (R->I_in_Element == NULL)
								I_NonRepeated = true;
							else
								I_NonRepeated = false;

							if (I_NonRepeated) R->I_in_Element = I_Branches[I_Branches_counter];
							if (R->P_to_previous != node)
								to_next_node = R->P_to_previous;
							else
								to_next_node = R->P_to_next;
							break;
						}
						case 'E':
						{
							Element *E = next_element;
							if (E->I_in_Element == NULL)
								I_NonRepeated = true;
							else
								I_NonRepeated = false;

							if (I_NonRepeated) E->I_in_Element = I_Branches[I_Branches_counter];
							if (to_next_node->index == E->Pve_node)
								sum_E += E->V;
							else
								sum_E -= E->V;


							if (E->P_to_previous != node)
								to_next_node = E->P_to_previous;
							else
								to_next_node = E->P_to_next;
							break;
						}
						}
					} while (!to_next_node->Is_essential);

					if (I_NonRepeated)
					{
						I_Branches[I_Branches_counter]->sum_E = sum_E;
						I_Branches[I_Branches_counter]->sum_R = sum_R;
						I_Branches[I_Branches_counter]->next = to_next_node;
						I_Branches_counter++;
					}


					if (is_one_loop)
					{
						solve_one_loop(-1 * sum_E, sum_R);
						GetVoltage_NE(node, num_nodes, is_one_loop);
						return;
					}

					coefficients[node[i]->index_essential - 2][node[i]->index_essential - 2] += 1 / sum_R;
					if (to_next_node->index != index_refrence_node)
						coefficients[node[i]->index_essential - 2][to_next_node->index_essential - 2] -= 1 / sum_R;
					coefficients[node[i]->index_essential - 2][num_essential_nodes - 1] += sum_E / sum_R;
					break;
				}

				case 'J':
				{
					if (node[i]->index == node[i]->Branches[b]->Nve_node)
						coefficients[node[i]->index_essential - 2][num_essential_nodes - 1] -= node[i]->Branches[b]->I;
					else
						coefficients[node[i]->index_essential - 2][num_essential_nodes - 1] += node[i]->Branches[b]->I;
					break;
				}

				}
			}

		}

	}

	solve_equations(coefficients, num_nodes, num_essential_nodes, node);

	get_currents(I_Branches, I_Branches_counter);

	GetVoltage_NE(node, num_nodes, is_one_loop);

	return;
}


void solve_equations(double **coefficients, int num_nodes, int num_essential_nodes, Node** nodes)
{
	double x, y;
	int n, i, j, k;
	bool test1 = false;
	bool test2 = false;


	double** nodes_voltage = new double*[num_essential_nodes - 1];
	for (int i = 0; i <= num_essential_nodes - 1; i++)
		nodes_voltage[i] = new double[num_essential_nodes];


	for (int a = 0; a < num_essential_nodes - 1; a++)
	{
		for (int b = 0; b < num_essential_nodes; b++)
		{
			nodes_voltage[a][b] = coefficients[a][b];
		}
	}


	n = num_essential_nodes - 1;

	for (k = 0; k<n; k++)
	{
		if (nodes_voltage[k][k] != 0)
			x = nodes_voltage[k][k];
		else
		{
			for (j = 0; j <n + 1; j++)
				for (i = 0; i<n; i++)
					if (i != k)    nodes_voltage[k][j] += nodes_voltage[i][j];
			x = nodes_voltage[k][k];
		}
		for (j = 0; j < n + 1; j++)
			nodes_voltage[k][j] /= x;
		for (i = 0; i<n; i++)
			if (i != k)
			{
				y = nodes_voltage[i][k];
				for (j = 0; j <n + 1; j++)
					nodes_voltage[i][j] -= y*nodes_voltage[k][j];
			}
	}
	int q;
	for (q = 1; q <= num_nodes; q++)
	{
		if (nodes[q]->Is_essential)
		{
			break;
		}
	}

	int num_row = 0;
	for (int w = q + 1; w <= num_nodes; w++)
	{
		if (nodes[w]->Is_essential)
		{
			nodes[w]->voltage = nodes_voltage[num_row][num_essential_nodes - 1];
			num_row++;
		}
	}
}


void get_currents(Current **I_Branches, int I_Branches_counter)
{
	for (int i = 0; i < I_Branches_counter; i++)
	{
		I_Branches[i]->I = (I_Branches[i]->previous->voltage - I_Branches[i]->next->voltage - I_Branches[i]->sum_E) / I_Branches[i]->sum_R;
	}
}



void GetVoltage_NE(Node *nodes[], int num_nodes, bool is_one_loop)
{
	for (int node_index = 1; node_index <= num_nodes; node_index++) // traverse all essential
	{
		if (nodes[node_index]->Is_essential)
		{
			for (int branch_index = 1; branch_index <= nodes[node_index]->Branches_counter; branch_index++)
			{
				// special case: if the branch have current source, dont
				if (nodes[node_index]->Branches[branch_index]->type == 'J') continue;

				double I_branch = nodes[node_index]->Branches[branch_index]->I_in_Element->I;
				if (!is_one_loop && nodes[node_index] != nodes[node_index]->Branches[branch_index]->I_in_Element->previous)
					I_branch = -1 * I_branch;

				// first element in the essential branch
				Element *brnch_elmnt = nodes[node_index]->Branches[branch_index];
				// first non_essential node in the branch to be calculated
				Node *non_essential = brnch_elmnt->P_to_next != nodes[node_index] ? brnch_elmnt->P_to_next : brnch_elmnt->P_to_previous;
				if (non_essential->is_calculated)
					continue;

				// get right essential node
				Element *E_Right = brnch_elmnt;
				Node *Right_Essential_node = non_essential;

				// get all non essential node voltage in the branch
				GetVoltage_NE_(nodes[node_index], non_essential, brnch_elmnt, I_branch);
			}
		}
	}
}



void GetVoltage_NE_(Node *essential, Node *pre, Element *E, double I)
{
	if (pre->Is_essential)	return;

	double V = E->V;
	if (essential->index == E->Nve_node)
		V *= -1;

	if (E->type == 'R')
		V = E->R *I;

	pre->voltage = essential->voltage - V;
	pre->is_calculated = true;

	Element *next_element = pre->Branches[2];
	if (pre->Branches[1] != E)
		next_element = pre->Branches[1];

	Node *next_node = next_element->P_to_previous;
	if (next_element->P_to_next != pre)
		next_node = next_element->P_to_next;

	GetVoltage_NE_(pre, next_node, next_element, I);

}



void OhmLaw(bool is_one_loop)
{
	int i = 1;
	while (R_indexes[i] != NULL)
	{
		Calculate_For_R(i);
		i++;
	}
	i = 1;
	while (E_indexes[i] != NULL)
	{
		Calculate_For_E(i, is_one_loop);
		i++;
	}
	i = 1;
	while (J_indexes[i] != NULL)
	{
		Calculate_For_J(i);
		i++;
	}

}



void Calculate_For_R(int i)
{
	R_indexes[i]->V = (R_indexes[i]->P_to_previous->voltage) - (R_indexes[i]->P_to_next->voltage);
	if (R_indexes[i]->V < 0)
		R_indexes[i]->V*-1;
	R_indexes[i]->I = R_indexes[i]->V / R_indexes[i]->R;
	R_indexes[i]->P = R_indexes[i]->I * R_indexes[i]->V;
}



void Calculate_For_E(int i, bool is_one_loop)
{
	E_indexes[i]->R = 0;
	E_indexes[i]->I = abs(E_indexes[i]->I_in_Element->I);
	Node* pos_node = E_indexes[i]->P_to_next;

	if (pos_node->index != E_indexes[i]->Pve_node)
		pos_node = E_indexes[i]->P_to_previous;

	int index_nearest_essential = get_nearest_essential(E_indexes[i], pos_node, is_one_loop);

	if (is_one_loop && ((index_nearest_essential == 1 && E_indexes[i]->I_in_Element->I > 0) || (index_nearest_essential == 2 && E_indexes[i]->I_in_Element->I < 0)))
		E_indexes[i]->I *= -1;

	if (!is_one_loop && ((index_nearest_essential != E_indexes[i]->I_in_Element->next->index && E_indexes[i]->I_in_Element->I > 0) || (index_nearest_essential == E_indexes[i]->I_in_Element->next->index && E_indexes[i]->I_in_Element->I < 0)))
		E_indexes[i]->I *= -1;


	E_indexes[i]->P = E_indexes[i]->I * E_indexes[i]->V;
}



void Calculate_For_J(int i)
{
	J_indexes[i]->R = 0;
	Node* pos_node = J_indexes[i]->P_to_next;
	Node* neg_node = J_indexes[i]->P_to_previous;
	if (pos_node->index != J_indexes[i]->Pve_node)
	{
		pos_node = J_indexes[i]->P_to_previous;
		neg_node = J_indexes[i]->P_to_next;
	}
	J_indexes[i]->V = pos_node->voltage - neg_node->voltage;
	J_indexes[i]->P = J_indexes[i]->I * J_indexes[i]->V;
}



void Balancecheck(double &P_dissipated, double &P_supplied)
{
	P_dissipated = 0;
	P_supplied = 0;
	int i = 1;
	while (R_indexes[i] != NULL)
	{
		P_dissipated += R_indexes[i]->P;
		i++;
	}
	i = 1;
	while (E_indexes[i] != NULL)
	{
		if (E_indexes[i]->P > 0)
			P_supplied += E_indexes[i]->P;
		else
			P_dissipated -= E_indexes[i]->P;
		i++;
	}
	i = 1;
	while (J_indexes[i] != NULL)
	{
		if (J_indexes[i]->P > 0)
			P_supplied += J_indexes[i]->P;
		else
			P_dissipated -= J_indexes[i]->P;
		i++;
	}

}



double get_In(Node**nodes, Element* R, int num_essential, int num_nodes, bool is_one_loop)
{
	double R_original = R->R;
	R->R = 0.0000000001;

	for (int i = 1; i <= num_nodes; i++)
	{
		nodes[i]->is_calculated = false;
	}

	for (int i = 1; R_indexes[i] != NULL; i++)
	{
		R_indexes[i]->I_in_Element = NULL;
	}
	for (int i = 1; E_indexes[i] != NULL; i++)
	{
		E_indexes[i]->I_in_Element = NULL;
	}
	for (int i = 1; J_indexes[i] != NULL; i++)
	{
		J_indexes[i]->I_in_Element = NULL;
	}

	Make_Equations(nodes, num_essential, num_nodes, is_one_loop);
	R->R = R_original;
	double Inor = R->I_in_Element->I;
	Make_Equations(nodes, num_essential, num_nodes, is_one_loop);
	OhmLaw(is_one_loop);
	return abs(Inor);
}



int get_nearest_essential(Element* element, Node* node, bool is_one_loop)
{
	if (node->Is_essential)
	{
		if (is_one_loop)
			if (node->Branches[2] == element)
				return 2;
			else
				return 1;

		return node->index;
	}

	Element* next_element = node->Branches[1];
	if (node->Branches[1] == element)
		next_element = node->Branches[2];
	Node* next_node = next_element->P_to_next;
	if (next_element->P_to_next == node)
		next_node = next_element->P_to_previous;
	return get_nearest_essential(next_element, next_node, is_one_loop);
}



void disable_sources(Element*onlysource)
{
	for (int i = 1; R_indexes[i] != NULL; i++)
	{
		R_indexes[i]->I_in_Element = NULL;
	}
	int i = 1;
	while (E_indexes[i] != NULL)
	{
		E_indexes[i]->original = E_indexes[i]->V;
		if (E_indexes[i] != onlysource)
			E_indexes[i]->V = 0;
		E_indexes[i]->I_in_Element = NULL;
		i++;
	}
	i = 1;
	while (J_indexes[i] != NULL)
	{
		J_indexes[i]->original = J_indexes[i]->I;
		if (J_indexes[i] != onlysource)
			J_indexes[i]->I = 0;
		J_indexes[i]->I_in_Element = NULL;
		i++;
	}

}



void enable_sources()
{
	int i = 1;
	while (E_indexes[i] != NULL)
	{
		E_indexes[i]->V = E_indexes[i]->original;
		i++;
	}
	i = 1;
	while (J_indexes[i] != NULL)
	{
		J_indexes[i]->I = J_indexes[i]->original;
		i++;
	}
}



double get_Rth(Node **node, int num_essential_nodes, int num_nodes, Element* need_Rth, bool is_one_loop)
{

	Node* new_node = new Node;
	new_node->Is_essential = false;
	num_nodes++;
	new_node->index = num_nodes;
	node[num_nodes] = new_node;

	Element *V = new Element;
	V->P_to_next = need_Rth->P_to_next;
	V->P_to_previous = new_node;
	V->I_in_Element = NULL;
	V->V = 1;
	V->type = 'E';
	V->Pve_node = num_nodes;
	V->Nve_node = need_Rth->P_to_next->index;
	int i;
	for (i = 1; i <= need_Rth->P_to_next->Branches_counter; i++)
	{
		if (need_Rth->P_to_next->Branches[i] == need_Rth)
		{
			need_Rth->P_to_next->Branches[i] = V;
			break;
		}
	}

	int ind_new_E = 1;
	while (E_indexes[ind_new_E] != NULL)
		ind_new_E++;
	E_indexes[ind_new_E] = V;
	V->index = ind_new_E;

	new_node->Branches[1] = need_Rth;
	new_node->Branches[2] = V;
	new_node->Branches_counter = 2;

	need_Rth->original = need_Rth->R;
	need_Rth->P_to_next = new_node;
	need_Rth->R = 1;

	disable_sources(V);
	for (int i = 1; i <= num_nodes; i++)
	{
		node[i]->is_calculated = false;
	}
	Make_Equations(node, num_essential_nodes, num_nodes, is_one_loop);
	double Rth_needded = (1 / abs(need_Rth->I_in_Element->I)) - 1;
	enable_sources();

	need_Rth->R = need_Rth->original;
	need_Rth->P_to_next = V->P_to_next;
	need_Rth->P_to_next->Branches[i] = need_Rth;
	E_indexes[ind_new_E] = NULL;

	delete V;
	delete new_node;
	return Rth_needded;
}



double P_max(double In, double Rth)
{
	return In*In*Rth / 4;
}



void solve_one_loop(double sum_E, double sum_R)
{
	double I = sum_E / sum_R;
	int i = 1;
	while (R_indexes[i] != NULL)
	{
		R_indexes[i]->I_in_Element->I = I;
		i++;
	}
	i = 1;
	while (E_indexes[i] != NULL)
	{
		E_indexes[i]->I_in_Element->I = I;
		i++;
	}
}



double super_position(char required, Element* element, Element* source, Node** nodes, int num_essential_nodes, int num_nodes, bool is_one_loop, int V1, int V2)
{
	disable_sources(source);
	Make_Equations(nodes, num_essential_nodes, num_nodes, is_one_loop);
	OhmLaw(is_one_loop);
	double returned;
	switch (required)
	{
	case 'I':
		if (element->type == 'J')
			returned = element->I;
		else
			returned = element->I_in_Element->I;
		break;
	case 'V':
		returned = nodes[V1]->voltage - nodes[V2]->voltage;
		break;
	}
	enable_sources();
	Make_Equations(nodes, num_essential_nodes, num_nodes, is_one_loop);
	OhmLaw(is_one_loop);
	return returned;
}



void check_all_elements()
{

	int x = 100;
	for (int i = 1; i<100; i++)
	{
		if ((R_indexes[i] != NULL) && (R_indexes[i]->P_to_next == NULL || R_indexes[i]->P_to_previous == NULL))
		{
			cout << "Error: You entered R" << i << " only once\n";
			Exit_prog = true;
		}
		if (R_indexes[i] == NULL && x == 100)
			x = i;
		if (R_indexes[i] != NULL && i>x)
		{
			cout << "Error: You entered R" << i << " without enter R" << x << endl;
			Exit_prog = true;
			break;
		}
	}

	x = 100;
	for (int i = 1; i< 100; i++)
	{
		if ((E_indexes[i] != NULL) && (E_indexes[i]->P_to_next == NULL || E_indexes[i]->P_to_previous == NULL))
		{
			cout << "Error: You entered E" << i << " only once\n";
			Exit_prog = true;
		}
		if (E_indexes[i] == NULL && x == 100)
			x = i;
		if (E_indexes[i] != NULL && i>x)
		{
			cout << "Error: You entered E" << i << " without enter E" << x << endl;
			Exit_prog = true;
			break;
		}
		if ((E_indexes[i] != NULL) && (E_indexes[i]->P_to_next->Is_essential && E_indexes[i]->P_to_previous->Is_essential))
		{
			cout << "E" << i << " is an ideal source\n\n";
			Exit_prog = true;
			break;
		}
	}

	x = 100;
	for (int i = 1; i<100; i++)
	{
		if ((J_indexes[i] != NULL) && (J_indexes[i]->P_to_next == NULL || J_indexes[i]->P_to_previous == NULL))
		{
			cout << "Error: You entered J" << i << " only once\n";
			Exit_prog = true;
		}
		if (J_indexes[i] == NULL && x == 100)
			x = i;
		if (J_indexes[i] != NULL && i>x)
		{
			cout << "Error: You entered J" << i << " without enter J" << x << endl;
			Exit_prog = true;
			break;
		}
		if ((J_indexes[i] != NULL) && (!J_indexes[i]->P_to_next->Is_essential || !J_indexes[i]->P_to_previous->Is_essential))
		{
			cout << "J" << i << " is an ideal source\n\n";
			Exit_prog = true;
			break;
		}
	}

}



void print_outputs(Node ** nodes, int num_essential_nodes, int num_nodes, bool is_one_loop)
{
	string line;
	cin.ignore();
	while (1)
	{
		cout << "Enter the requrements or press k to exit\n\n";
		getline(cin, line);
		if (line.length() == 0)
			continue;
		if (line.length() == 1)
			break;

		int k = 0;
		while (k < line.length()-1)
		{
			if (line[k] == ' ' && line[k + 1] == ' ')
				line.erase(k+1, 1);
			else
				k++;
		}

		if (line[k] == ' ')
			line.erase(k, 1);
		if (line[0] == ' ')
			line.erase(0, 1);

		string new_line = line;
		int i = 0;
		while ( i < new_line.length())
		{
			if (new_line[i] == ' ')
				new_line.erase(i, 1);
			else
				i++;
		}

		bool for_super_position = false;

		for  (i = 1; i < new_line.length(); i++)
		{
			if (new_line[i] == 'E' || new_line[i] == 'J')
			{
				for_super_position = true;
				break;
			}
		}

		if ((new_line.length() < 5) || (new_line[0] == 'V' && !for_super_position))
		{
			if (line[0] == 'I')
			{

				int i = 3;
				while (line[i] == ' ')
					i++;
				string sindex = line.substr(i, line.length() - i);
				int index = stoi(sindex);
				switch (line[2])
				{
				case 'R':
				{
					if (R_indexes[index] == NULL)
					{
						cout << "You have entered an index that isn't valid\n\n";
						continue;
					}
					if (abs(R_indexes[index]->I) < 1E-15)
						R_indexes[index]->I = 0;
					cout << "I in R" << index << " = " << abs(R_indexes[index]->I) << endl << endl;
					break;
				}
				case 'E':
				{
					if (E_indexes[index] == NULL)
					{
						cout << "You have entered an index that isn't valid\n\n";
						continue;
					}
					if (abs(E_indexes[index]->I) < 1E-15)
						E_indexes[index]->I = 0;
					cout << "I in E" << index << " = " << abs(E_indexes[index]->I) << endl << endl;
					break;
				}
				case 'J':
				{
					if (J_indexes[index] == NULL)
					{
						cout << "You have entered an index that isn't valid\n\n";
						continue;
					}
					cout << "I in J" << index << " = " << abs(J_indexes[index]->I) << endl << endl;
					break;
				}

				}
			}
			else if (line[0] == 'V')
			{
				i = 2;
				while (line[i] != ' ' &&  i < line.length())
					i++;

				int index1 = stoi(line.substr(2, i - 2));

				int j = i + 1;
				while (line[j] != ' ' && j < line.length())
					j++;

				int index2 = stoi(line.substr(i + 1, j - i - 1));

				if (nodes[index1] == NULL || nodes[index2] == NULL)
				{
					cout << "You have entered an index that is not valid\n\n";
					continue;
				}
				if (abs(nodes[index1]->voltage - nodes[index2]->voltage) < 1E-15)
					nodes[index1]->voltage = nodes[index2]->voltage;
				cout << "V between " << index1 << " and " << index2 << " = " << nodes[index1]->voltage - nodes[index2]->voltage << endl << endl;
			}
			else if (line[0] == 'P')
			{
				int i = 3;
				while (line[i] == ' ')
					i++;
				string sindex = line.substr(i, line.length() - i);
				int index = stoi(sindex);

				switch (line[2])
				{
				case 'R':
				{
					if (R_indexes[index] == NULL)
					{
						cout << "You have entered an index that isn't valid\n\n";
						continue;
					}
					if (abs(R_indexes[index]->P) < 1E-15)
						R_indexes[index]->P = 0;
					cout << "P in R" << index << " = " << abs(R_indexes[index]->P) << endl << endl;
					break;
				}
				case 'E':
				{
					if (E_indexes[index] == NULL)
					{
						cout << "You have entered an index that isn't valid\n\n";
						continue;
					}
					if (abs(E_indexes[index]->P) < 1E-15)
						E_indexes[index]->P = 0;
					cout << "P in E" << index << " = " << E_indexes[index]->P << endl << endl;
					break;
				}
				case 'J':
				{
					if (J_indexes[index] == NULL)
					{
						cout << "You have entered an index that isn't valid\n\n";
						continue;
					}
					if (abs(J_indexes[index]->P) < 1E-15)
						J_indexes[index]->P = 0;
					cout << "P in J" << index << " = " << J_indexes[index]->P << endl << endl;
					break;
				}


				}
			}
			else
			{
				int i = 3;
				while (line[i] == ' ')
					i++;
				string sindex = line.substr(i, line.length() - i);
				int index = stoi(sindex);

				if (R_indexes[index] == NULL)
				{
					cout << "You have entered an index that isn't valid\n\n";
					continue;
				}

				double Rmax = get_Rth(nodes, num_essential_nodes, num_nodes, R_indexes[index], is_one_loop);
				double In = get_In(nodes, R_indexes[index], num_essential_nodes, num_nodes, is_one_loop);
				double Pmax = P_max(In, Rmax);
				if (abs(Pmax) < 1E-15)
					Pmax = 0;
				cout << "R max = " << Rmax << "\tP max = " << Pmax << endl << endl;
			}
		}
		else
		{
			switch (line[0])
			{
			case 'I':
			{
				Element* element = NULL;
				string sources;
				double I = 0;

				int i = 4;
				while (line[i] != 'E' && line[i] != 'J')
					i++;

				sources = line.substr(i, line.length() - i);
				line.erase(i-1, line.length() - i+1);

				i = 3;
				while (line[i] == ' ')
					i++;
				string sindex = line.substr(i, line.length() - i);
				int index = stoi(sindex);

				switch (line[2])
				{
				case 'R':
				{
					if (R_indexes[index] == NULL)
					{
						cout << "You have entered an index that isn't valid\n\n";
						continue;
					}
					element = R_indexes[index];
					break;
				}
				case 'E':
				{
					if (E_indexes[index] == NULL)
					{
						cout << "You have entered an index that isn't valid\n\n";
						continue;
					}
					element = E_indexes[index];
					break;
				}
				case 'J':
				{
					if (J_indexes[index] == NULL)
					{
						cout << "You have entered an index that isn't valid\n\n";
						continue;
					}
					element = J_indexes[index];
					break;
				}

				}
				bool not_valid = false;

				while (sources.length() != 0)
				{
					char source_char = sources.at(0);

					i = 1;
					while (sources[i] == ' ' && i< sources.length())
						i++;

					int j = i;
					while (sources[j] != ' ' && j< sources.length())
						j++;


					int index = stoi(sources.substr(i, j - i));
					Element* source;

					
					if (source_char == 'E')
					{ 
						if (E_indexes[index] == NULL)
						{
							cout << "You have entered an index that isn't valid\n\n";
							not_valid = true;
							break;
						}
						source = E_indexes[index];
					}
					else
					{
						if (J_indexes[index] == NULL)
						{
							cout << "You have entered an index that isn't valid\n\n";
							not_valid = true;
							break;
						}
						source = J_indexes[index];
					}

					

					I += super_position('I', element, source, nodes, num_essential_nodes, num_nodes, is_one_loop, 0, 0);


					sources.erase(0, j+1);
				}

				if (not_valid)
					continue;

				if (abs(I) < 1E-15)
					I = 0;

				cout << "I = " << abs(I) << endl << endl;

				break;
			}
			case 'V':
			{
				Element* element = NULL;
				string sources;
				double V = 0;

				int i = 4;
				while (line[i] != 'E' && line[i] != 'J')
					i++;

				sources = line.substr(i, line.length() - i);
				line.erase(i - 1, line.length() - i + 1);

				i = 2;
				while (line[i] != ' ' &&  i < line.length())
					i++;

				int index1 = stoi(line.substr(2, i-2));

				int j = i+1;
				while (line[j] != ' ' && j < line.length())
					j++;

				int index2 = stoi(line.substr(i+1, j-i-1));

				if (nodes[index1] == NULL || nodes[index2] == NULL)
				{
					cout << "You have entered an index that isn't valid\n\n";
					continue;
				}

				bool not_valid = false;

				while (sources.length() != 0 )
				{
					char source_char = sources.at(0);

					i = 1;
					while (sources[i] == ' ')
						i++;

					int j = i;
					while (sources[j] != ' ' && j< sources.length())
						j++;


					int index = stoi(sources.substr(i, j - i));
					Element* source;

					
					if (source_char == 'E')
					{
						if (E_indexes[index] == NULL)
						{
							cout << "You have entered an index that isn't valid\n\n";
							not_valid = true;
							break;
						}
						source = E_indexes[index];
					}
					else
					{
						if (J_indexes[index] == NULL)
						{
							cout << "You have entered an index that isn't valid\n\n";
							not_valid = true;
							break;
						}
						source = J_indexes[index];
					}

					

					V += super_position('V', element, source, nodes, num_essential_nodes, num_nodes, is_one_loop, index1, index2);
					sources.erase(0, j + 1);
				}

				if (not_valid)
					continue;

				if (abs(V) < 1E-15)
					V = 0;

				cout << "V = " << V << endl << endl;

				break;
			}

			}
		}
	}
	
}



void frist_print()
{
	HANDLE hStdout = GetStdHandle(STD_OUTPUT_HANDLE);
	SetConsoleTextAttribute(hStdout, FOREGROUND_RED | FOREGROUND_GREEN |/* FOREGROUND_BLUE |*/ FOREGROUND_INTENSITY);
	cout << "\t \t \t Team 3 CMP2020 \nplease enter the your Circuits \nthe inputs should be by SI units \n";
	cout << "please enter the elements -bettry,resetance,current source- for each node \nbetween each node press k \nthe inputs are: \n\n";
	SetConsoleTextAttribute(hStdout, FOREGROUND_RED | FOREGROUND_GREEN | FOREGROUND_BLUE /*| FOREGROUND_INTENSITY*/);

}



void second_print()
{
	HANDLE hStdout = GetStdHandle(STD_OUTPUT_HANDLE);
	SetConsoleTextAttribute(hStdout, FOREGROUND_RED | FOREGROUND_GREEN |/* FOREGROUND_BLUE |*/ FOREGROUND_INTENSITY);
	cout << "the circuit is balanced\nIf you want current through an element write I then a space then write the element\nIf you want voltage betwen two nodes write V then space then the index of the two nodes separated by a space\n";
	cout << "If you want max resistance such that receive max power write M then a space then the resistance\nFor superposition just after the below add the wanted source or sources\n\n";
	SetConsoleTextAttribute(hStdout, FOREGROUND_RED | FOREGROUND_GREEN | FOREGROUND_BLUE /*| FOREGROUND_INTENSITY*/);

}

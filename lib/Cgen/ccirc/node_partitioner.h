/*--------------------------------------------------------------------------*
 * Copyright 2002 by Paul D. Kundarewich, Michael Hutton, Jonathan Rose     *
 * and the University of Toronto. 											*
 * Use is permitted, provided that this attribution is retained  			*
 * and no part of the code is re-distributed or included in any commercial	*
 * product except by written agreement with the above parties.              *
 *                                                                          *
 * For more information, contact us directly:                               *
 *	  Paul D. Kundarewich (paul.kundarewich@utoronto.ca)					*
 *    Jonathan Rose  (jayar@eecg.toronto.edu)                               *
 *    Mike Hutton  (mdhutton@cs.toronto.edu, mdhutton@eecg.toronto.edu)     *
 *    Department of Electrical and Computer Engineering                     *
 *    University of Toronto, 10 King's College Rd.,                         *
 *    Toronto, Ontario, CANADA M5S 1A4                                      *
 *    Phone: (416) 978-6992  Fax: (416) 971-2286                            *
 *--------------------------------------------------------------------------*/



#ifndef node_partitioner_H
#define node_partitioner_H

#include "circ.h"
#include "circuit.h"
#include <fstream>
#include <map>

typedef vector<int> PARTITION_WEIGHTS_TYPE_ALSO;
typedef map<PORT *, int> PORT_PTR_to_INT_MAP;

//
// Class_name NODE_PARTITIONER
//
// Description
//
//	Partitions the nodes into clusters.
//

class NODE_PARTITIONER
{
public:
	enum PARTITION_TYPE {CIRCUIT_TYPE, CLUSTER_TYPE};
	NODE_PARTITIONER();
	NODE_PARTITIONER(const NODE_PARTITIONER & another_node_generator);
	NODE_PARTITIONER & operator=(const NODE_PARTITIONER & another_node_generator);
	~NODE_PARTITIONER();

	COST_TYPE 	partition_circuit(CIRCUIT * circuit);
	COST_TYPE 	partition_cluster(CLUSTER * circuit);
private:
	CIRCUIT * 				m_circuit;
	CLUSTER *				m_cluster;
	PORT_PTR_to_INT_MAP		m_port_index_map;
	PARTITION_TYPE			m_partition_type;

		
	fstream					m_output_file;

	COST_TYPE	partition_graph();
	void fill_in_edge_info(int * start_of_edges_array, int * edge_index);
	void create_port_ptr_to_index_map();
	COST_TYPE get_scaled_cost(const int & number_of_edges, const int & number_partitions,
								const int & number_of_nodes, int * edges_info, int * results);
								
	PARTITION_WEIGHTS_TYPE_ALSO get_partition_weights(int * results, const int & number_of_nodes, const int & number_of_partitions);

	void label_nodes_with_the_partition_they_belong_to(int * results, const int & number_of_nodes,
				    									const int & number_of_partitions);
	int get_size();
	int get_nEdges();
	void construct_clusters(const NUM_ELEMENTS & number_partitions);
	void construct_single_cluster();
};


#endif

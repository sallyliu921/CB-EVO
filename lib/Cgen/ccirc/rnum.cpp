/*-------------------------------------------------------------------------*
 * Copyright 1996 by Michael Hutton, Jonathan Rose and the University of   *
 * Toronto.  Use is permitted, provided that this attribution is retained  *
 * and no part of the code is re-distributed or included in any commercial *
 * product except by written agreement with the above parties.             *
 *                                                                         *
 * For more information, contact us directly:                              *
 *    Mike Hutton  (mdhutton@cs.toronto.edu, mdhutton@eecg.toronto.edu)    *
 *    Jonathan Rose  (jayar@eecg.toronto.edu)                              *
 *    Department of Electrical and Computer Engineering                    *
 *    University of Toronto, 10 King's College Rd.,                        *
 *    Toronto, Ontario, CANADA M5S 1A4                                     *
 *    Phone: (416) 978-6992  Fax: (416) 971-2286                           *
 *-------------------------------------------------------------------------*/
 /* $Revision: by GUANGLEI (Gabriel) ZHOU, July 2021 $ */
/*

(12/7/2021): Remove all the Non-sparse matrix setting (CONFIRM = 0) 
Roadmap: Implement the rnum_of_node -> rnum_entire_graph   
rnum_of_node: doing the quick_calculate , log2(fanin) every node that is the outcone of a node

Some background about quick/full: 

\item [rnum] Output the reconvergence number, R, of the circuit as a whole, and 
		the minimum and maximum reconvergence number of any PI.  The
		complete definition of reconvergence is beyond the scope
		of this document, but it is a quantification of the 
		``amount" of reconvergence in the circuit, as a real number
		between 0 and the $\log_2($k$)$ where $k$ is the maximum
		fanin (usually LUT-size) of any node.  If neither
		rnum=quick or rnum=full is specified, both will be 
		calculated.   Note that quick and full are equivalent for
		combinational circuits, so the full calculation will
		never be done on a circuit without DFF nodes.

\item [rnum=quick] Perform the quicker (O(nLOG * nPI)), but less accurate,
		reconvergence calculation.  

\item [rnum=full] Perform the slow (O($\mbox{nLOG}^2$ * nPI)), 
		but more accurate, reconvergence calculation.  NOTE: this 
		can take a long time, and a lot of memory, as it requires	
		several matrix determinant calculations on matrices
		of size equal to the number of out-cones of an input.  
		For this reason, the `all' option specifies rnum=quick, not rnum=full.


For the interest of my project , we would be primarily targetting the combinational circuit
Only quick mode will be implemented for this integration as quick and full are equivalent for combinational circuits.  

The next step: 
[] look into the ccirc to see the graph definition, and replace the graph definition in here(circ)  

*/






/* $Revision: 3.0 $ */

/*
 *  Calculate the reconvergence number and related values for graph.
 *
 *  M. Hutton, November 1994.
 * 
 *  The reconvergence number is the weighted average, over all PIs x
 *  of log_2(det(K(x)))/conesize(x), where K(x) is the Kirchoff matrix
 *  of G induced by the recursive fanout of x.  The Kirchoff matrix
 *  of a graph is defined by defined by:
 *        K[i,i] = num_fanins(i)
 *        K[i,j] = -1 if (i,j) is an edge in G
 *        K[i,j] = 0 otherwise.
 *  and conesize is the size of the fanout-cone from PI x.
 *
 *  The logdet is the log_2 of the number of spanning out-trees from the
 *  relevant PI.  See the external doc'n for the theory behind this.
 *
 *  Before calulating determinant, need to take minor with respect to
 *  root node, which has fanin 0 and would make the determinant 0.
 *
 *  For connected single-source (x) graphs with no cycles, this can 
 *  be evaluated simply as product_i(log_2(num_fanins(x_i))).  For
 *  graphs with cycles, we have to actually calculate the n^3 size
 *  Kirchoff matrix.  This requires, for any reasonable size graphs, 
 *  a sparse matrix approach to the determinant (or a lot of waiting).
 *  
 *  Option "quick" causes this faster method to be used even if there
 *  are back-edges.
 *  
 *  For more details, read the external documentation. 
 *
 *  In order to test the sparse code, have CONFIRM to use both sparse
 *  and full matrix calculations.  Should normally be off.
 *
 *  The COLLAPSE feature is new -- I re-wrote the marking routines to
 *  identify nodes with a single marked fanin to that fanin, recursively.
 *  This should speed things up significantly.  Will take out the ifdefs
 *  after I have tested it further and done the formal proof.
 */





#include "rnum.h"


/*
 *  Calculate rnum. 
 *  return the reconvergence + its max/min
 */
void
rnum(CIRCUIT * circuit,double *m_R0,double *m_R0max,double *m_R0min)
{
    PORTS PI;
    PORTS::iterator port_iter;
    PORT * port;
    NODE * PI_node;
    double R0min, R0max;
    double R0sum, R0;
    int    R0num;
    double n0;		/* numerators for rnum calc */
    int    d0;		/* numerators for rnum calc */

    R0sum = 0.0;
    R0min = 9999999;
    R0max = 0.0;
    R0num = 0;
	PI = circuit->get_PI();

	Log("Start Calculating the rnum" );

    /* Calculate rnum of each PI; update the max/min R values                */
	// for each PI break cycles connected to it
	for (port_iter = PI.begin(); port_iter != PI.end(); port_iter++)
	{
		port = *port_iter;
        if (port->get_io_direction() != PORT::CLOCK){ /* we're ignoring clocks */
	    /* n0 are logdets,  d0 are conesizes */
	    /* Together they are the numerator, denominators for rnum */       
        PI_node = port->get_edge()->get_sink()->get_my_node(); 
		if(PI_node){
			_rnum_of_node(circuit, PI_node, &n0, &d0);
			if (d0 > 0) {
				R0min = MIN(R0min, n0/d0); 
				R0max = MAX(R0max, n0/d0);
			}
			R0sum += n0; 
			R0num += d0;
			}
		}
	}
    /* Calculate the R values of the entire graph */
    R0 = R0num > 0 ? R0sum / R0num : 0.0;

    *m_R0 = R0; *m_R0min = R0min; *m_R0max = R0max; 

}


/*
 *  Calculate reconvergence from a specified root, specifying
 *  the logdet value and the conesize.
 *
 *  Array marked is owned, throughout the life of rnum(), by rnum itself.
 *  To take advantage of sparsity, we keep marked as a list of marked nodes,
 *  which we sort to do topological traversal of all nodes.  This is 
 *  a log n penalty if the conesize is actually O(graph->size), but it
 *  saves us a pile of time with cones that are really small.
 *
 *  After calculating R for node, have to re-zero marks for all marked nodes.
 */
static void _rnum_of_node(CIRCUIT * circuit, NODE * node,double *n0, int *d0)
{

    *n0 = 0.0;
    *d0 = 0;

    _reset_marks(circuit);
    _mark_outcone(circuit, node);
    _quick_count(circuit, n0, d0);
}



/*
 *  Re-set all marks to 0
 */
static void
_reset_marks(CIRCUIT * circuit)
{
	assert(circuit);

	NODES & nodes = circuit->get_nodes();
	NODES::iterator node_iter;
	NODE * node = 0;

	for (node_iter = nodes.begin(); node_iter != nodes.end(); node_iter++)
	{
		node = *node_iter;
		assert(node);
		node->set_colour(NODE::UNMARKED);
	}
}


/*
 *   Mark the recursive fanout of the original node.
 */
static void _mark_outcone(CIRCUIT * circuit, NODE * node)
{
    int i;
    EDGES m_edges;
    EDGES::const_iterator edge_iter;	
	EDGE * edge = 0;
	NODE * sink_node = 0;

    /* quit on visited already, dummy node or is PO */
    //what is dummy node ?
    if (node->get_colour() == NODE::MARKED_OUTCONE|| node->get_output_port()->get_type() == PORT::PO ){
		return;			/* been here already; */
	}

    node->set_colour(NODE::MARKED_OUTCONE); 		/* mark self before recurse to avoid cycles */
    m_edges = node->get_output_edges();

    for (edge_iter = m_edges.begin(); edge_iter != m_edges.end(); edge_iter++)
	{
		edge = *edge_iter;
		assert(edge);
		sink_node = edge->get_sink_node();

		if (sink_node->get_type() != NODE::SEQ)
		{
			_mark_outcone(circuit, sink_node);
		}
	}
}


/*
 *  For quick rnum calculation, just have to get sum of the log-det of
 *  each node, and the number of nodes.  Can do this by just counting 
 *  marked fanins of marked nodes in a single pass.
 */
static void 
_quick_count(CIRCUIT * circuit, double *n0, int *d0)
{
    int    num, num_fanins;
    double sum;

    sum = 0.0;
    num = 0;
	NODES & nodes = circuit->get_nodes();
	NODES::iterator node_iter;
	NODE * node = 0;

	for (node_iter = nodes.begin(); node_iter != nodes.end(); node_iter++)
	{
        node = *node_iter;
		assert(node);
        if (node->get_type() == NODE::COMB) 
        {
                if (node->get_colour() == NODE::MARKED_OUTCONE)
				{
                num += 1;
                num_fanins = _count_marked_fanin(node);
                sum += num_fanins ? log2( (double) num_fanins) : 0.0;
                }
        }
    }
    *n0 = sum;
    *d0 = num;
}

/*
 *  For quick rnum calculation, just counting marked fanins of a node in a single pass.
 *  A fanin should be only counted if it is coming from the outcone(x)     
 *  Return the marked fanin number of a node  
 */
static int _count_marked_fanin(NODE *node)
{
    EDGES m_edges;
    EDGE * edge = 0;
    NODE * sink_node = 0;
    EDGES::const_iterator edge_iter;	
    int num;
	m_edges = node->get_output_edges();

    num = 0;
	for (edge_iter = m_edges.begin(); edge_iter != m_edges.end(); edge_iter++)
	{
		edge = *edge_iter;
		assert(edge);
		sink_node = edge->get_sink_node();

		if (sink_node->get_colour() == NODE::MARKED_OUTCONE)
		{
			++num;
		}
	}
    return num;
}

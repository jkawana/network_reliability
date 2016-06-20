#include <mpi.h>
#include <climits>
#include <algorithm>
#include "graph.h"
#include <stdlib.h>
#include <time.h>
#include <set>
#include "graph.cpp"

bool Dijkstra(int source, int destination, Graph G, int diameterConstraint);

const double totalTrials = 10000000  ;

int main(int argc, char *argv[])
{

        clock_t startTime = clock(); //do we use this?

	MPI_Init(&argc, &argv);
	int my_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	int num_procs;
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
	
	//to test monte carlo w. dijkstra

	//	double expectedR = .234521; //change accordingly. 
	double expectedR = .997224; //change accordingly. 

	double diameterConstraint = 8; // change accordingly
	double totalNumOfVertices = 25; //change accordingly. 
	int source = 0; 

	//	double destination = 18; //change accordingly
	double destination = 12; //change accordingly

	//To do: read expectedR, diameter, |V|, destination in from file
       
	double reliability;
	double totalSuccess = 0;
	srand(time(NULL)*my_rank);
	double randomNum;
	double sum = 0;
	double start, end, totaltime, outputtime;
	Graph G(totalNumOfVertices);
	G.create();
	
	MPI_Barrier(MPI_COMM_WORLD);
	start = MPI_Wtime();
	
	//looping n number of times
	for ( int i = 0 ; i < totalTrials; i++)
	{
		for ( int j = 0; j < G.getTotalEdges(); j++ )
		{
		        randomNum = rand() % 100 + 1;
			if ( randomNum > G.edge[j].successRate * 100)
				G.edge[j].determined = 0;
			else
				G.edge[j].determined = 1;
		}
		if ( Dijkstra(source, destination, G, diameterConstraint) ) 
		{	totalSuccess++;  
}
	}
	
       	printf("Processor: %d     totalSuccess: %f\n",my_rank,totalSuccess); 

	end = MPI_Wtime();
	totaltime = end - start;
	
	MPI_Reduce(&totaltime, &outputtime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	MPI_Reduce(&totalSuccess, &sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	if (my_rank == 0)
	{
		cout << "Total number of proccessors = " << num_procs << endl;
		cout << "Total number of hits = " << sum << endl;
		cout << "Each Proccessor performs " << totalTrials << "  Trials" << endl;
		cout << "Total number of trials = " << (totalTrials * num_procs) << endl;
		
		reliability = sum / (totalTrials * num_procs);
		
		cout << "Monte Carlo Result: " << endl;
		cout << "5x5grid" << endl; //change
		cout << "Diameter Constraint: " << diameterConstraint << endl;
		cout << "Source Node: " << source << endl;
		cout << "Destination Node: " << destination << endl;
		cout << "Reliability: " << reliability << endl;
		cout << "Expected Reliability: " << expectedR << endl;
		cout << "Difference: " << (reliability - expectedR) << endl;
		cout << "Time: " << outputtime << endl;
	}

	MPI_Barrier(MPI_COMM_WORLD); //what does this do?
	MPI_Finalize();

	return 0;
}

bool Dijkstra(int source, int destination, Graph G, int diameterConstraint)
{
  	
    set<int> T; //the set of verticies used to make sure all verticies distances were
    
	int *parent = new int[G.getTotalVertices()];
	int *lambda = new int[G.getTotalVertices()];
    
    //set up the vector to visit all nodes, T, the parent to have no parents and the min index to be the max
	for (int i = 0; i < G.getTotalVertices(); i++)
	{
        T.insert(i);
		parent[i] = -1; // at the beginning every vertex's parent is -1;
		lambda[i] = INT_MAX; // every vertex is assumed to be far away from s.
	}
    
	lambda[source] = 0; // s is at distance 0 of itself.

    
    int min_index, min_value; //adjNode varibles to find the smallest λ
    node * adjNode;
    
    //continue until T is empty (meanging all distances between the source and each node is visited)
    while (!T.empty()){
        
        min_index = *(T.begin()), min_value = lambda[min_index];
        
        
        //gets the smallest λ
        for (set<int>::iterator i = T.begin(); i != T.end(); ++i)
        {
            if (lambda[*i] < min_value)
            {
                min_value = lambda[*i];
                min_index = *i;
            }
        }
        
        
        //gets the address of the lowest lambda
        adjNode = &G.headnodes[min_index];
        
        //itterate over all adjacent vertecies to the min index
        for (adjNode = adjNode->next; adjNode; adjNode = adjNode->next)
        {
            
            if (abs(G.edge[adjNode->edgeID].determined) == 1 &&
                lambda[min_index] != INT_MAX && T.count(adjNode->vertex) > 0 &&
                (lambda[adjNode->vertex] > lambda[min_index] + adjNode->edgeWeight))
                
            {
                lambda[adjNode->vertex] = lambda[min_index] + adjNode->edgeWeight;
                parent[adjNode->vertex] = min_index;
            }
            
        }
        
        T.erase(min_index); //no longer consider this node for the adjacent list
        
    }

    
	/* Dijkstra Results: */
    //if the graph does connect and the distance between the between the source and destination is in the
    bool works = (parent[destination] != -1 && lambda[destination] <= diameterConstraint);
    
    delete[] parent;
    delete[] lambda;
    return works;
    



}

#include <mpi.h>
#include <climits>
#include <algorithm>
#include "graph.h"
#include <stdlib.h>
#include <time.h>
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
	double expectedR = .234521; //change accordingly. 
	double diameterConstraint = 8; // change accordingly
	double totalNumOfVertices = 25; //change accordingly. 
	int source = 0; 
	double destination = 18; //change accordingly

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
		cout << "Reliability = " << reliability << endl;
		cout << "Expected Reliability" << expectedR << endl;
		cout << "Difference: " << (reliability - expectedR) << endl;
		cout << "Time: " << outputtime << endl;
	}

	MPI_Barrier(MPI_COMM_WORLD); //what does this do?
	MPI_Finalize();

	return 0;
}
bool Dijkstra(int source, int destination, Graph G, int diameterConstraint)
{
	int s = source;
	int u = s; // u is the vertex with smallest Lambda, which at beginning is s.
	bool uIsSet = true; //true because u is set to s.
	int totalVertices = G.getTotalVertices();
	int numOfVerticesToVisit = totalVertices;
	bool *T = new bool[totalVertices];
	int *father = new int[totalVertices];
	int *Lambda = new int[totalVertices];
	for (int i = 0; i < totalVertices; i++)
	{
		T[i] = true; // all the vertices are eligible to be visited;
		father[i] = -1; // at the beginning every vertex's father is -1;
		Lambda[i] = INT_MAX; // every vertex is assumed to be far away from s.
	}
	Lambda[s] = 0; // s is at distance 0 of itself.

	//To do: Change loop structure
	while (1)
	{
		/* set u to the vertex with the smallest Lambda. */
		for (int i = 0; i < totalVertices; i++)
		{
			if (!uIsSet && T[i])
			{
				u = i;	//to set u to the first available vertex.
				uIsSet = true;
			}
			else if (uIsSet && T[i] && Lambda[i] < Lambda[u])
				u = i;	//to set u to the smallest Lambda (aka distance from source).
		}

		node* adjnode = G.headnodes[u].next;
		while (adjnode != NULL)
		{
			if (abs(G.edge[adjnode->edgeID].determined) == 1 &&
				Lambda[u] != INT_MAX && T[adjnode->vertex] &&
				(Lambda[adjnode->vertex] > Lambda[u] + adjnode->edgeWeight))
			{
				Lambda[adjnode->vertex] = Lambda[u] + adjnode->edgeWeight;
				father[adjnode->vertex] = u;
			}
			adjnode = adjnode->next;
		}
		uIsSet = false; // to indicate u is not set.
		T[u] = false; //u is not available to be visited (aka. to delete u from T).
		numOfVerticesToVisit--;
		if (numOfVerticesToVisit == 0)
			break;
	}

	/* Dijkstra Results: */
	int next = father[destination];
	if (next != -1 && Lambda[destination] <= diameterConstraint)
	{
		/* delete dynamically allocated variables. */
		delete[] T;
		delete[] father;
		delete[] Lambda;
		return true;
	}
	else
	{
		/* delete dynamically allocated variables. */
		delete[] T;
		delete[] father;
		delete[] Lambda;
		return false;
	}

}

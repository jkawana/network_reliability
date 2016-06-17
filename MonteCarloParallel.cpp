#include <mpi.h>
#include <iostream>
#include <climits>
#include <algorithm>
#include "Graph.h"
#include "Graph.cpp"
#include <stdlib.h>
#include <time.h>


using namespace std;

bool Dijkstra(int source, int destination, Graph G, int diameterContraint);

const double totalTrails = 1000000000;	//change accordingly.

int main(int argc, char *argv[])
{

	clock_t startTime = clock();

	MPI_Init(&argc, &argv);
	int world_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	/* To test Monte Carlo with Dijkstra: */
	double expectedR = .234521; //change accordingly.
	int diameterContraint = 8;	//change this accordingly.
	int totalNumOfVertices = 25;	//change this accordingly.
	int source = 0;				//change this accordingly.
	int destination = 18;		//change this accordingly.
	double reliability;
	double totalSuccess = 0;
	srand(time(NULL));
	double randomNum;
	double sum = 0;

	double start, end, totaltime, outputtime, outputtimemin;

	Graph G(totalNumOfVertices);
	G.create();

	//Start timing here
	MPI_Barrier(MPI_COMM_WORLD);
	start = MPI_Wtime();

	/* looping n number of times. */
	for (int i = 0; i < totalTrails; i++)
	{

		/* set the determined value. */
		for (int j = 0; j < G.getTotalEdges(); j++)
		{
			randomNum = rand() % 100 + 1;
			if (randomNum > G.edge[j].successRate * 100)
				G.edge[j].determined = 0;
			else
				G.edge[j].determined = 1;
		}
		/* use dijkstra to get totalSuccess aka increment. */
		if (Dijkstra(source, destination, G, diameterContraint))
			totalSuccess++;
	} 	//end loop
	double num = 1;
	double totalNum = 0;
	end = MPI_Wtime();
	totaltime = end - start;

	MPI_Reduce(&totaltime, &outputtime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	MPI_Reduce(&totalSuccess, &sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&num, &totalNum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	if (world_rank == 0)
	{
		cout << "Total number of proccessors = " << totalNum << endl;
		cout << "Each Proccessor performs " << totalTrails << "  Trials" << endl;
		cout << "Total number of trials = " << (totalTrails * totalNum) << endl;
		reliability = sum / (totalTrails * world_size);
		cout << "Monte Carlo Result: " << endl;
		cout << "5x5grid" << endl; //change
		cout << "Diameter Constraint: " << diameterContraint << endl;
		cout << "Source Node: " << source << endl;
		cout << "Destination Node: " << destination << endl;
		cout << "Reliability = " << reliability << endl;
		cout << "Expected Reliability" << expectedR << endl;
		cout << "Difference: " << (reliability - expectedR) << endl;
		cout << "Time: " << outputtime << endl;
	}

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();

	return 0;
}

bool Dijkstra(int source, int destination, Graph G, int diameterContraint)
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
	while (1)
	{
		/* set u to the vertex with the smallest Lambda. */
		for (int i = 0; i < totalVertices; i++)
		{
			if (!uIsSet && T[i])
			{
				u = i;	//to set u to the first availible vertex.
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
		T[u] = false; //u is not availible to be visited (aka. to delete u from T).
		numOfVerticesToVisit--;
		if (numOfVerticesToVisit == 0)
			break;
	}

	/* Dijkstra Results: */
	int next = father[destination];
	if (next != -1 && Lambda[destination] <= diameterContraint)
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


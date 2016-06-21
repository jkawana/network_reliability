#include "graph.h"

Graph::Graph() // default construtor.
{
	totalVertices = 0;
	totalEdges = 0;
	headnodes = NULL;
	edge = NULL;
}
Graph::Graph(int nodes) // construtor.
{
	setVertices(nodes);
}
void Graph::setVertices(int nodes)
{
	totalEdges = 0;
	totalVertices = nodes;
	headnodes = new node[totalVertices]; // ptr to an array of headnodes.
	edge = NULL;
	for (int i = 0; i < totalVertices; i++)
	{
		headnodes[i].vertex = i;
		headnodes[i].next = NULL;
		headnodes[i].edgeWeight = 0;	//assume weight = 0 for this edge to itself.
		headnodes[i].edgeID = -1;	//all edgeID to itself is default -1.
	}
}
void Graph::create()
{
	int weight = 1;	//This assumes all edges weight = 1. Change accordingly.
	int edgeV1, edgeV2;
	double successRate;
	char comma;
	string miscString;
	node *tempPtr;

	ifstream fin;	//to read from a file.
	fin.open("numbers.txt");

	fin >> totalEdges;

	edge = new Edge[totalEdges];	//to setup the edge array.

	if (totalEdges > 0)
	{
		for (int i = 0; i < totalEdges; i++)
		{
			fin >> edgeV1 >> edgeV2 >> successRate;
			
			edge[i].determined = -1;	//-1 means this edge is not determined yet.
			edge[i].successRate = successRate;
			
			tempPtr = &headnodes[edgeV1];	//adding to headnode[edgeV1] a link to edgeV2.
			while (tempPtr->next != NULL)
				tempPtr = tempPtr->next;
			tempPtr->next = new node;
			tempPtr->next->vertex = edgeV2;
			tempPtr->next->next = NULL;
			tempPtr->next->edgeWeight = weight;
			tempPtr->next->edgeID = i;

			tempPtr = &headnodes[edgeV2];	//adding to headnode[edgeV2] a link to edgeV1.
			while (tempPtr->next != NULL)
				tempPtr = tempPtr->next;
			tempPtr->next = new node;
			tempPtr->next->vertex = edgeV1;
			tempPtr->next->next = NULL;
			tempPtr->next->edgeWeight = weight;
			tempPtr->next->edgeID = i;			
		}
	}
	
	fin.close();	//close the file.
}
void Graph::create(const char *fileName)
{
	int weight = 1;	//This assumes all edges weight = 1. Change accordingly.
	int edgeV1, edgeV2;
	double successRate;
	char comma;
	string miscString;
	node *tempPtr;

	ifstream fin;	//to read from a file.
	fin.open(fileName);

	fin >> totalEdges;

	edge = new Edge[totalEdges];	//to setup the edge array.

	if (totalEdges > 0)
	{
		for (int i = 0; i < totalEdges; i++)
		{
			fin >> edgeV1 >> edgeV2 >> successRate;

			edge[i].determined = -1;	//-1 means this edge is not determined yet.
			edge[i].successRate = successRate;

			tempPtr = &headnodes[edgeV1];	//adding to headnode[edgeV1] a link to edgeV2.
			while (tempPtr->next != NULL)
				tempPtr = tempPtr->next;
			tempPtr->next = new node;
			tempPtr->next->vertex = edgeV2;
			tempPtr->next->next = NULL;
			tempPtr->next->edgeWeight = weight;
			tempPtr->next->edgeID = i;

			tempPtr = &headnodes[edgeV2];	//adding to headnode[edgeV2] a link to edgeV1.
			while (tempPtr->next != NULL)
				tempPtr = tempPtr->next;
			tempPtr->next = new node;
			tempPtr->next->vertex = edgeV1;
			tempPtr->next->next = NULL;
			tempPtr->next->edgeWeight = weight;
			tempPtr->next->edgeID = i;
		}
	}

	fin.close();	//close the file.
}
int Graph::getTotalVertices()
{
	return totalVertices;
}
int Graph::getTotalEdges()
{
	return totalEdges;
}

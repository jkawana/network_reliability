#include <mpi.h>
#include <climits>
#include <algorithm>
#include "graph.h"
#include <stdlib.h>
#include <time.h>
#include <set>
#include "graph.cpp"

bool Dijkstra(int source, int destination, Graph G, int diameter);

int main(int argc, char *argv[])
{
    
    double expectedR = .997224, //0.339515,
           calcR;
    
    int    diameter = 8,
           numVertices = 25,
           source = 0,
           destination = 12,
           totalHits = 0,
           myRank,
           numProcs;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    
    double randomNum, start, end, totaltime, outputtime;
    
    const int localTrials = 1000000,
    totalTrials = localTrials * numProcs;
    
    int hits = 0;
    
    Graph G(numVertices);
    G.create();
    
    MPI_Barrier(MPI_COMM_WORLD);
    start = MPI_Wtime();
    
    srand(time(NULL)*myRank);
    
    //looping over number of times
    for ( int i = 0 ; i < localTrials; i++)
    {
        for ( int j = 0; j < G.getTotalEdges(); j++ )
        {
            randomNum = rand() % 100 + 1;
            if ( randomNum > G.edge[j].successRate * 100)
                G.edge[j].determined = 0;
            else
                G.edge[j].determined = 1;
        }
        if ( Dijkstra(source, destination, G, diameter) )
        {
            hits++;
        }
    }
    
    printf("Processor: %d     hits: %d\n",myRank,hits);
    
    end = MPI_Wtime();
    totaltime = end - start;
    
    MPI_Reduce(&totaltime, &outputtime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&hits, &totalHits, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    
    if (myRank == 0)
    {
        calcR = (double) totalHits / totalTrials;
        
        cout << "Total number of proccessors = " << numProcs << endl;
        cout << "Total number of hits = " << totalHits << endl;
        cout << "Each Proccessor performs " << localTrials << "  Trials" << endl;
        cout << "Total number of trials = " << totalTrials << endl;
        
        cout << "Monte Carlo Result " << endl;
        cout << "5x5grid" << endl; //change
        cout << "Diameter Constraint = " << diameter << endl;
        cout << "Source Node = " << source << endl;
        cout << "Destination Node = " << destination << endl;
        cout << "Calculated Reliability = " << calcR << endl;
        cout << "Expected Reliability = " << expectedR << endl;
        cout << "Difference = " << (calcR - expectedR) << endl;
        cout << "Time = " << outputtime << endl;
    }
    
    MPI_Finalize();
    
    return 0;
}

bool Dijkstra(int source, int destination, Graph G, int diameter)
{
    
    set<int> T; //the set of vertices used to make sure all verticies distances were
    
    int totalVertices = G.getTotalVertices();
    int *parent = new int[totalVertices];
    int *lambda = new int[totalVertices];
    
    //set up the vector to visit all nodes, T, the parent to have no parents and the min index to be the max
    for (int i = 0; i < totalVertices; i++)
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
    bool works = (parent[destination] != -1 && lambda[destination] <= diameter);
    
    delete[] parent;
    delete[] lambda;
    return works;
}

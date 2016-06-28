#include <mpi.h>
#include "graph.h"
#include "graph.cpp"
#include <stdlib.h>
#include <time.h>

bool Dijkstra(int source, int destination, Graph G, int diameter);

int main(int argc, char *argv[])
{
    
    double      expectedR = 0.234521,
                calcR;
        
    int         diameter = 8,
                numVertices = 25,
                source = 0,
                destination = 18,
                totalHits = 0,
                myRank,
                numProcs;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    
    double      randomNum,
                start,
                end,
                totaltime,
                outputtime;
    
    const int   localTrials = 1000000,
                totalTrials = localTrials * numProcs;
    
    int         hits = 0;
    
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
            
            
            //maybe just make this a function to determine it and send it the randnumber, so it can account for the number of edges currently in the graph which would help with finding the ri and ru stuff in
            if ( randomNum > G.edge[j].successRate * 100 )
                G.edge[j].determined = 0;
            else
                G.edge[j].determined = 1;
            
            
            
        }
        if ( Dijkstra(source, destination, G, diameter) )
        {
            hits++;
        }
    }
    
    printf("\nProcessor: %d     hits: %d\n\n", myRank, hits);
    
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
        
        cout << "\nMonte Carlo Result " << endl;
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
    int totalVertices = G.getTotalVertices();
    bool *T = new bool[totalVertices];
    int *father = new int[totalVertices];
    int *Lambda = new int[totalVertices];
    
    int minIndex= source; // u is the vertex with smallest Lambda, which at beginning is s.
    bool minIndexIsSet = true; //used to make sure that the minIndex is set to the first avail value each time
    
    for (int i = 0; i < totalVertices; i++)
    {
        T[i] = true; // all the vertices are eligible to be visited;
        father[i] = -1; // at the beginning every vertex's father is -1;
        Lambda[i] = INT_MAX; // every vertex is assumed to be far away from s.
    }
    Lambda[source] = 0; // s is at distance 0 of itself.
    
    //can this loop be while numOfVerticesToVisit !=0 or > 0!?

    for (int numVerticesToVisit = totalVertices; numVerticesToVisit>0; numVerticesToVisit--)
    {
        /* set u to the vertex with the smallest Lambda. */
        for (int i = 0; i < totalVertices; i++)
        {
            //cant this be done with something else?
            if (!minIndexIsSet && T[i])
            {
                minIndex= i;	//to set minIndexto the first available vertex.
                minIndexIsSet = true;
            }
            else if (minIndexIsSet && T[i] && Lambda[i] < Lambda[minIndex])
                minIndex= i;	//to set minIndexto the smallest Lambda (aka distance from source).
        }
        
        //itterate over all adjacent nodes
        for( node* adjnode = G.headnodes[minIndex].next; adjnode; adjnode = adjnode->next)
        {
            if (abs(G.edge[adjnode->edgeID].determined) == 1 &&
                Lambda[minIndex] != INT_MAX && T[adjnode->vertex] &&
                (Lambda[adjnode->vertex] > Lambda[minIndex] + adjnode->edgeWeight))
            {
                Lambda[adjnode->vertex] = Lambda[minIndex] + adjnode->edgeWeight;
                father[adjnode->vertex] = minIndex;
            }
        }
        minIndexIsSet = false; // to indicate minIndexis not set, I.E TAKE VIRST AVAIL INDEX IN FOR LOOP ABOVE
        T[minIndex] = false; //minIndexis not available to be visited (aka. to delete minIndex from T).
    }
    
    /* Dijkstra Results: */
    
    bool results = (father[destination] != -1 && Lambda[destination] <= diameter);
    /* delete dynamically allocated variables. */
    delete[] T;
    delete[] father;
    delete[] Lambda;
    return results;

}

#include <mpi.h>
#include "graph.h"
#include "graph.cpp"
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <climits>

bool    Dijkstra(int source, int destination, Graph G, int diameter);

double  factorial(int x);

double  combination(int x, int y);

double  F(int lo, int hi, int m, double p);

void comboSetUp(&dobule arr, int size, string file){
	ifstream inFile(file.c_str());
	for (int i =0; i<size; i++){
		infill>>arr[i];
	}
}

double combo40[41];

int main(int argc, char *argv[])
{

	comboSetUp(combo40, 41, "5x5combo.txt");
    
    double      expectedR = 0.997224,
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
                otherRandom,
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
        G.resetAliveEdges();
        for ( int j = 0; j < G.getTotalEdges(); j++ )
        {
            randomNum = rand() % 100 + 1; 
            G.setSuccRate(randomNum, j);
        }
        
        //need to add more edges from the dead (MAKE ALIVE)
        if (G.edgesAlive.size() < G.minPath)    //when t= 18, min cut is 6
        {
            randomNum = rand() % (G.maxAlive-G.edgesAlive.size()) + (G.minPath-G.edgesAlive.size());
            
            for (int i = 0 ; i < randomNum; i++){
                otherRandom = rand() % (G.edgesDead.size()-1);
                G.edge[ G.edgesDead[otherRandom] ].determined = 1;
                G.edgesDead.erase(G.edgesDead.begin() + otherRandom);
            }
        }
        
        //need to remove edges from the alive (KILL)
        else if (G.edgesAlive.size() > G.maxAlive) //max alive is 38
        {
            randomNum = rand() % (G.edgesAlive.size()-G.minPath) + (G.edgesAlive.size()-G.maxAlive);

            for (int i = 0 ; i < randomNum; i++){
                otherRandom = rand() % (G.edgesAlive.size()-1);
                G.edge[ G.edgesAlive[otherRandom] ].determined = 0;
                G.edgesAlive.erase(G.edgesAlive.begin() + otherRandom);
            }
        }
        
        if ( Dijkstra(source, destination, G, diameter) )
        {
            hits++;
        }
    }
    
    end = MPI_Wtime();
    totaltime = end - start;
    
    MPI_Reduce(&totaltime, &outputtime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&hits, &totalHits, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    
    if (myRank == 0)
    {

        double Rl = F(G.totalEdges-G.minCut+1, G.totalEdges, G.totalEdges, 0.95);
        double Ru = 1 - F(0, G.minPath-1, G.totalEdges, 0.95);
	 
        calcR = (double) totalHits / totalTrials;
        
        
        double finalResult = Rl + (Ru-Rl)*calcR ;
        
        
		cout << "Rl, Ru = " << Rl << ' ' << Ru << endl;
		cout << "MC calculated reliability = " << calcR << endl;

        cout << "\nTotal number of proccessors = " << numProcs << endl;
        cout << "Total number of hits = " << totalHits << endl;
        cout << "Each Proccessor performs " << localTrials << "  Trials" << endl;
        cout << "Total number of trials = " << totalTrials << endl;
        
        cout << "\nMonte Carlo Result " << endl;
        cout << "5x5grid" << endl; 
        cout << "Diameter Constraint = " << diameter << endl;
        cout << "Source Node = " << source << endl;
        cout << "Destination Node = " << destination << endl;
        cout << "Calculated Reliability = " << finalResult << endl;
        cout << "Expected Reliability = " << expectedR << endl;
        cout << "Difference = " << (finalResult - expectedR) << endl;
        cout << "Time = " << outputtime << endl;
    }
    
    MPI_Finalize();
    
    return 0;
}

bool Dijkstra(int source, int destination, Graph G, int diameter)
{
    int totalVertices = G.getTotalVertices();
    
    bool *T = new bool[totalVertices];
    int  *father = new int[totalVertices];
    int  *Lambda = new int[totalVertices];
    
    int minIndex= source; // u is the vertex with smallest Lambda, which at beginning is s.
    bool minIndexIsSet = true; //used to make sure that the minIndex is set to the first avail value each time
    
    for (int i = 0; i < totalVertices; i++)
    {
        T[i] = true; // all the vertices are eligible to be visited;
        father[i] = -1; // at the beginning every vertex's father is -1;
        Lambda[i] = INT_MAX; // every vertex is assumed to be far away from s.
    }
    Lambda[source] = 0; // s is at distance 0 of itself.
    
    for (int numVerticesToVisit = totalVertices; numVerticesToVisit > 0; numVerticesToVisit--)
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
        minIndexIsSet = false; // to indicate minIndexis not set, I.E TAKE FIRST AVAIL INDEX IN FOR LOOP ABOVE
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

double factorial(int x)
//this or the factorial choose was causing the problem before.
{
    double sum=1;
    for (int i = 1; i<=x; i++) sum*=i;
    
    return sum;
}

double combination(int x, int y){
    return factorial(x) / (factorial(y)*factorial(x-y));
}


double F(int lo, int hi, int m, double p)
{
    double sum = 0;
    for ( int j = lo; j <= hi; j++ )
    {
        sum += (combination(m,j) * pow(p,j) * pow((1-p),(m-j)));
    }
    return sum;
}





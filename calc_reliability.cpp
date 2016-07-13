#include <mpi.h>
#include "graph.h"
#include "graph.cpp"
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <climits>
#include <algorithm>

bool    Dijkstra(int source, int destination, Graph G, int diameter);

double  factorial(int x);

double  combination(int x, int y);

double  F(int lo, int hi, int m, double p);

double * combo;

int * trialAmounts;

double O(int totalEdges, int edgesAlive, double prob);

double S(int lo, int hi, int totalEdges, int subgraphSize, double prob);

int main(int argc, char *argv[])
{

    //Get rid of edgesAlive and edgesDead

    ifstream enviVar("environmentVariables.txt");
	
    double      expectedR,
                calcR,
                edgeRel;
    
    int         diameter,
                numVertices,  //possibly move to graph file
                numEdges,
                source,
                destination,
                totalHits = 0,
                myRank,
                numProcs,
                minPath,
                minCut;
    

    string graphFile, combinationFile;
    
    double      randomNum,
                otherRandom,
                start,
                end,
                totaltime,
                outputtime;
    
    const int   localTrials = 1000000,
                totalTrials = localTrials * numProcs;
    
    int         hits = 0, indices[numEdges];

    //so now that we know that everything gets done four times, maybe we should move diameter, vertecies, edges, edgeRel? to the graph file

    //but i think keep the min cut, source, destination in the file theyre in now because even if the graph is the same, 
    //different source and destination for the same graph, which would result in different mincut / min path

    //so after doing all of this, create the graph so then we have acces to those varibles

    //should we put local trials in the file as well?

    enviVar >> expectedR >> edgeRel >> diameter >> numVertices >> numEdges >> minPath;
    enviVar >> minCut >> source >> destination >> graphFile >> combinationFile;

    enviVar.close();

    combo = new double[numEdges+1];
    ifstream comboReader(combinationFile.c_str()); 
    for (int i = 0 ; i < numEdges + 1; i++)
    {
        //combo[i]=combination(totalEdges, i);
        comboReader >> combo[i];
    }

    comboReader.close();

    int monkey=0;
    trialAmounts = new int [numEdges - minCut - minPath + 1];
    for (int i = 0 ; i < numEdges - minCut - minPath + 1; i++)
    {
        trialAmounts[i] = localTrials * S(minPath, numEdges - minCut, numEdges, minPath + i, edgeRel);
        monkey += trialAmounts[i];
	    cout << "i = " << i+6 << "  " << trialAmounts[i] << endl;
    }
    cout<<monkey<<endl;


    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    

    Graph G(numVertices);
    G.create();
    
    MPI_Barrier(MPI_COMM_WORLD);
    start = MPI_Wtime();
    
    srand(time(NULL)*myRank);
    

    for (int i = 0; i < numEdges; i++)
    {
        indices[i]=i;
    }


    for (int i = 0 ; i < numEdges - minCut - minPath + 1; i++)
    {
        //should this be up to trialAmounts[i] to get that many trials?
        for (int j = 0; j < i; j++)
        {
            //note: need + 1 to get the last index to be swapped
            random_shuffle(&indices[0], &indices[numEdges - minCut - minPath+1]);   //our rand vs their rand???????
        }
        G.resetEdges(); //here, right?...
        for (int k = 0; k < minPath + i; k++)
        {
            G.edge[k].determined = 1;
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
        double Rl = F(numEdges-minCut+1, numEdges, numEdges, edgeRel);
        double Ru = 1 - F(0, minPath-1, numEdges, edgeRel);
	 
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
        sum += combo[j] * pow(p,j) * pow((1-p),(m-j));
    }
    return sum;
}

double O(int totalEdges, int subgraphSize, double prob)
{
	//O(M) =  C(40,M) * (.95)^M * (.05) ^ (40-M).
	return combo[subgraphSize] * pow(prob, subgraphSize) * pow((1-prob), totalEdges - subgraphSize);
}

double S(int lo, int hi, int totalEdges, int subgraphSize, double prob)
{

	double divisor=0;

	for (int i = lo; i <= hi; i++) divisor += O(totalEdges, i, prob);

	return O(totalEdges, subgraphSize, prob) / divisor;

}

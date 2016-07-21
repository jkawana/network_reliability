#ifndef graph_h
#define graph_h
#include <iostream>
#include <fstream>


using namespace std;

struct Edge
{
    int determined;                     //-1 is not determined, 0 is determined this edge failed, 1 is determined this edge exist.
};

struct node
{
    int vertex;                         //vertex number aka name.
    node * next;                        //ptr to the next adjacent vertex of the headnode.
    int edgeWeight;                     //weight of this edge relative to the headnode, which is assumed to be 1. (changed accordingly in the implementation file.)
    int edgeID;                         //real edgeID values range from 0 to totalEdge-1, max.
};

class Graph
{
private:
    int totalVertices;                  // total number of vertices.
    int totalEdges;                     // total number of edges.
    
public:
    Edge *edge;                         //ptr to array of the Edges and the indexes are the edge names.
    node *headnodes;                    //ptr to the headnodes.
    
    Graph();                            //set all data members to zero or null.
    Graph(int nodes, const char *fileName);                   //sets the totalVertices and the headnodes accordingly.
    
    void setVertices(int nodes);        //sets the totalVertices and the headnodes accordingly.
    void create();                      //will setup the edge array according to the file "numbers.txt".
    void create(const char *fileName);	//will setup the edge array according to the parameterized fileName.
    int getTotalVertices();             //returns the totalVertices.
    int getTotalEdges();                //returns the totalEdges.    
    void resetEdges();


};
#endif

/*******************************************************************************
//  2TSP_DFS_BFS.cpp                Author: Ian Nobile
//  Section: 50                     Due Date: 30 August 2021
//
//  This memory-leak- and bug-free program implements Breadth- and Depth-First 
//  Search (BFS and DFS) algorithms to find the shortest path to the goal city 
//  of the given dataset (11PointDFSBFS.tsp), reports the optimum paths and the 
//  total minimum distances and compares the time each algorithm took to arrive 
//  at each.
//
*******************************************************************************/

#include <chrono>   // time the speed of the program
#include <iostream> // print to console
#include <fstream>  // read files
#include <string>   // header buffer for seeking inside a file
#include <vector>   // easy arraying
#include <cmath>    // distance formula
#include <queue>    // queue for BFS
#include <stack>    // stack for DFS

using namespace std;
using namespace std::chrono;


// class declarations:
class Node {
public:
    int num;
    float x;
    float y;
    float dist;
    Node* parent;   // for bfs
    bool visited;   // for bfs
    vector<Node*> children;

    // constructor:
    Node() {
        num = 0;
        x = 0.0;
        y = 0.0;
        dist = FLT_MAX;
        parent = NULL;
        visited = false;
    }
};

class Graph {
public:
    int dimension;
    vector<Node> nodes;
};


// function prototypes:
Graph buildGraph(char*);
float distCheck(float, float, float, float);
void edger(Graph&);
void bfs(Graph&);
void dfs(Node*, stack<Node*>*, float*, stack<Node*>*, float*);
void dfs(Graph&);


//------------------------------------------------------------------------------
//  Main Function
//------------------------------------------------------------------------------
int main(int argc, char* argv[]) {
    // check if path was passed as arg:
    if (argc == 1) {
        cout << "Please pass the path to the .TSP file as a command line argument" << endl;
        return 1;
    }
    
    // begin with a friendly greeting:
    cout << "Hello and welcome to the (BFS/DFS) Travelling Salesperson Problem Solver" << endl;

    Graph graph = buildGraph(argv[1]);

    // bfs starts here
    auto start = high_resolution_clock::now();  // start timer
    bfs(graph);
    auto stop = high_resolution_clock::now();   // stop timer
    auto bfsDuration = duration_cast<microseconds>(stop - start);  // calculate elapsed time

    // quick reset of node distances for dfs
    for (Node node : graph.nodes) { node.dist = FLT_MAX; }
    
    // dfs starts here
    start = high_resolution_clock::now();
    dfs(graph);
    stop = high_resolution_clock::now();
    auto dfsDuration = duration_cast<microseconds>(stop - start);
    
    // print speeds to console:
    cout << "Breadth-First Search execution took " << bfsDuration.count() / 1000000.0 << "s in total" << endl;
    cout << "Depth-First Search execution took " << dfsDuration.count() / 1000000.0 << "s in total" << endl << endl;

    // compare time results:
    if (bfsDuration.count() < dfsDuration.count()) {
        cout << "BFS was the faster algorithm, beating DFS by a full " << (dfsDuration.count() - bfsDuration.count()) / 1000000.0 << " seconds!" << endl << endl;
    }
    else {
        cout << "DFS was the faster algorithm, beating BFS by a full " << (bfsDuration.count() - dfsDuration.count()) / 1000000.0 << " seconds!" << endl << endl;
    }
    return 0;
}


//  function definitions:
//------------------------------------------------------------------------------
//  Creates edges between each node and their children; currently, just a 
//  hardcoding of the table found in the project prompt
//------------------------------------------------------------------------------
void edger(Graph& graph) {
    /*  The map:
        0:1,2,3
        1:2
        2:3,4
        3:4,5,6
        4:6,7
        5:7
        6:8,9
        7:8,9,10
        8:10
        9:10
        10:X    */

    // surely there's a more intelligent way to do this...
    graph.nodes[0].children.push_back(&graph.nodes[1]);
    graph.nodes[0].children.push_back(&graph.nodes[2]);
    graph.nodes[0].children.push_back(&graph.nodes[3]);

    graph.nodes[1].children.push_back(&graph.nodes[2]);

    graph.nodes[2].children.push_back(&graph.nodes[3]);
    graph.nodes[2].children.push_back(&graph.nodes[4]);

    graph.nodes[3].children.push_back(&graph.nodes[4]);
    graph.nodes[3].children.push_back(&graph.nodes[5]);
    graph.nodes[3].children.push_back(&graph.nodes[6]);

    graph.nodes[4].children.push_back(&graph.nodes[6]);
    graph.nodes[4].children.push_back(&graph.nodes[7]);

    graph.nodes[5].children.push_back(&graph.nodes[7]);

    graph.nodes[6].children.push_back(&graph.nodes[8]);
    graph.nodes[6].children.push_back(&graph.nodes[9]);

    graph.nodes[7].children.push_back(&graph.nodes[8]);
    graph.nodes[7].children.push_back(&graph.nodes[9]);
    graph.nodes[7].children.push_back(&graph.nodes[10]);

    graph.nodes[8].children.push_back(&graph.nodes[10]);

    graph.nodes[9].children.push_back(&graph.nodes[10]);
}

//------------------------------------------------------------------------------
//  Read .TSP file, creates nodes from coordinates and combines all in an
//  undirected graph object
//------------------------------------------------------------------------------
Graph buildGraph(char* argv) {
    // open .TSP file in read-only mode:
    ifstream tspfile;
    tspfile.open(argv, ios::in);
    // ensure file exists:
    if (!tspfile.is_open()) {
        Graph graph;
        graph.dimension = 0;
        return graph;
    }
    Graph graph;
    graph.dimension = 0;
    string heading = "";
    // advance buffer to the dimension section:
    while (heading.compare("DIMENSION:") != 0) {
        tspfile >> heading;
    }
    tspfile >> graph.dimension;
    // advance buffer to the coordinates section:
    while (heading.compare("NODE_COORD_SECTION") != 0) {
        tspfile >> heading;
    }
    // create nodes and push to graph vector
    Node newNode = Node();
    for (int i = 0;i < graph.dimension;i++) {
        tspfile >> newNode.num;
        tspfile >> newNode.x;
        tspfile >> newNode.y;
        newNode.dist = FLT_MAX;
        newNode.parent = NULL;
        newNode.visited = false;
        graph.nodes.push_back(newNode);
    }
    
    // The graph is now created, and we are finished with the .TSP file
    tspfile.close();

    // connect edges
    edger(graph);

    return graph;
}

//------------------------------------------------------------------------------
//  Returns the distance between two graph nodes using the pythagoran formula: 
//  dist = sqrt((x2 - x1)^2 + (y2 - y1)^2)
//------------------------------------------------------------------------------
float distCheck(float x1, float y1, float x2, float y2) {
    return sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2));
}

//------------------------------------------------------------------------------
//  This shoddy implementation of Dijkstra's Algorithm starts at the root node 
//  and checks the distance between it and each of its children, pushing each to
//  a queue in order of distance before iterating to the first therein and 
//  repeating until the goal node is encountered, at which point the optimum 
//  distance and all parent nodes will have already been recorded
//------------------------------------------------------------------------------
void bfs(Graph& graph) {
    Node* node = &graph.nodes[0]; // start at node 1
    node->dist = 0.0;
    queue<Node*> bfsq;
    bfsq.push(node);
    float distComp = 0.0;

    while (!bfsq.empty()) {
        node = bfsq.front();
        bfsq.pop();
        if (node->num == 11) { break; }
        for (Node* child : node->children) {
            distComp = node->dist + distCheck(node->x, node->y, child->x, child->y);
            //set new distance and parent if better than others
            if (distComp < child->dist) {
                child->dist = distComp;
                child->parent = node;
            }
        }
        // now, find minimum distance, set visited and push to q, then repeat until all visited
        for (int i = 0;i < node->children.size();i++) {
            distComp = FLT_MAX;  // reset minimum
            for (Node* child : node->children) {
                if (!child->visited && child->dist < distComp) {
                    distComp = child->dist;
                }
            }
            for (Node* child : node->children) {
                if (child->dist == distComp && !child->visited) {
                    child->visited = true;
                    bfsq.push(child);
                }
            }
        }
    } // end while
    
    // print winning route
    cout << "According to the BFS algorithm, The shortest path between 1 and 11 is:" << endl << endl << "\t";
    distComp = node->dist;
    stack<Node*> bfstack;
    bfstack.push(node);
    while (node->parent != NULL) {
        node = node->parent;
        bfstack.push(node);
    }    
    while (!bfstack.empty()) {
        cout << bfstack.top()->num;
        bfstack.pop();
        if (bfstack.size() == 0) { break; }
        cout << " -> ";
    }
    cout << endl << endl << "at a total distance of " << distComp << endl << endl << endl;
}

//------------------------------------------------------------------------------
//  This is the recursive portion of the depth-first search algorithm 
//  implementation. It too starts at the root but, unlike bfs, continues down
//  via the first child node until it reaches the goal. It then backtracks and 
//  examines parallel paths to see if their total distances are less; 
//  computationally resource-heavy, eventually exploring every possible 
//  alternative; also naive in that it too assumes all paths will  lead to the 
//  goal node
//------------------------------------------------------------------------------
void dfs(Node* node, stack<Node*>* minPath, float* minDist, stack<Node*>* path, float* currDist){
    if (node->num==11) {    //base case
        path->push(node);
        //update minPath and minDist if better:
        if (*currDist < *minDist) { 
            *minDist = *currDist;
            *minPath = *path;
        }
        path->pop();
        return;
    }
    for (Node* child : node->children) {
        *currDist += distCheck(node->x, node->y, child->x, child->y);
        path->push(node);
        dfs(child, minPath, minDist, path, currDist);
        *currDist -= distCheck(node->x, node->y, child->x, child->y);
        path->pop();
    }
}

//------------------------------------------------------------------------------
//  This is a wrapper function for this implementation of a depth-first search 
//  algorithm; declutters the main function significantly, while also handling 
//  the results-printing for the recursive function above
//------------------------------------------------------------------------------
void dfs(Graph& graph) {
    Node* node = &graph.nodes[0];
    node->dist = 0.0;
    stack<Node*> minPath;
    float minDist = FLT_MAX;
    stack<Node*> path;
    float currDist = 0.0;

    // begin recursive bit
    dfs(node, &minPath, &minDist, &path, &currDist);

    // print winning route
    while (!minPath.empty()) {
        path.push(minPath.top());
        minPath.pop();
    }
    cout << "According to the DFS algorithm, The shortest path between 1 and 11 is" << endl << endl << "\t";
    while (!path.empty()) {
        node = path.top();
        path.pop();
        cout << node->num;
        if (node->num == 11) { break; }
        cout << " -> ";
    }
    cout << endl << endl << "at a total distance of " << minDist << endl << endl << endl;
}


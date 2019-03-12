#define INT_MAX 1789654
#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <string>
#include <math.h>

Eigen::MatrixXd V;
Eigen::MatrixXi F;

int minDistance(int dist[], bool sptSet[], int V_m)
{
   // Initialize min value
   int min = 0, min_index;

   for (int v = 0; v < V_m; v++)
     if (sptSet[v] == false && dist[v] <= min)
         min = dist[v], min_index = v;

   return min_index;
}

void dijkstra(Eigen::MatrixXd graph, int src)
{

     V_m = graph.rows();
     int dist[V_m];     // The output array.  dist[i] will hold the shortest
                      // distance from src to i

     bool sptSet[V_m]; // sptSet[i] will be true if vertex i is included in shortest
                     // path tree or shortest distance from src to i is finalized

     // Initialize all distances as INFINITE and stpSet[] as false
     for (int i = 0; i < V_m; i++)
        dist[i] = INT_MAX, sptSet[i] = false;

     // Distance of source vertex from itself is always 0
     dist[src] = 0;

     // Find shortest path for all vertices
     for (int count = 0; count < V_m-1; count++)
     {
       // Pick the minimum distance vertex from the set of vertices not
       // yet processed. u is always equal to src in the first iteration.
       int u = minDistance(dist, sptSet, V_m);

       // Mark the picked vertex as processed
       sptSet[u] = true;

       // Update dist value of the adjacent vertices of the picked vertex.
       for (int v = 0; v < V_m; v++)

         // Update dist[v] only if is not in sptSet, there is an edge from
         // u to v, and total weight of path from src to  v through u is
         // smaller than current value of dist[v]
         if (!sptSet[v] && graph(u,v) && dist[u] != INT_MAX
                                       && dist[u]+graph(u,v) < dist[v])
            dist[v] = dist[u] + graph(u,v);
     }

     // print the constructed distance array
     // printSolution(dist, V_m);
}

int main(int argc, char *argv[])
{
  // Load a mesh in OFF format
  igl::readOFF(argv[1], V, F);

  std::string start, finish;
  start = argv[2],
  finish = argv[3];
  int st = 0;
  int fin = 0;
  int m_v = V.rows();
  int n_v = V.cols();
  int m_f = F.rows();
  int n_f = F.cols();
  Eigen::MatrixXd V_cost(m_v, m_v);

  for(size_t i=0; i<m_v; i++)
  {
    for(size_t j=0; j<m_v; j++)
    {
      double length = sqrt(pow(V(i,0) - V(j,0), 2) + pow(V(i,1) - V(j,1), 2) +
        pow(V(i,2) - V(j,2), 2));
        V_cost(i,j) = length;
    }
  }

  for(size_t i=0; i<10; i++)
  {
    std::cout << "length is " << V_cost(i,0) << std::endl;
  }

  try
  {
	  st = std::stoi(start);
	  fin = std::stoi(finish);
  }
  catch(...)
  {
	  std::cout << "Input vertices cannot be converted to int" <<
		  std::endl;
  }


  std::cout << "The start is " << st << " and finish is " <<
	  fin << std::endl;


  // for(size_t i=0; i<10; i++)
 //  {
// 	  std::cout << "V is " << V(i,0) << " and F is " << F(i,0) <<
//		 std::endl;
//  }



  // Plot the mesh
  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(V, F);
  viewer.launch();
}

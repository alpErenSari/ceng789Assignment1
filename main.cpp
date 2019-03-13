#define INT_MAX 1000
#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/vertex_triangle_adjacency.h>
#include <igl/adjacency_list.h>
#include <string>
#include <math.h>

Eigen::MatrixXd V;
Eigen::MatrixXi F;

int minDistance(double dist[], bool sptSet[], int V_m)
{
   // Initialize min value
   double min = INT_MAX, min_index;

   for (int v = 0; v < V_m; v++)
     if (sptSet[v] == false && dist[v] <= min)
         min = dist[v], min_index = v;

   return min_index;
}

double l2_norm(Eigen::MatrixXd &V, int i, int j)
{
  return sqrt(pow(V(i,0) - V(j,0), 2) + pow(V(i,1) - V(j,1), 2) +
    pow(V(i,2) - V(j,2), 2));
}

struct dij_out
{
  double min_dist;
  std::vector<int> route;
};

dij_out dijkstra(Eigen::MatrixXd graph, int src, int des)
{

      int V_m = graph.rows();
      std::cout << "V_m is " << V_m << std::endl;
     double dist[V_m];     // The output array.  dist[i] will hold the shortest
     int prev[V_m];
                      // distance from src to i
     std::vector<int> parent;

     bool sptSet[V_m]; // sptSet[i] will be true if vertex i is included in shortest
                     // path tree or shortest distance from src to i is finalized

     // Initialize all distances as INFINITE and stpSet[] as false
     for (int i = 0; i < V_m; i++)
        dist[i] = INT_MAX, sptSet[i] = false, prev[i] = -1;

     // Distance of source vertex from itself is always 0
     dist[src] = 0;

     // Find shortest path for all vertices
     for (int count = 0; count < V_m-1; count++)
     {
       // std::cout << "count number " << count << std::endl;
       // Pick the minimum distance vertex from the set of vertices not
       // yet processed. u is always equal to src in the first iteration.
       int u = minDistance(dist, sptSet, V_m);
       // std::cout << "min vertex is " << u << std::endl;

       // Mark the picked vertex as processed
       sptSet[u] = true;

       // Update dist value of the adjacent vertices of the picked vertex.
       for (int v = 0; v < V_m; v++)

         // Update dist[v] only if is not in sptSet, there is an edge from
         // u to v, and total weight of path from src to  v through u is
         // smaller than current value of dist[v]
         if (!sptSet[v] && graph(u,v) && dist[u] != INT_MAX
                                       && dist[u]+graph(u,v) < dist[v])
                                       {
                                         dist[v] = dist[u] + graph(u,v);
                                         prev[v] = u;
                                       }

     }


     int u = des;
     if(prev[u]>=0 || u==src)
     {
       while(u>=0)
       {
         parent.push_back(u);
         u = prev[u];
       }
     }

     dij_out retval;
     retval.min_dist = dist[des];
     retval.route = parent;

     return retval;

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
  Eigen::MatrixXi V_adj(m_v, m_v);
  V_cost = Eigen::MatrixXd::Zero(m_v, m_v);
  V_adj = Eigen::MatrixXi::Zero(m_v, m_v);
  // std::vector<std::vector<int> > VF;
  // std::vector<std::vector<int> > VFi;

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
  // for(size_t i=0; i<m_v; i++)
  // {
  //   for(size_t j=0; j<m_v; j++)
  //   {
  //     double length = sqrt(pow(V(i,0) - V(j,0), 2) + pow(V(i,1) - V(j,1), 2) +
  //       pow(V(i,2) - V(j,2), 2));
  //       V_cost(i,j) = length;
  //   }
  // }

  std::vector<std::vector<double>> A;
  igl::adjacency_list(F,A);


    for(size_t j=0; j<A[0].size(); j++)
    {
      std::cout << A[st][j] << " ";
    }
    std::cout << std::endl;



  for(size_t i=0; i<A.size(); i++)
  {
    for(size_t j=0; j<A[0].size(); j++)
    {
      int k = A[i][j];
      V_cost(i, k) = l2_norm(V, i, k);
    }
  }

  for(size_t i=0; i<10; i++)
  {
    // std::cout << "The results are \n";
    for(size_t j=0; j<10; j++)
    {
      std::cout << V_cost(i,j) << " ";
    }
    std::cout << std::endl;
  }

  double min_dist;
  std::vector<int> route;
  dij_out my_out;
  my_out = dijkstra(V_cost, st, fin);
  std::cout << "The minimum distance is " << my_out.min_dist << std::endl;


  for(size_t i=0; i<my_out.route.size(); i++)
  {
	  std::cout << my_out.route[i] << " ";
 }
 std::cout << std::endl;



  // Plot the mesh
  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(V, F);
  viewer.launch();
}

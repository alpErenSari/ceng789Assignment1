#define INT_MAX 1000
#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/vertex_triangle_adjacency.h>
#include <igl/adjacency_list.h>
#include "fiboheap.h"
#include <string>
#include <queue>
#include <algorithm>
#include <math.h>
#include <ctime>
#include <fstream>

Eigen::MatrixXd V;
Eigen::MatrixXi F;
Eigen::MatrixXd C;

typedef std::pair<double, int> dij_Pair;

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



dij_out dijkstra(std::vector<std::vector<dij_Pair>> adj_vec, int src, int des)
{


    int V_m = adj_vec.size();
    // std::cout << "V_m is " << V_m << std::endl;
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
      if(u == des)
        break;
      // std::cout << "min vertex is " << u << std::endl;

      // Mark the picked vertex as processed
      sptSet[u] = true;

      // Update dist value of the adjacent vertices of the picked vertex.
      std::vector<dij_Pair>::iterator i;
      for (i = adj_vec[u].begin(); i != adj_vec[u].end(); i++){

      // Update dist[v] only if is not in sptSet, there is an edge from
      // u to v, and total weight of path from src to  v through u is
      // smaller than current value of dist[v]
      int v = (*i).second;
      int weight = (*i).first;
      if (!sptSet[v] && dist[u] != INT_MAX
      && dist[u]+weight < dist[v])
      {
        dist[v] = dist[u] + weight;
        prev[v] = u;
      }
      }
    }


    int k = des;
    if(prev[k]>=0 || k==src)
    {
      while(k>=0)
      {
        parent.push_back(k);
        k = prev[k];
      }
    }

    dij_out retval;
    retval.min_dist = dist[des];
    retval.route = parent;

    return retval;
}


dij_out dijkstra_min_heap(std::vector<std::vector<dij_Pair>> adj_vec, int src, int des)
{
      int V_m = adj_vec.size();
      // std::cout << "V_m is " << V_m << std::endl;
     double dist[V_m];     // The output array.  dist[i] will hold the shortest
     int prev[V_m];
                      // distance from src to i
     std::vector<int> parent;
     std::priority_queue <dij_Pair, std::vector<dij_Pair>, std::greater<dij_Pair> > pq;

     bool sptSet[V_m]; // sptSet[i] will be true if vertex i is included in shortest
                     // path tree or shortest distance from src to i is finalized
     // Distance of source vertex from itself is always 0
     dist[src] = 0;
     pq.push(std::make_pair(0, src));

     // Initialize all distances as INFINITE and stpSet[] as false
     for (int i = 0; i < V_m; i++)
     {
       if(i != src)
       {
         pq.push(std::make_pair(INT_MAX, i)), dist[i] = INT_MAX;
       }
       prev[i] = -1;
     }

     // Find shortest path for all vertices
     while (!pq.empty())
     {
       // std::cout << "count number " << count << std::endl;
       // Pick the minimum distance vertex from the set of vertices not
       // yet processed. u is always equal to src in the first iteration.
       int u = pq.top().second;
       pq.pop();
       // std::cout << "min vertex is " << u << std::endl;
       if(u == des)
         break;

       std::vector<dij_Pair>::iterator i;
       // Update dist value of the adjacent vertices of the picked vertex.
       for (i = adj_vec[u].begin(); i != adj_vec[u].end(); i++)
       {
         // Update dist[v] only if is not in sptSet, there is an edge from
         // u to v, and total weight of path from src to  v through u is
         // smaller than current value of dist[v]
         int v = (*i).second;
         double weight = (*i).first;

         // If there is shorted path to v through u.
         if (dist[v] > dist[u] + weight)
         {
             // Updating distance of v
             dist[v] = dist[u] + weight;
             pq.push(std::make_pair(dist[v], v));
             prev[v] = u;
         }
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

dij_out dijkstra_fibo_heap(std::vector<std::vector<dij_Pair>> adj_vec, int src, int des)
{
      int V_m = adj_vec.size();
      // std::cout << "V_m is " << V_m << std::endl;
     double dist[V_m];     // The output array.  dist[i] will hold the shortest
     int prev[V_m];
                      // distance from src to i
     std::vector<int> parent;
     FibHeap <dij_Pair> pq;

     bool sptSet[V_m]; // sptSet[i] will be true if vertex i is included in shortest
                     // path tree or shortest distance from src to i is finalized
     // Distance of source vertex from itself is always 0

     // Initialize all distances as INFINITE and stpSet[] as false
     for (int i = 0; i < V_m; i++)
     {
       if(i == src)
       {
         pq.push(std::make_pair(0, i));
         dist[i] = 0;
       }
       else
       {
         pq.push(std::make_pair(INT_MAX, i)), dist[i] = INT_MAX;
       }
       prev[i] = -1;
     }

     // Find shortest path for all vertices
     while (!pq.empty())
     {
       // std::cout << "count number " << count << std::endl;
       // Pick the minimum distance vertex from the set of vertices not
       // yet processed. u is always equal to src in the first iteration.
       int u = pq.top().second;
       pq.pop();
       // std::cout << "min vertex is " << u << std::endl;
       if(u == des)
         break;

       std::vector<dij_Pair>::iterator i;
       // Update dist value of the adjacent vertices of the picked vertex.
       for (i = adj_vec[u].begin(); i != adj_vec[u].end(); i++)
       {
         // Update dist[v] only if is not in sptSet, there is an edge from
         // u to v, and total weight of path from src to  v through u is
         // smaller than current value of dist[v]
         int v = (*i).second;
         double weight = (*i).first;

         // If there is shorted path to v through u.
         if (dist[v] > dist[u] + weight)
         {
             // Updating distance of v
             dist[v] = dist[u] + weight;
             pq.push(std::make_pair(dist[v], v));
             prev[v] = u;
         }
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
}


int main(int argc, char *argv[])
{
  // Load a mesh in OFF format
  igl::readOFF(argv[1], V, F);

  std::string start, finish, fps_opt;
  start = argv[2],
  finish = argv[3];
  int st = 0;
  int fin = 0;
  int m_v = V.rows();
  int n_v = V.cols();
  int m_f = F.rows();
  int n_f = F.cols();
  Eigen::MatrixXd V_geodesic(m_v, m_v);
  // Eigen::MatrixXd V_cost(m_v, m_v);
  // Eigen::MatrixXi V_adj(m_v, m_v);
  // V_cost = Eigen::MatrixXd::Zero(m_v, m_v);
  // V_adj = Eigen::MatrixXi::Zero(m_v, m_v);
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
  igl::adjacency_list(F,A, true);
  std::vector<std::vector<dij_Pair>> my_pair(A.size());

    // for(size_t j=0; j<A[0].size(); j++)
    // {
    //   std::cout << A[st][j] << " ";
    // }
    // std::cout << std::endl;

  for(size_t i=0; i<A.size(); i++)
  {
    for(size_t j=0; j<A[i].size(); j++)
    {
      int k = A[i][j];
      double weight = l2_norm(V, i, k);
      // V_cost(i, k) = weight;
      my_pair[i].push_back(std::make_pair(weight, k));
    }
  }

  // for(size_t i=0; i<10; i++)
  // {
  //   // std::cout << "The results are \n";
  //   for(size_t j=0; j<10; j++)
  //   {
  //     std::cout << V_cost(i,j) << " ";
  //   }
  //   std::cout << std::endl;
  // }

  double min_dist;
  std::vector<int> route;

  std::clock_t start_time = std::clock();
  double duration;
  dij_out my_out;

  int bar_number = m_v/100;
  double max_geodesic = 0;
  for (size_t i = 0; i < m_v; i++) {
    for (size_t j = 0; j < m_v; j++) {
      dij_out my_out_temp = dijkstra_fibo_heap(my_pair, i, j);
      V_geodesic(i,j) = my_out_temp.min_dist;
      if(my_out_temp.min_dist > max_geodesic)
        max_geodesic = my_out_temp.min_dist;
      if(i == st && j == fin)
        my_out = my_out_temp;
    }
    // if(i%bar_number==0)
    //   std::cout << "Geodesic distance matrix progress: " << i/bar_number << "% \n";
  }

  duration = ( std::clock() - start_time ) / (double) CLOCKS_PER_SEC;
  std::cout << "Operation took "<< duration << "seconds for array dij" << std::endl;
  // my_out = dijkstra(my_pair, st, fin);
  // my_out = dijkstra_min_heap(my_pair, st, fin);
  // my_out = dijkstra_fibo_heap(my_pair, st, fin);
  std::cout << "The minimum distance is " << my_out.min_dist << std::endl;
  std::cout << "The maximum distance is " << max_geodesic << std::endl;

  std::ofstream myfile;
  myfile.open ("first_geodesic.txt");
  for (size_t i = 0; i < m_v; i++) {
    for (size_t j = 0; j < m_v; j++) {
      myfile << V_geodesic(i,j) << " ";
    }
    myfile << "\n";
  }
  myfile.close();


  for(size_t i=0; i<my_out.route.size(); i++)
  {
	  std::cout << my_out.route[i] << " ";
 }
 std::cout << std::endl;

// put the coordinates of the vertices in matrices to draw the line
Eigen::MatrixXd P1(my_out.route.size()-1, 3);
Eigen::MatrixXd P2(my_out.route.size()-1, 3);
for(size_t i=0; i<my_out.route.size()-1; i++)
{
  int vertex1 = my_out.route[i];
  int vertex2 = my_out.route[i+1];
  // adding 0.1 to coordinates to make the line visible
  P1(i,0) = V(vertex1,0) + 0.1, P1(i,1) = V(vertex1,1) + 0.1, P1(i,2) = V(vertex1,2) + 0.1;
  P2(i,0) = V(vertex2,0) + 0.1, P2(i,1) = V(vertex2,1) + 0.1, P2(i,2) = V(vertex2,2) + 0.1;
}

// this the FPS part
// int N_samples = 20;
// bool isVertexUsed[m_v];
// for (size_t i = 0; i < m_v; i++) {
//   isVertexUsed[i] = false;
// }
// isVertexUsed[0] = true;
// std::vector<int> fps(N_samples);
// fps[0] = 0;
// for (size_t i = 1; i < N_samples; i++) {
//   // int len_try = m_v-i;
//   double max = 0;
//   int max_place = 0;
//   // std::vector<double> distances(len_try);
//   for (size_t j = 0; j < m_v; j++) {
//     if(isVertexUsed[j])
//       continue;
//     double min_inner = INT_MAX;
//     for (size_t k = 0; k < i+1; k++) {
//       // dij_out my_dij = dijkstra(V_cost, fps[k], j);
//       dij_out my_dij = dijkstra_fibo_heap(my_pair, fps[k], j);
//       if(my_dij.min_dist < min_inner)
//         min_inner = my_dij.min_dist;
//     }
//     if(max < min_inner)
//     {
//       max = min_inner;
//       max_place = j;
//     }
//
//   }
//   fps[i] = max_place;
//   isVertexUsed[max_place] = true;
//   std::cout << "The iteration number is " << i << '\n';
//   std::cout << "The selected vertex is " << max_place << '\n';
// }
// // put the fps points into to P for displaying them
// Eigen::MatrixXd P(N_samples, 3);
// for(size_t i=0; i<N_samples; i++)
// {
//   int vertex = fps[i];
//   P(i,0) = V(vertex,0), P(i,1) = V(vertex,1), P(i,2) = V(vertex,2);
// }

  // start the geodesic iso-curve signature
  int k_iso = 20;
  double k_unit = max_geodesic/k_iso;
  std::vector<double> iso_curve_bin(k_iso, 0);
  for (size_t i = 0; i < m_f; i++) {
    int a = F(i,0), b = F(i,1), c = F(i,2);
    double k_a = V_geodesic(st, a);
    double k_b = V_geodesic(st, b);
    double k_c = V_geodesic(st, c);
    int radii_number = round((k_a + k_b + k_c)/(3*k_unit));
    double r_i = k_unit*radii_number;
    if((k_a < r_i && k_b > r_i && k_c > r_i) || (k_a > r_i && k_b < r_i && k_c < r_i))
    {
      double a1, a2;
      Eigen::VectorXd p1, p2;
      a1 = abs(r_i - k_a)/abs(k_b - k_a);
      a2 = abs(r_i - k_a)/abs(k_c - k_a);
      p1 = (1 - a1)*V.row(a) + a1*V.row(b);
      p2 = (1 - a2)*V.row(a) + a2*V.row(c);
      iso_curve_bin[radii_number] += (p1 - p2).norm();
      std::cout << "The norm is " << p2(1) << " " << p2(2) << '\n';
    }
    else if((k_b < r_i && k_a > r_i && k_c > r_i) || (k_b > r_i && k_a < r_i && k_c < r_i))
    {
      double a1, a2;
      Eigen::VectorXd p1, p2;
      a1 = abs(r_i - k_b)/abs(k_a - k_b);
      a2 = abs(r_i - k_b)/abs(k_c - k_b);
      p1 = (1 - a1)*V.row(b) + a1*V.row(a);
      p2 = (1 - a2)*V.row(b) + a2*V.row(c);
      iso_curve_bin[radii_number] += (p1 - p2).norm();
    }
    else if((k_c < r_i && k_a > r_i && k_b > r_i) || (k_c > r_i && k_a < r_i && k_b < r_i))
    {
      double a1, a2;
      Eigen::VectorXd p1, p2;
      a1 = abs(r_i - k_c)/abs(k_a - k_c);
      a2 = abs(r_i - k_c)/abs(k_b - k_c);
      p1 = (1 - a1)*V.row(c) + a1*V.row(a);
      p2 = (1 - a2)*V.row(c) + a2*V.row(b);
      iso_curve_bin[radii_number] += (p1 - p2).norm();
    }

  }

  std::cout << "ISO-Curve Signature: " << '\n';
  for (size_t i = 0; i < iso_curve_bin.size(); i++) {
    std::cout << iso_curve_bin[i] << " ";
  }
  std::cout << '\n';

  // Initialize the bilateral_hist with zeros
  std::vector<int> bilateral_hist(k_iso, 0);

  std::vector<std::vector<int> > iso_k_classes(k_iso);
  std::vector<int> roi_elementes;
  std::vector<bool> is_roi_element(m_v, false);
  double distance_threshold = max_geodesic/5;
  for (size_t i = 0; i < m_v; i++) {
    double min_geo_to_path = INT_MAX;
    for (size_t j = 0; j < my_out.route.size(); j++) {
      if(V_geodesic(i,j) < min_geo_to_path && i!=j)
        min_geo_to_path = V_geodesic(i,j);
    }
    if(min_geo_to_path < distance_threshold &&
      V_geodesic(i,st) < distance_threshold && V_geodesic(i,fin) < distance_threshold)
    {
      roi_elementes.push_back(i);
      is_roi_element[i] = true;
      int k_class = round(min_geo_to_path/k_unit);
      iso_k_classes[k_class].push_back(i);
      bilateral_hist[k_class]++;
    }
  }
  // iterate over the faces and compute bilateral histogram
  Eigen::VectorXd Z = Eigen::VectorXd::Zero(m_f);
  for (size_t i = 0; i < m_f; i++) {
    int a = F(i,0), b = F(i,1), c = F(i,2);
    if (is_roi_element[a] && is_roi_element[b] && is_roi_element[c]) {
      int k_class = round(distance_threshold/V_geodesic(st, a));
      Z(i) = distance_threshold/V_geodesic(st, a);
    }
  }

  igl::jet(Z,true,C);

  std::cout << "ISO-K classes" << '\n';
  for (size_t i = 0; i < 20; i++) {
    std::cout << bilateral_hist[i] << " ";
  }
  std::cout << '\n';
  // put the roi point into matrix P for visualization
  Eigen::MatrixXd P(roi_elementes.size(), 3);
  for(size_t i=0; i<roi_elementes.size(); i++)
  {
    int vertex = roi_elementes[i];
    P(i,0) = V(vertex,0), P(i,1) = V(vertex,1), P(i,2) = V(vertex,2);
  }


  // Plot the mesh
  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(V, F);
  viewer.data().add_edges(P1, P2, Eigen::RowVector3d(1,0,0));
  // viewer.data().add_points(P, Eigen::RowVector3d(0,0,1));
  viewer.data().set_colors(C);
  viewer.launch();
}

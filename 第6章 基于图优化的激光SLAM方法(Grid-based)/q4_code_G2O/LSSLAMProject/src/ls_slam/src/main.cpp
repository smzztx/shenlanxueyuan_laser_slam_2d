#include <gaussian_newton.h>
#include <readfile.h>

#include <ros/ros.h>
#include <visualization_msgs/MarkerArray.h>
#include <chrono>
#include <iostream>

#include <g2o/core/block_solver.h>
#include <g2o/core/optimization_algorithm_levenberg.h>
#include <g2o/core/optimization_algorithm_gauss_newton.h>
#include <g2o/core/optimization_algorithm_dogleg.h>
// #include <g2o/solvers/dense/linear_solver_dense.h>
// #include <g2o/solvers/eigen/linear_solver_eigen.h>
#include <g2o/solvers/csparse/linear_solver_csparse.h>
#include <g2o/types/slam2d/vertex_se2.h>
#include <g2o/types/slam2d/edge_se2.h>
#include <g2o/core/sparse_optimizer_terminate_action.h>

using namespace std;
using namespace g2o;
//for visual
void PublishGraphForVisulization(ros::Publisher* pub,
                                 std::vector<Eigen::Vector3d>& Vertexs,
                                 std::vector<Edge>& Edges,
                                 int color = 0)
{
    visualization_msgs::MarkerArray marray;

    //point--red
    visualization_msgs::Marker m;
    m.header.frame_id = "map";
    m.header.stamp = ros::Time::now();
    m.id = 0;
    m.ns = "ls-slam";
    m.type = visualization_msgs::Marker::SPHERE;
    m.pose.position.x = 0.0;
    m.pose.position.y = 0.0;
    m.pose.position.z = 0.0;
    m.scale.x = 0.1;
    m.scale.y = 0.1;
    m.scale.z = 0.1;

    if(color == 0)
    {
        m.color.r = 1.0;
        m.color.g = 0.0;
        m.color.b = 0.0;
    }
    else
    {
        m.color.r = 0.0;
        m.color.g = 1.0;
        m.color.b = 0.0;
    }

    m.color.a = 1.0;
    m.lifetime = ros::Duration(0);

    //linear--blue
    visualization_msgs::Marker edge;
    edge.header.frame_id = "map";
    edge.header.stamp = ros::Time::now();
    edge.action = visualization_msgs::Marker::ADD;
    edge.ns = "karto";
    edge.id = 0;
    edge.type = visualization_msgs::Marker::LINE_STRIP;
    edge.scale.x = 0.1;
    edge.scale.y = 0.1;
    edge.scale.z = 0.1;

    if(color == 0)
    {
        edge.color.r = 0.0;
        edge.color.g = 0.0;
        edge.color.b = 1.0;
    }
    else
    {
        edge.color.r = 1.0;
        edge.color.g = 0.0;
        edge.color.b = 1.0;
    }
    edge.color.a = 1.0;

    m.action = visualization_msgs::Marker::ADD;
    uint id = 0;

    //加入节点
    for (uint i=0; i<Vertexs.size(); i++)
    {
        m.id = id;
        m.pose.position.x = Vertexs[i](0);
        m.pose.position.y = Vertexs[i](1);
        marray.markers.push_back(visualization_msgs::Marker(m));
        id++;
    }

    //加入边
    for(int i = 0; i < Edges.size();i++)
    {
        Edge tmpEdge = Edges[i];
        edge.points.clear();

        geometry_msgs::Point p;
        p.x = Vertexs[tmpEdge.xi](0);
        p.y = Vertexs[tmpEdge.xi](1);
        edge.points.push_back(p);

        p.x = Vertexs[tmpEdge.xj](0);
        p.y = Vertexs[tmpEdge.xj](1);
        edge.points.push_back(p);
        edge.id = id;

        marray.markers.push_back(visualization_msgs::Marker(edge));
        id++;
    }

    pub->publish(marray);
}

int main(int argc, char **argv)
{
    ros::init(argc, argv, "ls_slam");

    ros::NodeHandle nodeHandle;

    // beforeGraph
    ros::Publisher beforeGraphPub,afterGraphPub;
    beforeGraphPub = nodeHandle.advertise<visualization_msgs::MarkerArray>("beforePoseGraph",1,true);
    afterGraphPub  = nodeHandle.advertise<visualization_msgs::MarkerArray>("afterPoseGraph",1,true);


    // std::string VertexPath = "/home/txcom-ubuntu64/LSSLAMProject/src/ls_slam/data/test_quadrat-v.dat";
    // std::string EdgePath = "/home/txcom-ubuntu64/LSSLAMProject/src/ls_slam/data/test_quadrat-e.dat";

   std::string VertexPath = "/home/txcom-ubuntu64/LSSLAMProject/src/ls_slam/data/intel-v.dat";
   std::string EdgePath = "/home/txcom-ubuntu64/LSSLAMProject/src/ls_slam/data/intel-e.dat";

//    std::string VertexPath = "/home/txcom-ubuntu64/LSSLAMProject/src/ls_slam/data/killian-v.dat";
//    std::string EdgePath = "/home/txcom-ubuntu64/LSSLAMProject/src/ls_slam/data/killian-e.dat";


    std::vector<Eigen::Vector3d> Vertexs;
    std::vector<Edge> Edges;

    ReadVertexInformation(VertexPath,Vertexs);
    ReadEdgesInformation(EdgePath,Edges);

    PublishGraphForVisulization(&beforeGraphPub,
                                Vertexs,
                                Edges);

    double initError = ComputeError(Vertexs,Edges);
    std::cout <<"initError:"<<initError<<std::endl;

    auto last_time = std::chrono::high_resolution_clock::now();
    typedef g2o::BlockSolver<g2o::BlockSolverTraits<3, 3>> Block; // 每个误差项优化变量维度为3，误差值维度为3

    // Block::LinearSolverType* linearSolver = new g2o::LinearSolverDense<Block::PoseMatrixType>();
    // Block::LinearSolverType* linearSolver = new g2o::LinearSolverEigen<Block::PoseMatrixType>();
    Block::LinearSolverType* linearSolver = new g2o::LinearSolverCSparse<Block::PoseMatrixType>();

    Block* solver_ptr = new Block(linearSolver);
    // g2o::OptimizationAlgorithmLevenberg* solver = new g2o::OptimizationAlgorithmLevenberg(solver_ptr);
    g2o::OptimizationAlgorithmGaussNewton* solver = new g2o::OptimizationAlgorithmGaussNewton( solver_ptr );
    // g2o::OptimizationAlgorithmDogleg* solver = new g2o::OptimizationAlgorithmDogleg( solver_ptr );
    g2o::SparseOptimizer optimizer;
    optimizer.setAlgorithm(solver);

    for (size_t i = 0; i < Vertexs.size(); i++) {
        VertexSE2* v = new VertexSE2();
        v->setEstimate(Vertexs[i]);
        v->setId(i);
        if (i == 0) {
            v->setFixed(true);
        }
        optimizer.addVertex(v);
    }

    for (size_t i = 0; i < Edges.size(); i++) {
        EdgeSE2* edge = new EdgeSE2();

        Edge tmpEdge = Edges[i];

        edge->setId(i);
        edge->setVertex(0, optimizer.vertices()[tmpEdge.xi]);
        edge->setVertex(1, optimizer.vertices()[tmpEdge.xj]);

        edge->setMeasurement(tmpEdge.measurement);
        edge->setInformation(tmpEdge.infoMatrix);
        optimizer.addEdge(edge);
    }

    optimizer.setVerbose(true);
    optimizer.initializeOptimization();
    SparseOptimizerTerminateAction* terminateAction = new SparseOptimizerTerminateAction;
    terminateAction->setGainThreshold(1e-4);
    optimizer.addPostIterationAction(terminateAction);
    optimizer.optimize(100);

    for (size_t i = 0; i < Vertexs.size(); i++) {
        VertexSE2* v = static_cast<VertexSE2*>(optimizer.vertices()[i]);
        Vertexs[i] = v->estimate().toVector();
    }
    auto current_time = std::chrono::high_resolution_clock::now();
    std::cout  << "all time(ms): " << std::chrono::duration_cast<std::chrono::milliseconds>(current_time - last_time).count() << std::endl;

    double finalError  = ComputeError(Vertexs,Edges);

    std::cout <<"FinalError:"<<finalError<<std::endl;

    PublishGraphForVisulization(&afterGraphPub,
                                Vertexs,
                                Edges,1);

    ros::spin();

    return 0;
}





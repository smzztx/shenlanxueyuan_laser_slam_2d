#include <gaussian_newton.h>
#include <readfile.h>

#include <ros/ros.h>
#include <visualization_msgs/MarkerArray.h>
#include <chrono>
#include <iostream>

#include <ceres/ceres.h>

using namespace std;
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

// Normalizes the angle in radians between [-pi and pi).
template <typename T>
inline T NormalizeAngle(const T& angle_radians) {
  // Use ceres::floor because it is specialized for double and Jet types.
  T two_pi(2.0 * M_PI);
  return angle_radians -
         two_pi * ceres::floor((angle_radians + T(M_PI)) / two_pi);
}

template <typename T>
Eigen::Matrix<T, 2, 2> RotationMatrix2D(T yaw_radians) {
  const T cos_yaw = ceres::cos(yaw_radians);
  const T sin_yaw = ceres::sin(yaw_radians);

  Eigen::Matrix<T, 2, 2> rotation;
  rotation << cos_yaw, -sin_yaw, sin_yaw, cos_yaw;
  return rotation;
}

// Computes the error term for two poses that have a relative pose measurement
// between them. Let the hat variables be the measurement.
//
// residual =  information^{1/2} * [  r_a^T * (p_b - p_a) - \hat{p_ab}   ]
//                                 [ Normalize(yaw_b - yaw_a - \hat{yaw_ab}) ]
//
// where r_a is the rotation matrix that rotates a vector represented in frame A
// into the global frame, and Normalize(*) ensures the angles are in the range
// [-pi, pi).
class PoseGraph2dErrorTerm {
 public:
  PoseGraph2dErrorTerm(double x_ab,
                       double y_ab,
                       double yaw_ab_radians,
                       const Eigen::Matrix3d& sqrt_information)
      : p_ab_(x_ab, y_ab),
        yaw_ab_radians_(yaw_ab_radians),
        sqrt_information_(sqrt_information) {}

  template <typename T>
  bool operator()(const T* const x_a,
                  const T* const y_a,
                  const T* const yaw_a,
                  const T* const x_b,
                  const T* const y_b,
                  const T* const yaw_b,
                  T* residuals_ptr) const {
    const Eigen::Matrix<T, 2, 1> p_a(*x_a, *y_a);
    const Eigen::Matrix<T, 2, 1> p_b(*x_b, *y_b);

    Eigen::Map<Eigen::Matrix<T, 3, 1>> residuals_map(residuals_ptr);

    residuals_map.template head<2>() =
        RotationMatrix2D(*yaw_a).transpose() * (p_b - p_a) - p_ab_.cast<T>();
    residuals_map(2) = NormalizeAngle(
        (*yaw_b - *yaw_a) - static_cast<T>(yaw_ab_radians_));

    // Scale the residuals by the square root information matrix to account for
    // the measurement uncertainty.
    residuals_map = sqrt_information_.template cast<T>() * residuals_map;

    return true;
  }

  static ceres::CostFunction* Create(double x_ab,
                                     double y_ab,
                                     double yaw_ab_radians,
                                     const Eigen::Matrix3d& sqrt_information) {
    return (new ceres::
                AutoDiffCostFunction<PoseGraph2dErrorTerm, 3, 1, 1, 1, 1, 1, 1>(
                    new PoseGraph2dErrorTerm(
                        x_ab, y_ab, yaw_ab_radians, sqrt_information)));
  }

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

 private:
  // The position of B relative to A in the A frame.
  const Eigen::Vector2d p_ab_;
  // The orientation of frame B relative to frame A.
  const double yaw_ab_radians_;
  // The inverse square root of the measurement covariance matrix.
  const Eigen::Matrix3d sqrt_information_;
};

// Defines a local parameterization for updating the angle to be constrained in
// [-pi to pi).
class AngleLocalParameterization {
 public:
  template <typename T>
  bool operator()(const T* theta_radians,
                  const T* delta_theta_radians,
                  T* theta_radians_plus_delta) const {
    *theta_radians_plus_delta =
        NormalizeAngle(*theta_radians + *delta_theta_radians);

    return true;
  }

  static ceres::LocalParameterization* Create() {
    return (new ceres::AutoDiffLocalParameterization<AngleLocalParameterization,
                                                     1,
                                                     1>);
  }
};

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

    // 构建最小二乘问题
    ceres::Problem problem;
    ceres::LossFunction* loss_function = NULL;
    ceres::LocalParameterization* angle_local_parameterization =
        AngleLocalParameterization::Create();
    for (int i=0; i < Edges.size(); ++i)
    {
        Edge tmpEdge = Edges[i];


        const Eigen::Matrix3d sqrt_information =
            tmpEdge.infoMatrix.llt().matrixL();
        // Ceres will take ownership of the pointer.
        ceres::CostFunction* cost_function = PoseGraph2dErrorTerm::Create(
            tmpEdge.measurement(0), tmpEdge.measurement(1), tmpEdge.measurement(2), sqrt_information);
        problem.AddResidualBlock(cost_function,
                                loss_function,
                                &Vertexs[tmpEdge.xi](0),
                                &Vertexs[tmpEdge.xi](1),
                                &Vertexs[tmpEdge.xi](2),
                                &Vertexs[tmpEdge.xj](0),
                                &Vertexs[tmpEdge.xj](1),
                                &Vertexs[tmpEdge.xj](2));

        problem.SetParameterization(&Vertexs[tmpEdge.xi](2),
                                    angle_local_parameterization);
        problem.SetParameterization(&Vertexs[tmpEdge.xj](2),
                                    angle_local_parameterization);
    }
    problem.SetParameterBlockConstant(&Vertexs[0](0));
    problem.SetParameterBlockConstant(&Vertexs[0](1));
    problem.SetParameterBlockConstant(&Vertexs[0](2));

    ceres::Solver::Options options;
    options.max_num_iterations = 100;
    options.linear_solver_type = ceres::SPARSE_NORMAL_CHOLESKY;
    options.gradient_tolerance = 10e-4;
    // options.function_tolerance = 10e-4;
    // options.parameter_tolerance = 10e-4;
    options.trust_region_strategy_type = ceres::DOGLEG;

    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);

    std::cout << summary.FullReport() << '\n';


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





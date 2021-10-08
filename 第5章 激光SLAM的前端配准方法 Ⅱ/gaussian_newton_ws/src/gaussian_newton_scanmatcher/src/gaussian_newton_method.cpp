#include <map.h>
#include "gaussian_newton_method.h"

const double GN_PI = 3.1415926;

//进行角度正则化．
double GN_NormalizationAngle(double angle)
{
    if(angle > GN_PI)
        angle -= 2*GN_PI;
    else if(angle < -GN_PI)
        angle += 2*GN_PI;

    return angle;
}

Eigen::Matrix3d GN_V2T(Eigen::Vector3d vec)
{
    Eigen::Matrix3d T;
    T  << cos(vec(2)),-sin(vec(2)),vec(0),
            sin(vec(2)), cos(vec(2)),vec(1),
            0,           0,     1;

    return T;
}

//对某一个点进行转换．
Eigen::Vector2d GN_TransPoint(Eigen::Vector2d pt,Eigen::Matrix3d T)
{
    Eigen::Vector3d tmp_pt(pt(0),pt(1),1);
    tmp_pt = T * tmp_pt;
    return Eigen::Vector2d(tmp_pt(0),tmp_pt(1));
}



//用激光雷达数据创建势场．
map_t* CreateMapFromLaserPoints(Eigen::Vector3d map_origin_pt,
                                std::vector<Eigen::Vector2d> laser_pts,
                                double resolution)
{
    map_t* map = map_alloc();

    map->origin_x = map_origin_pt(0);
    map->origin_y = map_origin_pt(1);
    map->resolution = resolution;

    //固定大小的地图，必要时可以扩大．
    map->size_x = 10000;
    map->size_y = 10000;

    map->cells = (map_cell_t*)malloc(sizeof(map_cell_t)*map->size_x*map->size_y);

    //高斯平滑的sigma－－固定死
    map->likelihood_sigma = 0.5;

    Eigen::Matrix3d Trans = GN_V2T(map_origin_pt);

    //设置障碍物
    for(int i = 0; i < laser_pts.size();i++)
    {
        Eigen::Vector2d tmp_pt = GN_TransPoint(laser_pts[i],Trans);

        int cell_x,cell_y;
        cell_x = MAP_GXWX(map,tmp_pt(0));
        cell_y = MAP_GYWY(map,tmp_pt(1));

        map->cells[MAP_INDEX(map,cell_x,cell_y)].occ_state = CELL_STATUS_OCC;
    }

    //进行障碍物的膨胀--最大距离固定死．
    map_update_cspace(map,0.5);

    return map;
}


/**
 * @brief InterpMapValueWithDerivatives
 * 在地图上的进行插值，得到coords处的势场值和对应的关于位置的梯度．
 * 返回值为Eigen::Vector3d ans
 * ans(0)表示市场值
 * ans(1:2)表示梯度
 * @param map
 * @param coords
 * @return
 */
Eigen::Vector3d InterpMapValueWithDerivatives(map_t* map,Eigen::Vector2d& coords)
{
    Eigen::Vector3d ans;
    //TODO
    // std::cout << "coords: " << coords << std::endl;
    // std::cout << "map: " << map->origin_x << std::endl;
    // std::cout << "map: " << map->origin_y << std::endl;
    // int map_index_x = MAP_GXWX(map, coords[0]);
    int map_index_x = floor((coords[0] - map->origin_x) / map->resolution) + map->size_x / 2;
    // std::cout << "map_index_x: " << map_index_x << std::endl;
    // int map_index_y = MAP_GYWY(map, coords[1]);
    int map_index_y = floor((coords[1] - map->origin_y) / map->resolution) + map->size_y / 2;
    // std::cout << "map_index_y: " << map_index_y << std::endl;
    if(MAP_VALID(map, map_index_x, map_index_y) && MAP_VALID(map, map_index_x + 1, map_index_y) && MAP_VALID(map, map_index_x + 1, map_index_y + 1) && MAP_VALID(map, map_index_x, map_index_y + 1))
    {
        int index1 = MAP_INDEX(map, map_index_x, map_index_y);
        double score1 = (map->cells+index1)->score;
        int index2 = MAP_INDEX(map, map_index_x + 1, map_index_y);
        double score2 = (map->cells+index2)->score;
        int index3 = MAP_INDEX(map, map_index_x + 1, map_index_y + 1);
        double score3 = (map->cells+index3)->score;
        int index4 = MAP_INDEX(map, map_index_x, map_index_y + 1);
        double score4 = (map->cells+index4)->score;

        double world_coords_x0 = MAP_WXGX(map, map_index_x);
        double world_coords_y0 = MAP_WYGY(map, map_index_y);
        double world_coords_x1 = MAP_WXGX(map, map_index_x + 1);
        double world_coords_y1 = MAP_WYGY(map, map_index_y + 1);
        // std::cout << "score: " << score1 << " " << score2 << " " << score3 << " " << score4 << std::endl;
        // std::cout << "world_coords: " << world_coords_x0 << " " << world_coords_x1 << " " << world_coords_y0 << " " << world_coords_y1 << std::endl;
        ans[0] = (coords[0] - world_coords_x1)/(world_coords_x0 - world_coords_x1) * (coords[1] - world_coords_y1)/(world_coords_y0 - world_coords_y1) * score1 + (coords[0] - world_coords_x0)/(world_coords_x1 - world_coords_x0) * (coords[1] - world_coords_y1)/(world_coords_y0 - world_coords_y1) * score2 + (coords[0] - world_coords_x0)/(world_coords_x1 - world_coords_x0) * (coords[1] - world_coords_y0)/(world_coords_y1 - world_coords_y0) * score3 + (coords[0] - world_coords_x1)/(world_coords_x0 - world_coords_x1) * (coords[1] - world_coords_y0)/(world_coords_y1 - world_coords_y0) * score4;
        ans[1] = 1/(world_coords_x0 - world_coords_x1) * (coords[1] - world_coords_y1)/(world_coords_y0 - world_coords_y1) * score1 + 1/(world_coords_x1 - world_coords_x0) * (coords[1] - world_coords_y1)/(world_coords_y0 - world_coords_y1) * score2 + 1/(world_coords_x1 - world_coords_x0) * (coords[1] - world_coords_y0)/(world_coords_y1 - world_coords_y0) * score3 + 1/(world_coords_x0 - world_coords_x1) * (coords[1] - world_coords_y0)/(world_coords_y1 - world_coords_y0) * score4;
        ans[2] = (coords[0] - world_coords_x1)/(world_coords_x0 - world_coords_x1) * 1/(world_coords_y0 - world_coords_y1) * score1 + (coords[0] - world_coords_x0)/(world_coords_x1 - world_coords_x0) * 1/(world_coords_y0 - world_coords_y1) * score2 + (coords[0] - world_coords_x0)/(world_coords_x1 - world_coords_x0) * 1/(world_coords_y1 - world_coords_y0) * score3 + (coords[0] - world_coords_x1)/(world_coords_x0 - world_coords_x1) * 1/(world_coords_y1 - world_coords_y0) * score4;
    }else{
        std::cerr << "coords is invalid, the given map coords lie without the absolute map bounds!!!" << std::endl;
    }
    //END OF TODO

    return ans;
}


/**
 * @brief ComputeCompleteHessianAndb
 * 计算H*dx = b中的H和b
 * @param map
 * @param now_pose
 * @param laser_pts
 * @param H
 * @param b
 */
void ComputeHessianAndb(map_t* map, Eigen::Vector3d now_pose,
                        std::vector<Eigen::Vector2d>& laser_pts,
                        Eigen::Matrix3d& H, Eigen::Vector3d& b)
{
    H = Eigen::Matrix3d::Zero();
    b = Eigen::Vector3d::Zero();

    //TODO
    for(auto point : laser_pts)
    {
        // std::cout << "ComputeHessianAndb start" << std::endl;
        Eigen::Matrix3d Trans = GN_V2T(now_pose);
        Eigen::Vector2d world_pt = GN_TransPoint(point, Trans);
        Eigen::Vector3d res = InterpMapValueWithDerivatives(map, world_pt);
        // std::cout << "ComputeHessianAndb end" << std::endl;
        // std::cout << "res: " << res << std::endl;
        Eigen::MatrixXd nabla_M(1, 2);
        nabla_M << res[1], res[2];
        Eigen::MatrixXd diff_S_T(2, 3);
        diff_S_T << 1, 0, -sin(now_pose[2])*point[0] - cos(now_pose[2])*point[1],
                    0, 1, cos(now_pose[2])*point[0] - sin(now_pose[2])*point[1];
        Eigen::MatrixXd J_tmp(1, 3);
        J_tmp = nabla_M * diff_S_T;
        H += J_tmp.transpose() * J_tmp;
        b += J_tmp.transpose() * (1 - res[0]);
    }
    //END OF TODO
}


/**
 * @brief GaussianNewtonOptimization
 * 进行高斯牛顿优化．
 * @param map
 * @param init_pose
 * @param laser_pts
 */
void GaussianNewtonOptimization(map_t*map,Eigen::Vector3d& init_pose,std::vector<Eigen::Vector2d>& laser_pts)
{
    int maxIteration = 20;
    Eigen::Vector3d now_pose = init_pose;

    for(int i = 0; i < maxIteration;i++)
    {
        //TODO
        Eigen::Matrix3d H;
        Eigen::Vector3d b;
        // std::cout << "GaussianNewtonOptimization start" << std::endl;
        ComputeHessianAndb(map, now_pose, laser_pts, H, b);
        // std::cout << "GaussianNewtonOptimization end" << std::endl;
        // Verify that H isn't singular
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(H);
        double cond = svd.singularValues()(0) / svd.singularValues()(svd.singularValues().size() - 1);
        if (cond <= 1000)
        {
            Eigen::Vector3d delta_T;
            delta_T = H.colPivHouseholderQr().solve(b);
            if(std::fabs(delta_T(2)) >= 0.17)
            {
                std::cout << "delta_T[2] is too large!!!" << std::endl;
                break;
            }
            now_pose += delta_T;
            //迭代条件是否满足
            if(std::sqrt( std::pow(delta_T(0),2) + std::pow(delta_T(1),2)) < 0.001 && delta_T(2) < (0.01/57.295))
            {
                std::cout << "delta_T is too small, break" << std::endl;
                break;
            }
        }else{
            std::cout << "Matrix H is almost singular." << "cond: " << cond << std::endl;
            break;
        }
        
        // if(H.determinant() != 0)
        // {
        //     // std::cout << "b: " << b << std::endl;
        //     // std::cout << "H: " << H << std::endl;
        //     Eigen::Vector3d delta_T;
        //     delta_T = H.colPivHouseholderQr().solve(b);
        //     // delta_T = H.inverse() * b;
        //     // Eigen::Matrix3d R;
        //     // double theta = now_pose(2);
        //     // R << cos(theta), -sin(theta), 0, 
        //     //     sin(theta),  cos(theta), 0,
        //     //     0,          0,      1;
        //     // // std::cout << "now_pose: " << now_pose << std::endl;
        //     // // std::cout << "delta_T: " << delta_T << std::endl;
        //     // now_pose += R * delta_T;
        //     now_pose += delta_T;
        // }else{
        //     std::cout << "the matrix H is singular!!!" << std::endl;
        // }
        
        
        //END OF TODO
    }
    init_pose = now_pose;

}

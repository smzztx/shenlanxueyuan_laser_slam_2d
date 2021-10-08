#include "occupany_mapping.h"
#include "nav_msgs/GetMap.h"
#include "sensor_msgs/PointCloud.h"
#include "sensor_msgs/PointCloud2.h"
#include "geometry_msgs/Point32.h"

/**
 * Increments all the grid cells from (x0, y0) to (x1, y1);
 * //不包含(x1,y1)
 * 2D画线算法　来进行计算两个点之间的grid cell
 * @param x0
 * @param y0
 * @param x1
 * @param y1
 */
std::vector<GridIndex> TraceLine(int x0, int y0, int x1, int y1)
{
    GridIndex tmpIndex;
    std::vector<GridIndex> gridIndexVector;

    bool steep = abs(y1 - y0) > abs(x1 - x0);
    if (steep)
    {
        std::swap(x0, y0);
        std::swap(x1, y1);
    }
    if (x0 > x1)
    {
        std::swap(x0, x1);
        std::swap(y0, y1);
    }

    int deltaX = x1 - x0;
    int deltaY = abs(y1 - y0);
    int error = 0;
    int ystep;
    int y = y0;

    if (y0 < y1)
    {
        ystep = 1;
    }
    else
    {
        ystep = -1;
    }

    int pointX;
    int pointY;
    for (int x = x0; x <= x1; x++)
    {
        if (steep)
        {
            pointX = y;
            pointY = x;
        }
        else
        {
            pointX = x;
            pointY = y;
        }

        error += deltaY;

        if (2 * error >= deltaX)
        {
            y += ystep;
            error -= deltaX;
        }

        //不包含最后一个点．
        if (pointX == x1 && pointY == y1)
            continue;

        //保存所有的点
        tmpIndex.SetIndex(pointX, pointY);

        gridIndexVector.push_back(tmpIndex);
    }

    return gridIndexVector;
}

void SetMapParams(void)
{
    mapParams.width = 2000;
    mapParams.height = 2000;
    mapParams.resolution = 0.05;

    //每次被击中的log变化值，覆盖栅格建图算法需要的参数
    mapParams.log_free = -1;
    mapParams.log_occ = 2;

    //每个栅格的最大最小值．
    mapParams.log_max = 100.0;
    mapParams.log_min = 0.0;

    mapParams.origin_x = 0.0;
    mapParams.origin_y = 0.0;

    //地图的原点，在地图的正中间
    mapParams.offset_x = 1000;
    mapParams.offset_y = 1000;

    pMap = new unsigned char[mapParams.width * mapParams.height];

    //计数建图算法需要的参数
    //每个栅格被激光击中的次数
    pMapHits = new unsigned long[mapParams.width * mapParams.height];
    //每个栅格被激光通过的次数
    pMapMisses = new unsigned long[mapParams.width * mapParams.height];

    //TSDF建图算法需要的参数
    pMapW = new unsigned long[mapParams.width * mapParams.height];
    pMapTSDF = new double[mapParams.width * mapParams.height];

    //初始化
    for (int i = 0; i < mapParams.width * mapParams.height; i++)
    {
        pMap[i] = 50;
        pMapHits[i] = 0;
        pMapMisses[i] = 0;
        pMapW[i] = 0;
        pMapTSDF[i] = -1;
    }
}

//从世界坐标系转换到栅格坐标系
GridIndex ConvertWorld2GridIndex(double x, double y)
{
    GridIndex index;

    index.x = std::ceil((x - mapParams.origin_x) / mapParams.resolution) + mapParams.offset_x;
    index.y = std::ceil((y - mapParams.origin_y) / mapParams.resolution) + mapParams.offset_y;

    return index;
}

int GridIndexToLinearIndex(GridIndex index)
{
    int linear_index;
    linear_index = index.y + index.x * mapParams.width;
}

//判断index是否有效
bool isValidGridIndex(GridIndex index)
{
    if (index.x >= 0 && index.x < mapParams.width && index.y >= 0 && index.y < mapParams.height)
        return true;

    return false;
}

void DestoryMap()
{
    if (pMap != NULL)
        delete pMap;
}

//
void OccupanyMapping(std::vector<GeneralLaserScan> &scans, std::vector<Eigen::Vector3d> &robot_poses)
{
    std::cout << "开始建图，请稍后..." << std::endl;
    //枚举所有的激光雷达数据
    for (int i = 0; i < scans.size(); i++)
    {
        GeneralLaserScan scan = scans[i];
        Eigen::Vector3d robotPose = robot_poses[i];

        //机器人的下标
        GridIndex robotIndex = ConvertWorld2GridIndex(robotPose(0), robotPose(1));

        for (int id = 0; id < scan.range_readings.size(); id++)
        {
            double dist = scan.range_readings[id];
            double angle = -scan.angle_readings[id]; // 激光雷达逆时针转，角度取反

            if (std::isinf(dist) || std::isnan(dist))
                continue;

            //计算得到该激光点的世界坐标系的坐标
            double theta = -robotPose(2); // 激光雷达逆时针转，角度取反
            double laser_x = dist * cos(angle);
            double laser_y = dist * sin(angle);

            double world_x = cos(theta) * laser_x - sin(theta) * laser_y + robotPose(0);
            double world_y = sin(theta) * laser_x + cos(theta) * laser_y + robotPose(1);

            //start of TODO 对对应的map的cell信息进行更新．（1,2,3题内容）
            #define QUESTION_3  //通过宏定义切换问题
            //1
            #ifdef QUESTION_1
                GridIndex beamPointIndex = ConvertWorld2GridIndex(world_x, world_y);
                std::vector<GridIndex> beamTraceindexes = TraceLine(robotIndex.x, robotIndex.y, beamPointIndex.x, beamPointIndex.y);
                for(auto index : beamTraceindexes)
                {
                    if(isValidGridIndex(index))
                    {
                        int tmpLinearIndex = GridIndexToLinearIndex(index);
                        if(pMap[tmpLinearIndex] == 0) continue;
                        pMap[tmpLinearIndex] += mapParams.log_free;
                    }else{
                        std::cerr << "index is invalid!!!" << std::endl;
                    }
                }
                if(isValidGridIndex(beamPointIndex))
                {
                    int tmpLinearIndex = GridIndexToLinearIndex(beamPointIndex);
                    
                    pMap[tmpLinearIndex] += mapParams.log_occ;
                    if(pMap[tmpLinearIndex] >= 100) pMap[tmpLinearIndex] = 100;
                }else{
                    std::cerr << "beamPointIndex is invalid!!!" << std::endl;
                }
            #endif

            //2
            #ifdef QUESTION_2
                GridIndex beamPointIndex = ConvertWorld2GridIndex(world_x, world_y);
                std::vector<GridIndex> beamTraceindexes = TraceLine(robotIndex.x, robotIndex.y, beamPointIndex.x, beamPointIndex.y);
                for(auto index : beamTraceindexes)
                {
                    if(isValidGridIndex(index))
                    {
                        int tmpLinearIndex = GridIndexToLinearIndex(index);
                        ++pMapMisses[tmpLinearIndex];
                    }else{
                        std::cerr << "index is invalid!!!" << std::endl;
                    }
                }
                if(isValidGridIndex(beamPointIndex))
                {
                    int tmpLinearIndex = GridIndexToLinearIndex(beamPointIndex);
                    ++pMapHits[tmpLinearIndex];
                }else{
                    std::cerr << "beamPointIndex is invalid!!!" << std::endl;
                }
            #endif

            //3
            #ifdef QUESTION_3
                double far_dist = dist + 0.142; //0.05*2*sqrt(2)
                //计算得到该激光点的世界坐标系的坐标
                double far_laser_x = far_dist * cos(angle);
                double far_laser_y = far_dist * sin(angle);

                double far_world_x = cos(theta) * far_laser_x - sin(theta) * far_laser_y + robotPose(0);
                double far_world_y = sin(theta) * far_laser_x + cos(theta) * far_laser_y + robotPose(1);
                GridIndex farBeamPointIndex = ConvertWorld2GridIndex(far_world_x, far_world_y);
                std::vector<GridIndex> farBeamTraceindexes = TraceLine(robotIndex.x, robotIndex.y, farBeamPointIndex.x, farBeamPointIndex.y);
                for(auto index : farBeamTraceindexes)
                {
                    if(isValidGridIndex(index))
                    {
                        //栅格坐标系转换到世界坐标系
                        double x = (index.x - mapParams.offset_x) * mapParams.resolution + mapParams.origin_x;
                        double y = (index.y - mapParams.offset_y) * mapParams.resolution + mapParams.origin_y;
                        double d = std::sqrt((x-robotPose(0))*(x-robotPose(0)) + (y-robotPose(1))*(y-robotPose(1)));
                        double sdf = dist - d;
                        double tsdf = std::max(-1.0, std::min(1.0, sdf/0.1));
                        int tmpLinearIndex = GridIndexToLinearIndex(index);
                        pMapTSDF[tmpLinearIndex] = (pMapW[tmpLinearIndex]*pMapTSDF[tmpLinearIndex] + tsdf) / (pMapW[tmpLinearIndex] + 1);
                        pMapW[tmpLinearIndex] = pMapW[tmpLinearIndex] + 1;
                    }else{
                        std::cerr << "index is invalid!!!" << std::endl;
                    }
                }
            #endif
            //end of TODO
        }
    }
    //start of TODO 通过计数建图算法或TSDF算法对栅格进行更新（2,3题内容）
    //2
    #ifdef QUESTION_2
        for (int i = 0; i < mapParams.width * mapParams.height; ++i)
        {
            if((pMapHits[i] + pMapMisses[i]) != 0 )
            {
                pMap[i] = (double)pMapHits[i]/(pMapHits[i] + pMapMisses[i]) * 100;
                if(pMap[i] >=35) pMap[i] = 100;
                // if(pMapHits[i] != 0)
                //     std::cout << pMapHits[i] << " " << pMapMisses[i] << " " << (int)pMap[i] << std::endl;
            }
        }
    #endif

    //3
    #ifdef QUESTION_3
        //test code
        // for (int i = 0; i < mapParams.width * mapParams.height; ++i)
        // {
        //     pMap[i] = pMapTSDF[i] * 100;
        // }

        // for (int i = 0; i < mapParams.width; ++i)
        // {
        //     for(int j = 0; j < mapParams.height-1; ++j)
        //     {
        //         GridIndex tmp1Index;
        //         tmp1Index.SetIndex(i, j);
        //         int tmp1LinearIndex = GridIndexToLinearIndex(tmp1Index);

        //         GridIndex tmp2Index;
        //         tmp2Index.SetIndex(i, j+1);
        //         int tmp2LinearIndex = GridIndexToLinearIndex(tmp2Index);

        //         double tmp1TSDF = pMapTSDF[tmp1LinearIndex];
        //         double tmp2TSDF = pMapTSDF[tmp2LinearIndex];
        //         if((tmp1TSDF<0 && tmp2TSDF>0) || (tmp1TSDF>0 && tmp2TSDF<0))
        //         {
        //             if(std::fabs(tmp1TSDF) < std::fabs(tmp2TSDF))
        //             {
        //                 pMap[tmp1LinearIndex] = 100;
        //             }else{
        //                 pMap[tmp2LinearIndex] = 100;
        //             }
        //         }
        //     }
        // }
        // for (int i = 0; i < mapParams.height; ++i)
        // {
        //     for(int j = 0; j < mapParams.width-1; ++j)
        //     {
        //         GridIndex tmp1Index;
        //         tmp1Index.SetIndex(j, i);
        //         int tmp1LinearIndex = GridIndexToLinearIndex(tmp1Index);

        //         GridIndex tmp2Index;
        //         tmp2Index.SetIndex(j+1, i);
        //         int tmp2LinearIndex = GridIndexToLinearIndex(tmp2Index);

        //         double tmp1TSDF = pMapTSDF[tmp1LinearIndex];
        //         double tmp2TSDF = pMapTSDF[tmp2LinearIndex];
        //         if((tmp1TSDF<0 && tmp2TSDF>0) || (tmp1TSDF>0 && tmp2TSDF<0))
        //         {
        //             if(std::fabs(tmp1TSDF) < std::fabs(tmp2TSDF))
        //             {
        //                 pMap[tmp1LinearIndex] = 100;
        //             }else{
        //                 pMap[tmp2LinearIndex] = 100;
        //             }
        //         }
        //     }
        // }
        for (int i = 0; i < mapParams.width-1; ++i) //x
        {
            for(int j = 0; j < mapParams.height-1; ++j) //y
            {
                GridIndex tmpOrgIndex;
                tmpOrgIndex.SetIndex(i, j);
                int tmpOrgLinearIndex = GridIndexToLinearIndex(tmpOrgIndex);
                double tmpOrgTSDF = pMapTSDF[tmpOrgLinearIndex];
                if(tmpOrgTSDF==1 || tmpOrgTSDF==-1) continue;   //去除未击中点

                GridIndex tmpUpIndex;
                tmpUpIndex.SetIndex(i, j+1);
                int tmpUpLinearIndex = GridIndexToLinearIndex(tmpUpIndex);
                double tmpUpTSDF = pMapTSDF[tmpUpLinearIndex];
                // if(tmpUpTSDF==1 || tmpUpTSDF==-1) continue;   //去除未击中点

                GridIndex tmpRightIndex;
                tmpRightIndex.SetIndex(i+1, j);
                int tmpRightLinearIndex = GridIndexToLinearIndex(tmpRightIndex);
                double tmpRightTSDF = pMapTSDF[tmpRightLinearIndex];
                // if(tmpRightTSDF==1 || tmpRightTSDF==-1) continue;   //去除未击中点
                // if(((tmpOrgTSDF<0 && tmpUpTSDF>0) || (tmpOrgTSDF>0 && tmpUpTSDF<0)) && ((tmpOrgTSDF<0 && tmpRightTSDF>0) || (tmpOrgTSDF>0 && tmpRightTSDF<0)))
                // {
                //     double x, y;
                //     if(std::fabs(tmpOrgTSDF) < std::fabs(tmpUpTSDF))
                //     {
                //         y = j;
                //     }else{
                //         y = j+1;
                //     }
                //     if(std::fabs(tmpOrgTSDF) < std::fabs(tmpRightTSDF))
                //     {
                //         x = i;
                //     }else{
                //         x = i + 1;
                //     }
                //     GridIndex tmpIndex;
                //     tmpIndex.SetIndex(x, y);
                //     int tmpLinearIndex = GridIndexToLinearIndex(tmpIndex);
                //     pMap[tmpLinearIndex] = 100;
                //     continue;
                // }
                if((tmpOrgTSDF<0 && tmpUpTSDF>0) || (tmpOrgTSDF>0 && tmpUpTSDF<0))
                {
                    if(std::fabs(tmpOrgTSDF) < std::fabs(tmpUpTSDF))
                    {
                        pMap[tmpOrgLinearIndex] = 100;
                    }else{
                        pMap[tmpUpLinearIndex] = 100;
                    }
                }
                if((tmpOrgTSDF<0 && tmpRightTSDF>0) || (tmpOrgTSDF>0 && tmpRightTSDF<0))
                {
                    if(std::fabs(tmpOrgTSDF) < std::fabs(tmpRightTSDF))
                    {
                        pMap[tmpOrgLinearIndex] = 100;
                    }else{
                        pMap[tmpRightLinearIndex] = 100;
                    }
                }
            }
        }
    #endif
    //end of TODO
    std::cout << "建图完毕" << std::endl;
}

//发布地图．
void PublishMap(ros::Publisher &map_pub)
{
    nav_msgs::OccupancyGrid rosMap;

    rosMap.info.resolution = mapParams.resolution;
    rosMap.info.origin.position.x = 0.0;
    rosMap.info.origin.position.y = 0.0;
    rosMap.info.origin.position.z = 0.0;
    rosMap.info.origin.orientation.x = 0.0;
    rosMap.info.origin.orientation.y = 0.0;
    rosMap.info.origin.orientation.z = 0.0;
    rosMap.info.origin.orientation.w = 1.0;

    rosMap.info.origin.position.x = mapParams.origin_x;
    rosMap.info.origin.position.y = mapParams.origin_y;
    rosMap.info.width = mapParams.width;
    rosMap.info.height = mapParams.height;
    rosMap.data.resize(rosMap.info.width * rosMap.info.height);

    //0~100
    int cnt0, cnt1, cnt2;
    cnt0 = cnt1 = cnt2 = 100;
    for (int i = 0; i < mapParams.width * mapParams.height; i++)
    {
        if (pMap[i] == 50)
        {
            rosMap.data[i] = -1.0;
        }
        else
        {

            rosMap.data[i] = pMap[i];
        }
    }

    rosMap.header.stamp = ros::Time::now();
    rosMap.header.frame_id = "map";

    map_pub.publish(rosMap);
}

int main(int argc, char **argv)
{
    ros::init(argc, argv, "OccupanyMapping");

    ros::NodeHandle nodeHandler;

    ros::Publisher mapPub = nodeHandler.advertise<nav_msgs::OccupancyGrid>("laser_map", 1, true);

    std::vector<Eigen::Vector3d> robotPoses;
    std::vector<GeneralLaserScan> generalLaserScans;

    std::string basePath = "/home/txcom-ubuntu64/OccupanyMappingProject/src/data";

    std::string posePath = basePath + "/pose.txt";
    std::string anglePath = basePath + "/scanAngles.txt";
    std::string scanPath = basePath + "/ranges.txt";

    //读取数据
    ReadPoseInformation(posePath, robotPoses);

    ReadLaserScanInformation(anglePath,
                             scanPath,
                             generalLaserScans);

    //设置地图信息
    SetMapParams();

    OccupanyMapping(generalLaserScans, robotPoses);

    PublishMap(mapPub);

    ros::spin();

    DestoryMap();

    std::cout << "Release Memory!!" << std::endl;
}

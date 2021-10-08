### install
安装自定义消息到系统，如果你的ROS版本不是kinetic，请将kinetic改成你的ROS版本名

   `sudo su`

   `source /opt/ros/kinetic/setup.bash`

   `catkin_make -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/opt/ros/kinetic install`

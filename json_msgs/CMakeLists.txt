cmake_minimum_required(VERSION 2.8.3)
project(json_msgs)

find_package(catkin REQUIRED COMPONENTS
  std_msgs
  actionlib_msgs
  message_generation
)

add_message_files(FILES JsonMsg.msg)
add_service_files(FILES JsonSrv.srv)
add_action_files(FILES Json.action)

generate_messages(DEPENDENCIES actionlib_msgs std_msgs)

catkin_package(CATKIN_DEPENDS message_runtime)

include_directories(${catkin_INCLUDE_DIRS})
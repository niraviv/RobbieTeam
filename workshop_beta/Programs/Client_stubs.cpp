#include "stdafx.h"
#include "Programs\Client_stubs.h"

#include "Utils\Communication_utils\socket.h"
#include "Utils\Communication_utils\client_utils.h"

#include "Mms_example.h"
#include "Path_planning\Motion_sequence.h"

#include <boost/thread/thread.hpp>

//global typedefs used in this file
typedef Environment<>                     Env;
typedef mms::Mms_path_planner_example<>   Planner;
typedef Motion_sequence<Planner::K>       Motion;
typedef Env::Reference_point              Ref_p;

//globals for the sake of the example - bad programing!
Socket_client*            socket_client_ptr = NULL;
bool                      finished_game = false;

void dbg_log(char* function_name, char* action)
{
  std::cout << "client, thread id " <<boost::this_thread::get_id()<<" " 
            << "in function " << function_name
            << " " <<action
            <<std::endl;
}

void plan(double remaining_time)
{
  dbg_log("plan", "entering");
  //here you should implement your planner
  //instead we just sleep for the remaining time

  if (finished_game)  //update this flag when you finished all queries
    return;

  boost::posix_time::seconds sleep_time(remaining_time);
  boost::this_thread::sleep(sleep_time);

  dbg_log("plan", "exiting");
  return;
}
void construct_motion(double remaining_time, Motion& motion_sequence)
{
  dbg_log("construct_motion", "entering");
  //construct a motion, 
  //make sure that you have enough time to complete the motion 
  //or else no write approval will granted

  if (finished_game)
    return;

  //we just sleep this time ...
  boost::posix_time::seconds sleep_time(remaining_time);
  boost::this_thread::sleep(sleep_time);

  dbg_log("construct_motion", "exiting");
  return;
}

void move(double remaining_time)
{
  dbg_log("move", "entering");
  
  if (finished_game)
    return;

  //here you should construct a motion
  Motion motion_sequence;
  construct_motion(remaining_time, motion_sequence);
  
  double motion_length(motion_sequence.motion_time( configuration.get_translational_speed(),
                                                    configuration.get_rotational_speed()));
  if (motion_length == 0) //no motion
    return;
  std::string path_filename;

  dbg_log("move", "requesting to write");
  if (request_to_write(*socket_client_ptr, motion_length, path_filename))
  {
    dbg_log("move", "request granted");
    //request granted
    std::ofstream out(path_filename.c_str(), std::ios::app); //it is essential to append
    motion_sequence.write(out);    
  }
  else
  {
    dbg_log("move", "request denied");
  }
  
  dbg_log("move", "exiting");
  return ;
}
void moveable_planner(double remaining_time)
{
  dbg_log("moveable_planner", "entering");
  
  CGAL::Timer timer;
  timer.start();
  while (timer.time() < remaining_time)
  {
    //plan and then move (if time remains)
    plan(remaining_time - timer.time());
    if (remaining_time > timer.time())
      move(remaining_time - timer.time());
  }

  dbg_log("moveable_planner", "exiting");
  return;
}

void static_planner(double remaining_time)
{
  dbg_log("static_planner", "entering");

  CGAL::Timer timer;
  timer.start();
  while (timer.time() < remaining_time)
  {
    //plan
    plan(remaining_time - timer.time());
  }

  dbg_log("static_planner", "exiting");
  return;
}

void client_stubs_main(int argc, char* argv[])
{
  ////////////////////////////////////////////////////////////
  
  Env env(argc,argv);
  Input_reader input_reader;
  

  ////////////////////////////////////////////////////////////
  //loading scene from environment
  //should be used in stub functions:
  //  plan, move, construct_motion
  Env::Polygon_vector&      workspace(env.get_workspace());
  Env::Extended_polygon&    robot(env.get_robot_a());
  Env::Extended_polygon&    dynamic_obstacle(env.get_robot_b());

  Env::Reference_point      source_ref_p(env.get_source_configiration_a());
  Env::Reference_point_vec* target_configurations_ptr = &(env.get_target_configirations());

  Env::Reference_point      dynamic_obstacle_configiration(env.get_source_configiration_b());
    
  ////////////////////////////////////////////////////////////
  //initialize connection to server
  socket_client_ptr = new Socket_client(configuration.get_host_name(), configuration.get_host_port());

  
  ////////////////////////////////////////////////////////////
  //begin planning
  bool read_additional_configurations = true;
  while (finished_game == false)
  {    
    //before doing any planning, see if something changed
    dbg_log("main client", "getting scene status");
    Scene_status scene_status = get_scene_status(*socket_client_ptr); 
    dynamic_obstacle_configiration = 
      input_reader.read_reference_point<Rational_kernel>(scene_status.quasi_dynamic_obs_location_filename);
    if (scene_status.updated_target_configurations && (read_additional_configurations) )
    {
      input_reader.read_reference_points<Rational_kernel>(scene_status.updated_target_configurations_filename,
                                                          std::back_inserter(*target_configurations_ptr));
      read_additional_configurations = false;
    }   
    
    //now plan
    dbg_log("main client", "planning");
    Time_frame_status time_frame_status = get_time_frame_status(*socket_client_ptr); 
    if (time_frame_status.is_moveable)
      moveable_planner(time_frame_status.remaining_time);
    else //(time_frame_status.is_moveable == false)
      static_planner(time_frame_status.remaining_time);

    std::cout <<std::endl<<std::endl;
  }

  std::cout <<"FINISHED..."<<std::endl;

  terminate_connection(*socket_client_ptr);
  delete (socket_client_ptr);

  return;
}


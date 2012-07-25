#include "stdafx.h"
#include "Programs\Path_planner.h"

#include "Utils\Rotation_utils\Rotation.h"
#include "Utils\ReferencePoint.h"

#include "Utils\UI_utils\Environment.h"
#include "Utils\UI_utils\TimeManager.h"
#include "Utils\UI_utils\Graphics.h"

#include "Utils\Number_utils\AK_conversions_1.h"
#include "Utils\Polygon_utils\PolygonIO.h"
#include "Path_planning\Motion_sequence.h"
#include "Path_planning\Motion_sequence_gui_converter.h"
#include "Mms_example.h"


void display_scene(int argc, char* argv[])
{
  //typedefs
  typedef mms::Mms_path_planner_example<>   Planner;

  //loading files from configuration file
  Time_manager tm;
  tm.write_time_log(std::string("start"));
  
  Environment<> env(argc,argv);
  tm.write_time_log(std::string("set environment"));

  //loading scene from environment
  Planner::Polygon_vec&       workspace(env.get_workspace());
  Planner::Extended_polygon   my_robot(env.get_robot_a());	

  //load query
  Planner::Reference_point q_s (env.get_source_configiration_a());
  Planner::Reference_point q_t (env.get_target_configirations().front());

  global_graphics->draw_polygons<Planner::K> (workspace, BLUE);
  my_robot.move_absolute(q_s);
  global_graphics->draw_polygon<Planner::K> (my_robot.get_absolute_polygon(), GREEN);
  my_robot.move_absolute(q_t);
  global_graphics->draw_polygon<Planner::K> (my_robot.get_absolute_polygon(), RED);
  global_graphics->display();
  return;
}
void single_robot_planner_example(int argc, char* argv[])
{
  //typedefs
  typedef mms::Mms_path_planner_example<>             Planner;
  typedef Motion_sequence_gui_converter<Planner::K>   Converter;

  //loading files from configuration file
  Time_manager tm;
  tm.write_time_log(std::string("start"));
  
  Environment<> env(argc,argv);
  tm.write_time_log(std::string("set environment"));

  //loading scene from environment
  Planner::Polygon_vec&       workspace(env.get_workspace());
  Planner::Extended_polygon   my_robot(env.get_robot_a());	
  
  //initialize the planner and preprocess 
  Planner planner(workspace, my_robot);    
  planner.preprocess(); 
   
  //load query
  Planner::Reference_point source(env.get_source_configiration_a());

  //NIR: now we store all the target configurations instead of just the front,
  //and perform a multiple-target query.
  std::vector<Planner::Reference_point> targets = env.get_target_configirations();
  int number_of_targets = targets.size();
  std::vector<Planner::Motion_sequence> motion_sequences(number_of_targets);
  if (!planner.query(source, targets, motion_sequences)) {
    std::cout << "no paths found :-(" << std::endl;
	return;
  }
  std::cout << "at least one path was found :-)" << std::endl;

  //Find the first nonempty motion sequence:
  Planner::Motion_sequence shortest_motion_sequence = motion_sequences[0];
  int j(0);
  while (shortest_motion_sequence.get_sequence().empty()) {
	  shortest_motion_sequence = motion_sequences[++j];
  }
  double translational_speed = configuration.get_translational_speed(),
	  rotational_speed = configuration.get_rotational_speed();
  double shortest_motion_time = shortest_motion_sequence.motion_time(translational_speed, rotational_speed);
  Planner::Reference_point nearest_target = targets[j];

  //See if there are any shorter nonempty motion sequences:
  for (j++; j < number_of_targets; j++) {
	  Planner::Motion_sequence current_motion_sequence = motion_sequences[j];
	  if (current_motion_sequence.get_sequence().empty()) {
		  continue;
	  }
	  double current_motion_time = current_motion_sequence.motion_time(translational_speed, rotational_speed);
	  if (current_motion_time < shortest_motion_time) {
		  shortest_motion_time = current_motion_time;
		  shortest_motion_sequence = current_motion_sequence;
		  nearest_target = targets[j];
	  }
  }

  //example of how to create a path that can be loaded by the GUI
  Converter::Path_3d path;
  Converter converter;
  converter(shortest_motion_sequence, std::back_inserter(path));
  
  std::ofstream os("path.txt");
  os << path.size() << std::endl;
  for (unsigned int i(0); i < path.size(); ++i) {
	  if (i == 0) {
		  os  << path[i].x << " " << path[i].y << " " << path[i].t << " "
			  << CGAL::to_double(nearest_target.get_location().x()) << " " //we place the second robot at the target for ease of visualization
			  << CGAL::to_double(nearest_target.get_location().y()) << " " //we place the second robot at the target for ease of visualization
			  << nearest_target.get_rotation().to_angle()                  //we place the second robot at the target for ease of visualization
			  << std::endl;
	  } else {		//i >0
		  os  << path[i].x << " " << path[i].y << " " << path[i].t << " "  //first robot location
			  << 0 << " " << 0 << " " << 0 << " "                          //second robot is static
			  << 2                                                         //speed is arbitrarily chosen to be 2
			  << std::endl;
	  }
  }
  return;
}
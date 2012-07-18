#ifndef MMS_EXAMPLE_H
#define MMS_EXAMPLE_H

#include "Manifolds\MMSTypedefs.h"
#include "Manifolds\Fsc_indx.h"
#include "Manifolds\Fixed_angle\Fixed_angle_manifold_container.h"
#include "Manifolds\Fixed_point\Fixed_point_manifold_container.h"
#include "Manifolds\Intersect_manifolds.h"
#include "Path_planning\PathPlanningUtils.h"
#include "heuristic_utils.h"
#include "FSC.h"
#include "Graph\Graph.h"
#include "Fsc_path_planning.h"
#include "CGAL/intersections.h"
#include "CGAL\Cartesian\function_objects.h"

namespace mms{

	template <typename K_ = Rational_kernel, 
		typename AK_ = CGAL::Algebraic_kernel_d_1 <typename CGAL::CORE_arithmetic_kernel::Integer>, 
		typename AK_conversions = Algebraic_kernel_d_1_conversions_rational<typename CGAL::CORE_arithmetic_kernel> >
	class Mms_path_planner_example
	{
	public:
		typedef typename K_                             K;
		typedef typename AK_                            AK;

		typedef typename K::Point_2                     Point;
		typedef Rotation<typename K::FT>                Rotation;
		typedef typename Reference_point<K>             Reference_point;
		typedef Rotation_range_absolute<typename K::FT> Rotation_range;

		typedef typename K::Segment_2                   Segment;
		typedef typename K::Line_2                      Line;
		typedef typename K::Ray_2                       Ray;

		typedef CGAL::Polygon_2 <K>                     Polygon;
		typedef CGAL::Polygon_with_holes_2<K>           Polygon_with_holes;
		typedef typename Extended_polygon<K>            Extended_polygon;
		typedef typename Smart_polygon_with_holes<K>    Smart_polygon;

		typedef std::vector<typename Polygon>           Polygon_vec;
		typedef CGAL::Polygon_set_2<K>                  Polygon_set;
		typedef std::vector<Reference_point>      Reference_point_vec;

		typedef Motion_step_rotational<K>               Motion_step_rotational;
		typedef Motion_step_translational<K>            Motion_step_translational;
		typedef Motion_sequence<K>                      Motion_sequence;


	public:
		typedef Fsc_indx<K>                             Fsc_indx;
		typedef FSC<K, AK, AK_conversions>              Fsc;

		typedef Fixed_angle_manifold_container<K>       Layers;
		typedef typename Layers::Manifold               Layer;

		typedef Fixed_point_manifold_container<K, AK, AK_conversions> C_space_lines;
		typedef typename C_space_lines::Manifold                      C_space_line;

		typedef Graph<Fsc_indx, Less_than_fsc_indx<K> > Connectiivity_graph;
		typedef Random_utils<K>                         Random_utils;  
	private:
		Polygon_vec&            _workspace;
		Polygon_vec             _decomposed_workspace;
		CGAL::Bbox_2            _workspace_bbox;
		Extended_polygon&       _robot;

		Connectiivity_graph     _graph;

		std::vector<Rotation>   _rotations;
		Layers                  _layers;
		C_space_lines           _lines;

		Random_utils            _rand;
		AK                      _ak;
	public:
		//constructor
		Mms_path_planner_example  (Polygon_vec &workspace, Extended_polygon& robot)
			: _workspace (workspace), _robot(robot),
			_graph(0,true), _rand(time(NULL))
		{
			compute_workspace_bbox();
			//the minkowski sum algorithm works faster on polygons are convex
			//hence we decompose the workspace obstacles to convex polygons
			decompose_workspace_into_convex_polygons();
		}
		//preprocess
		void preprocess (Reference_point start_point, Reference_point_vec target_points, const unsigned int num_of_angles = configuration.get_slices_granularity())
		{
			generate_rotations(num_of_angles);

			BOOST_FOREACH (Rotation rotation, _rotations)
				add_layer(rotation);
			//layers for start and targets
			add_layer(start_point.get_rotation());
			for (int i=0;i<target_points.size();i++)
				add_layer(target_points[i].get_rotation());
			//till here
			global_tm.write_time_log(std::string("finished layers"));

			generate_connectors();    
			//connectors for start and targets
			generate_connector_at_point(start_point.get_location());
			for (int i=0;i<target_points.size();i++)
				generate_connector_at_point(target_points[i].get_location());
			//till here
			global_tm.write_time_log(std::string("finished connectors"));

			global_tm.write_time_log(std::string("finished preproccesing"));
			return;
		}
		//query
		bool query( const Reference_point& source, const Reference_point& target,
			Motion_sequence& motion_sequence) 
		{
			////////////////////////////////////
			//connect source and target to graph
			////////////////////////////////////
			Motion_sequence source_motion_sequence, target_motion_sequence;
			Reference_point perturbed_source = connect_to_graph(source, source_motion_sequence);
			Reference_point perturbed_target = connect_to_graph(target, target_motion_sequence);

			if (perturbed_source ==  Reference_point() || 
				perturbed_target ==  Reference_point())
			{
				std::cout <<"failed to connect to pre-processed configuration space"<<std::endl;
				return false;
			}

			target_motion_sequence.reverse_motion_sequence();

			////////////////////////////////////
			//find path of fscs(if exists)
			////////////////////////////////////
			Fsc_indx source_fsc_indx (get_containig_fsc(perturbed_source));
			CGAL_postcondition (source_fsc_indx != Fsc_indx());
			Fsc_indx target_fsc_indx (get_containig_fsc(perturbed_target));
			CGAL_postcondition (target_fsc_indx != Fsc_indx());

			std::list<Fsc_indx> fsc_indx_path;
			if (source_fsc_indx == target_fsc_indx)
				fsc_indx_path.push_back(source_fsc_indx);
			else
				_graph.find_path( source_fsc_indx, target_fsc_indx, fsc_indx_path);

			if (fsc_indx_path.empty())
				return false;

			////////////////////////////////////
			//construct motion sequence
			////////////////////////////////////
			//(1) add source motion
			motion_sequence.add_motion_sequence(source_motion_sequence);

			//(2) construct real motion sequence from connectivity graph
			int             curr_fsc_indx = 0;
			Reference_point curr_ref_p    = perturbed_source;

			std::list<Fsc_indx>::iterator curr, next;

			next = curr = fsc_indx_path.begin();
			++next;
			while (next != fsc_indx_path.end())
			{
				Reference_point  next_ref_p = get_intersection(*curr, *next);

				Fsc* fsc_ptr = get_fsc(*curr);
				CGAL_postcondition(fsc_ptr->contains(curr_ref_p));
				CGAL_postcondition(fsc_ptr->contains(next_ref_p));
				CGAL_precondition ( (motion_sequence.get_sequence().empty()) || 
					(motion_sequence.get_sequence().back()->target() == curr_ref_p) );
				plan_path(fsc_ptr, curr_ref_p, next_ref_p, motion_sequence);
				CGAL_precondition ( (motion_sequence.get_sequence().empty()) || 
					(motion_sequence.get_sequence().back()->target() == next_ref_p) );

				curr++;
				next++;
				curr_ref_p = next_ref_p;
				delete fsc_ptr;
			}

			Fsc* fsc_ptr = get_fsc(*curr);
			plan_path(fsc_ptr, curr_ref_p, perturbed_target, motion_sequence);
			delete fsc_ptr;

			//(3) add source motion
			motion_sequence.add_motion_sequence(target_motion_sequence);
			return true;
		}
	private: //layer methods
		void generate_rotations(const unsigned int num_of_angles)
		{
			std::vector<Rotation> tmp_rotations;
			double step ((double)360/num_of_angles);
			//Rotation init_rotation(_rand.get_random_rotation());
			Rotation init_rotation(to_rotation<K::FT>(0));

			for (double r(0); r<360; r+=step)
				tmp_rotations.push_back(init_rotation* to_rotation<K::FT>(r,DEG,FINE));

			//maybe some rotations are the same as approximation is not good enough, remove duplicates 
			Less_than_rotation<K::FT> less_than;
			std::sort(tmp_rotations.begin(),tmp_rotations.end(), less_than);
			std::vector<Rotation>::iterator new_end = std::unique(tmp_rotations.begin(), tmp_rotations.end());

			//insert unique rotations (except maybe last) to the vector
			_rotations.clear();
			_rotations.insert(_rotations.begin(), tmp_rotations.begin(), new_end);

			//maybe last rotation is zero like the first rotation
			if (_rotations.back() == _rotations.front())
				_rotations.pop_back();

			return;
		}

		void add_layer(const Rotation& rotation)
		{
			//create layer
			Layer* layer_ptr = new Layer (Layer::Constraint(rotation));
			layer_ptr->decompose(_robot, _decomposed_workspace);
			int layer_id = _layers.add_manifold(layer_ptr);

			update_connectivity_graph_vertices(*layer_ptr, layer_id);
			return;
		}
	private:
		void generate_connectors()
		{
			compute_connectors_3_2();
			//generate_connectors_random();
		}

		void compute_connectors_3_2_help(Polygon polygon,Point vertex,Ray ray){
			Segment o_edge;
			Polygon::Edge_const_iterator o_edge_iter;
			CGAL::Object intersection_result;
			Point near,p,connector_p;
			near=Point(INT_MAX,INT_MAX);
			BOOST_FOREACH (Polygon o_polygon, _workspace){
				if (o_polygon!=polygon){
					o_edge_iter=o_polygon.edges_begin();
					while (o_edge_iter!=o_polygon.edges_end()){
						const Point *p;
						const Segment *s;
						o_edge=*o_edge_iter;
						intersection_result=CGAL::intersection(ray,o_edge);
						if (p=CGAL::object_cast<Point>(&intersection_result)){
							if (CGAL::compare_distance_to_point(vertex,*p,near)==CGAL::SMALLER){
								if (*p!=vertex)
									near=*p;
							}
						}
						else{
							if (s=CGAL::object_cast<Segment>(&intersection_result)){
								if (CGAL::compare_distance_to_point(vertex,(*s).vertex(0),(*s).vertex(1))==CGAL::SMALLER)
									p=&((*s).vertex(0));
								else
									p=&((*s).vertex(1));
								if (CGAL::compare_distance_to_point(vertex,*p,near)==CGAL::SMALLER){
									if (*p!=vertex)
										near=*p;
								}
							}
						}
						o_edge_iter++;
					}
				}
			}
			std::cout<<"("<<near.x()<<","<<near.y()<<")"<<endl;
			std::cout<<"("<<vertex.x()<<","<<vertex.y()<<")"<<endl;
			std::cout<<std::endl;
			connector_p=vertex+((near-vertex)/2);
			try{generate_connector_at_point(connector_p);}catch(...){};
		}

		void compute_connectors_3_2(){
			//for each vertex:
			//create ray starting the vertex, outside oriented, in the edge direction for both edges
			//for each edge of obstacle find intersection with the ray, if exist
			//find the closest intersection point
			//find the point in the middle between this point and the vertex
			//create connector

			Polygon::Vertex_iterator ver_iter;
			Point connector_p,vertex0,vertex1,near;
			Polygon::Edge_const_iterator edge_iter,o_edge_iter;
			Segment edge1,edge2,o_edge;
			CGAL::Object intersection_result;
			Line bisec;
			Ray ray;
			BOOST_FOREACH (Polygon polygon, _workspace){
				edge_iter=polygon.edges_begin();
				while (edge_iter!=polygon.edges_end()){

					edge1=*edge_iter;
					edge_iter++;

					vertex0=edge1.vertex(0);
					vertex1=edge1.vertex(1);

					try{
						ray=Ray(vertex0,edge1.opposite().supporting_line());
					} catch(...){};
					compute_connectors_3_2_help(polygon,vertex0,ray);
					try{
						ray=Ray(vertex1,edge1.supporting_line());
					} catch(...){};
					compute_connectors_3_2_help(polygon,vertex1,ray);
					/*BOOST_FOREACH (Polygon o_polygon, _workspace){
					if (o_polygon!=polygon){
					o_edge_iter=o_polygon.edges_begin();
					while (o_edge_iter!=o_polygon.edges_end()){
					const Point *p;
					const Segment *s;
					o_edge=*o_edge_iter;
					intersection_result=CGAL::intersection(ray,o_edge);
					if (p=CGAL::object_cast<Point>(&intersection_result)){
					if (CGAL::compare_distance_to_point(vertex0,*p,near)==CGAL::SMALLER)
					near=*p;
					}
					else{
					if (s=CGAL::object_cast<Segment>(&intersection_result)){
					if (CGAL::compare_distance_to_point(vertex0,(*s).vertex(0),(*s).vertex(1))==CGAL::SMALLER)
					p=&((*s).vertex(0));
					else
					p=&((*s).vertex(1));
					if (CGAL::compare_distance_to_point(vertex0,*p,near)==CGAL::SMALLER)
					near=*p;
					}
					}
					o_edge_iter++;
					}
					}
					}
					connector_p=vertex0+((near-vertex0)/2);
					try{generate_connector_at_point(connector_p);}catch(...){};*/

					/*try{
					ray=Ray(vertex1,edge1.supporting_line());
					}catch(...){};
					BOOST_FOREACH (Polygon o_polygon, _workspace){
					if (o_polygon!=polygon){
					o_edge_iter=o_polygon.edges_begin();
					while (o_edge_iter!=o_polygon.edges_end()){
					const Point *p;
					const Segment *s;
					o_edge=*o_edge_iter;
					intersection_result=CGAL::intersection(ray,o_edge);
					if (p=CGAL::object_cast<Point>(&intersection_result)){
					if (CGAL::compare_distance_to_point(vertex1,*p,near)==CGAL::SMALLER)
					near=*p;
					}
					else{
					if (s=CGAL::object_cast<Segment>(&intersection_result)){
					if (CGAL::compare_distance_to_point(vertex1,(*s).vertex(0),(*s).vertex(1))==CGAL::SMALLER)
					p=&((*s).vertex(0));
					else
					p=&((*s).vertex(1));
					if (CGAL::compare_distance_to_point(vertex1,*p,near)==CGAL::SMALLER)
					near=*p;
					}
					}
					o_edge_iter++;
					}
					}
					}
					connector_p=vertex1+((near-vertex1)/2);
					try{generate_connector_at_point(connector_p);} catch(...){};*/
				}
			}
		}

		void compute_connectors_3(){
			//for each vertex:
			//find bisector
			//create ray starting the vertex, outside oriented, in the bisector direction
			//for each edge of obstacle find intersection with the ray, if exist
			//find the closest intersection point
			//find the point in the middle between this point and the vertex
			//create connector

			Polygon::Vertex_iterator ver_iter;
			Point connector_p,vertex,near;
			Polygon::Edge_const_iterator edge_iter,o_edge_iter;
			Segment edge1,edge2,o_edge;
			CGAL::Object intersection_result;
			Line bisec;
			Ray ray;
			BOOST_FOREACH (Polygon polygon, _workspace){
				edge_iter=polygon.edges_begin();
				while (edge_iter!=polygon.edges_end()){


					// const Line *line1,*line2;
					//const CGAL::Line_2<K>
					//  Line line1,line2;

					edge1=*edge_iter;
					edge_iter++;
					edge2=*edge_iter;
					vertex=edge1.vertex(1);
					/*
					line1=(edge1.supporting_line());
					line2=(edge1.supporting_line().opposite());
					//CGAL::bisector(*line1,*line2);
					//Construct_bisector_2 bis(edge1.supporting_line(),edge2.supporting_line().opposite());
					K::FT a,b,c;
					CGAL::bisector_of_linesC2(line1.a(),line1.b(),line1.c(),line2.a(),line2.b(),line2.c(),a,b,c);
					bisec=Line(a,b,c);
					bisec=CGAL::bisector((const Line &)(edge1.supporting_line()),(const Line &)(edge2.supporting_line().opposite()));
					//bisec=CGAL::bisector(line1,line2);
					// bisec=K::ConstructBisector_2(*line1,*line2);
					*/


					ray=Ray(vertex,bisec);
					BOOST_FOREACH (Polygon o_polygon, _workspace){
						if (o_polygon!=polygon){
							o_edge_iter=o_polygon.edges_begin();
							while (o_edge_iter!=o_polygon.edges_end()){
								const Point *p;
								const Segment *s;
								o_edge=*o_edge_iter;
								intersection_result=CGAL::intersection(ray,o_edge);
								if (p=CGAL::object_cast<Point>(&intersection_result)){
									if (CGAL::compare_distance_to_point(vertex,*p,near)==CGAL::SMALLER)
										near=*p;
								}
								else{
									if (s=CGAL::object_cast<Segment>(&intersection_result)){
										if (CGAL::compare_distance_to_point(vertex,(*s).vertex(0),(*s).vertex(1))==CGAL::SMALLER)
											p=&((*s).vertex(0));
										else
											p=&((*s).vertex(1));
										if (CGAL::compare_distance_to_point(vertex,*p,near)==CGAL::SMALLER)
											near=*p;
									}
								}
								o_edge_iter++;
							}
						}
					}
					connector_p=vertex+((near-vertex)/2);
					generate_connector_at_point(connector_p);
				}
			}
		}
		/*
		void compute_connectors_2_1(){
		Polygon::Edge_const_iterator iter;
		Polygon::Segment_2 edge1,edge2;
		Polygon::Vertex_iterator intersec;
		CGAL::Line_2<K> bisec;
		Point vertex,p,near;
		K::FT len;

		BOOST_FOREACH(Polygon polygon1, _workspace){
		iter=polygon1.edges_begin();
		while (iter!=polygon1.edges_end()){
		edge1=*iter;
		iter++;
		edge2=(*iter).opposite();
		vertex=edge1.vertex(1);
		bisec=CGAL::bisector(edge1.supporting_line(),edge2.supporting_line().opposite());
		//TO DO: initialize near
		BOOST_FOREACH(Polygon polygon2, _workspace){
		if (polygon2!=polygon1){
		p=get_closest_point(vertex,polygon2);
		if (CGAL::compare_distance_to_point(vertex,p,near)==CGAL::SMALLER)
		near=p;
		}
		}
		len=((vertex-near).squared_length())/4; //squared radious of free circle
		//p=(((bisec.to_vector())*len/(bisec.to_vector().squared_length()))).end();
		p=CGAL::intersection(edge1.supporting_line(),CGAL::Circle_2<K>(vertex,len),intersec);
		BOOST_FOREACH(Point p,intersec){
		try{
		generate_connector_at_point(p);
		}
		catch(...){
		std::cout<<"i"<<std::endl;
		};
		}
		p=CGAL::intersection(edge2.supporting_line(),CGAL::Circle_2<K>(vertex,len));
		try{
		generate_connector_at_point(p);
		}
		catch(...){
		std::cout<<"i"<<std::endl;
		};
		p=CGAL::intersection(bisec,CGAL::Circle_2<K>(vertex,len));
		try{
		generate_connector_at_point(p);
		}
		catch(...){
		std::cout<<"i"<<std::endl;
		};
		}
		}
		}*/
		/*
		void compute_connectors_2(){
		int i=0;
		Polygon::Edge_const_iterator iter1,iter2;
		Polygon::Segment_2 edge1,edge2;
		CGAL::Line_2 bisec;
		//Point vertex1,vertex2,p;
		Point p,vertex,near;

		BOOST_FOREACH(Polygon polygon1, _workspace){
		BOOST_FOREACH(Polygon polygon2, _workspace){
		iter1=polygon1.edges_begin();
		iter2=polygon2.edges_begin();
		while (iter1!=polygon1.edges_end()){
		edge1=*iter1;
		iter1++;
		edge2=*iter1;
		vertex=edge1.vertex(1);//intersection of edge1 and edge2
		bisec=CGAL::bisector(edge1.supporting_line(),edge2.supporting_line().opposite());
		while (iter2!=polygon2.edges_end()){
		p=CGAL::intersection(bisec,*iter2);
		if (CGAL::Compare(CGAL::squared_distance(vertex,p),len)==CGAL::SMALLER)
		len=squared_distance(vertex,p);

		//p=vertex1+((vertex2-vertex1)/2);
		//std::cout <<i<<std::endl;
		//std::cout <<"x: "<<vertex1.x()<<";   y: "<<vertex1.y()<<std::endl;
		//std::cout <<"x: "<<vertex2.x()<<";   y: "<<vertex2.y()<<std::endl;
		//std::cout <<"x: "<<p.x()<<";   y: "<<p.y()<<std::endl;
		//std::cout <<std::endl;
		try{
		generate_connector_at_point(p);
		}
		catch(...){
		std::cout<<i<<std::endl;
		};
		i++;
		iter2++;
		}
		iter2=polygon2.edges_begin();
		iter1++;
		}
		}
		}
		}*/

		void compute_connectors(){
			int i=0;
			Polygon::Vertex_iterator iter1,iter2;
			Point vertex1,vertex2,p;

			BOOST_FOREACH(Polygon polygon1, _workspace){
				BOOST_FOREACH(Polygon polygon2, _workspace){
					iter1=polygon1.vertices_begin();
					iter2=polygon2.vertices_begin();
					while (iter1!=polygon1.vertices_end()){
						vertex1=*iter1;
						while (iter2!=polygon2.vertices_end()){
							vertex2=*iter2;
							p=vertex1+((vertex2-vertex1)/2);
							//std::cout <<i<<std::endl;
							//std::cout <<"x: "<<vertex1.x()<<";   y: "<<vertex1.y()<<std::endl;
							//std::cout <<"x: "<<vertex2.x()<<";   y: "<<vertex2.y()<<std::endl;
							//std::cout <<"x: "<<p.x()<<";   y: "<<p.y()<<std::endl;
							//std::cout <<std::endl;
							try{
								generate_connector_at_point(p);
							}
							catch(...){
								std::cout<<i<<std::endl;
							};
							i++;
							iter2++;
						}
						iter2=polygon2.vertices_begin();
						iter1++;
					}
				}
			}
		}
		void generate_connectors_random()
		{
			for (int i(0); i < configuration.get_max_num_of_intra_connections(); ++i)
				generate_connector_random();
		}
		void generate_connector_random()
		{
			generate_connector_at_point(_rand.get_random_point_in_bbox(_workspace_bbox));
		}
		void generate_connector_at_point(Point p)
		{
			////////////////////////////////////////////////////////////////
			//get free point in the configuration space on one of the layers
			////////////////////////////////////////////////////////////////
			Rotation& r  = _rotations[_rand.get_random_int(0, _rotations.size()-1)]; //random layer
			int layer_id = _layers.get_containing_manifold_id(r);
			Layer* layer_ptr = _layers.get_manifold(layer_id);

			if (layer_ptr->is_free(p)==false)
				return;
			//while (layer_ptr->is_free(p) == false)
			//	p = _rand.get_random_point_in_bbox(_workspace_bbox);

			int cell_id  = layer_ptr->get_containing_cell(p); 
			if (cell_id == NO_ID)
				return;// this should NOT be NO_ID but there is a patch in the configuration space of the angle primitive due to a bug in polygon_est_2



			C_space_line* c_space_line_ptr; 
			C_space_line::Constraint constraint;
			////////////////////////////////////////////////////////////////
			//choose roi
			////////////////////////////////////////////////////////////////
			double cell_size_ratio = get_size_percentage(layer_ptr->get_fsc(cell_id).cell() );
			CGAL_precondition (cell_size_ratio >=0 && cell_size_ratio <=1);

			if (configuration.get_use_region_of_interest() &&
				cell_size_ratio < 1 )
			{
				double small_rotation(configuration.get_rotation_range()/2);
				CGAL_precondition(small_rotation < 180);
				double half_range_double = small_rotation + cell_size_ratio * (180 - small_rotation);
				Rotation half_range = to_rotation<K::FT>(half_range_double, DEG);
				Rotation_range range(r*(-half_range), r*half_range);      

				constraint = C_space_line::Constraint(p, range);
			}
			else
			{
				constraint = C_space_line::Constraint(p);
			}

			////////////////////////////////////////////////////////////////
			//attempt to filter
			////////////////////////////////////////////////////////////////
			if (filter_out(constraint))
				return;

			////////////////////////////////////////////////////////////////
			//create connector
			////////////////////////////////////////////////////////////////
			c_space_line_ptr = new C_space_line (constraint, _ak);
			c_space_line_ptr->decompose(_robot, _decomposed_workspace);
			int c_space_line_id = _lines.add_manifold(c_space_line_ptr);

			////////////////////////////////////////////////////////////////
			//update connectivity graph
			////////////////////////////////////////////////////////////////
			update_connectivity_graph(c_space_line_id);
			return;

		}
	private: //filtering methods
		bool filter_out(typename C_space_line::Constraint& constraint)
		{
			if (configuration.get_use_filtering() == false)
				return false;

			std::vector<int> intersecting_layer_ids;
			_layers.get_interseceting_manifolds(constraint.region_of_interest(), std::back_inserter(intersecting_layer_ids));

			int cc_id = NO_ID;
			BOOST_FOREACH(int layer_id, intersecting_layer_ids)
			{
				Layer* layer_ptr = _layers.get_manifold(layer_id);
				const Rotation& rotation = layer_ptr->constraint().restriction();

				//find point in layer;
				int fsc_id = layer_ptr->get_containing_cell(constraint.restriction());
				if (fsc_id == NO_ID)
					continue;
				Fsc_indx curr_fsc_index(FIXED_ANGLE, layer_id, fsc_id);

				int curr_cc_id = _graph.get_cc_id(curr_fsc_index);
				if (cc_id == NO_ID)
					cc_id = curr_cc_id;
				else if (cc_id != curr_cc_id)
					return false;
			}

			return true;
		}

	private: //Connectiivity_graph methods
		void update_connectivity_graph_vertices(Layer& layer, int layer_id)
		{
			for (int fsc_id(0); fsc_id<layer.num_of_fscs(); ++fsc_id)
			{
				Fsc_indx layer_fsc_indx(FIXED_ANGLE, layer_id, fsc_id);
				_graph.add_vertex(layer_fsc_indx);
			}
			return;
		}
		void update_connectivity_graph(int c_space_line_id)
		{
			CGAL_precondition (c_space_line_id != NO_ID);

			C_space_line* line_ptr = _lines.get_manifold(c_space_line_id);
			const Rotation_range& range = line_ptr->constraint().region_of_interest();
			std::vector<int> intersecting_layer_ids;

			_layers.get_interseceting_manifolds( range, std::back_inserter(intersecting_layer_ids));
			BOOST_FOREACH(int layer_id, intersecting_layer_ids)
			{
				Layer* layer_ptr = _layers.get_manifold(layer_id);

				Int_pair edges_ids;
				intersect<K> (*layer_ptr, *line_ptr, edges_ids);

				CGAL_postcondition (  ((edges_ids.first == NO_ID) && (edges_ids.second == NO_ID)) ||
					((edges_ids.first != NO_ID) && (edges_ids.second != NO_ID)) );
				if (edges_ids.first == NO_ID)
					continue;

				Fsc_indx layer_fsc_indx(FIXED_ANGLE, layer_id, edges_ids.first);
				Fsc_indx line_fsc_indx (FIXED_POINT, c_space_line_id, edges_ids.second);

				_graph.add_vertex(layer_fsc_indx);
				_graph.add_vertex(line_fsc_indx);
				_graph.add_edge(layer_fsc_indx, line_fsc_indx);
			}
			return;
		}  
	private: //query related methods
		Reference_point connect_to_graph( const Reference_point& ref_p,
			Motion_sequence& motion_sequence)
		{
			//(1) find closest base layer
			Point location(ref_p.get_location());
			Rotation rotation(ref_p.get_rotation());

			int closest_layer_id = _layers.get_closest_layer_id(rotation);

			//(2) get the layer
			Layer* layer_ptr = _layers.get_manifold(closest_layer_id);    
			Rotation closest_rotation(layer_ptr->constraint().restriction());
			if (rotation == closest_rotation)
			{
				//no motion to do
				return ref_p;
			}

			if (layer_ptr->is_free(location) == false)
				return Reference_point();

			//(3) construct point predicate toward that layer
			C_space_line::Constraint constraint(location);
			C_space_line* line_ptr = new C_space_line (constraint, _ak);
			line_ptr->decompose(_robot, _decomposed_workspace);

			//(4) find fsc containing source point and target point
			int source_fsc_id = line_ptr->free_space_location_hint(ref_p);
			int target_fsc_id = line_ptr->free_space_location_hint(Reference_point(location, closest_rotation));
			if (source_fsc_id != target_fsc_id)
				return Reference_point();

			//(5) plan path on line
			std::vector<K::FT> tau_path;
			plan_path_in_interval<K>( FixedPoint::get_parametrization_theta<K>(ref_p.get_rotation()),
				FixedPoint::get_parametrization_theta<K>(closest_rotation),
				line_ptr->get_fsc(source_fsc_id).cell(),
				std::back_inserter(tau_path));
			CGAL_postcondition(tau_path.size() == 3); //source, target, mid


			//(6) convert point path to Motion_step_rotational
			CGAL::Orientation orientation = get_orientation(tau_path[0], tau_path[1], tau_path[2]);  
			Motion_step_rotational* motion_step_ptr = new Motion_step_rotational( location, 
				ref_p.get_rotation(), closest_rotation, 
				orientation);
			motion_sequence.add_motion_step(motion_step_ptr);

			delete line_ptr;
			return Reference_point(location, closest_rotation);
		}
	private: //Fsc_indx related methods
		Fsc_indx get_containig_fsc(const Reference_point& ref_p)
		{
			int manifold_id = _layers.get_containing_manifold_id(ref_p.get_rotation());
			if (manifold_id == NO_ID)
				return Fsc_indx();

			int fsc_id = _layers.get_manifold(manifold_id)->get_fsc_id(ref_p);

			if (fsc_id == NO_ID)
				return Fsc_indx();

			return Fsc_indx(FIXED_ANGLE,      //_type
				manifold_id,      //_manifold_id
				fsc_id            //_fsc_id;
				);
		}
		Fsc* get_fsc(const Fsc_indx& fsc_indx)
		{
			Fsc* fsc_ptr;
			if (fsc_indx._type == FIXED_ANGLE)
			{
				Layer* layer_ptr = _layers.get_manifold(fsc_indx._manifold_id);
				fsc_ptr = new Fsc(  layer_ptr->get_fsc(fsc_indx._fsc_id),
					layer_ptr->constraint());
			}
			else if (fsc_indx._type == FIXED_POINT)
			{
				C_space_line* c_space_line_ptr = _lines.get_manifold(fsc_indx._manifold_id);
				fsc_ptr = new Fsc(  c_space_line_ptr->get_fsc(fsc_indx._fsc_id),
					c_space_line_ptr->constraint());
			}
			return fsc_ptr;
		}

		Reference_point  get_intersection(const Fsc_indx& fsc_indx_1, const Fsc_indx& fsc_indx_2)
		{
			CGAL_precondition(  ((fsc_indx_1._type == FIXED_ANGLE) && (fsc_indx_2._type == FIXED_POINT)) || 
				((fsc_indx_1._type == FIXED_POINT) && (fsc_indx_2._type == FIXED_ANGLE)) );

			Rotation r;
			Point p;
			if (fsc_indx_1._type == FIXED_ANGLE)
			{
				r = _layers.get_manifold(fsc_indx_1._manifold_id)->constraint().restriction();
				p = _lines .get_manifold(fsc_indx_2._manifold_id)->constraint().restriction();
			}
			else //(fsc_indx_2._type == FIXED_ANGLE)
			{
				r = _layers.get_manifold(fsc_indx_2._manifold_id)->constraint().restriction();
				p = _lines .get_manifold(fsc_indx_1._manifold_id)->constraint().restriction();
			}

			return Reference_point(p,r);
		}
	private:
		void compute_workspace_bbox()
		{
			_workspace_bbox = _workspace[0].bbox();
			for (unsigned int i(1); i < _workspace.size(); ++i)
				_workspace_bbox = _workspace_bbox + _workspace[i].bbox();
			return;
		}
		void decompose_workspace_into_convex_polygons()
		{
			BOOST_FOREACH(Polygon polygon, _workspace)
				decompose_into_convex_polygons(polygon, std::back_inserter(_decomposed_workspace) );
		}
		template <typename OutputIterator>
		void decompose_into_convex_polygons(const Polygon& polygon, OutputIterator& oi)
		{
			CGAL::Small_side_angle_bisector_decomposition_2<K> decomp;
			decomp( polygon, oi);
			return;
		}

	};  //Mms_path_planner_example

} //mms
#endif //MMS_EXAMPLE_H
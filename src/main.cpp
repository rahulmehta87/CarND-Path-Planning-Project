#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "spline.h"

using namespace std;

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}
int ClosestWaypoint(double x, double y, const vector<double> &maps_x, const vector<double> &maps_y)
{

	double closestLen = 100000; //large number
	int closestWaypoint = 0;

	for(int i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x,y,map_x,map_y);
		if(dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}

	}

	return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{

	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2((map_y-y),(map_x-x));

	double angle = fabs(theta-heading);
  angle = min(2*pi() - angle, angle);

  if(angle > pi()/4)
  {
    closestWaypoint++;
  if (closestWaypoint == maps_x.size())
  {
    closestWaypoint = 0;
  }
  }

  return closestWaypoint;
}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

	int prev_wp;
	prev_wp = next_wp-1;
	if(next_wp == 0)
	{
		prev_wp  = maps_x.size()-1;
	}

	double n_x = maps_x[next_wp]-maps_x[prev_wp];
	double n_y = maps_y[next_wp]-maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
	double proj_x = proj_norm*n_x;
	double proj_y = proj_norm*n_y;

	double frenet_d = distance(x_x,x_y,proj_x,proj_y);

	//see if d value is positive or negative by comparing it to a center point

	double center_x = 1000-maps_x[prev_wp];
	double center_y = 2000-maps_y[prev_wp];
	double centerToPos = distance(center_x,center_y,x_x,x_y);
	double centerToRef = distance(center_x,center_y,proj_x,proj_y);

	if(centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for(int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
	}

	frenet_s += distance(0,0,proj_x,proj_y);

	return {frenet_s,frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int prev_wp = -1;

	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
	{
		prev_wp++;
	}

	int wp2 = (prev_wp+1)%maps_x.size();

	double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
	// the x,y,s along the segment
	double seg_s = (s-maps_s[prev_wp]);

	double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
	double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

	double perp_heading = heading-pi()/2;

	double x = seg_x + d*cos(perp_heading);
	double y = seg_y + d*sin(perp_heading);

	return {x,y};

}

int main() {
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;

  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
  	istringstream iss(line);
  	double x;
  	double y;
  	float s;
  	float d_x;
  	float d_y;
  	iss >> x;
  	iss >> y;
  	iss >> s;
  	iss >> d_x;
  	iss >> d_y;
  	map_waypoints_x.push_back(x);
  	map_waypoints_y.push_back(y);
  	map_waypoints_s.push_back(s);
  	map_waypoints_dx.push_back(d_x);
  	map_waypoints_dy.push_back(d_y);
  }

  // Start in lane 1
  int lane = 1;

  // Define a threshold of 6.0 seconds to prevent quick successive lane changes
  // A second lane change should not take place before 6 seconds after another lane change
  int lane_change_threshold = ceil(6.0 / 0.02);
  int lane_change_counter = 0;

  // Define a reference velocity to target
  double ref_v = 0.0; // mph

  // Define the anchor way-point spacing
  double wp_spacing = 30.0; // meters

  h.onMessage([&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy,
               &lane,&ref_v,&wp_spacing,&lane_change_counter,&lane_change_threshold]
               (uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);
        
        string event = j[0].get<string>();
        
        if (event == "telemetry") {
          // j[1] is the data JSON object
          
        	// Main car's localization Data
          	double car_x = j[1]["x"];
          	double car_y = j[1]["y"];
          	double car_s = j[1]["s"];
          	double car_d = j[1]["d"];
          	double car_yaw = j[1]["yaw"];
          	double car_speed = j[1]["speed"];

          	// Previous path data given to the Planner
          	auto previous_path_x = j[1]["previous_path_x"];
          	auto previous_path_y = j[1]["previous_path_y"];
          	// Previous path's end s and d values 
          	double end_path_s = j[1]["end_path_s"];
          	double end_path_d = j[1]["end_path_d"];

          	// Sensor Fusion Data, a list of all other cars on the same side of the road.
          	auto sensor_fusion = j[1]["sensor_fusion"];

          	json msgJson;

          	vector<double> next_x_vals;
          	vector<double> next_y_vals;

          	int prev_size = previous_path_x.size();

          	if (prev_size > 0) {
          	  car_s = end_path_s;
          	}

            bool too_close = false;

            // Define flags which determine whether it is safe to change lanes
            bool shift_left_safe = lane > 0 && lane_change_counter >= lane_change_threshold;
            bool shift_right_safe = lane < 2 && lane_change_counter >= lane_change_threshold;

            // Define the safe distance to be maintained from the vehicle to the front; assumes a minimum gap of 2 seconds or 10 m
            double front_safe_distance = max(10.0, 2.0 * car_speed * 0.44704);

            // Initialize the cost of being in each lane and the room that is available in each lane
            double lane_cost[3] = {0.0, 0.0, 0.0};
            double lane_room[3] = {1000.0, 1000.0, 1000.0};

            // Initialize the speed of the next vehicle in each lane to 50 mph
            double lane_next_vehicle_speeds[3] = {50.0 * 0.44704, 50.0 * 0.44704, 50.0 * 0.44704};

            // Update lane shift safety and lase costs with information from each vehicle around
            for (int i = 0; i < sensor_fusion.size(); i++) {
              float d = sensor_fusion[i][6];
              int check_lane = (int) floor(d / 4.0);
              double vx = sensor_fusion[i][3];
              double vy = sensor_fusion[i][4];
              double check_speed = sqrt(vx * vx + vy * vy);

              // Define the safe distance to be maintained from a vehicle behind if shifting into its lane
              double rear_safe_distance = max(10.0, 1.0 * check_speed);

              double check_car_s = sensor_fusion[i][5];
              check_car_s += prev_size * 0.02 * check_speed;

              // Current position of the vehicle
              double start_s = j[1]["s"];

              // Current position of the vehicle being checked
              double start_check_s = sensor_fusion[i][5];

              // Estimated position of vehicle after lane change; assuming it takes 6 seconds
              double end_s = car_s + (6.0 - prev_size * 0.02) * ref_v * 0.44704;

              // Estimated position of the vehicle being checked after lane change
              double end_check_s = check_car_s + (6.0 - prev_size * 0.02) * check_speed;

              // Update the available room and the speed of the next vehicle in each lane
              if (check_car_s > car_s && check_car_s - car_s < lane_room[check_lane]) {
                lane_room[check_lane] = check_car_s - car_s;
                lane_next_vehicle_speeds[check_lane] = check_speed;
              }

              // If the vehicle being checked is in the same lane, decide if we are too close
              if (d < (2 + 4 * lane + 3) && d > (2 + 4 * lane - 3)) {
                if (check_car_s > car_s && check_car_s - car_s < front_safe_distance) {
                  too_close = true;
                }
              }

              // If the vehicle being checked is in the lane to the right, decide if it is okay to shift to the right lane
              if (d < (6 + 4 * lane + 3) && d > (6 + 4 * lane - 3)) {
                if ((end_check_s > end_s - rear_safe_distance && end_check_s < end_s + front_safe_distance)
                 || (start_check_s > start_s - rear_safe_distance && start_check_s < start_s + front_safe_distance)
                 || (end_check_s >= end_s && start_check_s <= start_s)
                 || (end_check_s <= end_s && start_check_s >= start_s)) {
                  shift_right_safe = false;
                }
              }

              // If the vehicle being checked is in the lane to the left, decide if it is okay to shift to the left lane
             if (d < (-2 + 4 * lane + 3) && d > (-2 + 4 * lane - 3)) {
                if ((end_check_s > end_s - rear_safe_distance && end_check_s < end_s + front_safe_distance)
                 || (start_check_s > start_s - rear_safe_distance && start_check_s < start_s + front_safe_distance)
                 || (end_check_s >= end_s && start_check_s <= start_s)
                 || (end_check_s <= end_s && start_check_s >= start_s)) {
                  shift_left_safe = false;
                }
              }
            }

            // Update the lane costs based on how much forward progress can be made in each lane in the next 8 seconds
            for (int i = 0; i < 3; i++) {
              lane_cost[i] = 200 - min(ref_v * 0.44704 * 8.0, lane_room[i] + lane_next_vehicle_speeds[i] * 8.0 - front_safe_distance);
            }

            // Shift lanes if needed
            if (too_close) {
              // If we are too close, reduce the speed
              ref_v -= 0.223 * min(1.0, front_safe_distance / lane_room[lane]);

              if (shift_left_safe && shift_right_safe) {
                // If it is safe to shift left or right, choose the lane with the lower cost
                if (lane_cost[lane + 1] < lane_cost[lane - 1]) {
                  lane += 1;
                  lane_change_counter = 0;
                } else {
                  lane -= 1;
                  lane_change_counter = 0;
                }
              } else if (shift_left_safe) {
                lane -= 1;
                lane_change_counter = 0;
              } else if (shift_right_safe) {
                lane += 1;
                lane_change_counter = 0;
              }
            } else {
              // If we are not too close, accelerate up to 49.5 mph
              if (ref_v < 49.5) {
                ref_v += 0.223;
              }

              // Switch lanes if it is safe and it helps reduce the lane cost
              if (shift_left_safe && shift_right_safe) {
                if (lane_cost[lane + 1] < lane_cost[lane - 1] && lane_cost[lane + 1] < lane_cost[lane]) {
                  lane += 1;
                  lane_change_counter = 0;
                } else if (lane_cost[lane - 1] < lane_cost[lane]) {
                  lane -= 1;
                  lane_change_counter = 0;
                }
              } else if (shift_left_safe && ((lane_cost[lane - 1] < lane_cost[lane]) || (lane == 2 && lane_cost[0] < lane_cost[2]))) {
                lane -= 1;
                lane_change_counter = 0;
              } else if (shift_right_safe && ((lane_cost[lane + 1] < lane_cost[lane]) || (lane == 0 && lane_cost[2] < lane_cost[0]))) {
                lane += 1;
                lane_change_counter = 0;
              }
            }

            lane_change_counter++;

          	// Create a list of widely spaced way-points, evenly spaced at wp_spacing.
          	// Later, we will interpolate these way-points with a spline and fill it in with more points that control
          	// the speed
          	vector<double> pts_x;
          	vector<double> pts_y;

          	// Reference x, y, yaw states
          	// This is either where the car is, or at the end point of the previous path.
          	double ref_x = car_x;
          	double ref_y = car_y;
          	double ref_yaw = deg2rad(car_yaw);

          	// If previous size is almost empty, use the car as starting reference.
          	if (prev_size < 2) {
          	  // Use two points that make the path tangent to the car's current path.
          	  double prev_car_x = car_x - cos(ref_yaw);
          	  double prev_car_y = car_y - sin(ref_yaw);

          	  pts_x.push_back(prev_car_x);
          	  pts_y.push_back(prev_car_y);

          	  pts_x.push_back(car_x);
          	  pts_y.push_back(car_y);
          	} else {
          	  // Use the previous end point as the starting reference.
          	  ref_x = previous_path_x[prev_size - 1];
          	  ref_y = previous_path_y[prev_size - 1];

              // Use two points that make the path tangent to the path at the reference
          	  double ref_x_prev = previous_path_x[prev_size - 2];
          	  double ref_y_prev = previous_path_y[prev_size - 2];
              ref_yaw = atan2(ref_y - ref_y_prev, ref_x - ref_x_prev);

          	  pts_x.push_back(ref_x_prev);
          	  pts_y.push_back(ref_y_prev);

          	  pts_x.push_back(ref_x);
          	  pts_y.push_back(ref_y);
          	}

          	vector<double> next_wp0 = getXY(car_s + wp_spacing,
          	                               (2 + 4 * lane),
          	                               map_waypoints_s,
          	                               map_waypoints_x,
          	                               map_waypoints_y);
          	vector<double> next_wp1 = getXY(car_s + 2 * wp_spacing,
          	                               (2 + 4 * lane),
          	                               map_waypoints_s,
          	                               map_waypoints_x,
          	                               map_waypoints_y);
          	vector<double> next_wp2 = getXY(car_s + 3 * wp_spacing,
          	                               (2 + 4 * lane),
          	                               map_waypoints_s,
          	                               map_waypoints_x,
          	                               map_waypoints_y);

          	pts_x.push_back(next_wp0[0]);
          	pts_y.push_back(next_wp0[1]);

          	pts_x.push_back(next_wp1[0]);
          	pts_y.push_back(next_wp1[1]);

          	pts_x.push_back(next_wp2[0]);
          	pts_y.push_back(next_wp2[1]);

            // Transform the anchor way-points to the reference coordinate system
            for (int i = 0; i < pts_x.size(); i++) {
              double shift_x = pts_x[i] - ref_x;
              double shift_y = pts_y[i] - ref_y;

              pts_x[i] = shift_x * cos(0 - ref_yaw) - shift_y * sin(0 - ref_yaw);
              pts_y[i] = shift_x * sin(0 - ref_yaw) + shift_y * cos(0 - ref_yaw);
            }

            // Create a spline
            tk::spline s;

            // Set the anchor way-points to the spline
            s.set_points(pts_x, pts_y);

            // Start with all of the previous path points from last time
            for (int i = 0; i < prev_size; i++) {
              next_x_vals.push_back(previous_path_x[i]);
              next_y_vals.push_back(previous_path_y[i]);
            }

            // Calculate how to break up spline points so that we travel at our desired reference velocity

            double target_x = wp_spacing;
            double target_y = s(target_x);
            double target_dist = sqrt(target_x * target_x + target_y * target_y);

            double x_add_on = 0;

            // Fill up the rest of our path planner
            for (int i = 0; i < 50 - prev_size; i++) {
              double N = target_dist / (0.02 * (ref_v / 2.23694));
              double x_point = x_add_on + target_x / N;
              double y_point = s(x_point);

              x_add_on = x_point;

              double x_ref = x_point;
              double y_ref = y_point;

              // Transform back to map coordinates

              x_point = x_ref * cos(ref_yaw) - y_ref * sin(ref_yaw);
              y_point = x_ref * sin(ref_yaw) + y_ref * cos(ref_yaw);

              x_point += ref_x;
              y_point += ref_y;

              next_x_vals.push_back(x_point);
              next_y_vals.push_back(y_point);
            }

          	msgJson["next_x"] = next_x_vals;
          	msgJson["next_y"] = next_y_vals;

          	auto msg = "42[\"control\","+ msgJson.dump()+"]";

          	//this_thread::sleep_for(chrono::milliseconds(1000));
          	ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
          
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}

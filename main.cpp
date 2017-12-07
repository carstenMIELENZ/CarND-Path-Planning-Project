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
// added spline CAM
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





// ------------------------------------------------------------------------
// New variables & constants
// + lane
// + velocity
// ------------------------------------------------------------------------

int     lane     = 1;    // start lane 
double  velocity = 0;    // start velocity 

// constants
const double VELOCITY_EMERGENCY  = 0.33; // speed of front car in terms percentance of ego velocity 
const double VELOCITY_MAX        = 49.5;
const double VELOCITY_STEP       = 0.336;
const double VELOCITY_DEC        = 0.85;
const int    DISTANCE_NUM        = 50;
const double POINTS_PER_SEC      = 0.02;
const double REF_DISTANCE        = 30.0;
const double BACK_DISTANCE       = 0.33;
const int PREVIOUS_SIZE_LIMIT    = 2;


// -------------------------------------------------------------------------------------
// function check_lane
// 
// + check if a lane change can be preformed based on lane gap and velocity
// + of the cars in the new lane. 
// 
// + RETURN: min. velocity of "front" cars n new lane if possible, else
//           negative, if new has no "front" cars then max. velocity 
// 
// + CONSTANTS:  
//    VELOCITY_MAX    = speed limit 
//    VELOCITY_DEC    = velocity tolerance of lane change 
//    POINTS_PER_SEC  = steps for velocity 
//    REF_DISTANCE    = reference distance from ego car to next car  
//    BACK_DISTANCE   = factor of reference distance of backwards cars   
//  
// + GAP definition 
//                            B                 F
//  Lane current 	  <--------> ego car <--------->  
//  Lane to change     car (backward)                  car (front)
//  
//  B = REF_DISTANCE * BACK_DISTANCE
//  F = REF_DISTANCE 
//  
// + LANE speed definition:  speed of new lane > speed of current lane * VELOCITY_DEC  
//  
// -------------------------------------------------------------------------------------

double check_lane (const vector<vector<double>> &fusion, double ref_s, int lane_ref, double speed_ref, int lane_off_set,int prev_size) { 
  
  // analyze cars on to switch lane
  bool    lane_ok    = true;
  double  lane_speed = VELOCITY_MAX;
  int     lane       = lane_ref+lane_off_set;

  // analyze 
  for (int a=0;a<fusion.size();++a) {

     // get line info
     float d = fusion[a][6];

     // check if inside lane
     if (d<(2+4*lane+2) && (d>2+4*lane-2)) {
      
       lane_ok = false;

       // get data 
       double vx           = fusion[a][3];                    
       double vy           = fusion[a][4];                    
       double check_speed  = sqrt(vx*vx+vy*vy);
       double check_car_s  = fusion[a][5];                    
              check_car_s += ((double)prev_size*POINTS_PER_SEC*check_speed);
       
       // position check for lane change
               
       // cars in front of vehicle
       if (check_car_s>ref_s) { 
        
         // no collision + speed is higher of speed of car infront current lane 
         if ( ((check_car_s-ref_s)>REF_DISTANCE) && (check_speed>(VELOCITY_DEC*speed_ref)) ) {
        
           lane_ok = true;
           // set smallest velocity of lane
           if (check_speed<lane_speed) lane_speed = check_speed;
         }
      } 
      // backwards check
      else if (check_car_s<=(ref_s-REF_DISTANCE*BACK_DISTANCE)) {
        
             lane_ok = true;
           }
    } // end inside lane


    // final check not ok => lane not ok  
    if (lane_ok == false) break;                           

  } // end list fusion 


  // return status
  if (lane_ok) return lane_speed;
  else         return -1.0; 

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

  h.onMessage([&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
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

                // vehicle data to simulator
          	vector<double> next_x_vals;
          	vector<double> next_y_vals;

          	// TODO: define a path made up of (x,y) points that the car will visit sequentially every .02 seconds
		
		/* ************************************************************************************
 		 * Adedd code Dec-03-2017 
 		 * +  
 		 * + 
		 * ************************************************************************************ */	

             	 // get from simulator
                 int prev_size = previous_path_x.size(); 

                 // -------------------------------------------------------------------------------------
                 // sensor fusion analysis
                 // -------------------------------------------------------------------------------------
                 
                 // check previous trajectory
                 if (prev_size>0) car_s = end_path_s;

                 // trigger for reducing 
		 bool too_close       = false;
		 bool emergency_break = false;

                 // analyze sensor data for each lane
                 for (int i=0;i<sensor_fusion.size();++i) {
                   
                   // check if other vehicles on my lane 
                   float d = sensor_fusion[i][6];

                   // check for vehicles inside my lane 
                   if (d<(2+4*lane+2) && (d>2+4*lane-2)) {
                    
		     double vx          = sensor_fusion[i][3];                    
		     double vy          = sensor_fusion[i][4];                    
		     double check_speed = sqrt(vx*vx+vy*vy);                    
		     double check_car_s = sensor_fusion[i][5];                    

		     // calculate position
     		     check_car_s += ((double)prev_size*POINTS_PER_SEC*check_speed);

                     // too close ?
		     if ((check_car_s>car_s) && ((check_car_s-car_s)<REF_DISTANCE)) {

                       // reduce velocity 
                       too_close = true;

                       // ------------------------------------------------------------------------------- 
                       // check for lane change or emergency break 
                       // ------------------------------------------------------------------------------- 

                       if (car_speed*VELOCITY_EMERGENCY > check_speed) emergency_break = true; // emergency break
                       else {

                       // lane change analysis
                       
                       double speed_lane   = -1.0;
                       double speed_lane_l = -1.0;
                                            
                       // ------------------------------------------------------------------------------- 
                       // lane analysis by checking gap for lane change & speed at new lane
                       // ------------------------------------------------------------------------------- 

                       //  change from left to middle lane

                       if (lane == 0) {
                         speed_lane = check_lane(sensor_fusion, car_s, lane, check_speed, 1,prev_size); 
                         if (speed_lane > 0) {
                           lane = 1;  
                         } 

                         // change from right to middle lane       
                       } 
                       else if (lane == 2) {
                              speed_lane = check_lane(sensor_fusion, car_s, lane, check_speed, -1,prev_size); 
                              if (speed_lane > 0) {
                                lane = 1;  
                              } 
                            } 

                            // change from middle to left or right lane       
                                     
                            else {
                                   speed_lane_l = check_lane(sensor_fusion, car_s, lane, check_speed, -1,prev_size);
                                   speed_lane   = check_lane(sensor_fusion, car_s, lane, check_speed,  1,prev_size);
                                   if (speed_lane_l > speed_lane) {     
                                     lane = 0;
                                   }
                                   else if (speed_lane_l < speed_lane) {
                                          lane = 2;
                                        }
                                        else if (speed_lane > 0) {
                                               lane = 0;
                                             }  
                                 } 
 
                       } // end change lane or emergency break 
 
                     } // end too close

                   } // end in lane 

                 } // end fusion analysis


                 // -------------------------------------------------------------------------------------
                 //  velocity adaption 
                 // -------------------------------------------------------------------------------------
                
                 // reduce velocity 
       		 if (too_close) { 
                 
                   velocity -= VELOCITY_STEP;
                
                   // emergncy break 
                   if (emergency_break) { 
                     velocity -= 2*VELOCITY_STEP;
                   }
                 }
                 // increase velocity
 		 else if (velocity < VELOCITY_MAX) velocity += VELOCITY_STEP;                                
               

                 // -------------------------------------------------------------------------------------
                 // Trajectory 
                 // -------------------------------------------------------------------------------------
 
                 // wavepoint list (x,y) 
          	 vector<double> ptsx;
          	 vector<double> ptsy;

                 // start reference of vehicle               
                 double ref_x    = car_x;
                 double ref_y    = car_y;
                 double ref_yaw  = deg2rad(car_yaw);
 
                 //  check if list size 
                 if (prev_size < PREVIOUS_SIZE_LIMIT) {
                  
                   // take 2 points for path - kind of linarization (car_x + previous_x calculated from past linear)
                   double prev_car_x = car_x-cos(ref_yaw);
                   double prev_car_y = car_y-sin(ref_yaw);
                   
                   // build list
                   ptsx.push_back(prev_car_x);
                   ptsx.push_back(car_x);
                   ptsy.push_back(prev_car_y);
                   ptsy.push_back(car_y);
                 }
                 else {
     
                   // take reference from last entry of list
                   ref_x = previous_path_x[prev_size-1];
                   ref_y = previous_path_y[prev_size-1];
                   // one more from past
                   double ref_x_prev = previous_path_x[prev_size-2];
                   double ref_y_prev = previous_path_y[prev_size-2];

                   // calculate yaw
                   ref_yaw = atan2(ref_y-ref_y_prev,ref_x-ref_x_prev);
                  
                   // build list
                   ptsx.push_back(ref_x_prev);
                   ptsx.push_back(ref_x);
                   ptsy.push_back(ref_y_prev);
                   ptsy.push_back(ref_y);
                }

                //
		// Add Fenet line of 3 x REF_DISTANCE - lane selects middle of lane
                //
                vector<double> next_wp0 = getXY(car_s+REF_DISTANCE  ,(2+4*lane),map_waypoints_s,map_waypoints_x,map_waypoints_y);
                vector<double> next_wp1 = getXY(car_s+REF_DISTANCE*2,(2+4*lane),map_waypoints_s,map_waypoints_x,map_waypoints_y);
                vector<double> next_wp2 = getXY(car_s+REF_DISTANCE*3,(2+4*lane),map_waypoints_s,map_waypoints_x,map_waypoints_y);
		
                // add to vector (after this vector size is 5 = 2 + 3 
          	ptsx.push_back(next_wp0[0]);
          	ptsx.push_back(next_wp1[0]);
          	ptsx.push_back(next_wp2[0]);
          	ptsy.push_back(next_wp0[1]);
          	ptsy.push_back(next_wp1[1]);
          	ptsy.push_back(next_wp2[1]);

                //
                // Transform to vehicles coordinate system - see MPC 
                //
                for (int i=0;i<ptsx.size();++i) {
                   
                   // shift by 0 degrees
                   double shift_x = ptsx[i]-ref_x;
                   double shift_y = ptsy[i]-ref_y;
                  
                   ptsx[i] =(shift_x*cos(0-ref_yaw)-shift_y*sin(0-ref_yaw));
                   ptsy[i] =(shift_x*sin(0-ref_yaw)+shift_y*cos(0-ref_yaw));
                }

                //
                // Create a spline wavepoints ptsx/y
                //
                tk::spline s;
                // set points
                s.set_points(ptsx,ptsy);
                
                //
                // Create/reuse path based on previous - it only unses points which have not used in simulator
                // org. list size = 50, simulator used (simulated) 5 points => previous list size = 45
                //
                for (int i=0;i<previous_path_x.size(); ++i) {
                
          	   next_x_vals.push_back(previous_path_x[i]);
          	   next_y_vals.push_back(previous_path_y[i]);
                }

                //
                // Calculate how to breakup spline for trajectory - triangle linearization
                //
                double target_x    = REF_DISTANCE; // x horizon
                double target_y    = s(target_x);  
                double target_dist = sqrt((target_x*target_x)+(target_y*target_y));   
                double x_add_on    = 0;
                // calc. of the following loop moved for perf.
                double N           = (target_dist/(POINTS_PER_SEC*velocity/2.24)); // note 2.24 factor for m/s vs. mph - may put out of loop
                double target_step = (target_x)/N;                         

                // 
                // Create Trajectory 
                // 
                for (int i=0;i<DISTANCE_NUM-previous_path_x.size();++i) {
                  
                   double x_point = x_add_on+target_step;
                   double y_point = s(x_point);  
                   // new start      
                   x_add_on = x_point;

                   // Transfrom from vehicle coord. to simulator coord.
                   double x_ref = x_point;                   
                   double y_ref = y_point;                   
                    
                   x_point  = (x_ref*cos(ref_yaw)-y_ref*sin(ref_yaw)); 
                   y_point  = (x_ref*sin(ref_yaw)+y_ref*cos(ref_yaw)); 
                   
                   x_point += ref_x; 
                   y_point += ref_y;

                   // add to vector 
          	   next_x_vals.push_back(x_point);
          	   next_y_vals.push_back(y_point);
                }
     
                // these variables has to be set for the simulator
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

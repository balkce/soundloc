#ifndef MultiTargetTracker_HPP
#define MultiTargetTracker_HPP

#include "opencv2/highgui/highgui.hpp"
#include "Track.hpp"
#include "scans.hpp"
#include "matrix.h"
#include "munkres.h"
#include <fstream>

class MultiTargetTracker
{    
public:
     MultiTargetTracker(double init_sample_period = 0.1, 
			uint init_max_miss_detections = 20, 
			uint init_min_continuous_detections = 3,
			double init_distance_threshold = 0.3, bool init_arg_conf = true)
	  {
	       kRadians = M_PI / 180;
	       kDegrees = 180 / M_PI;
	       min_doas_assigned = 2;
	       init(init_sample_period, init_max_miss_detections, init_min_continuous_detections, init_distance_threshold);
	       
	       if(!init_arg_conf){
		    //initializing arguments from file
			
		    //getting directory of this header file
		    size_t found;
		    std::string head_location(__FILE__);
		    found=head_location.find_last_of("/\\");
			
		    //assuming conf file is located in same dir as this header
		    //and named mcmcda.conf
		    std::string conf_location = head_location.substr(0,found) + "/kalman.conf";
			
		    //reading file
		    std::cout << "KALMAN: Reading conf file... " << conf_location <<"\n";
		    std::ifstream conf_file(conf_location.c_str());
		    std::string line;
			
		    //getting max_miss_detections
		    if( std::getline( conf_file, line ) ) {
			 max_miss_detections = atof(line.c_str());
			 std::cout << "KALMAN: max_miss_detections :  " << max_miss_detections << "\n";
		    }else{
			 exit(1);
		    }
		    //getting sample_period
		    if( std::getline( conf_file, line ) ) {
			 sample_period = atof(line.c_str());
			 std::cout << "KALMAN: sample_period :  " << sample_period << "\n";
		    }else{
			 exit(1);
		    }
		    //getting min_doas_assigned
		    if( std::getline( conf_file, line ) ) {
			 min_doas_assigned = atof(line.c_str());
			 std::cout << "KALMAN: min_doas_assigned :  " << min_doas_assigned << "\n";
		    }else{
			 exit(1);
		    }
		    //getting kDistThres
		    if( std::getline( conf_file, line ) ) {
			 kDistThres = atof(line.c_str());
			 std::cout << "KALMAN: kDistThres :  " << kDistThres << "\n";
		    }else{
			 exit(1);
		    }
			
		    conf_file.close();
			
	       }else{
		    std::cout << "KALMAN: Using default values... \n";
		    std::cout << "KALMAN: sample_period :  " << sample_period << "\n";
		    std::cout << "KALMAN: max_miss_detections :  " << max_miss_detections << "\n";
		    std::cout << "KALMAN: min_doas_assigned :  " << min_doas_assigned << "\n";
	       }
	       
	  }
     
     void init(double init_sample_period, 
	       uint init_max_miss_detections = 20, 
	       uint init_min_continuous_detections = 3,
	       double init_distance_threshold = 0.3)
	       {
		    kDistThres = init_distance_threshold;
		    tracked_doa_number = 0;
		    time = 0;
		    sample_period = init_sample_period;
		    max_miss_detections = init_max_miss_detections;
	       }
     
     void update(Scan& new_scan)
	  {
	       time++;
	       std::cout << std::endl;

	       std::cout << "============== [ " << time << " ] ====================" << std::endl;
	       std::cout << " Observations :: ";
	       if (new_scan.size() > 1) {
		    for (Scan::iterator it = new_scan.begin(); it != new_scan.end(); ++it)
			 std::cout << it->value << "   ";

		    for (int i = 0;  i < new_scan.size() - 1; ++i){
			 for (int j = i + 1;  j < new_scan.size(); ++j){
			      if (std::abs(new_scan[i].value - new_scan[j].value) < 5.0){
				   new_scan[i].value = (new_scan[i].value + new_scan[j].value) / 2;
				   new_scan.erase(new_scan.begin() + j);
				   j--;
			      }
			 }
		    }
	       }
	       else if (!new_scan.empty())
		    std::cout << new_scan[0].value << "   ";
	       std::cout << std::endl;
	       std::cout << " No Dup :: ";
	       for (Scan::iterator it = new_scan.begin(); it != new_scan.end(); ++it)
		    std::cout << it->value << "   ";
	       std::cout << std::endl;
		       
	       std::vector<std::pair<float,float>> doas;
	       for (Scan::iterator it = new_scan.begin(); it != new_scan.end(); ++it){
		    std::pair<float,float> new_doa;
		    new_doa.first = std::cos(it->value * kRadians);
		    new_doa.second = std::sin(it->value * kRadians);
		    doas.push_back(new_doa);
	       }

	       std::cout << " Cartesian :: ";
	       for (std::vector<std::pair<float,float>>::iterator it = doas.begin(); it != doas.end(); ++it)
		    std::cout << "(" << it->first << ", " << it->second << ")   ";
	       std::cout << std::endl;
	       
	       if(!doas.empty()){
		    if(!tracks.empty()) { 	
			 Matrix<double> cost_matrix = computeCostMatrix(tracks, doas);
			 Munkres munkres;
			 munkres.solve(cost_matrix);
			 std::vector<int> associations = getAssignments(cost_matrix);
			 std::vector<int> assigned(doas.size());
			 for(int i = 0; i < associations.size(); i++){
			      std::cout << i << " not assigned " << DOADistance(tracks[i], doas[associations[i]]) << std::endl;
			      if(associations[i] >= 0 && DOADistance(tracks[i], doas[associations[i]]) < kDistThres) {// valid association
				   tracks[i].addObservation(doas[associations[i]], time);
				   assigned[associations[i]] = 1;
			      }
			 }
			 
			 // create a new track for each unassigned observation
			 for(int i = 0; i < doas.size(); i++){
			      if(assigned[i] == 0){
				   kTrack new_track(doas[i], time, tracked_doa_number, sample_period, max_miss_detections, min_doas_assigned);
				   tracks.push_back(new_track);
				   tracked_doa_number++;
			      }
			 }
		    } else {
			 // no existing tracks, create one for each observation
			 for (std::vector<std::pair<float,float>>::iterator it = doas.begin(); it != doas.end(); ++it){
			      kTrack new_track(*it, time, tracked_doa_number,sample_period, max_miss_detections, min_doas_assigned);
			      tracks.push_back(new_track);
			      tracked_doa_number++;
			 }
		    }
	       } 
			
	       for (std::vector<kTrack>::iterator it = tracks.begin(); it != tracks.end(); ++it){
		    if (it->getCurrentTime() != time)			 
			 it->addObservation(time);
	       }

	       // delete invalid tracks
	       for (std::vector<kTrack>::iterator it = tracks.begin(); it != tracks.end();){ 
		    if (!it->valid())
			 it = tracks.erase(it);
		    else 
			 ++it;

		    if (tracks.empty())
			 break;				  
            
	       }


	       std::cout << " kTracks :: ";
	       for (std::vector<kTrack>::iterator it = tracks.begin(); it != tracks.end(); ++it){
		    std::pair<float,float> corrected_pos = it->getCorrection(); 
		    float angle_deg = std::atan2(corrected_pos.second, corrected_pos.first) * kDegrees;
		    std::cout << "(" << it->getTrackNumber() << ")" << angle_deg << "   ";
	       }
	       std::cout << std::endl;

	  }


     Scan getCurrentStates(void)
	  {
	       Scan current_states = Scan();

               // checks if there are valid tracks
	       if(!tracks.empty()){
		    for (std::vector<kTrack>::iterator it = tracks.begin(); it != tracks.end(); it++)
			 if (it->in_range()) {
			      Observation current_value;
			      if (it->getLastTime() == time){
				   std::pair<float,float> last_pos = it->getLastPos();
				   current_value.value = std::atan2(last_pos.second, last_pos.first) * kDegrees;
				   current_states.push_back(current_value);

			      }
			      else{
				   std::pair<float,float> corrected_pos = it->getCorrection();
				   current_value.value = std::atan2(corrected_pos.second, corrected_pos.first) * kDegrees;
				   current_states.push_back(current_value);
			      }
			 }
	       }

	       return current_states;
	  }
	  
	  double DOADistance(kTrack track, std::pair<float,float> doa)
	  {
	       std::pair<float,float> corrected_pos = track.getCorrection();

	       return std::sqrt((corrected_pos.first - doa.first) * (corrected_pos.first - doa.first) + (corrected_pos.second - doa.second) * (corrected_pos.second - doa.second));
	  }

	  Matrix<double> computeCostMatrix(std::vector<kTrack> tracks, std::vector< std::pair<float,float> > doas)
	  {
	       Matrix<double> cost_matrix(tracks.size(), doas.size());

	       for (int row = 0 ; row < tracks.size() ; row++) 
		    for (int col = 0 ; col < doas.size() ; col++)
			 cost_matrix(row,col) = (double) DOADistance(tracks[row], doas[col]);

	       return cost_matrix;
	  }

	  std::vector<int> getAssignments(Matrix<double> cost_matrix)
	  {
	       std::vector<int> assignments(cost_matrix.rows(), -1);
		  
	       for (int row = 0 ; row < cost_matrix.rows() ; row++) 
		    for (int col = 0 ; col < cost_matrix.columns() ; col++)
			 if (cost_matrix(row,col) == 0)
			      assignments[row] = col;
		  
	       return assignments;
	  }

private:
     int time;
     int tracked_doa_number;
     std::vector<kTrack> tracks;

     double sample_period;
     int max_miss_detections;
     int min_doas_assigned;
     
     double kDistThres;
     double kRadians;
     double kDegrees;
};
#endif

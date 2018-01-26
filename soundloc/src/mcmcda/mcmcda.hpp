/**
 * @brief Classes and structures for Markov Chain Monte Carlo Data Association.
 * 
 * @author Gibran Fuentes Pineda <gibranfp@turing.iimas.unam.mx>
 * @date 2013
 */
#ifndef MCMCDA_HPP
#define MCMCDA_HPP
#include <string>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <algorithm>
#include "partition.hpp"
#include "scans.hpp"
#include "partition.hpp"
#include "posterior.hpp"
#include "proposal.hpp"
#include "utils.hpp"

class MCMCDA
{    
public:
     /**
      * @brief Function to determine whether two observations are neighbors
      *
      * @param Observation1 First observation
      * @param Observation2 Second observation
      */
     bool isNeighbor(Observation observation1, Observation observation2) 
	  {
	       double d_x = std::cos(observation2.value * (M_PI / 180)) - std::cos(observation1.value * (M_PI / 180)); 
	       double d_y = std::sin(observation2.value * (M_PI / 180)) - std::sin(observation1.value * (M_PI / 180)); 

	       double distance = std::sqrt(d_x * d_x + d_y * d_y);
	       if (distance < max_neighbor_distance)
		    return true;
	       else
		    return false;
	  }

     /**
      * @brief MCMCDA class constructor.
      *
      * @param init_scample_size Number of samples for Metropolis-Hastings
      * @param init_max_miss_detections Maximum time a target is not observed
      * @param init_time_window Window of time where MCMCDA is going to be applied
      * @param init_max_target_speed Maximum speed of the target in a sampling time
      */
     MCMCDA(uint init_sample_size, uint init_max_miss_detections, uint init_time_window = 100,
	    uint init_max_target_speed = 1, bool init_posterior_fixed = true)
	  {
	       sample_size = init_sample_size;
	       max_miss_detections = init_max_miss_detections;
	       time_window = init_time_window;
	       max_target_speed = init_max_target_speed;
	       posterior_fixed = init_posterior_fixed;
		
	       if(!posterior_fixed){
		    //initializing posterior arguments if not fixed
			
		    //getting directory of this header file
		    size_t found;
		    std::string head_location(__FILE__);
		    found=head_location.find_last_of("/\\");
			
		    //assuming conf file is located in same dir as this header
		    //and named mcmcda.conf
		    std::string conf_location = head_location.substr(0,found) + "/mcmcda.conf";
			
		    //reading file
		    std::cout << "MCMCDA: Reading conf file... " << conf_location <<"\n";
		    std::ifstream conf_file(conf_location.c_str());
		    std::string line;
			
		    //getting max_miss_detections
		    if( std::getline( conf_file, line ) ) {
			 max_miss_detections = atof(line.c_str());
			 std::cout << "MCMCDA: max_miss_detections :  " << max_miss_detections << "\n";
		    }else{
			 exit(1);
		    }
		    //getting posterior_size_prior_a
		    if( std::getline( conf_file, line ) ) {
			 posterior_size_prior_a = atof(line.c_str());
			 std::cout << "MCMCDA: posterior_size_prior_a :  " << posterior_size_prior_a << "\n";
		    }else{
			 exit(1);
		    }
		    //getting posterior_size_prior_b
		    if( std::getline( conf_file, line ) ) {
			 posterior_size_prior_b = atof(line.c_str());
			 std::cout << "MCMCDA: posterior_size_prior_b :  " << posterior_size_prior_b << "\n";
		    }else{
			 exit(1);
		    }
		    //getting posterior_track_length_prior_a
		    if( std::getline( conf_file, line ) ) {
			 posterior_track_length_prior_a = atof(line.c_str());
			 std::cout << "MCMCDA: posterior_track_length_prior_a :  " << posterior_track_length_prior_a << "\n";
		    }else{
			 exit(1);
		    }
		    //getting posterior_track_length_prior_b
		    if( std::getline( conf_file, line ) ) {
			 posterior_track_length_prior_b = atof(line.c_str());
			 std::cout << "MCMCDA: posterior_track_length_prior_b :  " << posterior_track_length_prior_b << "\n";
		    }else{
			 exit(1);
		    }
			
		    conf_file.close();
			
	       }else{
		    max_miss_detections = 10;
		    posterior_size_prior_a = 0.01;
		    posterior_size_prior_b = -0.6;
		    posterior_track_length_prior_a = 0.3;
		    posterior_track_length_prior_b = 0.6;
		    std::cout << "MCMCDA: Using default values...\n";
		    std::cout << "MCMCDA: max_miss_detections :  " << max_miss_detections << "\n";
		    std::cout << "MCMCDA: posterior_size_prior_a :  " << posterior_size_prior_a << "\n";
		    std::cout << "MCMCDA: posterior_size_prior_b :  " << posterior_size_prior_b << "\n";
		    std::cout << "MCMCDA: posterior_track_length_prior_a :  " << posterior_track_length_prior_a << "\n";
		    std::cout << "MCMCDA: posterior_track_length_prior_b :  " << posterior_track_length_prior_b << "\n";
	       }

	       max_neighbor_distance = max_miss_detections * max_target_speed;
	       current_time = 0;
	       partition.tracks = TrackSet();
	       partition.false_positives = Track();
	  }

     /*
      * @brief Performs markov chain monte carlo data association
      *
      * @param init Initial partition      
      *
      * @return maxpost Partition with maximum posterior
      */
     Partition multiScanMCMCDA(Partition init)
	  {
	       std::random_device seed2{};
	       std::default_random_engine generator(seed2());
	       std::uniform_real_distribution<double> uniform_distribution(0.0,1.0);
	       
	       Proposal proposal(nforest, current_time, max_miss_detections);
	       Posterior posterior(scans);
	       if (!posterior_fixed){
		    posterior.set_prior_arguments(posterior_size_prior_a,posterior_size_prior_b,posterior_track_length_prior_a,posterior_track_length_prior_b);
	       }
	       
	       // initialize the markov chain
	       Partition map_partition = init; // maximum a posteriori partition

	       double posterior_map_partition = posterior(map_partition);

	       for (uint t = 0; t < sample_size; ++t){ // Metropolis-Hastings sampler

		    // propose partition based on proposed move and previous partition
		    Partition sample = proposal(map_partition, scans);	

		    if (!sample.tracks.empty()){
			 // compute posteriors
			 double posterior_sample = posterior(sample);

			 // compute transition probabilities
			 double fwd_posterior = posterior_sample * proposal.backward_prob();
			 double back_posterior = posterior_map_partition * proposal.forward_prob();

			 // compute acceptance probability
			 double acceptance = std::min(1.0, fwd_posterior / back_posterior);

			 // accept or reject sample
			 double u = uniform_distribution(generator);
			 if (u < acceptance){
			      // find new maximum posterior
			      if ((posterior_sample / posterior_map_partition) > 1){
				   map_partition = sample;
				   posterior_map_partition = posterior_sample;
				   scans = proposal.getScans();
			      }
			 }
		    }
	       }
		
	       return map_partition;
	  }

     /*
      * @brief Performs markov chain monte carlo data association from a set of observations
      */
     void sample(void)
	  {
	       partition = multiScanMCMCDA(partition);
	  }

     ScanSet getScans()
	  {
	       return scans;
	  }

     /*
      * @brief Adds new set of observations
      *
      * @param new_scan Set of observations
      */
     void addObservations(Scan new_scan)
	  {     
	       if (new_scan.size() > 1){
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

	       // adding new observations
	       scans.push_back(new_scan);

	       // Add new_scan to neighbor forest
	       NeighborScan neighbor_scan;
	       for (Scan::iterator ob = new_scan.begin();  ob != new_scan.end(); ++ob){
		    // Add new observations to false positives
		    Item item = {current_time - 1, (uint) (ob - new_scan.begin())};
		    partition.false_positives.push_back(item);
		    
		    // At first there won't be any neighbors as it is the last observations
		    NeighborObservation neighbor_obs = NeighborObservation();
		    for (uint i = 0; i < max_miss_detections; ++i)
			 neighbor_obs.push_back(NeighborTime());
		    
		    neighbor_scan.push_back(neighbor_obs);
	       }
	       nforest.push_back(neighbor_scan);
	       
	       // new time
	       current_time++;

	       // updating neighborhood tree of observations
	       if (current_time > 1){
		    uint time;
		    uint scan_number;

		    // There are less observations than neighbor window
		    if (current_time <= max_miss_detections)
			 time = 0;
		    else
			 time = current_time - max_miss_detections - 1;
		    
		    // Update past observation's neighbors with new scan
		    while (time < current_time - 1){
			 for (uint number = 0; number < scans[time].size(); ++number){
			      for (Scan::iterator ob = new_scan.begin();  ob != new_scan.end(); ++ob){
				   if (isNeighbor(*ob, scans[time][number])){
					Item item = {current_time - 1, (uint) (ob - new_scan.begin())};
					nforest[time][number][current_time - time - 2].push_back(item);
				   }
			      }
			 }
			 ++time;
		    }
	       }
	  }
     
     /*
      * @brief Adds a new set of observations to the model and finds the new partition
      *
      * @param new_scan New set of observations
      */
     void update(Scan& new_scan)
	  {     
	       addObservations(new_scan);	       

	       if (scans.size() > 1)
		    sample();
	  }

     /*
      * @brief Returns the current found tracks
      */
     uint getCurrentTime()
	  {
	       return current_time;
	  }

     /*
      * @brief Returns the current found tracks
      */
     TrackSet getTracks()
	  {
	       return partition.tracks;
	  }

     /*
      * @brief Returns the forest of observations
      */
     NForest getForest()
	  {
	       return nforest;
	  }

     /*
      * @brief Returns valid track's last observations
      */
     Scan getCurrentStates(TrackSet& tracks)
	  {
	       Scan current_states;
	       for (TrackSet::iterator it = tracks.begin(); it != tracks.end(); ++it){
		    if (it->size() > 2 && current_time - it->back().time < max_miss_detections){
			 if (current_time - 1 == it->back().time){
			      current_states.push_back(scans[it->back().time][it->back().number]);
			 }
			 else{
			      cv::KalmanFilter kalman_filter(4, 2, 0);
			      // transition matrix: angle and radial speed                1    sample_time
			      //                                                          0        1
			      kalman_filter.transitionMatrix = *(cv::Mat_<float>(4, 4) << 1,0,0.1,0,   0,1,0,0.1,   0,0,1,0,   0,0,0,1);
					
			      // measurement vector: angle
			      cv::Mat_<float> measurement(2,1); 
			      measurement.setTo(cv::Scalar(0));
					
			      // initilize all Kalman matrices
			      setIdentity(kalman_filter.measurementMatrix);
			      setIdentity(kalman_filter.processNoiseCov, cv::Scalar::all(1e-4));
			      setIdentity(kalman_filter.measurementNoiseCov, cv::Scalar::all(1e-1));
			      setIdentity(kalman_filter.errorCovPost, cv::Scalar::all(1e-1));

			      // First observation
			      kalman_filter.predict();
			      measurement(0) = std::cos(scans[it->at(0).time][it->at(0).number].value); 
			      measurement(1) = std::sin(scans[it->at(0).time][it->at(0).number].value); 
			      kalman_filter.correct(measurement);
					
			      // second observation
			      kalman_filter.predict();
			      measurement(0) = std::cos(scans[it->at(1).time][it->at(1).number].value); 
			      measurement(1) = std::sin(scans[it->at(1).time][it->at(1).number].value); 
			      kalman_filter.correct(measurement);

			      uint time = it->at(1).time;
			      cv::Mat predicted, corrected;
			      for (Track::iterator ob = it->begin() + 2; ob != it->end(); ++ob)
			      {
				   // Kalman's prediction phase
				   predicted = kalman_filter.predict();

				   if ((*ob).time == time + 1){// incorporating new measurement
					measurement(0) = std::cos(scans[ob->time][ob->number].value); 
					measurement(1) = std::sin(scans[ob->time][ob->number].value); 
				   }
				   else{ // no measurement, taking predicted
					measurement(0) = predicted.at<float>(0);
					measurement(1) = predicted.at<float>(1);
				   }
						
				   // Kalman correction phase
				   corrected = kalman_filter.correct(measurement);
						
				   ++time;
			      }
					
			      while(time < current_time - 1){
				   predicted = kalman_filter.predict();
				   measurement(0) = predicted.at<float>(0);
				   measurement(1) = predicted.at<float>(1);
				   corrected = kalman_filter.correct(measurement);
				   ++time;
			      }
				
			      predicted = kalman_filter.predict();
			      Observation current_value;
			      current_value.value = std::atan2(predicted.at<float>(1), predicted.at<float>(0)) * (180.0 / M_PI);
			      current_states.push_back(current_value);
			 }
		    }
	       }
		
	       return current_states;
	  }

	

private:
     uint sample_size;
     uint current_time;
     uint time_window;
     uint max_miss_detections;
     double max_target_speed;
     double max_neighbor_distance;
     ScanSet scans; 
     NForest nforest; 
     Partition partition;
     float posterior_size_prior_a;
     float posterior_size_prior_b;
     float posterior_track_length_prior_a;
     float posterior_track_length_prior_b;
     bool posterior_fixed;
};
#endif

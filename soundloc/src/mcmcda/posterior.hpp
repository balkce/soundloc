/**
 * @brief Classes and structures for Markov Chain Monte Carlo Data Association.
 * 
 * @author Gibran Fuentes Pineda <gibranfp@turing.iimas.unam.mx>
 * @date 2013
 */
#ifndef POSTERIOR_HPP
#define POSTERIOR_HPP

#include <iomanip>

#include <vector>
#include <opencv2/video/tracking.hpp>
#include "partition.hpp"
#include "scans.hpp"

const float kNorm = 180;

// functor for posterior distribution
class Posterior
{
public:
	/*
	 * @brief constructor
	 * 
	 * @param scans All measurements
	 *
	 */ 
	Posterior(ScanSet& init_scans): scans(init_scans) {
		posterior_size_prior_a = 0.01;
		posterior_size_prior_b = -0.6;
		posterior_track_length_prior_a = 0.3;
		posterior_track_length_prior_b = 0.6;
	}
	
	void set_prior_arguments(float sa, float sb, float la, float lb){
		posterior_size_prior_a = sa;
		posterior_size_prior_b = sb;
		posterior_track_length_prior_a = la;
		posterior_track_length_prior_b = lb;
	}
     
     
	/*
	 * @brief prior probability of a partition
	 * 
	 * @param partition Observation partition
	 *
	 * @return prior probability
	 */ 
    double prior(Partition& partition)
    {
        // negative exponential distribution to penalize large number of tracks
		double size_prior = posterior_size_prior_a * exp(posterior_size_prior_b * partition.tracks.size());
        
        // positive exponential distribution to favor longer tracks
        double track_length_prior = 1.0;
        for (TrackSet::iterator it = partition.tracks.begin(); it != partition.tracks.end(); ++it)
            track_length_prior *= (posterior_track_length_prior_a * exp(posterior_track_length_prior_b * (double)(it->size()))); 
        
        return size_prior * track_length_prior;
    }

	/*
	 * @brief Computes the likelihood of a partition. It uses a distribution that favors tracks with 
	 *        observations that are more consistent with a given dynamic model
	 * 
	 * @param partition Observation partition
	 *
	 * @return likelihood probability
	 */ 
	double likelihood(Partition& partition)
    {
        if (partition.tracks.empty()){
            return 0.0;
        }
       
	int width = 15;

	
	// std::cout.unsetf ( std::ios::floatfield );                // floatfield not set
	// std::cout.precision(2);
	// std::cout << "TNum" << std::setw(5) 
	// 	  << "Time" << std::setw(5) 
	// 	  << "ONum" << std::setw(5) 
	// 	  << " Observation" << std::setw(10) 
	// 	  << "        Predicted" << std::setw(width) 
	// 	  << "Estimated" << std::setw(width) 
	// 	  << "Error" << std::endl;
       double partition_likelihood = 1.0;
        for (TrackSet::iterator it = partition.tracks.begin(); it != partition.tracks.end(); ++it){
            double track_likelihood = 1.0;
            if (it->size() >= 2) {
                // construct kalman filter      
                cv::KalmanFilter kalman_filter(4, 2, 0);

                // transition matrix: angle and radial speed                1        0    sample_time        0
                //                                                          0        1       0            sample_time
		//                                                          0        0       1               0
		//                                                          0        0       0               1
                kalman_filter.transitionMatrix = *(cv::Mat_<float>(4, 4) << 1,0,0.1,0,   0,1,0,0.1,   0,0,1,0,   0,0,0,1);

                // measurement vector: angle
                cv::Mat_<float> measurement(2,1); 
		measurement.setTo(cv::Scalar(0));
		
                // initilize all Kalman matrices
                setIdentity(kalman_filter.measurementMatrix);
                setIdentity(kalman_filter.processNoiseCov, cv::Scalar::all(1e-4));
                setIdentity(kalman_filter.measurementNoiseCov, cv::Scalar::all(1e-4));
                setIdentity(kalman_filter.errorCovPost, cv::Scalar::all(1e-1));

		// First observation
		kalman_filter.predict();
		measurement(0) = std::cos(scans[it->at(0).time][it->at(0).number].value * (M_PI/180)); 
		measurement(1) = std::sin(scans[it->at(0).time][it->at(0).number].value * (M_PI/180)); 
		kalman_filter.correct(measurement);

		// second observation
		kalman_filter.predict();
		measurement(0) = std::cos(scans[it->at(1).time][it->at(1).number].value * (M_PI/180)); 
		measurement(1) = std::sin(scans[it->at(1).time][it->at(1).number].value * (M_PI/180)); 
		kalman_filter.correct(measurement);

                // compute likelihood for each track based on motion model (using Kalman filtering)
                Track::iterator ob = it->begin() + 2;
                uint time = it->at(1).time;
                while (ob != it->end()) {
		     // std::cout << "  " << (int)(it - partition.tracks.begin()) << std::setw(5);

		     // Kalman's prediction phase
                    cv::Mat predicted = kalman_filter.predict();
                    
                    if ((*ob).time == time + 1){// incorporating new measurement
			 measurement(0) = std::cos(scans[ob->time][ob->number].value * (M_PI/180)); 
			 measurement(1) = std::sin(scans[ob->time][ob->number].value * (M_PI/180)); 
			 
 			 // std::cout << ob->time << std::setw(5) 
			 // 	   << ob->number << std::setw(5); 

			 // std::cout << "[" << scans[ob->time][ob->number].value << "](" << measurement(0) << "," << measurement(1) << ")" << std::setw(10);
			 ++ob;
                    }
		    else{ // no measurement, taking predicted
			 measurement(0) = predicted.at<float>(0);
			 measurement(1) = predicted.at<float>(1);
			 // std::cout << time + 1 << std::setw(5) 
			 // 	   << "PR" << std::setw(5); 
			 // std::cout << "      none" << std::setw(10);
		    }

		    // Kalman correction phase
		    cv::Mat corrected = kalman_filter.correct(measurement);

		    // Innovation covariance
		    cv::Mat covariance = kalman_filter.temp3;

		    // error between measurement and predicted
		    cv::Mat error = kalman_filter.temp5;
		    
		    // Observation motion likelihood
		    double scale = 1 / std::sqrt(cv::determinant(2 * M_PI * covariance));
		    cv::Mat exponent = -0.5 * (error.t() * covariance.inv() * error);
		    double obs_likelihood  = scale * exp(exponent.at<float>(0));

		    // std::cout << "[" << std::atan2(predicted.at<float>(1), predicted.at<float>(0)) * (180.0 / M_PI)
			// 		  << "](" 
			// 		  << predicted.at<float>(0) << ","  
		    // 	      << predicted.at<float>(1) << ")"
			 //	    	      << predicted.at<float>(2) << ","
			 //<< predicted.at<float>(3) << ")" 
			      // << std::setw(width) 
		    	  //     << "[" << std::atan2(corrected.at<float>(1), corrected.at<float>(0)) * (180.0 / M_PI) << "](" << corrected.at<float>(0) << ","  
		    	  //     << corrected.at<float>(1) << ")"
		    	  //     // << corrected.at<float>(2) << ","
		    	  //     // << corrected.at<float>(3) << ")" 
			      // << std::setw(width) 
		    	  //     << "(" << error.at<float>(0) << "," << error.at<float>(1) << ")" 
			      // << std::endl << std::endl;
		    

		    // Track motion likelihood
		    track_likelihood *= obs_likelihood;		    
		    
                    ++time;
                }

		// Partition likelihood
		partition_likelihood *= track_likelihood;
		
		//std::cout << " Likelihood = " << partition_likelihood << " track = " << track_likelihood <<  std::endl;
            }
        }

        return partition_likelihood;
    }

	/*
	 * @brief Computes the likelihood of a partition. It uses a distribution that favors tracks with 
	 *        observations that are more consistent with a given dynamic model
	 * 
`	 * @param partition Observation partition
	 *
	 * @return likelihood probability
	 */ 
	double operator() (Partition& partition)
	{
	     double prior_val = 1.0;//prior(partition);
	     double likelihood_val = likelihood(partition);

	     return prior_val * likelihood_val;
	}
     
private:
	ScanSet scans;
	float posterior_size_prior_a;
	float posterior_size_prior_b;
	float posterior_track_length_prior_a;
	float posterior_track_length_prior_b;
};
#endif

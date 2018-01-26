#include <opencv2/video/tracking.hpp>
#include <iostream>
#include <stdio.h>
#include <time.h>
#include <string>

class kTrack
{
public:
     // Track constructor
     kTrack(std::pair<float,float> init_pos, 
	    uint init_time, 
	    uint init_track_number = 0,
	    double init_sample_period = 0.1, 
	    uint init_max_miss_detections = 10,
	    uint init_min_doas_assigned = 2)
	  { 
	       init(init_pos, init_time, init_track_number, init_sample_period, init_max_miss_detections, init_min_doas_assigned); 
	  }
     
     void init(std::pair<float,float> init_pos, 
	       uint init_time, 
	       uint init_track_number = 0,
	       double init_sample_period = 0.1, 
	       uint init_max_miss_detections = 10,
	       uint init_min_doas_assigned = 2)
	  {
	       first_pos = last_pos = init_pos;
	       sample_period = init_sample_period;
	       max_miss_detections = init_max_miss_detections;
	       track_number = init_track_number;
	       first_time = current_time = last_detection_time = init_time;
	       min_doas_assigned = init_min_doas_assigned;
	       doas_assigned = 0;

	       // initialize dynamic model parameters
	       kalman_filter = cv::KalmanFilter(4, 2, 0);
	       processNoise = cv::Mat(4, 1, CV_32F);
	       measurement = cv::Mat_<float>(2, 1);
	       measurement.setTo(cv::Scalar(0));
	       
	       kalman_filter.transitionMatrix = *(cv::Mat_<float>(4, 4) <<  
						  1, 0, sample_period,   0,            //x pixel position
						  0, 1,   0,   sample_period,	        //y pixel position
						  0, 0,   1,       0,		//vx velocity in x
						  0, 0,   0,       1);	       //vy velocity in y
	       
	       setIdentity(kalman_filter.measurementMatrix);
	       setIdentity(kalman_filter.processNoiseCov, cv::Scalar::all(1e-3)); 
	       setIdentity(kalman_filter.measurementNoiseCov, cv::Scalar::all(1e-3)); 
	       setIdentity(kalman_filter.errorCovPost, cv::Scalar::all(.1));

	       kalman_filter.statePost.at<float>(0) = init_pos.first;
	       kalman_filter.statePost.at<float>(1) = init_pos.second;
	       kalman_filter.statePost.at<float>(2) = 0.1;
	       kalman_filter.statePost.at<float>(3) = 0.1;

	       // prediction phase
	       prediction = kalman_filter.predict();
	       predicted_pos.first = prediction.at<float>(0);
	       predicted_pos.second = prediction.at<float>(1);

	       // add new measurements to the dynamic model
	       measurement(0) = init_pos.first;
	       measurement(1) = init_pos.second;

	       // correction phase
	       correction = kalman_filter.correct(measurement);
	       corrected_pos.first = correction.at<float>(0);
	       corrected_pos.second = correction.at<float>(1);
	  }
          
	void addObservation(std::pair<float,float> pos, int time)
	  {
	       // add new measurements
	       measurement(0) = pos.first;
	       measurement(1) = pos.second;

	       last_pos.first = pos.first;
	       last_pos.second = pos.second;


	       // Kalman filter's prediction phase
	       prediction = kalman_filter.predict();
	       predicted_pos.first = prediction.at<float>(0);
	       predicted_pos.second = prediction.at<float>(1);
	       

	       // Kalman filter's correction phase
	       correction = kalman_filter.correct(measurement);
	       corrected_pos.first =  correction.at<float>(0);
	       corrected_pos.second = correction.at<float>(1);

	       current_time = time;
	       last_detection_time = time;
	       doas_assigned++;
	  }

     void addObservation(int time)
	  {	       
	       // Kalman filter's prediction phase
	       prediction = kalman_filter.predict();
	       predicted_pos.first = prediction.at<float>(0);
	       predicted_pos.second = prediction.at<float>(1);
	       
	       // take predicted as measurement
	       measurement(0) = predicted_pos.first;
	       measurement(1) = predicted_pos.second;

	       // Kalman filter's correction phase
	       correction = kalman_filter.correct(measurement);
	       corrected_pos.first =  correction.at<float>(0);
	       corrected_pos.second = correction.at<float>(1);
	       
	       current_time = time;
	  }

     bool in_range(void)
	  {
	       return current_time - last_detection_time < max_miss_detections && doas_assigned >= min_doas_assigned;
	  }

     bool valid(void)
	  {
	       return (current_time - last_detection_time) <= max_miss_detections;
          }

     std::pair<float,float> getPrediction()
	  {
	       return predicted_pos;
	  }
     
     std::pair<float,float> getCorrection()
	  {
	       return corrected_pos;
	  }

     std::pair<float,float> getLastPos()
	  {
	       return last_pos;
	  }

     std::pair<float,float> getFirstPos()
	  {
	       return first_pos;
	  }

	int getCurrentTime()
	  {
		  return current_time;
	  }

	int getLastTime()
	  {
		  return last_detection_time;
	  }

     int getTrackNumber()
	  {
		  return track_number;
	  }

private:
     cv::KalmanFilter kalman_filter;
     
     cv::Mat prediction;
     cv::Mat processNoise;
     cv::Mat correction;
     
     cv::Mat_<float> state;
     cv::Mat_<float> measurement;
     
     std::pair<float,float> predicted_pos;
     std::pair<float,float> corrected_pos;

     std::pair<float,float> first_pos;
     std::pair<float,float> last_pos;

     uint first_time;
     uint last_detection_time;
     uint current_time;

     double sample_period;
     uint max_miss_detections;
     uint doas_assigned;
     
     uint track_number;
     uint min_doas_assigned;
     
};

/*
	Source coded by Caleb Rascon, 2010
	IIMAS, UNAM
	México
	
	Example that uses libmultisoundloc.a
*/

#include "multisoundloc.h"

#include <iostream>
#include <signal.h>
#include <string.h>
#include <stdio.h>

//ROS libs
#include "ros/ros.h"
#include "std_srvs/Empty.h"
#include <jsoncpp/json/json.h>
#include <json_msgs/JsonSrv.h>

//Glib, for threading (GThread)
#include <glib.h>

void finalize(int signum){
	std::cout << "\nSoundLoc: Caught signal " << signum << std::endl;
	std::cout << "SoundLoc: Exiting gracefully.\n";
	
	//Telling all threads that we are shutting down.
	
	soundloc_finalize();
	
	//gplot_stream.cmd("set terminal postscript");
	//gplot_stream.cmd("set output 'streams.eps'");
	//gplot_stream.cmd("replot");
	
	std::cout << "SoundLoc: Bye.\n";
	
	ROS_INFO("Closing ROS node.");
	ros::shutdown();
}

bool ros_reset_soundloc(json_msgs::JsonSrv::Request  &req, json_msgs::JsonSrv::Response &res){
	ROS_INFO("resetting soundloc");
	soundloc_clear();
	return true;
}

bool ros_get_sound_directions(json_msgs::JsonSrv::Request  &req, json_msgs::JsonSrv::Response &res){
	ROS_INFO("providing current sound directions");fflush(stdout);
	Json::FastWriter writer;

	g_mutex_lock (&mutex_sources);
	if(sources.size() > 0){
		Json::Value response;
		Json::Value json_sound_direction(Json::arrayValue);	
		bool first_angle = TRUE;
		for(int i = 0; i <sources.size(); i++){
			if(sources[i].real_source){
				Json::Value object_json;
				object_json["doa"] = sources[i].last_doa;
				json_sound_direction.append(object_json);
			}
		}
	    	response["doas"] = json_sound_direction;
	    	res.json = writer.write(response);
	}else{
		Json::Value root;
		Json::Reader reader;
		reader.parse("{}", root);
		assert(root != Json::nullValue);
		res.json = writer.write(root);
		assert(res.json == "{}\n");
	}
	g_mutex_unlock (&mutex_sources);
	
	return true;
}

bool ros_get_sound_directions_pre(json_msgs::JsonSrv::Request  &req, json_msgs::JsonSrv::Response &res){
	ROS_INFO("providing current sound directions");
	Json::FastWriter writer;

	g_mutex_lock (&mutex_sources);
	if(sources.size() > 0){
		Json::Value response;
		Json::Value json_sound_direction(Json::arrayValue);	
		bool first_angle = TRUE;
		for(int i = 0; i <sources.size(); i++){
			for(int j = 0; j <sources[i].doas.size(); j++){
				Json::Value object_json;
				object_json["doa"] = sources[i].doas[j].value;
				json_sound_direction.append(object_json);
			}
		}
	    	response["doas"] = json_sound_direction;
	    	res.json = writer.write(response);
	}else{
		Json::Value root;
		Json::Reader reader;
		reader.parse("{}", root);
		assert(root != Json::nullValue);
		res.json = writer.write(root);
		assert(res.json == "{}\n");
	}
	g_mutex_unlock (&mutex_sources);
	
	return true;
}

bool ros_get_vad_results(json_msgs::JsonSrv::Request  &req, json_msgs::JsonSrv::Response &res){
	ROS_INFO("providing current sound directions");
	Json::FastWriter writer;
	Json::Value json_vad;
	
	Json::Value json_sound_direction(Json::arrayValue);
	if(vad_res.doas.size() > 0){
		for(int i = 0; i <vad_res.doas.size(); i++){
			Json::Value object_json;
			object_json["doa"] =  vad_res.doas[i];
			json_sound_direction.append(object_json);
		}
		json_vad["doas"] = json_sound_direction;
		json_vad["file"] = vad_res.audio_file_path;
	}else{
		json_vad["doas"] = json_sound_direction;
		json_vad["file"] = "";
	}
	
	Json::Value response;
	response["doas"] = json_vad;
	res.json = writer.write(response);
	
	return true;
}

int main( int argc, char *argv[] )
{
	//if ( argc != 8 && argc != 9 && argc != 10) usage();
	
	//Assigning a CTRL+C handler
	signal(SIGINT, finalize);
	
	//ROS initialization
	std::cout << "SoundLoc: Starting up ROS connection...\n";
	
	ros::init(argc, argv, "soundloc");
	ros::NodeHandle ros_nh("~");
	ros::ServiceServer reset_soundloc = ros_nh.advertiseService("reset_soundloc", ros_reset_soundloc);
	ros::ServiceServer get_sound_directions = ros_nh.advertiseService("get_sound_directions", ros_get_sound_directions);
	ros::ServiceServer get_sound_directions_pre = ros_nh.advertiseService("get_sound_directions_pre", ros_get_sound_directions_pre);
	ros::ServiceServer get_vad_results = ros_nh.advertiseService("get_vad_results", ros_get_vad_results);
	
	// Obtaining parameters from ROS parameter server
	double distance_between_mics;
	if (ros_nh.getParam("distance_between_mics",distance_between_mics)){
		ROS_INFO("Distance between mics: %f",distance_between_mics);
	}else{
		distance_between_mics = 0.21;
		ROS_WARN("Distance between mics argument not found in ROS param server, using default value(%f).",distance_between_mics);
	}
	
	int number_of_samples_of_noise;
	if (ros_nh.getParam("number_of_samples_of_noise",number_of_samples_of_noise)){
		ROS_INFO("Number of samples of noise: %d",number_of_samples_of_noise);
	}else{
		number_of_samples_of_noise = 10;
		ROS_WARN("Number of samples of noise argument not found in ROS param server, using default value(%d).",number_of_samples_of_noise);
	}
	
	double noise_threshold;
	if (ros_nh.getParam("noise_threshold",noise_threshold)){
		ROS_INFO("Noise threshold: %f",noise_threshold);
	}else{
		noise_threshold = 0.0015;
		ROS_WARN("Noise threshold argument not found in ROS param server, using default value(%f).",noise_threshold);
	}
	
	double noise_peak_change;
	if (ros_nh.getParam("noise_peak_change",noise_peak_change)){
		ROS_INFO("Noise peak change: %f",noise_peak_change);
	}else{
		noise_peak_change = 0.002;
		ROS_WARN("Noise peak change argument not found in ROS param server, using default value(%f).",noise_peak_change);
	}
	
	double coherence_threshold;
	if (ros_nh.getParam("coherence_threshold",coherence_threshold)){
		ROS_INFO("Coherence threshold: %f",coherence_threshold);
	}else{
		coherence_threshold = 30;
		ROS_WARN("Coherence threshold argument not found in ROS param server, using default value(%f).",coherence_threshold);
	}
	
	bool reverse_angle;
	if (ros_nh.getParam("reverse_angle",reverse_angle)){
		ROS_INFO("Reverse angle: %d",reverse_angle);
	}else{
		reverse_angle = true;
		ROS_WARN("Reverse angle argument not found in ROS param server, using default value(%d).",reverse_angle);
	}
	
	bool write_vad_wavs;
	if (ros_nh.getParam("write_vad_wavs",write_vad_wavs)){
		ROS_INFO("Write VAD WAVs: %d",write_vad_wavs);
	}else{
		write_vad_wavs = false;
		ROS_WARN("Write VAD WAVs argument not found in ROS param server, using default value(%d).",write_vad_wavs);
	}
	
	bool graph_out;
	if (ros_nh.getParam("graph_out",graph_out)){
		ROS_INFO("Graph out: %d",graph_out);
	}else{
		graph_out = true;
		ROS_WARN("Graph out argument not found in ROS param server, using default value(%d).",graph_out);
	}
	
	bool connect_ports;
	if (ros_nh.getParam("connect_ports",connect_ports)){
		ROS_INFO("Connect JACK ports: %d",connect_ports);
	}else{
		connect_ports = true;
		ROS_WARN("Connect JACK ports argument not found in ROS param server, using default value(%d).",connect_ports);
	}
	
	int tde_method;
	if (ros_nh.getParam("tde_method",tde_method)){
		ROS_INFO("TDE method: %d",tde_method);
	}else{
		tde_method = 2;
		ROS_WARN("TDE method argument not found in ROS param server, using default value(%d).",tde_method);
	}
	
	int tracking_method;
	if (ros_nh.getParam("tracking_method",tracking_method)){
		ROS_INFO("Tracking method: %d",tracking_method);
	}else{
		tracking_method = 2;
		ROS_WARN("Tracking method argument not found in ROS param server, using default value(%d).",tracking_method);
	}
	
	
	std::cout << "\nSoundLoc: Stand by while everything is initialized.\n\n";
	
	soundloc_init(distance_between_mics, number_of_samples_of_noise, noise_threshold, noise_peak_change, coherence_threshold, tde_method, tracking_method, reverse_angle, graph_out, connect_ports, write_vad_wavs);
	//soundloc_init(distance_between_mics, number_of_samples_of_noise, noise_threshold, noise_peak_change, coherence_threshold, TDE_GCCPHAT, TRACK_MCMCDA, reverse_angle, graph_out, connect_ports, write_vad_wavs);
	//soundloc_init(distance_between_mics, number_of_samples_of_noise, noise_threshold, noise_peak_change, coherence_threshold, TDE_GCCPHAT, TRACK_KALMAN, reverse_angle, graph_out, connect_ports, write_vad_wavs);
	//soundloc_init(distance_between_mics, number_of_samples_of_noise, noise_threshold, noise_peak_change, coherence_threshold, TDE_GCCPHAT, TRACK_CLUSTER, reverse_angle, graph_out, connect_ports, write_vad_wavs);
	
	// For soundloc_init ALL ARGUMENTS ARE OPTIONAL EXCEPT FOR distance_between_mics
	//soundloc_init(distance_between_mics);
	
	
	ros::spin();
	return 0;
}

/*
	Source coded by Caleb Rascon, 2013
	IIMAS, UNAM
	México
	
	This library, using JACK, obtains the angle of different sources in the range of [-180 180].
	
	It requires three channels of independent directional microphones on top of the robot, set at 0 dB gain.
	The following is a simple diagram of their assumed setup:
	
														front
														/\
													left right
	
	
	The microphones should be set up in a triangle, with an equal distance between their ends and
	a 60 degree angle between left-front, front-right, and right-left.
	
	The angle of one source detected from each pair is compared with the other two to decide which
	provides the best resolution (e.g. the one most in front of a pair, e.g. the smallest angle of all).
	Then, a redundancy check is carried out to see if all the pairs are detecting the source coming from
	the same direction.
	
	The angle given starts from the line in which the front mic lies. A negative angle is to the
	left; a positive angle is to the right, all the way to the back of the front mic.
	
	If more than source is active, because the sample size is small, it is likely that one sample is able
	to get through redundancy check and provide an adequate angle of one source, and the next redundant
	angle of another source. For this purpose, an expectation system is built on top of the DOA estimator
	to mantain in memory the number of sources and each of their angles, which provides a basic CASA system.
	
	When running, it will output to standard out its estimate of the angle, and the current angle of the robot
	simulation until it reaches a good enough angle.
	
	Requirements:
		- A soundcard with at least three channels of input, and set to be the system's default.
		- A running JACK Server.
	
	
	ToDo's (2010-02-11):
	- Some sort of filtering. As of now, it doesn't filter the captured data, so any sound is
		considered a source. This causes the program to be very sensitive to noise.
		(2010-02-23) -> This point has been addressed by using a bandpass Iir filter from 1 to 4 kHz.
						A basic V.A.D. training system is also used. However, we need to:
	- (2010-02-23) Get a more sophisticated Voice Activity Detection system. As of now, only a
		comparison of mean magnitudes is used and seems artificial.
		(2012-05-23) -> Settled on an in-house V.A.D. system that is adaptive to ambient noise.
	- (2012-05-23) Build a multi-user DOA estimator using built-in expectations.
	- (2013-09-11) Started the port out of using STK.
	- (2013-10-21) Ported algorithm to library form: libmultisoundloc.a
*/

#include <cstdlib>
#include <stdio.h>
#include <cmath>
#include <string.h>
#include <iostream>
#include <unistd.h>
#include <algorithm>
#include <sys/time.h>
#include <sys/times.h>
#include <vector>


//JACK
#include <jack/jack.h>

//To write audio files
#include <sndfile.h>

//To convert between sample rates
#include <samplerate.h>

//PlPlot stuff
#include "plstream.h"

// Include FFTW header
#include <complex>
#include <fftw3.h>

// Include everything else we coded
#include "multisoundloc.h"

std::complex<double> *x_fft, *x_time, *y_fft, *y_time, *c_fft, *c_time;
fftw_plan x_forward, y_forward, c_inverse;

#define INVALID 999999
#define PI 3.141592653589793238462643383279502884
jack_port_t **input_ports;
jack_client_t *client;
double **frame_filt;
double **frame;
unsigned int frame_k = 0;
double **frame_res;
unsigned int frame_res_k = 0;

double *frame_wav;
unsigned int frame_wav_k = 0;

const int sound_speed = 343; //in m/s
int TDE = 0;
int TRACK = 0;
int ENDING = 0;
int ENDED = 0;
int JACK_FILLED_FRAME = 0;
int SOUNDLOC_CAPTURING = 0;
double refresh_angle_time = 10;

bool reverse_angle = 1;
bool write_vad_wavs = 1;
bool graph_out = 0;

double distance_between_mics; //empirically found: 0.1089
unsigned int number_of_samples;
int number_of_samples_of_noise;
double noise_threshold;
double noise_peak_change;
double coherence_threshold;
int number_of_doas_to_be_real_source = 2;
double min_resolution;
double source_diff;
double *w;
int maxdelay;
int *delays;
int delays_size;

std::vector < source > sources;
std::vector< doa > doas;
vad_results vad_res;

//GLib Thread stuff
GMutex mutex_sources;
GMutex mutex_doas;
GMutex mutex_sample_scan;
GMutex mutex_frame;

//sndfile stuff
SNDFILE * audio_file;
SF_INFO audio_file_info;
unsigned int audio_file_position = 0;
unsigned int audio_file_num = 0;
char audio_file_num_str[100];

//samplerate stuff
#define DEFAULT_CONVERTER SRC_SINC_MEDIUM_QUALITY
float * samplerate_buff_in;
SRC_STATE * samplerate_conv;
SRC_DATA samplerate_data;

//MCMCDA stuff
static const int MCMCMDA_population = 100;
static const int MCMCMDA_MaxMissDetections = 10;
MCMCDA mcmcda(MCMCMDA_population, MCMCMDA_MaxMissDetections, 100,5,false); 
MultiTargetTracker kalman_track(0.1,20,3,0.3,false); 
Scan sample_scan = Scan();

ofstream history_file;

unsigned int channels = 3;
unsigned int channel_offset = 0;
double sampleRate = 48000;

bool adding_stream = false;

double new_doa = INVALID;
unsigned int sample_number = 0;

plstream *pls;
PLFLT *plplot_x  = new PLFLT[1];
PLFLT *plplot_y  = new PLFLT[1];
PLFLT dtr = M_PI / 180.0;
std::vector < PLFLT > plplot_last_sources_x;
std::vector < PLFLT > plplot_last_sources_y;
std::vector < PLFLT > plplot_last_doas_x;
std::vector < PLFLT > plplot_last_doas_y;

//Sleep function in milliseconds
void millisleep(int milli){
	struct timespec st = {0};
	st.tv_sec = milli/1000;
	st.tv_nsec = (milli%1000)*1000000L;
	nanosleep(&st, NULL);
}

double diff_angle(double a, double b){
	double d = fabs(a-b);
	
	if (d > 180)
		d = 360-d;
	
	return d;
}

bool compare_observation (Observation i, Observation j) { return (i.value > j.value); }
bool compare_sources (source i, source j) { return (i.confidence > j.confidence); }

void write_to_history(){
	int i,j;
	stringstream ss;
	ss << setw(5) << setfill('0') << sample_number;
	
	history_file << endl << "Sample: " << ss.str() << endl;
	history_file << endl;
	
	//Writing Stream Info
	history_file << "DOAS: " << endl;
	g_mutex_lock (&mutex_sources);
	for(i = 0; i < sources.size(); i++){
		if(sources[i].real_source == 1){
			history_file << sources[i].main_doa << endl;
		}
	}
	g_mutex_unlock (&mutex_sources);
	history_file << endl;
	
	history_file << endl << "---" << endl;
}

double shift_doa(double d){
	double d_shifted = 0;
	d_shifted = d;
	d_shifted = -1*d_shifted;
	
	d_shifted += 90;
	if (d_shifted > 180)
		d_shifted -=360;
	
	return d_shifted;
}

void reload_doa_plot(){
	int i;
	double tmp_doa_value;
	//clearing page
	//pls->clear();
	
	// Draw doas
	//erasing last doas
	pls->col0( 0 );
	for ( i = 0; i < plplot_last_doas_x.size(); i++ ){
		plplot_x[0] = plplot_last_doas_x[i];
		plplot_y[0] = plplot_last_doas_y[i];
		pls->poin(1, plplot_x, plplot_y,1);
	}
	plplot_last_doas_x.clear();
	plplot_last_doas_y.clear();
	pls->col0( 2 );
	g_mutex_lock (&mutex_doas);
	for ( i = 0; i < doas.size(); i++ ){
		tmp_doa_value = doas[i].value;
		plplot_x[0] = cos( dtr * shift_doa(tmp_doa_value) );
		plplot_y[0] = sin( dtr * shift_doa(tmp_doa_value) );
		
		pls->poin(1, plplot_x, plplot_y,1);
		
		plplot_last_doas_x.push_back(plplot_x[0]);
		plplot_last_doas_y.push_back(plplot_y[0]);
	}
	g_mutex_unlock (&mutex_doas);
	
	// Draw sources
	//erasing last sources
	pls->col0( 0 );
	for ( i = 0; i < plplot_last_sources_x.size(); i++ ){
		plplot_x[0] = plplot_last_sources_x[i];
		plplot_y[0] = plplot_last_sources_y[i];
		pls->poin(1, plplot_x, plplot_y,4);
	}
	plplot_last_sources_x.clear();
	plplot_last_sources_y.clear();
	pls->col0( 11 );
	g_mutex_lock (&mutex_sources);
	for ( i = 0; i < sources.size(); i++ ){
		tmp_doa_value = sources[i].main_doa;
		plplot_x[0] = cos( dtr * shift_doa(tmp_doa_value) );
		plplot_y[0] = sin( dtr * shift_doa(tmp_doa_value) );
		
		pls->poin(1, plplot_x, plplot_y,4);
		
		plplot_last_sources_x.push_back(plplot_x[0]);
		plplot_last_sources_y.push_back(plplot_y[0]);
	}
	g_mutex_unlock (&mutex_sources);
	
	//pls->replot();
}

static gpointer plot_streams(gpointer data){
	std::cout << "SoundLoc.PlotSourceStreams: Starting plot_streams thread.\n";fflush(stdout);
	while (SOUNDLOC_CAPTURING == 0)
	        millisleep(200);
	time_t curr_time = time(NULL);
	
	
	std::cout << "SoundLoc.PlotSourceStreams: Setting up PlPlot...";fflush(stdout);
	// plplot initialization
	pls = new plstream(0,0,0,0,0,"xwin","");
	
	// Set device name
	pls->setopt("device",":0.0");
	pls->setopt("geometry","500x500");
	
	// Initialize PLplot.
	pls->init();
	
	// Set up viewport and window, but do not draw box.
	pls->env( -1.4, 1.4, -1.4, 1.4, 0, -2 );
	//pls->diplt(0,1,0,1);
	
	std::cout << "done\nSoundLoc.PlotSourceStreams: drawing base.\n";fflush(stdout);
	
	int   i, int_theta;
	char  text[5];
	PLFLT theta, dx, dy, r, offset;
	
	// Draw outer circle for polar grid
	pls->col0( 1);
	pls->arc( 0.0, 0.0, 0.8, 0.8, 0.0, 360.0, 0.0, 0 );
	for ( i = 0; i <= 11; i++ ){
		theta = 30.0 * i;
		dx = 0.8*cos( dtr * shift_doa(theta) );
		dy = 0.8*sin( dtr * shift_doa(theta) );
		
		// Draw radial spokes for polar grid.
		pls->join( 0.0, 0.0, dx, dy );
		int_theta = (int) ROUND( theta );
		if (int_theta > 180){
			int_theta -= 360;
		}
		sprintf( text, "%d", int_theta );
		
		dx = cos( dtr * shift_doa(theta) );
		dy = sin( dtr * shift_doa(theta) );
		// Write labels for angle.
		if ( theta < 9.99 ){
			offset = 0.45;
		}else if ( theta < 99.9 ){
			offset = 0.30;
		}else{
			offset = 0.15;
		}
		
		//Slightly off zero to avoid floating point logic flips at 90 and 270 deg.
		if ( dx >= -0.00001 )
			pls->ptex( dx, dy, dx, dy, -offset, text );
		else
			pls->ptex( dx, dy, -dx, -dy, 1. + offset, text );
	}
	
	//Draw title
	//pls->col0( 4 );
	//pls->mtex( "t", 2.0, 0.5, 0.5, "#frPLplot Example 3 - r(#gh)=sin 5#gh" );
	
	std::cout << "SoundLoc.PlotSourceStreams: starting to plot...\n";fflush(stdout);
	while(ENDING == 0){
		reload_doa_plot();
		millisleep(200);
	}
	
	std::cout << "SoundLoc.PlotSourceStreams: thread ended.\n";
}

double doa_mean(std::vector < doa > v){
	if (v.size() > 0){
		double this_mean = 0;
		
		for(int i = 0; i < v.size(); i++)
			this_mean += v[i].value;
		
		return this_mean/v.size();
	}else
		return 0;
}

void print_sources(){
	g_mutex_lock (&mutex_sources);
	if (sources.size() > 0){
		std::cout << "\n\n------------\n";
		for(int i = 0; i < sources.size(); i++){
			if (sources[i].real_source){
			std::cout << "\t------------\n";
				std::cout << "\tSource " << i << " ("<< sources[i].confidence <<"):\n";
				std::cout << "\t\t main doa : " << sources[i].main_doa << "\n";
				std::cout << "\t\t # of doas: " << sources[i].doas.size() << "\n";
			std::cout << "\t------------\n";
			}
		}
		std::cout << "------------\n\n\n";
		fflush(stdout);
	}
	g_mutex_unlock (&mutex_sources);
}

void remove_angle_from_source(source * s, int a){
	g_mutex_lock (&mutex_doas);
	if(doas.size() > 0)
		doas.erase(doas.begin());	//it is safe to assume that everytime this function is called
						//the oldest doa of the stream is the one being erased
	g_mutex_unlock (&mutex_doas);

	s->doas.erase(s->doas.begin()+a);
	
	s->main_doa = doa_mean(s->doas);
	
	if(s->doas.size() <= number_of_doas_to_be_real_source && s->real_source == 1){
		s->real_source = 0;
	}
	
	s->main_doa_time = time(NULL);
}

void add_angle_to_doa_list(double a){
	doa this_doa;
	this_doa.value = a;
	this_doa.reg_time = time(NULL) ;
	
	g_mutex_lock (&mutex_doas);
	doas.push_back(this_doa);
	g_mutex_unlock (&mutex_doas);
}

void add_angle_to_source(source * s, double a){
	doa this_doa;
	this_doa.value = a;
	this_doa.reg_time = time(NULL) ;
	
	s->doas.push_back(this_doa);
	
	s->main_doa = doa_mean(s->doas);
	s->main_doa_time = this_doa.reg_time;
	
	if(s->doas.size() > number_of_doas_to_be_real_source && s->real_source == 0){
		s->real_source = 1;
	}
	
	if (s->real_source)
		s->last_doa = a;
}

static gpointer constant_write_history(gpointer data){
	std::cout << "SoundLoc.WriteHistory: Starting constant_write_history thread.\n";fflush(stdout);
	while (SOUNDLOC_CAPTURING == 0)
		millisleep(225);
	time_t curr_time = time(NULL);
	while(ENDING == 0){
		curr_time = time(NULL);
		
		//Writing to history
		sample_number++;
		write_to_history();
		
		millisleep(100);
	}
	std::cout << "SoundLoc.WriteHistory: constant_write_history thread ended.\n";fflush(stdout);
}

static gpointer track_refresh(gpointer data){
	std::cout << "SoundLoc.TrackRefresh: Starting track_refresh thread.\n";fflush(stdout);
	while (SOUNDLOC_CAPTURING == 0)
		millisleep(200);
	int j;
	while(ENDING == 0){
		//Updating Sources
		g_mutex_lock (&mutex_sample_scan);
		if (!sample_scan.empty()){
			std::sort(sample_scan.begin(), sample_scan.end(), compare_observation);
			
			//Updating tracker
			if (TRACK == TRACK_MCMCDA){
				mcmcda.update(sample_scan);
				
				//Obtaining new scans "current_states"
				ScanSet scans_result = mcmcda.getScans();
				TrackSet tracks_result = mcmcda.getTracks();
				print_tracks( tracks_result , scans_result );
				Scan current_states = mcmcda.getCurrentStates(tracks_result);
				
				//Recreating "sources" variable
				g_mutex_lock (&mutex_sources);
				sources.clear();
				std::cout << "MCMCDA Current States:: ";
				for (Scan::iterator it = current_states.begin(); it != current_states.end(); ++it){
					std::cout << it->value << ", ";
					
					source this_source;
					this_source.real_source = 1;
					add_angle_to_source(&this_source,it->value);
					sources.push_back(this_source);
				}
				g_mutex_unlock (&mutex_sources);
				std::cout<< std::endl;
				sample_scan.clear();
			}else if (TRACK == TRACK_KALMAN){
				kalman_track.update(sample_scan);
				
				//Obtaining new scans "current_states"
				Scan current_states = kalman_track.getCurrentStates();
				
				//Recreating "sources" variable
				g_mutex_lock (&mutex_sources);
				sources.clear();
				std::cout << "Kalman Current States:: ";fflush(stdout);
				for (Scan::iterator it = current_states.begin(); it != current_states.end(); ++it){
					std::cout << it->value << ", ";fflush(stdout);
					
					source this_source;
					this_source.real_source = 1;
					add_angle_to_source(&this_source,it->value);
					this_source.confidence = 0;
					
					g_mutex_lock (&mutex_doas);
			        	for(j = 0; j < doas.size(); j++){
				        	if(diff_angle(doas[j].value, it->value) < 5){
						        this_source.confidence = this_source.confidence +1;
				        	}
				        }
				        g_mutex_unlock (&mutex_doas);
				        
					sources.push_back(this_source);
				}
				g_mutex_unlock (&mutex_sources);
				std::cout<< std::endl;fflush(stdout);
				
				
				sample_scan.clear();
			}
			std::cout << "Ordering sources:: ";fflush(stdout);
			g_mutex_lock (&mutex_sources);
			std::sort(sources.begin(),sources.end(),compare_sources);
			g_mutex_unlock (&mutex_sources);
			std::cout << std::endl;fflush(stdout);
		}
		g_mutex_unlock (&mutex_sample_scan);

		millisleep(100);
	}
	std::cout << "SoundLoc.TrackRefresh: track_refresh ended.\n";fflush(stdout);
}

static gpointer check_old_angles(gpointer data){
	std::cout << "SoundLoc.CheckOldAngles: Starting check_old_angles thread.\n";fflush(stdout);
	while (SOUNDLOC_CAPTURING == 0)
		millisleep(200);
	time_t curr_time = time(NULL);
	
	int i,j,k;
	
	while(ENDING == 0){
		curr_time = time(NULL);
		
		//purge out old registered doas
		if (TRACK == TRACK_CLUSTER){
			for(i = 0; i < sources.size(); i++){
				for(j = 0; j < sources[i].doas.size(); j++){
					if(curr_time - sources[i].doas[j].reg_time > refresh_angle_time){
						std::cout << "SoundLoc.CheckOldAngles: Removing doa " << sources[i].doas[j].value << " from -> " << i << ":" << sources[i].main_doa << "° ("<< sources[i].doas.size() <<")\n";
						remove_angle_from_source(&sources[i],j);
						std::cout << "\t new main_doa: "<< sources[i].main_doa << "°\n";
						fflush(stdout);
						
						if(sources[i].doas.size() == 0){
							std::cout << "SoundLoc.CheckOldAngles: Removing source -> " << i <<"\n";
							fflush(stdout);
							g_mutex_lock (&mutex_sources);
							sources.erase(sources.begin()+i);
							g_mutex_unlock (&mutex_sources);
							i--;
							print_sources();
							break;
						}else{
							print_sources();
						}
					}
				}
			}
		}else if (TRACK == TRACK_KALMAN || TRACK == TRACK_MCMCDA){
		    //purging source list
			for(i = 0; i < sources.size(); i++){
				if(curr_time - sources[i].doas[0].reg_time > refresh_angle_time){
					std::cout << "SoundLoc.CheckOldAngles: Purging source -> " << i << ":" << sources[i].main_doa << "°...";
					fflush(stdout);
					remove_angle_from_source(&sources[i],0);
					
					if(sources[i].doas.size() == 0){
						g_mutex_lock (&mutex_sources);
						sources.erase(sources.begin()+i);
						g_mutex_unlock (&mutex_sources);
						i--;
    					std::cout << "and removed.\n";
						fflush(stdout);
						print_sources();
						break;
					}else{
    					std::cout << "weird, it couldn't be removed.\n";
						fflush(stdout);
						print_sources();
					}
				}
			}

            //purgins doa list
			for(j = 0; j < doas.size(); j++){
				if(curr_time - doas[j].reg_time > refresh_angle_time){
					std::cout << "SoundLoc.CheckOldAngles: Removing doa " << doas[j].value << "° ("<< doas.size() <<")\n";
					g_mutex_lock (&mutex_doas);
					doas.erase(doas.begin()+j);
					g_mutex_unlock (&mutex_doas);
					break;
				}
			}
		}else{
			for(j = 0; j < doas.size(); j++){
				if(curr_time - doas[j].reg_time > refresh_angle_time){
					std::cout << "SoundLoc.CheckOldAngles: Removing doa " << doas[j].value << "° ("<< doas.size() <<")\n";
					g_mutex_lock (&mutex_doas);
					doas.erase(doas.begin()+j);
					g_mutex_unlock (&mutex_doas);
					break;
				}
			}
		}
		
		millisleep(200);
	}
	std::cout << "SoundLoc.CheckOldAngles: check_old_angles thread ended.\n";
}

void soundloc_finalize(){
	std::cout << "SoundLoc: Telling all threads to end.\n";
	
	ENDING = 1;
	do{millisleep(200);}while(ENDED == 0);
	
	history_file.close();
	jack_client_close (client);
}

void soundloc_clear(){
	g_mutex_lock (&mutex_sources);
	sources.clear();
	g_mutex_unlock (&mutex_sources);
	
	g_mutex_lock (&mutex_doas);
	doas.clear();
	g_mutex_unlock (&mutex_doas);
}

double sign(double x){
	double s = 1;
	if (x < 0)
		s = -1;
	return s;
}

double average_magnitude (double x[], int size){
	int i;
	double mag = 0;
	
	for (i=0; i<size; i++){
		mag += fabs(x[i]);
	}
	
	mag /= size;
	
	return mag;
}

std::vector<int> get_max_cross_correlation_location (double x[], double y[],bool v=false){
	//x and y are vectors of size n
	
	std::vector<int> this_delays;
	int i,j,delay;
	double mx,my,sx,sy,sxy,denom,r;
	
	double min_corr = INVALID, max_corr = -1, corr_mean = 0;
	int max_corr_location = -maxdelay;
	
	/* Calculate the mean of the two series x[], y[] */
	mx = 0;
	my = 0;   
	for (i=0;i<number_of_samples;i++) {
		mx += x[i];
		my += y[i];
	}
	mx /= number_of_samples;
	my /= number_of_samples;
	
	/* Calculate the denominator */
	sx = 0;
	sy = 0;
	for (i=0;i<number_of_samples;i++) {
		sx += (x[i] - mx) * (x[i] - mx);
		sy += (y[i] - my) * (y[i] - my);
	}
	denom = sqrt(sx*sy);
	
	if (denom == 0)
		return this_delays;
	
	/* Calculate the correlation series */
	for (delay=-maxdelay;delay<=maxdelay;delay++) {
		sxy = 0;
		for (i=0;i<number_of_samples;i++) {
			j = i + delay;
			
			if (j < 0 || j >= number_of_samples)
				//continue;
				sxy += (x[i] - mx) * (-my);
			else
				sxy += (x[i] - mx) * (y[j] - my);
		}
		r = sxy / denom;
		
		corr_mean += r;
		
		if (r > max_corr){
			max_corr = r;
			max_corr_location = delay;
		}
		
		if (r < min_corr){
		    min_corr = r;
		}
	}
	
	corr_mean /= (maxdelay*2)+1;
	
	if (v)
		std::cout << "SoundLoc: (delay, corr, corr_mean): " << max_corr_location << ", " << max_corr << ", " << corr_mean << "\n";
	
	if (max_corr > fabs(corr_mean)*5)
		this_delays.push_back(max_corr_location);
	
	return this_delays;
}

double hamming(unsigned int buffer_i, unsigned int buffer_size){
	// PI 3.14159265
	return 0.54 - 0.46*cos(2*3.14159265*buffer_i/(buffer_size-1));
}

double hann(unsigned int buffer_i, unsigned int buffer_size){
	// PI 3.14159265
	return 0.5 - 0.5*cos(2*3.14159265*buffer_i/(buffer_size-1));
}

std::vector<int> find_peaks(double *m, int num_peaks = 1){
	std::vector<int> peaks;
	std::vector<double> peaks_corr;
	double min_corr = INVALID, max_corr = -1, corr_mean = 0;
	int max_corr_location = -maxdelay;
	int i,j;
	
	for (i = 0; i < delays_size; i++) {
		corr_mean += m[i];
		
		if (m[i] > max_corr){
			max_corr = m[i];
			max_corr_location = delays[i];
		}
		
	}
	corr_mean /= delays_size;
	
	double lowest = fabs(corr_mean)+175;
	//std::cout << "SoundLoc: (delay, corr, corr_mean): " << max_corr_location << ", " << max_corr << ", " << corr_mean << "\n";
	
	if (max_corr <= lowest)
		return peaks;
	
	peaks.push_back(max_corr_location);
	peaks_corr.push_back(max_corr);
	for(i=1;i<num_peaks;i++){
		max_corr=-1;
		for(j=0;j<delays_size;j++){
			if(m[j]>max_corr&&m[j]<peaks_corr[i-1]){
				max_corr=m[j];
				max_corr_location = delays[j];
			}
		}
		if (max_corr > lowest){
			peaks.push_back(max_corr_location);
			peaks_corr.push_back(max_corr);
		}else{
			break;
		}
	}
	
	return peaks;
}

std::vector<int> get_max_cross_correlation_location_gcc_phat (double x[], double y[],bool v=false){
	int i,j;
	
	for(i = 0; i < number_of_samples; i++){
		x_time[i] = x[i]*hann(i, number_of_samples);// + 0.0*I;
		y_time[i] = y[i]*hann(i, number_of_samples);// + 0.0*I;
	}
	
	fftw_execute(x_forward);
	fftw_execute(y_forward);
	
	std::complex<double> tmp, this_x, this_y;
	double tmpnorm;
	for(i = 0;i < number_of_samples/2;i++){
		if(w[i] >= 1000 && w[i] <= 4000){
			this_x = conj(x_fft[i]);
			this_y = y_fft[i];
			tmp = this_y*this_x;
			tmpnorm = abs(tmp);
			c_fft[i] = tmp;
			
			if(tmpnorm > 0.000001)
				c_fft[i] /= tmpnorm;
			else
				c_fft[i] /= 0.000001;
			
			this_x = conj(x_fft[number_of_samples-1-i]);
			this_y = y_fft[number_of_samples-1-i];
			tmp = this_y*this_x;
			tmpnorm = abs(tmp);
			c_fft[number_of_samples-1-i] = tmp;
			
			if(tmpnorm > 0.000001)
				c_fft[number_of_samples-1-i] /= tmpnorm;
			else
				c_fft[number_of_samples-1-i] /= 0.000001;
		}else{
			c_fft[i] = std::complex<double>(0,0);//0.0 + 0.0*I;
			c_fft[number_of_samples-1-i] = std::complex<double>(0,0); //0.0 + 0.0*I;
		}
		
		//c_fft[i] *= hann(i,number_of_samples);
		//c_fft[number_of_samples-1-i] *= hann(number_of_samples-1-i,number_of_samples);
	}
	
	// getting the filtered result in the time domain
	fftw_execute(c_inverse);
	
	double *c = new double[delays_size];
	for(i = 0; i < maxdelay; i++){
		c[maxdelay+i] = real(c_time[i]);
		c[maxdelay-1-i] = real(c_time[number_of_samples-1-i]);
	}
	c[maxdelay*2] = real(c_time[maxdelay]);
	
	std::vector<int> this_delays = find_peaks(c);
	
	if (v){
		std::vector<double> this_c;
		std::vector<double> this_d;
		
		for(i = 0; i < delays_size; i++){
			this_c.push_back(c[i]);
			this_d.push_back(delays[i]);
		}
/*
		for(i = 0; i < number_of_samples/4; i++){
			this_c.push_back(x[i]);
			this_d.push_back(i);
		}
*/
	}
	
	return this_delays;
}

std::vector<double> get90angle(double left[], double right[], int index, unsigned int size){
	std::vector<int> delay;
	std::vector<double> angles;
	double this_angle;
	
	//using cross_correlation to obtain delay
	if (TDE == TDE_GCCPHAT){
		delay = get_max_cross_correlation_location_gcc_phat(right,left);
	}else if (TDE == TDE_CC){
		delay = get_max_cross_correlation_location(right,left);
	}
	
	if (delay.size() > 0){
		for (int i = 0; i <delay.size(); i++){
			//using cop-out model to obtain angle from delay (-90, 90)
			angles.push_back(asin(((delay[i]/sampleRate)*sound_speed)/distance_between_mics)*180/PI);
			
			//std::cout << "SoundLoc: angle, delay: \t" << this_angle << ",\t" << delay  << " (" << maxdelay <<")\n";
		}
	}
	return angles;
}

// Compute the angle of rotation and its direction (positive to left, negative to right; Player's weird)
void compare_angles(double *rotation, double *direction, double target_theta, double curr_theta){
	if (curr_theta >= 0 && target_theta <= 0){
		*rotation = curr_theta + fabs(target_theta);
		*direction = -1;   
	}else if (curr_theta >= 0 && target_theta >= 0){
		*rotation = fabs(target_theta - curr_theta);
		if (target_theta >= curr_theta) {
			*direction = 1;
		}else{
			*direction = -1;
		}
	}else if (curr_theta <= 0 && target_theta <= 0){
		*rotation = fabs(target_theta - curr_theta);
		if (target_theta <= curr_theta) {
			*direction = -1;
		}else{
			*direction = 1;
		}  
	}else if (curr_theta <= 0 && target_theta >= 0){
		*rotation = fabs(curr_theta) + target_theta;
		*direction = 1;   
	} 
	
	if (*rotation > 180) {
		*rotation = 360 - *rotation;  
		*direction = -(*direction); 
	}
	
	if (rotation < 0)
		*rotation = fabs(*rotation);
}

int jack_callback (jack_nframes_t nframes, void *arg){
	jack_default_audio_sample_t **buf_in;
	
	buf_in = (jack_default_audio_sample_t**)malloc(channels*sizeof(jack_default_audio_sample_t *));
	
	for(int i = 0;i <channels;i++){
		buf_in[i] = (jack_default_audio_sample_t*)jack_port_get_buffer (input_ports[i], nframes);
	}
	
	int k;
	
	//If there is still some residue, putting it in the current frame
	if(!JACK_FILLED_FRAME && frame_res_k != 0){
		for (k = 0; k < frame_res_k;k++){
			g_mutex_lock (&mutex_frame);
			for (int i = 0; i < channels;i++){
				frame[i][frame_k] = frame_res[i][k];
			}
			g_mutex_unlock (&mutex_frame);
			frame_k++;
		
			if (frame_k >= number_of_samples){
				JACK_FILLED_FRAME = 1;
				frame_k = 0;
				frame_res_k = 0;
				break;
			}
		}
		
		if (k == frame_res_k)
			frame_res_k = 0;
	}
	
	if(!JACK_FILLED_FRAME && frame_res_k == 0){
		for (k = 0; k < nframes;k++){
			g_mutex_lock (&mutex_frame);
			for (int i = 0; i < channels;i++){
				frame[i][frame_k] = buf_in[i][k];
			}
			g_mutex_unlock (&mutex_frame);
			frame_k++;
			
			if (frame_k >= number_of_samples){
				JACK_FILLED_FRAME = 1;
				frame_res_k = 0;
				frame_k = 0;
				break;
			}
		}
		
		if(JACK_FILLED_FRAME && k < nframes){
			for (int k2 = k; k2 < nframes;k2++){
				for (int i = 0; i < channels;i++){
					frame_res[i][frame_res_k] = buf_in[i][k2];
				}
				frame_res_k++;
			}
		}
	}else{
		for (k = 0; k < nframes;k++){
			for (int i = 0; i < channels;i++){
				frame_res[i][frame_res_k] = buf_in[i][k];
			}
			frame_res_k++;
		}
	}

/*
	printf("wrote %d bytes, frame size: %d\n",nframes,number_of_samples);
	printf("\tframe_k at %d\n",frame_k);
	printf("\tframe_res_k at %d\n",frame_res_k);
	printf("\tJACK_FILLED_FRAME at %d\n",JACK_FILLED_FRAME);
	fflush(stdout);
*/
	
	return 0;
}

void filter_frame(){
/*
//Set up filter (Band Pass: 1 kHz to 4 kHz, speech)

b
std::vector<StkFloat> numerator;
numerator.push_back(0.0348);
numerator.push_back(0);
numerator.push_back(-0.0696);
numerator.push_back(0);
numerator.push_back(0.0348);

a
std::vector<StkFloat> denominator;
denominator.push_back(1);
denominator.push_back(-3.2680);
denominator.push_back(4.1247);
denominator.push_back(-2.3984);
denominator.push_back(0.5466);
*/

	double b [] = {0.0348,0,-0.0696,0,0.0348};
	double a [] = {1,-3.2680,4.1247,-2.3984,0.5466};
	int poly_order = 5;
	
	int i,j,c;
	
	for (c = 0; c<channels; c++){
		for (i = 0; i<poly_order; i++){
			//frame_filt[c][i] = frame[c][i];
			frame_filt[c][i] = 0;
		}
		for (i = poly_order; i<number_of_samples; i++){
			frame_filt[c][i] = 0;
			for (j = 0; j<poly_order; j++){
				g_mutex_lock (&mutex_frame);
				frame_filt[c][i] += b[j]*frame[c][i-j];
				g_mutex_unlock (&mutex_frame);
				if (j > 0)
					frame_filt[c][i] -= a[j]*frame_filt[c][i-j];
			}
	    }
	}
}

std::vector<double *> find_candidates(std::vector<double> angleRL, std::vector<double> angleLF, std::vector<double> angleFR){
	std::vector<double *> angle_opt;
	double angleRL_opt[2], angleLF_opt[2], angleFR_opt[2];
	double angleRL_t = 0, angleLF_t = 0, angleFR_t = 0;
	double angle_coherence[8];
	unsigned int min_index;
	
	for (int rl = 0; rl < angleRL.size(); rl++){
		for (int lf = 0; lf < angleLF.size(); lf++){
			for (int fr = 0; fr < angleFR.size(); fr++){
				angleRL_opt[0] = (angleRL[rl] == 0) ? 0 : -angleRL[rl];
				angleRL_opt[1] = (angleRL[rl] > 0) ? -180+angleRL[rl] : 180+angleRL[rl];
				
				angleLF_opt[0] = 120-angleLF[lf];
				if (angleLF_opt[0] > 180)
				angleLF_opt[0] -= 360;
				angleLF_opt[1] = angleLF[lf]-60;
				
				angleFR_opt[0] = -120-angleFR[fr];
				if (angleFR_opt[0] <= -180)
				angleFR_opt[0] += 360;
				angleFR_opt[1] = angleFR[fr]+60;
				
/*
				std::cout << "angleRL        = " << angleRL[rl] << "\n";
				std::cout << "angleRL_opt[0] = " << angleRL_opt[0] << "\n";
				std::cout << "angleRL_opt[1] = " << angleRL_opt[1] << "\n";
				std::cout << "angleLF        = " << angleLF[lf] << "\n";
				std::cout << "angleLF_opt[0] = " << angleLF_opt[0] << "\n";
				std::cout << "angleLF_opt[1] = " << angleLF_opt[1] << "\n";
				std::cout << "angleFR        = " << angleFR[fr] << "\n";
				std::cout << "angleFR_opt[0] = " << angleFR_opt[0] << "\n";
				std::cout << "angleFR_opt[1] = " << angleFR_opt[1] << "\n\n";
*/
				
				angle_coherence[0] = (diff_angle(angleRL_opt[0],angleLF_opt[0])+diff_angle(angleRL_opt[0],angleFR_opt[0])+diff_angle(angleLF_opt[0],angleFR_opt[0]))/3; //000
				angle_coherence[1] = (diff_angle(angleRL_opt[0],angleLF_opt[0])+diff_angle(angleRL_opt[0],angleFR_opt[1])+diff_angle(angleLF_opt[0],angleFR_opt[1]))/3; //001
				angle_coherence[2] = (diff_angle(angleRL_opt[0],angleLF_opt[1])+diff_angle(angleRL_opt[0],angleFR_opt[0])+diff_angle(angleLF_opt[1],angleFR_opt[0]))/3; //010
				angle_coherence[3] = (diff_angle(angleRL_opt[0],angleLF_opt[1])+diff_angle(angleRL_opt[0],angleFR_opt[1])+diff_angle(angleLF_opt[1],angleFR_opt[1]))/3; //011
				angle_coherence[4] = (diff_angle(angleRL_opt[1],angleLF_opt[0])+diff_angle(angleRL_opt[1],angleFR_opt[0])+diff_angle(angleLF_opt[0],angleFR_opt[0]))/3; //100
				angle_coherence[5] = (diff_angle(angleRL_opt[1],angleLF_opt[0])+diff_angle(angleRL_opt[1],angleFR_opt[1])+diff_angle(angleLF_opt[0],angleFR_opt[1]))/3; //101
				angle_coherence[6] = (diff_angle(angleRL_opt[1],angleLF_opt[1])+diff_angle(angleRL_opt[1],angleFR_opt[0])+diff_angle(angleLF_opt[1],angleFR_opt[0]))/3; //110
				angle_coherence[7] = (diff_angle(angleRL_opt[1],angleLF_opt[1])+diff_angle(angleRL_opt[1],angleFR_opt[1])+diff_angle(angleLF_opt[1],angleFR_opt[1]))/3; //111
				
				min_index = 0;
				for (int i = 1; i < 8; i++){
					if(angle_coherence[i] < angle_coherence[min_index])
						min_index = i;
				}
				
				//std::cout << "SoundLoc.Capture: Coherence at " << angle_coherence[min_index];
				if (angle_coherence[min_index] < (coherence_threshold+min_resolution)){
					angleRL_t = angleRL_opt[(int)(floor(min_index/4))%2];
					angleLF_t = angleLF_opt[(int)(floor(min_index/2))%2];
					angleFR_t = angleFR_opt[min_index%2];
					
/*
					std::cout << "\nangleRL_t       = " << angleRL_t << "\n";
					std::cout <<   "angleLF_t       = " << angleLF_t << "\n";
					std::cout <<   "angleFR_t       = " << angleFR_t << "\n";
*/
					double *this_angle_opt = new double[6];
					this_angle_opt[0] = angleRL[rl];
					this_angle_opt[1] = angleLF[lf];
					this_angle_opt[2] = angleFR[fr];
					this_angle_opt[3] = angleRL_t;
					this_angle_opt[4] = angleLF_t;
					this_angle_opt[5] = angleFR_t;
					
					angle_opt.push_back(this_angle_opt);
					
				}
			}
		}
	}
	
	return angle_opt;
}

double get_angle_from_candidate(double * candidate){
	double angle = INVALID;
	double angleRL = candidate[0], angleLF = candidate[1], angleFR = candidate[2];
	double angleRL_t = candidate[3], angleLF_t = candidate[4], angleFR_t = candidate[5];
	
	if (fabs(angleRL) < fabs(angleLF) && fabs(angleRL) < fabs(angleFR)){
		if (fabs(30-fabs(angleRL)) <= 2*min_resolution){
			if (fabs(30-fabs(angleLF)) <= 2*min_resolution){
				angle = (angleRL_t+angleLF_t)/2;
			}else if (fabs(30-fabs(angleFR)) <= 2*min_resolution){
				angle = (angleRL_t+angleFR_t)/2;
			}else{
				angle = angleRL_t;
			}
		}else{
			angle = angleRL_t;
		}
		//std::cout << " -> New Direction Detected\n";
	}else if(fabs(angleLF) < fabs(angleRL) && fabs(angleLF) < fabs(angleFR)){
		if (fabs(30-fabs(angleLF)) <= 2*min_resolution){
			if (fabs(30-fabs(angleRL)) <= 2*min_resolution){
				angle = (angleLF_t+angleRL_t)/2;
			}else if (fabs(30-fabs(angleFR)) <= 2*min_resolution){
				angle = (angleLF_t+angleFR_t)/2;
			}else{
				angle = angleLF_t;
			}
		}else{
			angle = angleLF_t;
		}
		//std::cout << " -> New Direction Detected\n";
	}else if(fabs(angleFR) < fabs(angleLF) && fabs(angleFR) < fabs(angleRL)){
		if (fabs(30-fabs(angleFR)) <= 2*min_resolution){
			if (fabs(30-fabs(angleLF)) <= 2*min_resolution){
				angle = (angleFR_t+angleLF_t)/2;
			}else if (fabs(30-fabs(angleRL)) <= 2*min_resolution){
				angle = (angleFR_t+angleRL_t)/2;
			}else{
				angle = angleFR_t;
			}
		}else{
			angle = angleFR_t;
		}
		//std::cout << " -> New Direction Detected\n";
	}else{
		//std::cout << " -> Inapropriate Angle\n";
	}
	
	return angle;
}

void write_frame_to_audio_file(double *data_frames, int nframes, int end_of_file){
	int i,error;
	
	samplerate_data.input_frames = nframes;
	for (i = 0; i < nframes; i++)
		samplerate_data.data_in[i] = (float)(data_frames[i]);
	samplerate_data.end_of_input = end_of_file;
	
	if ((error = src_process (samplerate_conv, &samplerate_data))){
		printf ("\nSoundLoc: samplerate error -> %s\n", src_strerror (error)) ;
		exit (1) ;
	}
	
	//Output to file
	long written_frames = sf_writef_float(audio_file, samplerate_data.data_out, samplerate_data.output_frames_gen);
	if(written_frames != samplerate_data.output_frames_gen){
		printf("SoundLoc: samplerate error -> %s\n",sf_strerror(audio_file));
		exit(1);
	}
}

static gpointer capture_sound(gpointer data){
	std::cout << "SoundLoc.Capture: Starting capture_sound thread.\n";fflush(stdout);

	long samples, i,j,k;
	
	std::vector<double> angleRL, angleLF, angleFR;
	std::vector<double *> angle_opts;
	double anglev = 0, angle_t = 0;
	
	double angle;
	
	//arrays used to hold values to pass onto analytical functions
	double arr_left[number_of_samples];
	double arr_right[number_of_samples];
	double arr_front[number_of_samples];
	
	double mag_left = 0, mag_right = 0, mag_front = 0;
	double mag_main[number_of_samples_of_noise];
	double mag_avg = 0;
	double noise_mag = 0;
	double this_mag = 0;
	int frames_not_passed = 0;
	bool part_of_source = 0;
	bool in_silence = 1;
	
	min_resolution = asin(((1/sampleRate)*sound_speed)/distance_between_mics)*180/PI;
	std::cout << "SoundLoc.Capture: Using min resolution -> " << min_resolution << "\n";fflush(stdout);
	std::cout << "SoundLoc.Capture: Coherence threshold  -> " << (coherence_threshold+min_resolution) << "\n";fflush(stdout);
	
	if (TRACK == TRACK_CLUSTER){
		double source_diff_init = 10; //this value has not changed in tests, so hardcoding it
		
		std::cout << "SoundLoc.Capture: (TRACK_CLUSTER) Using source_diff    -> " << (double)((int)(source_diff_init/min_resolution)+1)*min_resolution << "\n";
		
		source_diff = (double)((int)(source_diff_init/min_resolution)+1)*min_resolution;
	}
	
	//learning magnitude of environmental noise
	std::cout << "TestVAD.Capture: Getting initial data for VAD -> ";
	fflush(stdout);
	
	for(j=0; j<number_of_samples_of_noise; j++){
		
		//wait for jack fill the filter frame
		while(!JACK_FILLED_FRAME){millisleep(5);}
		
		std::cout << "."; fflush(stdout);
		
		if (TDE == TDE_CC)
			filter_frame();
		
		mag_left = 0;
		mag_right = 0;
		mag_front = 0;
		
		for (i = 0; i<number_of_samples; i++){
			if(TDE == TDE_CC){
				arr_left[i] = frame_filt[0][i];
				arr_right[i] = frame_filt[1][i];
				arr_front[i] = frame_filt[2][i];
			}else if (TDE == TDE_GCCPHAT){
				g_mutex_lock (&mutex_frame);
				arr_left[i] = frame[0][i];
				arr_right[i] = frame[1][i];
				arr_front[i] = frame[2][i];
				g_mutex_unlock (&mutex_frame);
			}
			
			mag_left += fabs(arr_left[i]);
			mag_right += fabs(arr_right[i]);
			mag_front += fabs(arr_front[i]);
		}
		
		//tell jack to fill next frame
		JACK_FILLED_FRAME = 0;
		
		mag_main[j]= ((mag_left+mag_right+mag_front)/number_of_samples)/(channels);
	}
	
	mag_avg = 0;
	for(j=0; j<number_of_samples_of_noise; j++)
		mag_avg += mag_main[j];
	mag_avg /= number_of_samples_of_noise;
	
	noise_mag = mag_avg;
	
	std::cout << "done.\n";
	fflush(stdout);
	
	// Here's the runtime loop
	SOUNDLOC_CAPTURING = 1;
	
	//std::cout <<"SoundLoc.Capture: Runtime loop:\n";
	while ( ENDING == 0 ) {
		
		//wait for jack fill the filter frame
		while(!JACK_FILLED_FRAME){millisleep(5);}
		
		if (TDE == TDE_CC)
			filter_frame();
		
		mag_left = 0;
		mag_right = 0;
		mag_front = 0;
		
		for (i = 0; i<number_of_samples; i++){
			if(TDE == TDE_CC){
				arr_left[i] = frame_filt[0][i];
				arr_right[i] = frame_filt[1][i];
				arr_front[i] = frame_filt[2][i];
			}else if (TDE == TDE_GCCPHAT){
				g_mutex_lock (&mutex_frame);
				arr_left[i] = frame[0][i];
				arr_right[i] = frame[1][i];
				arr_front[i] = frame[2][i];
				g_mutex_unlock (&mutex_frame);
			}
			
			mag_left += fabs(arr_left[i]);
			mag_right += fabs(arr_right[i]);
			mag_front += fabs(arr_front[i]);
		}
		
		//tell jack to fill next frame
		JACK_FILLED_FRAME = 0;
		
		this_mag = ((mag_left+mag_right+mag_front)/number_of_samples)/(channels);
		//std::cout << "TestVAD.Capture: (" << in_silence << ") this_mag -> " << this_mag << " (noise_mag: " << noise_mag << ")" << " <> " << noise_mag+noise_threshold;
		
		//if ((thismag-mag_avg) > (mag_avg-min_mag)*mult_noise_threshold){
		if (!in_silence && this_mag > noise_mag+noise_threshold){
			//std::cout << " << PASSED\n"; fflush(stdout);
			frames_not_passed = 0;
			
			angleRL = get90angle(arr_right, arr_left, 0, number_of_samples);
			angleLF = get90angle(arr_left, arr_front, 1, number_of_samples);
			angleFR = get90angle(arr_front, arr_right, 2, number_of_samples);
			
			if (angleRL.size() > 0 && angleLF.size() > 0 && angleFR.size() > 0){
				
				angle_opts = find_candidates(angleRL,angleLF,angleFR);
				
				if (angle_opts.size() > 0){
					for (k = 0; k < angle_opts.size(); k++){
						angle = get_angle_from_candidate(angle_opts[k]);
						
						if(angle != INVALID){
							if (reverse_angle){
								std::cout << "\nangle before rev = " << angle << "\n"; fflush(stdout);
								if (angle > 0)
									angle = -179+angle;
								else
									angle = 180+angle;
							}
							
							if (angle <= 180 && angle >= -180){
								new_doa = angle;
								
								std::cout << "SoundLoc: Detected sound from " << angle  << "°\n"; fflush(stdout);
								add_angle_to_doa_list(angle);
								
								/* ********** START TRACKING PART ************* */
								
								if (TRACK == TRACK_CLUSTER){
									//CLUSTERING TRACKING
									//Check if there are no sources
									if (sources.size() == 0){
										//if so, create the first source
										source this_source;
										this_source.real_source = 0;
										add_angle_to_source(&this_source,angle);
										g_mutex_lock (&mutex_sources);
										sources.push_back(this_source);
										g_mutex_unlock (&mutex_sources);
										
										std::cout << "SoundLoc: First source -> 0:" << this_source.main_doa << "° ("<< this_source.doas.size() <<")\n"; fflush(stdout);
										fflush(stdout);
									}else{
										//if not, fist check if it's part of an existing source (close enough to its main_doa)
										for(j = 0; j < sources.size(); j++){
											//if the source has 1 element and the diff is less than 5 degress, or
											//if the source has more than 1 element and the diff is less than 10 degrees
											//the current angle is part of the source
											if((sources[j].doas.size() > number_of_doas_to_be_real_source && fabs(sources[j].main_doa - angle) < source_diff) || (sources[j].doas.size() <= number_of_doas_to_be_real_source && fabs(sources[j].main_doa - angle) < (source_diff/2))){
												std::cout << "SoundLoc: Part of source -> " << j << ":" << sources[j].main_doa << "° ("<< sources[j].doas.size() <<")\n";
												add_angle_to_source(&sources[j],angle);
												std::cout << "\t new main_doa: "<< sources[j].main_doa << "°\n";
												part_of_source = 1;
												fflush(stdout);
												break;
											}
										}
										
										//if it isn't part of a source, create a new one
										if(!part_of_source){
											source this_source;
											this_source.real_source = 0;
											add_angle_to_source(&this_source,angle);
											g_mutex_lock (&mutex_sources);
											sources.push_back(this_source);
											g_mutex_unlock (&mutex_sources);
											
											std::cout << "SoundLoc: New source -> " << sources.size()-1 << ":" << this_source.main_doa << "° ("<< this_source.doas.size() <<")\n";
											fflush(stdout);
										}
										
										part_of_source = 0;
									}
									g_mutex_lock (&mutex_sources);
									for(j = 0; j < sources.size(); j++){
									    sources[j].confidence = (float)sources[j].doas.size();
									}
									g_mutex_unlock (&mutex_sources);
									print_sources();
								}else if(TRACK == TRACK_MCMCDA || TRACK == TRACK_KALMAN){
									//MCMCMDA TRACKING
									Observation sample_observation;
									sample_observation.value = angle;
									sample_observation.tracked = false;

									g_mutex_lock (&mutex_sample_scan);
									sample_scan.push_back(sample_observation);
									g_mutex_unlock (&mutex_sample_scan);
								}
								
								/* ********** END TRACKING PART ************* */
								
							}else{
								std::cout << "SoundLoc.Player: Detected invalid sound from " << angle  << "°. Ignoring.\n";
								fflush(stdout);
								new_doa = INVALID;
							}
						}
					}
				}else{
					new_doa = INVALID;
				}
			}else{
				//std::cout << "SoundLoc: Invalid correlations.\n";
			}
		}
		else{
			frames_not_passed++;
			//std::cout << "\n";
		}
		
		//checking for change in activity
		mag_avg = 0;
		for(j=0; j<number_of_samples_of_noise; j++)
			mag_avg += mag_main[j];
		mag_avg /= number_of_samples_of_noise;
		
		if(in_silence && this_mag > mag_avg+noise_peak_change){
			in_silence = 0;
			noise_mag = mag_avg;
			
			//resetting magnitude frames
			for(j=0; j<number_of_samples_of_noise; j++){
				mag_main[j] = mag_avg;
			}
			
			if(write_vad_wavs){
				//open audio file
				audio_file_num++;
				sprintf(audio_file_num_str,"vadsegment_%010d.wav",audio_file_num);
				printf("SoundLoc: VAD activated. Writing data to %s\n",audio_file_num_str);
				audio_file = sf_open (audio_file_num_str,SFM_WRITE,&audio_file_info);
				if(audio_file == NULL){
					printf("SoundLoc: samplerate error -> %s\n",sf_strerror(NULL));
					exit(1);
				}
				src_reset (samplerate_conv);
				
				//fill audio file with this arr_left (mic1)
				write_frame_to_audio_file(arr_left, number_of_samples, SF_FALSE);
			}
		}else if(!in_silence && (this_mag < mag_avg-noise_peak_change || frames_not_passed > 5)){
			frames_not_passed = 0;
			in_silence = 1;
			
			//resetting magnitude frames
			for(j=0; j<number_of_samples_of_noise; j++){
				mag_main[j] = noise_mag;
			}
			
			if(write_vad_wavs){
				//fill audio file with this arr_left (mic1)
				write_frame_to_audio_file(arr_left,number_of_samples, SF_TRUE);
				
				//close audio file
				printf("SoundLoc: VAD deactivated. Closing %s\n",audio_file_num_str);
				sf_close(audio_file);
				
				//present VAD results
				vad_res.reg_time = time(NULL);
				strcpy(vad_res.audio_file_path,audio_file_num_str);
				vad_res.doas.clear();
				for(j = 0; j <sources.size(); j++){
					if(sources[j].real_source){
						vad_res.doas.push_back(sources[j].last_doa);
					}
				}
			}
		}else{
			//storing this magnitude for next frames
			for(j=1; j<number_of_samples_of_noise; j++){
				mag_main[j-1] = mag_main[j];
			}
			mag_main[number_of_samples_of_noise-1] = this_mag;
			
			if (!in_silence && write_vad_wavs){
				//fill audio file with this arr_left (mic1)
				write_frame_to_audio_file(arr_left,number_of_samples, SF_FALSE);
			}
		}
	}
	std::cout << "SoundLoc.Capture: capture_sound thread ended.\n";fflush(stdout);
	ENDED = 1;
}

void jack_shutdown (void *arg){
	exit (1);
}

void soundloc_init(double distance_between_mics_in, int number_of_samples_of_noise_in, double noise_threshold_in, double noise_peak_change_in, double coherence_threshold_in, int tde_method, int track_method, bool reverse_angle_in, bool graph_out_in, bool connect_ports, bool write_vad_wavs_in){
	const char *client_name = "soundloc_multi";
	jack_options_t options = JackNoStartServer;
	jack_status_t status;
	
	printf ("SoundLoc: initializing mutex locks... ");fflush(stdout);
	//if(!g_thread_supported()) g_thread_init(NULL);
	g_mutex_init(&mutex_sources);
	g_mutex_init(&mutex_doas);
	g_mutex_init(&mutex_sample_scan);
	g_mutex_init(&mutex_frame);
	printf ("done.\n");fflush(stdout);

	g_mutex_lock (&mutex_sources);
	millisleep(100);
	g_mutex_unlock (&mutex_sources);

	
	/* open a client connection to the JACK server */
	client = jack_client_open (client_name, options, &status);
	if (client == NULL){
		/* if connection failed, say why */
		printf ("SoundLoc: jack_client_open() failed, status = 0x%2.0x\n", status);
		if (status & JackServerFailed) {
			printf ("SoundLoc: Unable to connect to JACK server.\n");
		}
		exit (1);
	}
	if (status & JackServerStarted){
		printf ("SoundLoc: JACK server started\n");
	}
	if (status & JackNameNotUnique){
		client_name = jack_get_client_name(client);
		printf ("SoundLoc: unique name `%s' assigned\n", client_name);
	}
	
	jack_set_process_callback (client, jack_callback, 0);
	jack_on_shutdown (client, jack_shutdown, 0);
	sampleRate = jack_get_sample_rate (client);
	printf ("SoundLoc: engine sample rate: %.0f\n", sampleRate);
	fflush(stdout);
	
	input_ports = (jack_port_t **)malloc(channels*sizeof(jack_port_t *));
	char port_name[50];
	for(int i = 0;i <channels;i++){
		sprintf(port_name,"input_%d",i+1);
		input_ports[i] = (jack_port_t *)jack_port_register (client, port_name, JACK_DEFAULT_AUDIO_TYPE,JackPortIsInput, 0);
		if ((input_ports[i] == NULL)) {
			printf("SoundLoc: no more JACK ports available after %d\n",i);
			exit (1);
		}
	}
	
	distance_between_mics = distance_between_mics_in;
	number_of_samples_of_noise = number_of_samples_of_noise_in;
	noise_threshold = noise_threshold_in;
	noise_peak_change = noise_peak_change_in;
	coherence_threshold = coherence_threshold_in;
	TDE = tde_method;
	TRACK = track_method;
	
	if (TDE != TDE_CC && TDE != TDE_GCCPHAT){
		std::cout << "SoundLoc: Invalid Time Delay Estimation method (either TDE_CC or TDE_GCCPHAT). Exiting.\n";
		exit(1);
	}
	
	if (TRACK != TRACK_CLUSTER && TRACK != TRACK_MCMCDA && TRACK != TRACK_KALMAN){
		std::cout << "SoundLoc: Invalid Tracking method (can be TRACK_CLUSTER, TRACK_MCMCDA or TRACK_KALMAN). Exiting.\n";
		exit(1);
	}
	
	if(noise_threshold < 0){
		std::cout << "\nSoundLoc: Noise threshold must be greater than zero.\n";
		exit(1);
	}
	
	std::cout << "SoundLoc: Distance btw mics  -> " << distance_between_mics << "\n";
	std::cout << "SoundLoc: # of noise samples -> " << number_of_samples_of_noise << "\n";
	std::cout << "SoundLoc: Noise threshold    -> " << noise_threshold << "\n";
	std::cout << "SoundLoc: Noise peak change  -> " << noise_peak_change << "\n";
	std::cout << "SoundLoc: Max. Coherence     -> " << coherence_threshold << "\n";
	
	number_of_samples = (unsigned int)(sampleRate/10);
	std::cout << "SoundLoc: Sample size (base) -> " << number_of_samples << "\n";
	int number_of_fft_samples_curr = 0;
	int number_of_fft_samples_past = 0;
	do{
		number_of_fft_samples_past = number_of_fft_samples_curr;
		number_of_fft_samples_curr++;
		while(number_of_fft_samples_curr & (number_of_fft_samples_curr - 1)){ //searching for the next power of 2
			number_of_fft_samples_curr++;
		}
	}while(number_of_fft_samples_curr - (signed int)number_of_samples < 0);
	
	number_of_samples = number_of_fft_samples_past;
	std::cout << "SoundLoc: Sample size (end)  -> " << number_of_samples << "\n";
	
	//Allocating space for frame buffer
	frame = (double **)malloc(channels*sizeof(double *));
	frame_res = (double **)malloc(channels*sizeof(double *));
	frame_filt = (double **)malloc(channels*sizeof(double *));
	for(int i = 0;i <channels;i++){
		frame[i] = (double *)malloc(number_of_samples*sizeof(double));
		frame_res[i] = (double *)malloc(10*number_of_samples*sizeof(double));
		frame_filt[i] = (double *)malloc(number_of_samples*sizeof(double));
	}
	
	x_fft = (complex<double>*) fftw_malloc(sizeof(complex<double>) * number_of_samples);//(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * number_of_samples);
	x_time = (complex<double>*) fftw_malloc(sizeof(complex<double>) * number_of_samples);
	y_fft = (complex<double>*) fftw_malloc(sizeof(complex<double>) * number_of_samples);
	y_time = (complex<double>*) fftw_malloc(sizeof(complex<double>) * number_of_samples);
	c_fft = (complex<double>*) fftw_malloc(sizeof(complex<double>) * number_of_samples);
	c_time = (complex<double>*) fftw_malloc(sizeof(complex<double>) * number_of_samples);
	
	x_forward = fftw_plan_dft_1d(number_of_samples, reinterpret_cast<fftw_complex*>(x_time), reinterpret_cast<fftw_complex*>(x_fft), FFTW_FORWARD, FFTW_MEASURE);
	y_forward = fftw_plan_dft_1d(number_of_samples, reinterpret_cast<fftw_complex*>(y_time), reinterpret_cast<fftw_complex*>(y_fft), FFTW_FORWARD, FFTW_MEASURE);
	c_inverse = fftw_plan_dft_1d(number_of_samples, reinterpret_cast<fftw_complex*>(c_fft), reinterpret_cast<fftw_complex*>(c_time), FFTW_BACKWARD, FFTW_MEASURE);
	
	//Creating frequency array
	w = (double *)malloc(number_of_samples/2*sizeof(double));
	for (int i = 0; i < number_of_samples/2; i++){
		w[i] = ((sampleRate/2)/number_of_samples/2)*(i+1);
	}
	
	maxdelay = (int)((distance_between_mics/sound_speed)*sampleRate);
	delays_size = (2*maxdelay)+1;
	delays = (int *)malloc(delays_size*sizeof(int));
	int d = -maxdelay;
	printf("Soundloc: Using delays: ");
	for (int i = 0; i <delays_size; i++, d++){
		delays[i] = d;
		printf("(%d -> %d)  ",i,delays[i]);
	}
	printf("\n");fflush(stdout);
	
	//initializing stuff for libsndfile
	audio_file_info.samplerate = 16000;
	audio_file_info.channels = 1;
	int audio_file_format = SF_FORMAT_WAV | SF_FORMAT_PCM_16;
	audio_file_info.format = audio_file_format;
	
	printf("Soundloc: Setting audio file info to: \n");
	printf("\tSample Rate: %d\n",audio_file_info.samplerate);
	printf("\tChannels: %d\n",audio_file_info.channels);
	SF_FORMAT_INFO audio_format_info;
	audio_format_info.format = audio_file_format;
	sf_command(NULL,SFC_GET_FORMAT_INFO,&audio_format_info, sizeof (SF_FORMAT_INFO));
	printf("\tFormat: %s\n",audio_format_info.name);
	
	//initializing stuff for libsamplerate
	int samplerate_error;
	samplerate_conv = src_new (DEFAULT_CONVERTER,1,&samplerate_error);
	if(samplerate_conv == NULL){
		printf("Soundloc: libsamplerate error -> %s\n",src_strerror (samplerate_error));
		exit(1);
	}
	
	samplerate_data.src_ratio = (double)(((double)audio_file_info.samplerate)/((double)sampleRate));
	printf("Soundloc: Using samplerate ratio: %f\n", samplerate_data.src_ratio);
	if (src_is_valid_ratio (samplerate_data.src_ratio) == 0){
		printf ("Soundloc: libsamplerate error -> Sample rate change out of valid range.\n") ;
		exit (1) ;
	}
	samplerate_buff_in = (float *)malloc(number_of_samples*sizeof(float)); //necessary to avoid overlapping buffers
	samplerate_data.data_in = samplerate_buff_in;
	samplerate_data.data_out = (float *)malloc(number_of_samples*sizeof(float));
	samplerate_data.input_frames = 0;
	samplerate_data.output_frames = number_of_samples;
	samplerate_data.end_of_input = 0;
	
	reverse_angle = reverse_angle_in;
	write_vad_wavs = write_vad_wavs_in;
	
	if(reverse_angle)
		std::cout << "\nSoundLoc: Reversing angle.\n";
	else
		std::cout << "\nSoundLoc: Original angle.\n";
	
	graph_out = graph_out_in;
	
	if(graph_out)
		std::cout << "\nSoundLoc: Output to graph.\n";
	
	history_file.open ("./output.mdoa");
	
	if (jack_activate (client)) {
		printf ("SoundLoc: cannot activate client.\n");
		exit (1);
	}
	printf ("SoundLoc: JACK client activated.\n");fflush(stdout);

	std::cout << "SoundLoc: Connecting input ports.\n";fflush(stdout);
	if(connect_ports){
		const char **port_names = jack_get_ports (client, NULL, NULL, JackPortIsPhysical|JackPortIsOutput);
		if (port_names == NULL) {
			printf("SoundLoc: no physical capture ports\n");
			exit (1);
		}
		for(int i = 0;i <channels;i++){
			//printf ("3.7.%d\n",i);fflush(stdout);
			printf("SoundLoc: connecting %s to %s\n",port_names[i+channel_offset],jack_port_name (input_ports[i]));fflush(stdout);
			if (jack_connect (client, port_names[i+channel_offset], jack_port_name (input_ports[i]))) {
				printf("SoundLoc: cannot connect ports: %s <> %s\n",port_names[i+channel_offset],jack_port_name (input_ports[i]));
			}
		}
		free (port_names);
	}

	std::cout << "\nSoundLoc: Stand by while everything is initialized.\n\n";
		
	//Setting Threading Stuff
	std::cout << "SoundLoc: Starting up Threads...\n";fflush(stdout);
	GThread *capture_sound_thread = NULL;
	GThread *check_old_angles_thread = NULL;
	GThread *track_refresh_thread = NULL;
	GThread *constant_write_history_thread = NULL;
	GThread *plotstreams = NULL;
	
	//if(!g_thread_supported()) g_thread_init(NULL); //initializing threading system if it isn't already
	
	capture_sound_thread = g_thread_new( "capture_sound", (GThreadFunc) capture_sound, NULL);
	
	while(SOUNDLOC_CAPTURING==0){millisleep(1);}
	
	check_old_angles_thread = g_thread_new("check_old_angles", (GThreadFunc) check_old_angles, NULL);
	
	if(TRACK == TRACK_MCMCDA || TRACK == TRACK_KALMAN){
		track_refresh_thread = g_thread_new("track_refresh", (GThreadFunc) track_refresh, NULL);
	}
	
	constant_write_history_thread = g_thread_new("constant_write_history", (GThreadFunc) constant_write_history, NULL);
	
	if(graph_out){
		plotstreams = g_thread_new("plot_streams", (GThreadFunc) plot_streams, NULL);
	}
}

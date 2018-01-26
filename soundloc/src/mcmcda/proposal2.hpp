/**
 * @Brief Methods for performing Markov Chain Monte Carlo Data Association.
 * 
 * @author Gibran Fuentes Pineda <gibranfp@turing.iimas.unam.mx>
 * @date 2013
 */
#ifndef PROPOSAL_HPP
#define PROPOSAL_HPP

#include <random>

#include "partition.hpp"
#include "scans.hpp"

// Some constant numbers
static const uint kNumberOfMoves = 8;
static const uint kMinTrackSize = 2;
static const uint kTrackEnding = 2;

class Proposal
{
public:
	/*
	 * @brief Randomly selects a candidate from a set
	 * 
	 * @returns Selected candidate
	 */
	template<typename T>
	T sampleCandidate(std::vector<T>& candidates)
	{
		std::random_device seed{};
		std::default_random_engine generator{seed()};
		uint end = candidates.size() - 1;
		std::uniform_int_distribution<int> distribution{0, end};
     
		return candidates[distribution(generator)];
	}

	/*
	 * @brief Randomly selects a time from a range
	 * 
	 * @returns Selected time
	 */
	uint sampleTime(uint start_time, uint end_time)
	{          
		std::uniform_int_distribution<int> distribution_time{start_time, end_time};

		return distribution_time(generator_time);
	}

	/*
	 * @brief Randomly selects one track from a set of candidates
	 * 
	 * @returns Selected track
	 */
	TrackSet::iterator sampleTrack(TrackSet& tracks)
	{
		uint end = tracks.size() - 1;
		std::uniform_int_distribution<int> distribution{0, end};
		uint track_num = distribution(generator_track);
	       
		return tracks.begin() + track_num;
	}
     
	/*
	 * @brief Randomly decide whether a track must be ended or not
	 * 
	 * @returns True or false
	 */
	bool sampleEnding(void)
	{	       
		std::bernoulli_distribution distribution{0.01};
		     
		return distribution(generator_ending);
	}

	/*
	 * @brief Randomly selects a move
	 *
	 * @return Selected move
	 */
	uint sampleMove(Partition& current)
	{
		std::random_device seed{};
		std::default_random_engine generator_moves(seed());

		if (current.tracks.empty())
		    return 0;
		else if (current.tracks.size() == 1){ 
		    // When there is only one track, merge and switch moves are not possible
		    std::vector<int> possible_moves = {0, 1, 2, 4, 5, 6};
		    std::uniform_int_distribution<int> distribution{0, kNumberOfMoves - 3};
		    return possible_moves.at(distribution(generator_moves));
		}
		else{
		    // All moves are possible
		    std::uniform_int_distribution<int> distribution{0, kNumberOfMoves - 2};
		    return distribution(generator_moves);
		}

		return 0;
	}

	/*
	 * @brief Returns forward transition probability
	 * 
	 * @returns Forward probability
	 */
	double forward_prob()
	{
		return forward;
	}

	/*
	 * @brief Returns backward transition probability
	 *
	 * @returns Backward probability
	 */
	double backward_prob()
	{
		return backward;
	}


	/*
	 * @brief Class constructor
	 * 
	 * @param obs_db Observations
	 * @param nf Observation neighbor tree
	 * @param current_time Current sampling time 
	 * @param init_max_miss_detections Maximum time a target is not observed
	 */
	Proposal(NForest& nf, uint& current_time, uint& init_max_miss_detections) :
		nforest(nf), 
		max_miss_detections(init_max_miss_detections), 
		last_time (current_time) {}
     
	/**
	 * @brief Computes the probability of a given move
	 *
	 * @param move_number Number of the move
	 *
	 * @return Probability of the move
	 */
	double movePDF(uint move_number, Partition& partition)
	{
		if (partition.tracks.empty()){
		    if (move_number == 0)
				return 1.0;
		    else
				return 0.0;
		}
	       
		if (partition.tracks.size() < 2){
		    if (move_number == 3 || move_number == 7)
				return 0.0;
		    else
				return 1.0 / (double)(kNumberOfMoves - 2);
		}

		return 1.0 / (double)kNumberOfMoves;
	}
     
	/**
	 * @brief Extends a given track
	 *
	 * @param track Track to extend
	 * @param track false_positives Unassigned observations
	 *
	 * @returns Extension probability
	 */
	double extendTrack(Track& track, Track& false_positives)
	{
		bool end_track = false;
		double extend_probability = 1.0;
	    
		// track's last observation
		uint prev_time = track.back().time;
		uint prev_number = track.back().number; 
		while (prev_time < last_time && !end_track)
		{
		    std::vector<Item> candidate_pairs;
		    // find neighbors of previous observation as candidates
		    for (NeighborTime::iterator ob = nforest[prev_time][prev_number].begin(); ob != nforest[prev_time][prev_number].end(); ++ob){
				if (scans[ob->time][ob->number].tracked == false){
					Item item = {ob->time, ob->number};
					candidate_pairs.push_back(item);
				}
		    }
		
		    // randomly select a candidate observation to add to the track
		    if (!candidate_pairs.empty()){
				// randomly select a candidate observation
				Item item = sampleCandidate(candidate_pairs);		 
				track.push_back(item);

				// mark selected observation as tracked
				scans[item.time][item.number].tracked = true;

				// compute probability of extension
				extend_probability *= (1.0 / (double) candidate_pairs.size());
				
				prev_time = track.back().time;
				prev_number = track.back().number;
		    }
		    
		    // track termination distribution
		    end_track = sampleEnding();
		}
	       
		return extend_probability;
	}
     
	/**
	 * @brief Proposes a new track in the current partition from the unassigned observations
	 *
	 * @param current Current partition
	 */
	Partition birthMove(Partition current)
	{       
		std::random_device seed1{};
		generator_time.seed(seed1());

		double birth_move_prob = movePDF(0, current);
	       
		// select the new track birth time uniformly at random
		uint birth = sampleTime(0, scans.size() - 2);
	       
		// add initial observation to new track from possible observations at birth time
		Track new_track = Track();
		double track_probability = 1.0;
		std::vector<Item> candidates_first_observation;
		for (uint number = 0; number < scans[birth].size(); ++number)
		{
		    for (NeighborTime::iterator ob = nforest[birth][number].begin(); ob != nforest[birth][number].end(); ++ob){
				if (scans[birth][number].tracked == false){
					Item item = {birth, number};
					candidates_first_observation.push_back(item);
					break;
				}
		    }
		}
	       
		// randomly select a candidate for the new track
		if (!candidates_first_observation.empty()){		  
		    // select observation uniformly at random from candidates
		    Item item1_to_add = sampleCandidate(candidates_first_observation);
		    
		    // add selected observation to the new track
		    new_track.push_back(item1_to_add);
		
			// delete from false positives
			partition

		    // mark added observation as tracked
		    scans[item1_to_add.time][item1_to_add.number].tracked = true;
		    
		    // probability of selecting first observation
		    track_probability *= (1.0 / (double) candidates_first_observation.size());

		    // find the neighbors of the first observation
		    std::vector<Item> candidates_second_observation;
		    near_scan = sampleTime(0, max_miss_detections - 1);
		    if (last_time - birth >= max_miss_detections){
				for (NeighborTime::iterator ob = nforest[birth][item1_to_add.number][near_scan].begin(); 
					 ob != nforest[birth][item1_to_add.number][near_scan].end(); ++ob){
					if (scans[ob->time][ob->number].tracked == false){
						Item item = {ob->time, ob->number};
						candidates_second_observation.push_back(item);
					}
				}
		    }

		    if (!candidates_second_observation.empty()){
				// probability of selecting second observation
				track_probability *= (1.0 / (double) candidates_second_observation.size());
			 
				// select observation uniformly at random from candidates
				Item item2_to_add = sampleCandidate(candidates_second_observation);
			 
				// add selected observation to the new track
				new_track.push_back(item2_to_add);
			 
				// mark added observation as tracked
				scans[item2_to_add.time][item2_to_add.number].tracked = true;
		    }
		}

		if (new_track.size() >= kMinTrackSize){
		    // add further observations to new track
		    track_probability *= extendTrack(new_track, current.false_positives);
 		    
		    // add new track to current current
		    current.tracks.push_back( new_track );
		
		    // calculate transition probability
		    double birth_time_prob = (1.0 / ((double)scans.size() - 2));
		    
		    // p_forward_transition = p(move) * p(birth_time) * p(new_track)
		    forward = (double) birth_move_prob * birth_time_prob * track_probability;
		    
		    // p_backward_transition = p(move) * p(track_number)
		    backward = (double) movePDF(1, current) * (1.0 / (double)current.tracks.size());
		}
		else{ // move can not be done
		    forward = 0;
		    backward = 0;
		}
	       
		return current;
	}

	/**
	 * @brief Eliminates a randomly selected track from the current partition
	 *
	 * @param current Current partition
	 */ 
	Partition deathMove(Partition current)
	{
		double death_move_prob = movePDF(1, current);
	       
		// select a track uniformly at random
		TrackSet::iterator track_to_kill = sampleTrack(current.tracks);
     
		// add selected track observations to empty track
		for (Track::iterator it = track_to_kill->begin(); it != track_to_kill->end(); ++it){
		    //current.false_positives.push_back(*it);
		    scans[it->time][it->number].tracked = false;
		}

		double track_probability = 1.0;
		uint birth = track_to_kill->front().time;
		uint candidate_size = 0;
		uint observation_count = 0;
		// computes backward probability of first observation
		for (uint number = 0; number < scans[birth].size(); ++number)
		{
		    bool found_neighbor = false;
		    // find observations at birth time with at least one neighbor
		    for (int i = 0; i < max_miss_detections && found_neighbor == false; ++i) {
				for (NeighborTime::iterator ob = nforest[birth][number][i].begin(); 
					 ob != nforest[birth][number][i].end(); ++ob){
					if (scans[birth][number].tracked == false){
						candidate_size++;
						found_neighbor = true;
						break;
					}
				}
		    }
		}
	       
		// probability of selecting the first observation
		track_probability *= (1.0 / (double) candidate_size);
	       
		uint prev_time = track_to_kill->front().time;
		uint prev_number = track_to_kill->front().number;
		// computes backward probability of the rest of observations
		for (Track::iterator it = track_to_kill->begin() + 1; it != track_to_kill->end(); ++it)
		{
		    candidate_size = 0;
		    // find neighbors of previous observations as candidates
		    for (NeighborTime::iterator ob = nforest[prev_time][prev_number][it->time - prev_time - 1].begin(); 
				 ob != nforest[prev_time][prev_number][it->time - prev_time - 1].end(); ++ob){
				if (scans[ob->time][ob->number].tracked == false)
					candidate_size++;
		    }
		    
		    // probability of selecting the observation
		    track_probability *= (1.0 / (double) candidate_size);
		    prev_time = it->time;
		    prev_number = it->number;
		}
	       	            
		// calculate transition probability
		forward = death_move_prob * (1.0 / (double)current.tracks.size());
	       
		// delete track from current
		current.tracks.erase(track_to_kill);
	       
		// calculate transition probability backward
		backward = movePDF(0, current) * (1.0 / ((double) scans.size() - 2.0)) * track_probability;
	       
		return current;
	}

	/**
	 * @brief Splits a randomly selected track from the current partition
	 *
	 * @param current Current partition
	 */ 
	Partition splitMove(Partition current)
	{
		double split_move_prob = movePDF(2, current);
	       
		// find tracks with at least 4 observations
		std::vector<uint> splitable_tracks;
		for (uint i = 0; i < current.tracks.size(); ++i)
		    if (current.tracks[i].size() > 3)
				splitable_tracks.push_back(i);

		// if there are tracks to split
		if (!splitable_tracks.empty()){
		    // select the new track to be splitted
		    uint rnd = sampleCandidate(splitable_tracks);
		    TrackSet::iterator track_to_split = current.tracks.begin() + rnd;
		    
		    // randomly select a time
		    uint split_time = sampleTime(1, track_to_split->size() - 3);
		    Track::iterator split_point = track_to_split->begin() + split_time + 1;
		    
		    
		    // create new track with observations from the beginning to split time
		    Track splitted_rest(std::make_move_iterator(split_point),
								std::make_move_iterator(track_to_split->end()));
		    
		    // forward transition probability
		    forward = split_move_prob * (1.0 / (double) splitable_tracks.size());
		    forward *= 1.0 / (double) (track_to_split->size() - 3);
		    
		    // remove track's observations after selected time
		    track_to_split->erase(split_point, track_to_split->end());
		    
		    // add splitted track to partition
		    current.tracks.push_back(splitted_rest);
		    
		    // Compute backward probability
		    // find pairs of tracks to merge: 
		    // (T_i,T_j) | such that T_i last observation is near to T_j first observation
		    uint mergeable_tracks = 0;
		    for (uint i = 0; i < current.tracks.size() - 1; ++i){
				Item Ti_first = current.tracks[i].front();
				Item Ti_last = current.tracks[i].back();
				for (uint j = i + 1; j < current.tracks.size(); ++j){
					Item Tj_first = current.tracks[j].front();
					Item Tj_last = current.tracks[j].back();
					uint near_scan_Tj_last = Ti_first.time - Tj_last.time - 1;
					uint near_scan_Ti_last = Tj_first.time - Ti_last.time - 1;
					if (Ti_first.time > Tj_last.time && near_scan_Tj_last < max_miss_detections){
						NeighborObservation::iterator near_obs = nforest[Tj_last.time][Tj_last.number].begin() + near_scan_Tj_last;
						for (NeighborTime::iterator nt = near_obs->begin(); nt != near_obs->end(); ++nt){
							if (nt->time == Ti_first.time && nt->number == Ti_first.number){
								mergeable_tracks++;
								break;
							}
						}
					}
					else if (Tj_first.time > Ti_last.time && near_scan_Ti_last < max_miss_detections){
						NeighborObservation::iterator near_obs = nforest[Ti_last.time][Ti_last.number].begin() + near_scan_Ti_last;
						for (NeighborTime::iterator nt = near_obs->begin(); nt != near_obs->end(); ++nt){
							if (nt->time == Tj_first.time && nt->number== Tj_first.number){
								mergeable_tracks++;
								break;
							}
						}
					}
				}
		    }

			// backward transition probability
			backward = movePDF(3, current) * (1.0 / (double) mergeable_tracks);
		} 
		else { // move can not be done
		    forward = 0;
		    backward = 0;
		}

		return current;
	}

	/**
	 * @brief Merges two tracks from the current  partition
	 *
	 * @param current Current partition
	 */ 
	Partition mergeMove(Partition current)
	{
		double merge_move_prob = movePDF(3, current);
	       
		// find pairs of tracks to merge: 
		// (T_i,T_j) | such that T_i last observation is near to T_j first observation
		std::vector<Item> mergeable_tracks;
		for (uint i = 0; i < current.tracks.size() - 1; ++i){
			Item Ti_first = current.tracks[i].front();
			Item Ti_last = current.tracks[i].back();
			for (uint j = i + 1; j < current.tracks.size(); ++j){
				Item Tj_first = current.tracks[j].front();
				Item Tj_last = current.tracks[j].back();
				uint near_scan_Tj_last = Ti_first.time - Tj_last.time - 1;
				uint near_scan_Ti_last = Tj_first.time - Ti_last.time - 1;
				if (Ti_first.time > Tj_last.time && near_scan_Tj_last < max_miss_detections){
					NeighborObservation::iterator near_obs = nforest[Tj_last.time][Tj_last.number].begin() + near_scan_Tj_last;
					for (NeighborTime::iterator nt = near_obs->begin(); nt != near_obs->end(); ++nt){
						if (nt->time == Ti_first.time && nt->number== Ti_first.number){
							Item item = {j, i};
							mergeable_tracks.push_back(item);
							break;
						}
					}
				}
				else if (Tj_first.time > Ti_last.time && near_scan_Ti_last < max_miss_detections){
					NeighborObservation::iterator near_obs = nforest[Ti_last.time][Ti_last.number].begin() + near_scan_Ti_last;
					for (NeighborTime::iterator nt = near_obs->begin(); nt != near_obs->end(); ++nt){
						if (nt->time == Tj_first.time && nt->number== Tj_first.number){
							Item item = {i, j};
							mergeable_tracks.push_back(item);
							break;
						}
					}
				}
		    }
		}
     
		// if there are mergeable tracks to merge
		if (!mergeable_tracks.empty()){
		    // randomly select a candidate pair to merge
		    Item rnd = sampleCandidate(mergeable_tracks);
		    TrackSet::iterator start = current.tracks.begin() + rnd.time;
		    TrackSet::iterator finish = current.tracks.begin() + rnd.number;
			
		    // forward transition probability
		    forward = merge_move_prob * (1.0 / (double) mergeable_tracks.size());
		    
		    // merge tracks into a single track
		    Track merged(start->size() + finish->size());
		    std::copy(start->begin(), start->end(), merged.begin());
		    std::copy(finish->begin(), finish->end(), merged.begin() + start->size());
	  
		    // remove tracks that were merged
		    current.tracks.erase(start);
		    current.tracks.erase(finish);
	  
		    // adds merged track
		    current.tracks.push_back(merged);
	  
		    // find tracks greater than 4 observatios
		    uint splitable_tracks = 0;
		    for (uint i = 0; i < current.tracks.size(); ++i)
				if (current.tracks[i].size() > 3)
					splitable_tracks++;;
	  
		    // backward transition probability
		    backward = movePDF(2, current) * (1.0 / (double) splitable_tracks);
		    backward *= 1.0 / ((double) merged.size() - 3);
		}
		else { // move can not be done
		    forward = 0.0;
		    backward = 0.0;
		}

		return current;
	}

	/**
	 * @brief Extends a randomly selected current
	 *
	 * @param current Current current
	 */ 
	Partition extensionMove(Partition current)
	{
		double extend_move_prob = movePDF(4, current);
	       
		// randomly select a track to extend
		TrackSet::iterator track_to_extend = sampleTrack(current.tracks);
		       
		// get track last observation time
		uint track_last_time = track_to_extend->back().time;
	     
		// extend track
		double track_probability = extendTrack(*track_to_extend, current.false_positives);
	       
		// calculate transition probability
		forward = extend_move_prob * (1.0 / (double)current.tracks.size()) * track_probability;
		backward = movePDF(5, current) * (1.0 / (double) current.tracks.size()) * (1.0 / (double)track_to_extend->size());
	       
		return current;
	}
     
	/**
	 * @brief Reduces a randomly selected current
	 *
	 * @param current Current current
	 */ 
	Partition reductionMove(Partition current)
	{
		double reduce_move_prob = movePDF(5, current);
	 
		// randomly select a track to reduce
		TrackSet::iterator track_to_reduce = sampleTrack(current.tracks);
		if (track_to_reduce->size() > 2){
		    // randomly select a time
		    uint reduction_time = sampleTime(1, track_to_reduce->size() - 2);
		    Track::iterator reduction_point = track_to_reduce->begin() + reduction_time + 1;			

		    // move track's observations after sampled time to false positives
		    for (Track::iterator it = reduction_point; it < track_to_reduce->end(); ++it)
				scans[it->time][it->number].tracked = false;
	      
		    // forward transition probabilities
		    forward = reduce_move_prob * (1.0 / (double)current.tracks.size());
		    forward *= 1.0 / (double) track_to_reduce->size();
	      
		    // compute backward probability
		    double track_probability = 1.0;
		    uint prev_time = (reduction_point - 1)->time;
		    uint prev_number = (reduction_point - 1)->number;
		    uint candidate_size = 0;
		    for (Track::iterator it = reduction_point; it < track_to_reduce->end(); ++it) {
				uint near_scan = (it->time - prev_time - 1);
				for (NeighborTime::iterator ob = nforest[prev_time][prev_number][near_scan].begin();
					 ob != nforest[prev_time][prev_number][near_scan].end(); ++ob){
					if (scans[ob->time][ob->number].tracked == false)
						candidate_size++;
				}
				track_probability *= (1.0 / (double)max_miss_detections) * (1.0 / (double) candidate_size);
				prev_time = it->time;
				prev_number = it->number;
		    }

		    // backward probability
		    backward = movePDF(4, current) * (1.0 / (double) current.tracks.size());
		    backward *= track_probability;

		    // reduce track
		    track_to_reduce->erase(reduction_point, track_to_reduce->end());
		}
		else{ // move is rejected
		    forward = 0.0;
		    backward = 0.0;
		}
	 
		return current;
	}
     
	/**
	 * @brief Updates a randomly selected track from the current partition
	 *
	 * @param current Current partition
	 *
	 * @return Proposed partition
	 */ 
	Partition updateMove(Partition current)
	{
		// select a track to update uniformly at random
		TrackSet::iterator track_to_update = sampleTrack(current.tracks);
	       
		// select update time uniformly a random from the selected track
		uint update_time = sampleTime(0, track_to_update->size() - 1) + 1;
		Track::iterator update_point = track_to_update->begin() + update_time;
	       
		// add selected track observations to empty track
		for (Track::iterator it = update_point; it != track_to_update->end(); ++it)
		    scans[it->time][it->number].tracked = false;

		// delete observations from selected time
		track_to_update->erase(update_point, track_to_update->end());
     
		// reassign observations to track from selected time
		double track_probability = extendTrack(*track_to_update, current.false_positives);
	       
		// track must have at least two observations
		if (track_to_update->size() >= 2){
		    // backward and forward transition probability are the same
		    forward = 1.0;
		    backward = 1.0;
		}
		else{
		    forward = 0.0;
		    backward = 0.0;
		}

		return current;
	}
     
	/**
	 * @brief Switches two tracks in the current partition
	 *
	 * @param current Current Partition
	 *
	 * @return Proposed partition
	 */ 
	Partition switchMove(Partition current)
	{
		// find pairs of tracks that can be switched
		std::vector<PairUI> switchable_tracks;
		std::vector<PairUI> switch_times;
		for (uint i = 0; i < current.tracks.size() - 1; ++i){
		    TrackSet::iterator track_i = current.tracks.begin() + i;
		    for (uint j = i + 1; j < current.tracks.size(); ++j){
				TrackSet::iterator track_j = current.tracks.begin() + j;
				for (uint p = 1; p < track_i->size() - 2; ++p){
					for (uint q = 1; q < track_j->size() - 2; ++q){
						if (track_i->at(p).time + max_miss_detections > track_j->at(q + 1).time){
							Item item_i_p = track_i->at(p);
							Item item_i_p1 = track_i->at(p + 1);
							Item item_j_q = track_i->at(q);
							Item item_j_q1 = track_i->at(q + 1);
							uint found1 = 0;
							for (int k = 0; k < max_miss_detections; ++k)
							{
								for (NeighborTime::iterator it = nforest[item_j_q.time][item_j_q.number][k].begin(); it != nforest[item_j_q.time][item_j_q.number][k].end(); ++it)
								{
									if (item_i_p1.time == it->time && item_i_p1.number == it->number)
									{
										found1 = 1;
										break;
									}
								}
								if (found1 = 1)
									break;
							}

							uint found2 = 0;
							for (int k = 0; k < max_miss_detections; ++k)
							{
								for (NeighborTime::iterator it = nforest[item_i_p.time][item_i_p.number][k].begin(); it != nforest[item_i_p.time][item_i_p.number][k].end(); ++it)
								{
									if (item_j_q1.time == it->time && item_j_q1.number == it->number)
									{
										found2 = 1;
										break;
									}
								}
				      
								if (found2 = 1)
									break;
							}

							if (found1 != 0 && found2 != 0)
							{
								switchable_tracks.push_back(std::make_pair(i,j));
								switch_times.push_back(std::make_pair(p,q));
							}
							if (track_j->at(q).time + max_miss_detections < track_i->at(p + 1).time)
								break;
						}
					}
				}
		    }
		}

		// if there are tracks that can be switched
		if (!switchable_tracks.empty()){
		    // randomly select a candidate pair to merge
		    PairUI sample_pair = sampleCandidate(switchable_tracks);
		    TrackSet::iterator start = current.tracks.begin() + sample_pair.first;
		    TrackSet::iterator finish = current.tracks.begin() + sample_pair.second;

		    // remove tracks that were merged
		    current.tracks.erase(start);
		    current.tracks.erase(finish);
	    
		    // backward and forward transition probability are the same
		    forward = 1.0;
		    backward = 1.0;
		}
		else { // move can not be done
		    forward = 0;
		    backward = 0;
		}

		return current;
	}
          
	ScanSet getScans(){
		return scans;
	}
	/**
	 * @brief Proposes a new current from a set of observations and the current partition
	 *
	 * @param current Current partition
	 */     
	Partition operator() (Partition& current, ScanSet&obs_db) 
	{
		scans = obs_db;
		std::random_device seed1{};
		generator_time.seed(seed1());
		Partition proposed_partition;

		// Select move
		//uint move = sampleMove(current); 
		uint move = 4;
		switch (move){
		case 0:
			proposed_partition = birthMove(current);
			break;
		case 1:
			proposed_partition = deathMove(current);
			break;     
		case 2:
			proposed_partition = splitMove(current);
			break;     
		case 3:
			proposed_partition = mergeMove(current);
			break;     
		case 4:
			proposed_partition = extensionMove(current);
			break;     
		case 5:
			proposed_partition = reductionMove(current);
			break;     
		case 6:
			proposed_partition = updateMove(current);
			break;     
		case 7:
			proposed_partition = switchMove(current);
			break;     
		}
	       
		// move is rejected
		if (forward != 0 && backward != 0)
		    return proposed_partition; 
		else 
			return current;
	}
private:
	uint last_time;
	uint max_miss_detections;
	ScanSet scans;
	NForest nforest;
	double forward;
	double backward;
	std::default_random_engine generator_time;     
	std::default_random_engine generator_track;
	std::default_random_engine generator_ending;
};
#endif

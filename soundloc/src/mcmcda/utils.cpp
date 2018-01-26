#include "utils.hpp"

/*
 * @brief Prints forest of neighboring observations
 *
 * @param nforest Forest to be printed
 */
void print_forest(NForest& nforest)
{
	std::cout << "Printing " << nforest.size() << " neighbor trees" << std::endl;
	for (NForest::iterator it = nforest.begin(); it != nforest.end(); ++it){
		std::cout << "t" << (int) (it - nforest.begin()) << " ================= " << std::endl;
		for (NeighborScan::iterator nh = it->begin(); nh != it->end(); ++nh){
			std::cout << "\to(" << (int) (nh - it->begin()) << ") = ";
			for (NeighborObservation::iterator nt = nh->begin(); nt != nh->end(); ++nt){
				std::cout << "d" << (int) (nt - nh->begin()) << " [";
				for (NeighborTime::iterator ns = nt->begin(); ns != nt->end(); ++ns)
					std::cout << "(" << ns->time << ", " << ns->number << ")";
				std::cout << "] ";
				if (nt + 1 != nh->end())
					std::cout << " -- ";
			}
			std::cout << std::endl;
		}
	}
}

/*
 * @brief Prints a set of tracks
 *
 * @param tracks Set of tracks to be printed
 */
void print_tracks(TrackSet& tracks, ScanSet& scans)
{
	std::cout << "Printing " << tracks.size() << " tracks" << std::endl;
	for (TrackSet::iterator it = tracks.begin(); it != tracks.end(); ++it){
		std::cout << (int) (it - tracks.begin()) << " ::: ";
		for (Track::iterator tr = it->begin(); tr != it->end(); ++tr)	
			std::cout << "[" << tr->time << ", " << tr->number << ", " 
					  << scans[tr->time][tr->number].value << ", " 
					  << scans[tr->time][tr->number].tracked << " ] ";
		std::cout << std::endl;
	}
}

/*
 * @brief Prints a set of scans 
 *
 * @param scans Set of scans to be printed
 */
void print_scans(ScanSet& scans)
{
	std::cout << "Printing " << scans.size() << " scans" << std::endl;
	for (ScanSet::iterator it = scans.begin(); it != scans.end(); ++it){
		std::cout << (int) (it - scans.begin()) << " ::: ";
		for (Scan::iterator sc = it->begin(); sc != it->end(); ++sc)	
			std::cout << "[" << sc->value << ", " << sc->tracked << "] ";
		std::cout << std::endl;
	}     
}

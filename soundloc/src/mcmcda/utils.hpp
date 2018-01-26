/**
 * @brief Classes and structures for Markov Chain Monte Carlo Data Association.
 * 
 * @author Gibran Fuentes Pineda <gibranfp@turing.iimas.unam.mx>
 * @date 2013
 */
#ifndef UTILS_HPP
#define UTILS_HPP
#include <iostream>
#include "partition.hpp"
#include "scans.hpp"

void print_forest(NForest&);
void print_tracks(TrackSet&, ScanSet&);
void print_scans(ScanSet&);
#endif

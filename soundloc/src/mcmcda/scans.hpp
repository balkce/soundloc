#ifndef SCANS_HPP
#define SCANS_HPP
#include <vector>

typedef struct Observation {
     double value;
     bool tracked;
} Observation;

typedef unsigned int uint;
typedef std::vector<Observation> Scan;
typedef std::vector<Scan> ScanSet;
#endif

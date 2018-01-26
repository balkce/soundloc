#ifndef PARTITION_HPP
#define PARTITION_HPP
#include <vector>

typedef unsigned int uint;
typedef std::pair<double,double> PairD;
typedef std::pair<uint,uint> PairUI;

typedef struct Item
{
     uint time;
     uint number;
} Item;

typedef std::vector<Item> Track;
typedef std::vector<Item>::iterator TrackIterator;
typedef std::vector<Track> TrackSet;
typedef std::vector<Track>::iterator TrackSetIterator;

typedef struct Partition
{
     TrackSet tracks;
     Track false_positives;
} Partition;

// create a typedef for the Graph type
typedef std::vector<Item> NeighborTime;
typedef std::vector<NeighborTime> NeighborObservation;
typedef std::vector<NeighborObservation> NeighborScan;
typedef std::vector<NeighborScan> NForest;
#endif

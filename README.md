# NEW PACKAGE UPDATE
There is an updated version of this package called [soundloc_kmeans](https://github.com/balkce/soundloc_kmeans) that is recommended.

# soundloc
Lightweight multiple sound source localization, based on a triangular microphone array.

Uses JACK for input/audio audio server.

Outputs localizations through the JsonMsg topic, described by the json_msgs node.

Configured by soundloc_config.yaml:
* distance between microphones
* number of windows of noise samples for initial voice activity detection calibration
* energy threshold to consider change of inactive and active states
* energy threshold to trigger activity
* coherence threshold (in degrees) between local microphone-pair estimations
* angle reversion (180 degree shift)
* graphical representation of localizations (requires PLplot, see dependencies)
* automatic connection to JACK inputs
* record inputs in WAV files
* choose between Cross-Correlation (1) or Phase Transform (2) for time delay estimation
* choose between clustering method (1), Markov Chain Monte Carlo Data Association (2), or Kalman filtering (3) for tracking methodology

## Dependencies
Packages that can be installed trough apt official repositories:
* libjack-jackd2-dev: JACK development libraries
* libfftw3-dev: a very fast FFT C/C++ implementation
* libsndfile1-dev, libsamplerate0-dev: for WAV file creation
* libplplot-c++11, libplplot-dev, plplot12-driver-xwin: for plotting results in real-time
* libjsoncpp-dev: for constructing JSON structures


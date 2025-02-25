# trieff 0.8.6

## Major Features

* Added segment bootstrapping method for more accurate confidence intervals with temporal data
  * This new approach preserves temporal patterns by sampling segments of continuous queue activity
  * Maintains the dependency structure between patients who were in the ED during the same time period
  * Produces theoretically more valid confidence intervals for queue-based data

* Enhanced bootstrap vignette with detailed comparison between standard and segment bootstrapping

# trieff 0.8.5

* Initial release

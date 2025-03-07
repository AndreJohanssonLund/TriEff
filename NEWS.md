# NEWS

## v 1.0.0 - 2025-02-28

* Increased the version number to signify that this is the version that is sent in to review togeather with the original article. No changes to code or API from previous version


## v 0.8.7 - 2025-02-25

* Reworked segmentation logic to reduce memory overhead and increase speed.


## v 0.8.6 - 2025-02-25

* Added segment bootstrapping method for more accurate confidence intervals with temporal data
  * This new approach preserves temporal patterns by sampling segments of continuous queue activity
  * Maintains the dependency structure between patients who were in the ED during the same time period
  * Produces theoretically more valid confidence intervals for queue-based data

* Enhanced bootstrap vignette with detailed comparison between standard and segment bootstrapping

## v 0.8.5 - 2025-02-08

* Initial release

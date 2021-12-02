# zalpha 0.3.0

* Improved the speed by editing the equal_vector function
* L_plus_R now returns 0 instead of NA when n<2 for choose(n,2) 
* Zalpha_all now always calculates the diversity statistics LR and L_plus_R
* If a pair of SNPs are further apart than the maximum bin in the LD profile, the values in the maximum bin will now be returned

# zalpha 0.2.0

* Can now accept data with missing values
* New function create_LDprofile for generating an LD profile
* Updated and improved vignette
* Fixed noLD issues

# zalpha 0.1.0

* First version on CRAN

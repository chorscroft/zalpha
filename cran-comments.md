## Changes from previous version

* Updated email address for maintainer as will lose access to old one soon
* Improved the speed by editing the equal_vector function
* L_plus_R now returns 0 instead of NA when n<2 for choose(n,2) 
* Zalpha_all now always calculates the diversity statistics LR and L_plus_R
* If a pair of SNPs are further apart than the maximum bin in the LD profile, the values in the maximum bin will now be returned

## Test environments

* Windows 10 x64
* Windows Server 2008 R2 SP1, R-devel, 32/64 bit
* Ubuntu Linux 20.04 LTS, R-release, GCC
* Fedora Linux, R-devel, clang, gfortran
* Debian Linux, R-release, GCC
* Debian Linux, R-devel, GCC, no long double

## Results

### NOTES:

   * Email address change:
   New maintainer:
      Clare Horscroft <chorscroft@aol.co.uk>
   Old maintainer(s):
      Clare Horscroft <c.horscroft@soton.ac.uk>

Email address has been changed as access to the old email address will be lost soon - email from the old address will be sent to verify change.

   * Found the following (possibly) invalid DOIs:
     DOI: 10.1534/genetics.115.185900
       From: DESCRIPTION
       Status: Forbidden
       Message: 403
       
DOI is correct, https://doi.org/10.1534/genetics.115.185900 goes to the correct reference. 


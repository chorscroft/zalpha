## Resubmission

This is a resubmission. In this version I have:

* added an immediate call to on.exit after changing options so as to not change the user's options outside of the function

Thanks to Martina Schmirl for the suggestion!

## Test environments
* Windows 10 x64
* Windows Server 2008 R2 SP1, R-devel, 32/64 bit
* Ubuntu Linux 16.04 LTS, R-release, GCC
* Fedora Linux, R-devel, clang, gfortran

## R CMD check results
There were no ERRORs or WARNINGs or NOTES. 

## rhub check results
   * Maintainer: 'Clare Horscroft <c.horscroft@soton.ac.uk>'
    New submission
   
   * Possibly mis-spelled words in DESCRIPTION:
     Kivisild (7:122)
     Sluckin (7:110)
   
   * Found the following (possibly) invalid DOIs:
     DOI: 10.1534/genetics.115.185900
       From: DESCRIPTION
       Status: Forbidden
       Message: 403
       
Mis-spelled words are author names. DOI is correct, https://doi.org/10.1534/genetics.115.185900 goes to the correct reference. 

## win_devel check results
 * Maintainer: 'Clare Horscroft     <c.horscroft@soton.ac.uk>'

    New submission

 * Possibly mis-spelled words in DESCRIPTION:
      Kivisild (7:122)
      Sluckin (7:110)

Mis-spelled words are author names.

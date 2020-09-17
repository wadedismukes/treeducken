## Resubmission

* Changed wording in DESCRIPTION to avoid spelling notes

* Altered example to increase runtime

* Fixed error in unit test causing it to fail occasionally

* Updated URL in README.md

* Fixed seg fault occuring rarely on testing

* Fixed runtime error found in GCC ASAN/UBSAN

## Test environments
* local macOS install, R 4.0.2
* local Debian 10, R 4.0.2
* macOS-latest (release) (on Github actions), 4.0.2
* ubuntu 20.04 (release) (on Github actions), R 4.0.2
* ubuntu 20.04 (devel) (on Github actions), R 4.0.2
* windows-latest (on Github actions), R 4.0.2
* win-builder (release)
* Windows Server 2008 R2 SP1, R-devel, 32/64 bit (using r-hub)
* Fedora Linux, R-devel, clang, gfortran) (using r-hub)
* Ubuntu Linux 16.04 LTS, R-release, GCC (using r-hub)
* Debian Linux, R-devel, GCC ASAN/UBSAN (using r-hub)

## R CMD check results
There were no ERRORs or WARNINGs.

There were 3 NOTES:

* checking for future file timestamps ... NOTE
  unable to verify current time

This appears to be an issue with the web service being used to perform this check not with the package.

* checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Wade Dismukes <wade.dismukes@gmail.com>'

This is a first time submission to CRAN for me (outside of the previous submission which was rejected).

* checking installed package size ... NOTE
    installed size is 12.6Mb
    sub-directories of 1Mb or more:
      libs  12.2Mb

This only occurs on Linux systems. Mac and Windows builds do not contain this note. The large file is the 
shared library (i.e. the .so file). 

## Downstream dependencies
There are currently no downstream dependencies for this package.
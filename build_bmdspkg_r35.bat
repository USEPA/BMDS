C:\Users\theam\Documents\RBMDS-master\RBMDS-master>more build_bmdspkg_r35.bat
@echo Building BMDS package for R 3.5.x
set path=C:\Rtools\bin;%path%
set topdir=c:/Users/theam/Documents/libs
echo %topdir%
R CMD build ToxicR
#R CMD INSTALL --build ToxicR_1.0.tar.gz
R CMD INSTALL ToxicR
#/c/Users/drxba/Documents
#xcopy BMDS_0.1.zip BMDS_0.1.r35.zip

@echo Building BMDS package for R 3.5.x
set path=C:\Rtools\bin;C:\Users\lolszyk\myApps\R\R-3.5.3\bin;%path%
set topdir=%cd%
R CMD build BMDS
R CMD INSTALL --build BMDS_0.1.tar.gz
xcopy BMDS_0.1.zip BMDS_0.1.r35.zip


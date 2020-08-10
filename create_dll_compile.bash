#!/bin/bash
#Script to copy over the correct R files and then
#build the correct make scripts using configure. 

####################################################
#ALWAYS NEEDS A PATH TO THE EIGEN DIRECTORY
#BY DEFAULT IT ASSUMES IT IS IN /usr/local/include/eigen3
SCRIPT_EIGEN_INCLUDE=/usr/local/include/eigen3


##################################################
echo -e "\e[1;31m ------------------------------------------------\e[0m"
echo -e "\e[1;31m ------------------------------------------------\e[0m"
echo -e "\e[1;31m Assuming Eigen linear algebra library located at \e[0m"
echo -e "\e[1;32m $SCRIPT_EIGEN_INCLUDE  \e[0m"
echo -e "\e[1;31m ------------------------------------------------\e[0m"
echo -e "\e[1;31m ------------------------------------------------\e[0m"
##################################################
echo -e "\e[1;31m ------------------------------------------------\e[0m"
echo -e "\e[1;31m ------------------------------------------------\e[0m"
echo -e "\e[1;31m Warning removing all information in the directory \e[0m"
echo -e "\e[1;32m BMD_DLL_COMPILE  \e[0m"
echo -e "\e[1;31m ------------------------------------------------\e[0m"
echo -e "\e[1;31m ------------------------------------------------\e[0m"

rm -R ../BMD_DLL_COMPILE/*
rmdir ../BMD_DLL_COMPILE
##############################################################
mkdir ../BMD_DLL_COMPILE
mkdir ../BMD_DLL_COMPILE/include
mkdir ../BMD_DLL_COMPILE/code_base
cp -a ./ToxicR/src/include/* ../BMD_DLL_COMPILE/include 
cp -a ./ToxicR/src/code_base/* ../BMD_DLL_COMPILE/code_base
cp ./configure.ac ../BMD_DLL_COMPILE
cp ./Makefile.am  ../BMD_DLL_COMPILE
cp ./version.c    ../BMD_DLL_COMPILE

cd ../BMD_DLL_COMPILE
libtoolize
aclocal
autoconf
automake --add-missing
./configure EIGEN_INCLUDE=$SCRIPT_EIGEN_INCLUDE

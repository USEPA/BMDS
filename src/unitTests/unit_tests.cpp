//unit_tests.cpp
#include "assert.h"
#include <vector>
#include "bmds_helper.h"
#include "unit_tests.h"

int run_all_unitTests(){

	std::cout<<"Running unit tests"<<std::endl;
	objfunc_test();
	return 0;
}

void objfunc_test(){
	std::vector<double> x{1.5, 2.0, 3.2};
        std::vector<double> tmp;
        //assert(objfunc_bmdl(x, tmp, NULL)==1.5);
        expect_true(objfunc_bmdl(x, tmp, NULL)==1.5);
}

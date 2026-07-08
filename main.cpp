/**
*/

#include "General.h"
#include "Obj_Manager.h"
#include "LBM_Manager.hpp"
#include "Lat_Manager.hpp"

using namespace std;
ofstream* Log_function::logfile = nullptr;
Obj_Manager* Obj_Manager::pointer_me;
LBM_Manager* LBM_Manager::pointer_me;

int main()
{
	Timer tmr;
	ensureOutputDir();
	// Create application log file
	std::ofstream logfile;
	logfile.open(outputPath("log"), ios::out);
	Log_function::logfile = &logfile;
	if (!logfile.is_open()) {
		std::cout << "Can't open logfile" << std::endl;
	};

	Obj_Manager obj_manager;
	Obj_Manager::pointer_me = &obj_manager;

	double t0 = tmr.elapsed();

	obj_manager.initial();

	obj_manager.output();

	double t1 = tmr.elapsed();
	std::cout << "Time for initial mesh: " << t1 - t0 << std::endl;

	LBM_Manager lbm_manager;
	lbm_manager.fluidSimulate();

	double t2 = tmr.elapsed();
	std::cout << "Time for fluid simulation: " << t2 - t1 << std::endl;

	// obj_manager.time_marching_management();

	// obj_manager.output();

	logfile.close();
	return 0;
}

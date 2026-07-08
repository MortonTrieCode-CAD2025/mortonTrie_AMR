/**
* @file
* @brief Include stardard libraries, headers and functions.
* @note 1) include libraries; 2) define macros; 3) define inline funciotns.
*/

#pragma once

// libraries
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <array>
#include <cstdlib>
#include <cmath>
#include <bitset>
#include <stdexcept>

#define DEPRECATED [[deprecated("This function is deprecated. Use a newer version or alternative.")]]

// All runtime-generated files (log, Morton dump, Tecplot mesh, flow-field snapshots)
// are written under this directory, keeping the source tree clean.
constexpr const char* C_OUTPUT_DIR = "output";

typedef double D_real;
typedef int D_int;
typedef unsigned long long D_uint;
typedef double D_Phy_real;

typedef double D_Phy_DDF;
typedef double D_Phy_Velocity;
typedef double D_Phy_Rho;

#define C_BIT 64                               ///< Bits to store morton code 
using D_morton = std::bitset<C_BIT>;         // Bitset to store morton code, due to the method used to find neighbours in Morton_Assist.h, the bit is limited to unsigned long long

// Geometry parameter
// #define NONSLIP  0
// #define MOVING   1
// #define SOLID    2
// #define EXTRAP   3
// #define FLUID    4
// #define BOUNDARY 5

// detect system to choose approriate libraries for mikdir
#ifdef _WIN32
#include <io.h>
#include <direct.h>
#include <Windows.h>
#elif __linux__
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#endif

inline void ensureOutputDir()
{
#ifdef _WIN32
    _mkdir(C_OUTPUT_DIR);
#else
    mkdir(C_OUTPUT_DIR, 0777);
#endif
}

inline std::string outputPath(const std::string& filename)
{
    return std::string(C_OUTPUT_DIR) + "/" + filename;
}

// timer
// #ifdef __unix__
// #include <ctime>
// class Timer
// {
// public:
// 	Timer() { clock_gettime(CLOCK_REALTIME, &beg_); }

// 	double elapsed() {
// 		clock_gettime(CLOCK_REALTIME, &end_);
// 		return end_.tv_sec - beg_.tv_sec +
// 			(end_.tv_nsec - beg_.tv_nsec) / 1000000000.;
// 	}

// 	void reset() { clock_gettime(CLOCK_REALTIME, &beg_); }

// private:
// 	timespec beg_, end_;
// };
// #else
#include <chrono>
class Timer
{
public:
	Timer() : beg_(clock_::now()) {}
	void reset() { beg_ = clock_::now(); }
	double elapsed() const {
		return std::chrono::duration_cast<second_>
			(clock_::now() - beg_).count();
	}

private:
	typedef std::chrono::high_resolution_clock clock_;
	typedef std::chrono::duration<double, std::ratio<1> > second_;
	std::chrono::time_point<clock_> beg_;
};
// #endif

// headers
#include "Constants.h"

// define expressions
#define SQ(x) ((x) * (x)) ///< Square
#define CU(x) ((x) * (x) * (x))
#define SQRT3 1.7320508075688772

// Lattice constants
constexpr D_real L_C = 1.;
constexpr D_real L_Cs = 1./SQRT3;     ///< Lattice sound speed
constexpr D_real L_CsSq = 1.0/3.0;    ///< Square lattice sound speed
constexpr D_real L_Rho = 1.;

constexpr D_real C_eps = 1e-10;            ///< Small number for comparing floats to zero
const D_real C_pi = 4.*atan(1.);       ///< Pi

/**
* @brief The structure store functions to store log file
*/
struct Log_function
{
public:
	static std::ofstream* logfile;			///< pointer point to log file
};

/**
* @brief function to write normal information to the logfile.
* @param[in]  msg        information write to the logfile.
* @param[in]  logfile   pointer to the logfile address.
*/
inline void log_infor(const std::string &msg, std::ofstream *logfile)
{

	// *logfile << "Infor: " << msg << std::endl;
	// std::cout << "Infor:" << msg << std::endl;
	*logfile << msg << std::endl;
};

/**
* @brief function to write warning information to the logfile.
* @param[in]  msg        information write to the logfile.
* @param[in]  logfile   pointer to the logfile address.
*/
inline  void log_warning(const std::string &msg, std::ofstream *logfile)
{

	*logfile << "Warning: " << msg << std::endl;
	// SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), FOREGROUND_INTENSITY | FOREGROUND_RED | FOREGROUND_GREEN);
	std::cout << "Warning:" << msg << std::endl;
	// SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), FOREGROUND_INTENSITY | FOREGROUND_RED | FOREGROUND_GREEN | FOREGROUND_BLUE);
};

/**
* @brief function to write error information to the logfile.
* @param[in]  msg        information write to the logfile.
* @param[in]  logfile   pointer to the logfile address.
*/
inline  void log_error(const std::string &msg, std::ofstream *logfile)
{

	*logfile << "Error: " << msg << std::endl;
	// SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), FOREGROUND_INTENSITY | FOREGROUND_RED);
	std::cout << "Error:" << msg << std::endl;
	// SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), FOREGROUND_INTENSITY | FOREGROUND_RED | FOREGROUND_GREEN | FOREGROUND_BLUE);
	exit(0);
};

/**
* @brief function to calculate 2 power N.
* @param[in]  n        power exponent.
*/
inline D_uint two_power_n(D_uint n)
{
	D_uint two = 1;
	two  <<= n;
	return two;
};

template<typename T>
inline T min_of_two(T a, T b)
{
	return (a < b) ? a : b;
}

template<typename T>
inline T max_of_two(T a, T b)
{
	return (a > b) ? a : b;
}

template<typename T>
inline T min_of_three(T a, T b, T c)
{
	T temp = (b < c) ? b : c;
	return a < temp ? a : temp;
}

template<typename T>
inline T max_of_three(T a, T b, T c)
{
	T temp = (b > c) ? b : c;
	return a > temp ? a : temp;
}


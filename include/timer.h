#pragma once

#include <map>
#include <chrono>

namespace dimp3{;

class Timer{
public:
	std::chrono::high_resolution_clock::time_point  last;
	
public:
	/// count time in [us]
	int CountUS();

	/// sleep [ms]
	static void Sleep(int ms);
};

}

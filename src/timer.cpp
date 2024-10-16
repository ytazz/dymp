#include <dymp/timer.h>

#ifdef _WIN32
# include <windows.h>
#else
# include <sys/time.h>
# include <unistd.h>
#endif

namespace dymp{;

void Timer::Sleep(int ms){
#ifdef _WIN32
	::Sleep(ms);
#else
	usleep(ms*1000);
#endif
}

int Timer::CountUS(){
	auto t = std::chrono::high_resolution_clock::now();
	auto duration = t - last;
	int us = std::chrono::duration_cast<std::chrono::microseconds>(duration).count();

	last = t;

	return us;
}

}

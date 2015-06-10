#pragma once

#include <fstream>
#include <iostream>
#include <sstream>
#include "util/timer.h"
#include <cmath>

enum ProgressCounterStyle {
    ETA,
    SIMPLE
};

static const std::string clear = "\r" + std::string(80,' ') + "\r";

class ProgressCounter {
    private:
        std::ofstream tty;
        util::WallTimer timer;
        std::string task;
        std::size_t max;
        std::size_t count;

    public:
        ProgressCounter(std::string const & _task, std::size_t max);
        template <ProgressCounterStyle T> void progress(void);
        void inc(void);
        void reset(std::string const & _task);
};

inline
ProgressCounter::ProgressCounter(std::string const & _task, std::size_t _max)
    : tty("/dev/tty", std::ios_base::out), timer(),
    task(_task), max(_max), count(0) {}

inline void
ProgressCounter::inc(void) {
    #pragma omp atomic
    ++count;

    if(count >= max) {
        std::stringstream ss;
        ss << clear << task << " 100%... done. (Took "
            << timer.get_elapsed_sec() << "s)";
        #pragma omp critical
        std::cout << ss.rdbuf() << std::endl;
    }
}

inline void
ProgressCounter::reset(std::string const & _task) {
    timer.reset();
    count = 0;
    task = _task;
}

template <ProgressCounterStyle T> void
ProgressCounter::progress(void) {
    if (max <= 100 || (max > 100 && count % (max / 100) == 0)) {
        float percent = static_cast<float>(count) / max;
        int ipercent = std::floor(percent * 100.0f + 0.5f);


        std::stringstream ss;
        ss << clear << task << " " << ipercent << "%...";

        if (T == ETA && ipercent > 3){
            std::size_t const elapsed = timer.get_elapsed();
            std::size_t eta = (elapsed / percent - elapsed) / 1000;
            ss << " eta ~ " << eta << " s";
        }

        #pragma omp critical
        tty << ss.rdbuf() << std::flush;
    }
}

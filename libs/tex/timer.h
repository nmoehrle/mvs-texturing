/*
 * Copyright (C) 2015, Nils Moehrle
 * TU Darmstadt - Graphics, Capture and Massively Parallel Computing
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the BSD 3-Clause license. See the LICENSE.txt file for details.
 */

#ifndef TEX_TIMER_HEADER
#define TEX_TIMER_HEADER

#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <cerrno>
#include <cstring>

#include "util/timer.h"
#include "util/exception.h"

class Timer {
    private:
        util::ClockTimer ctimer;
        util::WallTimer wtimer;

        std::string header;

        struct Event{
            std::string name;
            std::size_t abs_clocks;
            std::size_t abs_milliseconds;
            std::size_t rel_clocks;
            std::size_t rel_milliseconds;
        };

        std::vector<Event> events;

    public:
        Timer(void);
        Timer(std::string const & _header);

        void measure(std::string const & eventname);
        void write_to_file(std::string const & filename) const;
};

#endif /* TEX_TIMER_HEADER */

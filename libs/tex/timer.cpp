/*
 * Copyright (C) 2015, Nils Moehrle
 * TU Darmstadt - Graphics, Capture and Massively Parallel Computing
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the BSD 3-Clause license. See the LICENSE.txt file for details.
 */

#include "timer.h"

Timer::Timer(void) {
    header = "";
    ctimer = util::ClockTimer();
    wtimer = util::WallTimer();
}

Timer::Timer(std::string const & _header) {
    header = _header;
    ctimer = util::ClockTimer();
    wtimer = util::WallTimer();
}

void Timer::measure(std::string const & eventname) {
    std::size_t abs_clocks = ctimer.get_elapsed();
    std::size_t abs_milliseconds = wtimer.get_elapsed();

    if(events.empty()) {
        Event event = {eventname, abs_clocks, abs_milliseconds, abs_clocks, abs_milliseconds};
        events.push_back(event);
    } else {
        Event last_event = events.back();
        std::size_t rel_clocks = abs_clocks - last_event.abs_clocks;
        std::size_t rel_milliseconds = abs_milliseconds - last_event.abs_milliseconds;

        Event event = {eventname, abs_clocks, abs_milliseconds, rel_clocks, rel_milliseconds};
        events.push_back(event);
    }
}

void Timer::write_to_file(std::string const & filename) const {
    std::ofstream out(filename.c_str());
    if (!out.good())
        throw util::FileException(filename, std::strerror(errno));

    if(!header.empty())
        out << "#" << header << std::endl;

    out << "Event, Absolute clocks, Absolute milliseconds, Relative clocks, Relative milliseconds"
        << std::endl;
    for (std::size_t i = 0; i < events.size(); ++i) {
        Event event = events[i];
        out << event.name << ", "
            << event.abs_clocks << ", "
            << event.abs_milliseconds << ", "
            << event.rel_clocks << ", "
            << event.rel_milliseconds << std::endl;
    }
    out.close();
}

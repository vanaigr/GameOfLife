#pragma once
#include<string>
#include<iostream>
#include"Timer.h"

template<class Units = std::chrono::microseconds>
class AutoTimer
{
private:
    Timer<Units> timer;
    const std::string name;
    std::ostream& stream;

public:
    AutoTimer(const std::string &&name_, std::ostream &stream_ = std::cout) : timer(), name(name_), stream(stream_) {}
    ~AutoTimer() {
        stream << name << ':' << timer.elapsedTime() << std::endl;
    }

    AutoTimer(const AutoTimer&) = delete;
    AutoTimer& operator=(const AutoTimer&) = delete;
};


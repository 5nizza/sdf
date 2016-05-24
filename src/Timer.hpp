//
// Created by ayrat on 24/05/16.
//

#ifndef SDF_TIMER_H
#define SDF_TIMER_H


#include <ctime>


class Timer {
public:
    Timer() {
        last = std::clock();
    }

    /* (restarts the counter)
       returns the number of seconds since the last call to the function (or since the creation if just created).
    */
    clock_t sec_restart() {
        clock_t elapsed = (std::clock() - last) / CLOCKS_PER_SEC;
        last = clock();
        return elapsed;
    }

private:
    std::clock_t last;
};


#endif //SDF_TIMER_H

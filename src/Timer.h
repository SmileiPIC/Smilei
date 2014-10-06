#ifndef TIMER_H
#define TIMER_H

#include <string>

class SmileiMPI;

class Timer {
public:
    Timer();
    ~Timer();
    void init(std::string name);
    void init(SmileiMPI *smpi, std::string name);
    void update();
    void restart();
    double getTime(){return time_acc_;}
    void print(double tot);
private:
    double last_start_;
    double time_acc_;
    std::string name_;
    SmileiMPI* smpi_;
};

#endif


#ifndef GRID_STATUS_H_
#define GRID_STATUS_H_

#include <string>

using namespace std;

enum State
{
    OK,
    LOAD_ERROR,
    FILE_OPEN_ERROR,
    MEMORY_ALLOC_ERROR,
    GRID_GENERATE_ERROR,
    UNKNOWN_ERROR
};

struct Status
{
    State state = State::OK;
    string msg = "";
};

class GridStatus
{

private:
    Status status;

public:
    GridStatus() = default;
    
    GridStatus(const Status &status_) : status(status_) {}
    
    inline Status GetStatus() { return status; }
    inline string GetMsg() { return status.msg; }
    inline State GetState() { return status.state; }
    inline void SetStatus(const State& state, const string &msg) { status = {state, msg}; }

    ~GridStatus() = default;
};

#endif
#include <iostream>
#include <boost/progress.hpp>
#include "TaskManager.h"

using namespace boost;

int main() {
    progress_timer timer;
    srand(time(NULL));
    TaskManager taskManager;
    taskManager.execute();
    return 0;
}
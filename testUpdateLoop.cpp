#include "disc.hpp"

int main()
{
    vec v1 = {10,0,0};
    std::vector<vec> pMap = {{10,200,0},{400,200,0},{400 - DISC_RADIUS,200 + sqrt(3)*DISC_RADIUS,0},{400 - DISC_RADIUS,200 - sqrt(3)*DISC_RADIUS,0}};
    StateQ p;
    p.initialiseQueue(v1,NO_DISCS,pMap);
    p.computationThread.join();
    std::cout << "final debug!!!!!!!!!11\n";
    while(p.q.size() > 0)
    {
        p.q.front().printDebug();
        p.q.pop();
    }
    std::cout << std::endl;
    return 0;
}

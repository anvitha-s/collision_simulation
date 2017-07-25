#include "disc.hpp"

int main()
{
    double origin = 400;
    vec v1 = {10,0,0};
    vec d1 = {10,200,0};
    vec v2 = {0,0,0};
    //vec d2 = {400,200,0};
    vec d2 = {origin - DISC_RADIUS,200 + DISC_RADIUS*sqrt(3) ,0};
    //vec d3 = {420,200,0};
    vec d3 = {origin - DISC_RADIUS,200 - DISC_RADIUS*sqrt(3) ,0};
    vec d4 = {origin, 200,0};

    vec d5 = {origin + DISC_RADIUS,200 - DISC_RADIUS*sqrt(3) ,0};
    vec d6 = {origin + DISC_RADIUS,200 + DISC_RADIUS*sqrt(3) ,0};
    vec d7 = {origin + 2*DISC_RADIUS,200,0};
    vec d8 = {origin - 2*DISC_RADIUS,200,0};

    std::vector<vec> position_map = {d1,d4,d3,d2,d8,d5,d6,d7};
    StateQ p;
    p.initialiseQueue(v1,NO_DISCS,position_map);
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

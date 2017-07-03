#include "disc.hpp"

int main()
{
    vec v1 = {10,0,0};
    std::vector<vec> pMap = {{500,200,0},{400,200,0}};
    StateQ p;
    p.initialiseQueue(v1,NO_DISCS,pMap);
    p.computationThread.join();
    std::cout << "final debug!!!!!!!!!11\n";
    while(p.q.size() > 0)
    {
        p.q.front().printDebug();
        std::cout << "contact pairs : \n";
        for(int m = 0;m < NO_DISCS;m++)
        {
            std::cout << "\n" << m << "=> ";
            for(auto n:p.q.front().discs[m].contact_pairs)
                std::cout << n << ", ";
        }
        std::cout << "STOPPING TIME : " << p.q.front().stopping_time << std::endl;
        
        p.q.pop();
    }
    std::cout << std::endl;
    return 0;
}

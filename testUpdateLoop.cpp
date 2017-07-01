#include "disc.hpp"

int main()
{
    vec v1 = {1,0,0};
    std::vector<vec> pMap = {{10,200,0},{200,200,0},{190,200 + 10*sqrt(3),0},{190,200 - 10*sqrt(3),0}};
    StateQ p;
    p.initialiseQueue(v1,4,pMap);
    p.computationThread.join();
    std::cout << "final debug!!!!!!!!!11\n";
    while(p.q.size() > 0)
    {
        p.q.front().printDebug();
        std::cout << "contact pairs : \n";
        for(int m = 0;m < 4;m++)
        {
            std::cout << "\n" << m << "=> ";
            for(auto n:p.q.front().discs[m].contact_pairs)
                std::cout << n << ", ";
        }
        
        p.q.pop();
    }
    std::cout << std::endl;
    return 0;
}

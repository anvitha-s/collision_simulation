#include "disc.hpp"

int main()
{
    State a(3);
    vec v1 = {1,0,0},v2 = {0,0,0},v3 = {0,0,0}, p1 = {100,200,0}, p2 = {200,210,0},p3 = {200,190,0};
    a.discs[0].velocity = v1;
    a.discs[1].velocity = v2;
    a.discs[2].velocity = v3;
    a.discs[0].position = p1;
    a.discs[1].position = p2;
    a.discs[2].position = p3;
    int i = 0,j = 1;
    double b = a.predicted_collision.collision_instance = a.PICgenerator(i,j);
    std::cout << a.predicted_collision.collision_instance << "\n";
    j = 2;
    double d = a.predicted_collision.collision_instance = a.PICgenerator(i,j);
    std::cout << a.predicted_collision.collision_instance << "\n";
    if(b==d)
        std::cout << "equal\n";
/*    a.printDebug();
    a.impactUpdate(i,j);
    a.printDebug();
    a.velocityUpdate(); 
    a.printDebug();*/
    return 0;
}

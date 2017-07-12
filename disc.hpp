#ifndef _DISC_HPP_
#define _DISC_HPP_

#include <iostream>
#include <armadillo>
#include <yaml-cpp/yaml.h>
#include <algorithm>
#include <queue>
#include <thread>

using namespace std;
using namespace arma;

#define INIT_STATE -2
#define DISC_MASS 1
#define DISC_RADIUS 10
#define NO_DISCS 4

class Disc
{
    public:
    double radius,I,mass;
    vec velocity,angular_velocity,position;
    vector<int> contact_pairs;
    Disc():velocity({0,0,0}),angular_velocity({0,0,0}),mass(DISC_MASS),radius(DISC_RADIUS) {
    std::cout << "default constr\n";
        I = (mass*pow(radius,2))*0.5;
    }
    Disc(YAML::Node node_)
    {
        radius = node_["radius"].as<double>();
        mass = node_["mass"].as<double>();
        I = (mass*pow(radius,2))/(4*(node_["alpha"].as<double>()));
        YAML::Node node1_ = node_["linear_velocity"];
        velocity << node1_[0].as<double>() << node1_[1].as<double>() << 0;
        node1_ = node_["angular_velocity"];
        angular_velocity << 0 << 0 <<node1_.as<double>();
        
    }
};

struct PIC
{
    public:
    map<int,std::vector<int>> disc_indices;
    double collision_instance;
    
    void addPIC(int b_,std::vector<int> a_)
    {
        disc_indices[b_] = a_;
    }
    PIC(){}
    PIC(double t,int a, int b):collision_instance(t) 
    {
        disc_indices[a] = {b};
    };
};

class State
{
    public:
    int no_discs;
    PIC predicted_collision;
    double stopping_time;
    std::vector<Disc> discs;
    std::map<int,std::vector<vec>> impulses;
    vector<vector<int>> stackGenerator(std::pair<int,vector<int>> p);

    State(){};
    State(int noDiscs_):predicted_collision(-2,-1,-1),no_discs(noDiscs_)
    {
        std::cout << "resizing\n";
        discs.resize(no_discs);
    }
    PIC PICgenerator(const int& d1);
    double PICgenerator(const int& d1,const int& d2);
    void updateState();
    vector<vec> findPOC(const int& d1,const int& d2);
    void impactUpdate(const int& d1,const int& d2);
    void inContact(const int& d1,const int& d2);
    void velocityUpdate(int d1 = -1);
    void collisionUpdate();
    void removeContact(const int& d1,const int& d2);
    void computeStoppingTime(std::vector<int> activeIndices_);
//    void updateVelocities(const int& d1,const int& d2);
    void printDebug();
    void timeUpdate();
};

class StateQ
{
    public:
    std::queue<State> q;
    void updateQueue();
    void initialiseQueue(vec& initVel_,int noDiscs_,std::vector<vec> positionMap_);
    std::thread computationThread;
};
/* 
bool solveQuadratic(double a, double b, double c, double &root);
bool solveCubic(double a, double b, double c, double d, double &root);
bool solveQuartic(double a, double b, double c, double d, double e, double &root);
*/
#endif //_DISC_HPP_

#include <iostream>
#include <armadillo>
#include <yaml-cpp/yaml.h>
#include "solvePoly.hpp"

using namespace std;
using namespace arma;


class Disc
{
    public:
    double radius,I,mass;
    vec velocity,angular_velocity,position;
    vector<int> contact_pairs;
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
    
    void addPIC(std::vector<int> a_)
    {
        disc_indices.push_back(a_);
    } // or std::vector??
    PIC(double t,int a, int b):collision_instance(t),disc_indices({a,{b}}) {};
};

class State
{
    public:
    PIC predicted_collision;
    std::vector<Disc> discs;
    std::map<int,std::vector<double>> impulses;
    vector<vector<int>> stackGenerator(std::pair<int,vector<int>>& p);

    PIC PICgenerator(const int& d1);
    double PICgenerator(const int& d1,const int& d2);
    void updateState();
    vec findPOC(const int& d1,const int& d2);
    void impactUpdate(const int& d1,const int& d2);
    void inContact(const int& d1,const int& d2);
};

#include <iostream>
#include <armadillo>
#include "disc.hpp"
#include <cmath>
#include <yaml-cpp/yaml.h>

using namespace std;
using namespace arma;

vec State::findPOC(const int& d1,const int& d2)
{
   vec p2 = discs[d2].position + collision_instance*discs[d2].velocity + 0.5*n*g*collision_instance*collision_*normalise(discs[d2].velocity);
   vec p1 = discs[d1].position + collision_instance*discs[d1].velocity + 0.5*n*g*collision_instance*collision_*normalise(discs[d1].velocity);
   vec c1c2 = p2 - p1;
   return normalise(c1c2);
}

void State::impactUpdate(const int& d1,const int& d2)
{
    //position vector of 2 wrt 1
    vec k,j;
    k = findPOC(d1,d2);
    double phi = (3*n*(1+e))/(b+1);
    vec G = discs[d1].radius*cross(discs[d1].angular_velocity,k) + discs[d1].velocity - discs[d2].velocity - discs[d2].radius*cross(discs[d2].angular_velocity,k); 
    j = cross(cross(G,k),k);
    j = normalise(j);
    double A;
    vec diff_v = discs[d1].velocity - discs[d2].velocity;
    std::cout << diff_v.n_rows << diff_v.n_cols << " " << k.n_rows << k.n_cols <<std::endl;
    std::cout << as_scalar((diff_v.t())*k) << std::endl;
    std::cout << diff_v.n_rows << diff_v.n_cols << " " << k.n_rows << k.n_cols <<std::endl;
    A = as_scalar((diff_v.t())*k);//*(d1.mass*d2.mass/(d1.mass + d2.mass))*(-1)*(1+e);
    double B;
    //finding angle between G and k.
    double theta = acos(as_scalar((k.t())*G));
    if(theta > phi) //slipping collision
        B = (-1)*n*A;
    else            //sticking collision
        B = -(b+1)*(as_scalar((G.t())*j))/((1/discs[d1].mass) + (1/discs[d2].mass) + (pow(discs[d1].radius,2)/discs[d1].I) + (pow(discs[d2].radius,2)/discs[d2].I)); 
    vec J = A*k + B*j;
    discs[d1].velocity += J/discs[d1].mass;
    std::cout << k.size() << std::endl; 
    discs[d1].angular_velocity += cross(k,J)*(discs[d1].radius/discs[d1].I);
    discs[d2].velocity -= J/discs[d2].mass;
    discs[d2].angular_velocity += cross(k,J)*(discs[d2].radius/discs[d2].I);
}

int main(int argc, char **argv)
{
    YAML::Node node_ = YAML::LoadFile("/home/bs221b/Documents/summer_project/data.yaml");
    for(YAML::const_iterator it = node_.begin(); it != node_.end(); ++it)
    {
        std::cout << "i\n";
        std::cout << it->first.as<std::string>();
    }
    YAML::Node node1 = node_["disc_1"];
    YAML::Node node2 = node_["disc_2"];
    double e = node_["e"].as<double>(); //coefficient of restitution along normal
    double n = node_["n"].as<double>(); //coefficient of friction
    double b = node_["b"].as<double>(); //tangential coefficient of restitution
    YAML::Node node1 = node_["position2_1"];
    vec position2_1; 
    position2_1 << node1[0].as<double>() << node1[1].as<double>() << 0;
    //position vector of 2 wrt 1
    k = findPOC(position2_1,d1,d2);
    double phi = (3*n*(1+e))/(b+1);
    Disc d1(node1),d2(node2);
    vec k,j;
    std::cout << "input k : ";
    double i1,i2;
    std::cin >> i1 >> i2;
    //normalising to a unit vector.
    i1 = i1/(sqrt(pow(i1,2) + pow(i2,2)));
    i2 = i2/(sqrt(pow(i1,2) + pow(i2,2)));
    k << i1 << i2 << 0;
    vec G = d1.radius*cross(d1.angular_velocity,k) + d1.velocity - d2.velocity - d2.radius*cross(d2.angular_velocity,k); 
    j = (cross(G,k),k);
    i1 = j(0,0)/(sqrt(pow(j(0,0),2) + pow(j(1,0),2)));
    i2 = j(1,0)/(sqrt(pow(j(0,0),2) + pow(j(1,0),2)));
    j << i1 << i2 << 0;
    double A;
    vec diff_v = d1.velocity - d2.velocity;
    std::cout << diff_v.n_rows << diff_v.n_cols << " " << k.n_rows << k.n_cols <<std::endl;
    std::cout << as_scalar((diff_v.t())*k) << std::endl;
    std::cout << diff_v.n_rows << diff_v.n_cols << " " << k.n_rows << k.n_cols <<std::endl;
    A = as_scalar((diff_v.t())*k);//*(d1.mass*d2.mass/(d1.mass + d2.mass))*(-1)*(1+e);
    double B;
    //finding angle between G and k.
    double theta = acos(as_scalar((k.t())*G));
    if(theta > phi) //slipping collision
        B = (-1)*n*A;
    else            //sticking collision
        B = -(b+1)*(as_scalar((G.t())*j))/((1/d1.mass) + (1/d2.mass) + (pow(d1.radius,2)/d1.I) + (pow(d2.radius,2)/d2.I)); 
    vec J = A*k + B*j;
    d1.velocity += J/d1.mass;
    std::cout << k.size() << std::endl; 
    d1.angular_velocity += cross(k,J)*(d1.radius/d1.I);
    d2.velocity -= J/d2.mass;
    d2.angular_velocity += cross(k,J)*(d2.radius/d2.I);
    std::cout << "disc 1 :\n linear velocity : " << d1.velocity(0,0) << "," <<d1.velocity(1,0)<<  "," <<d1.velocity(2,0) << std::endl ;
    std::cout << "angular velocity : " << d1.angular_velocity(0,0) << "," <<d1.angular_velocity(1,0)<<  "," <<d1.angular_velocity(2,0) << std::endl ;
    std::cout << "disc 2 :\n linear velocity : " << d2.velocity(0,0) << "," <<d2.velocity(1,0)<<  "," <<d2.velocity(2,0) << std::endl;
    std::cout << "angular velocity : " << d2.angular_velocity(0,0) << "," <<d2.angular_velocity(1,0)<< "," << d2.angular_velocity(2,0) << std::endl;
}

#include "disc.hpp"
#include "solvePoly.hpp"
//coeff of friction
double n = 0.1;
//acceleration due to gravity
double g = 1;
//coeff of restitution in the normal direction
double e = 1;
//coeff of restitution in the tangential direction
double b = -1;
double board_length = 650;

map<int,vec> boundaries = {{-1,{-1,0,0}},{-2,{0,1,0}},{-3, {1,0,0}},{-4, {0,-1,0}}};

vector<vec> State::findPOC(const int& d1,const int& d2)
{
    vec k,G;
    if(d2 < 0)
    {
        k = boundaries[d2];
        G = discs[d1].radius*cross(discs[d1].angular_velocity,k) + discs[d1].velocity;
    }
    else
    {    
        vec c1c2 = discs[d2].position - discs[d1].position;
        k = normalise(c1c2);
        G = discs[d1].radius*cross(discs[d1].angular_velocity,k) + discs[d1].velocity - discs[d2].velocity - discs[d2].radius*cross(discs[d2].angular_velocity,k);
    }
    return {k,G};
}

bool State::isColliding(const int& d1, const int& d2)
{
    vec v1v2 = discs[d1].velocity - discs[d2].velocity;
    vec c12 = discs[d2].position - discs[d1].position;
    std::cout << "v1v2 : " << v1v2;
    std::cout << "c12 : " << c12;
    std::cout << "dot(v1v2,c12) : " << dot(v1v2,c12) << std::endl; 
    if(dot(v1v2,c12) > 0)
        return true;
    return false;
}

void State::orderVector(int a,std::vector<int>& b)
{
    int min;
    double c = -1;
    vec c1c2;
    double d1;
    std::cout << "b : " << b[1] << std::endl;
    for(int i = 0; i < b.size();i++)
    {
        c1c2 = normalise(discs[b[i]].position - discs[a].position);
        d1 = dot(normalise(discs[a].velocity),c1c2);
        if(d1 > c)
        {
            min = i;
            c = d1;
        }
        std::cout << "b[i] : " << b[i] << " : " << d1 << std::endl;  
    }
    if(min != 1)
    {
        int c = b[1];
        b[1] = b[min];
        b[min] = c;
    }
    std::cout << "b : " << b[1] << std::endl;
}

void State::impactUpdate(const int& d1,const int& d2)
{
    if(!isColliding(d1,d2))
        return;
    //std::cout << "start impact update\n";
    //position vector of 2 wrt 
    //std::cout << "iu1\n";
    vec j;
    vector<vec> m = findPOC(d1,d2);
    vec k = m[0];
    vec G = m[1];
    //updateVelocities(d1,d2);
    //std::cout << "k : {" << k[0] << ", "<< k[1] << ", "<< k[2] << "}\n";
    //    std::cout << "iu2\n";
    double phi = (3*n*(1+e))/(b+1);
    //    std::cout << "iu3\n";
    j = cross(cross(G,k),k);
    //    std::cout << "iu4\n";
    j = normalise(j);
    //    std::cout << "iu5\n";
    //std::cout << "k : " << k;
    //std::cout << "G : " << G;
    double theta = acos(as_scalar((k.t())*G));
    double A;
    double B;
    vec J;
    std::cout << "impact update d1 : " << d1 << "d2 : " << d2 << std::endl;
    //    std::cout << "iu6\n";
    if(d2 < 0)
    {
        //std::cout << "iu7\n";
        A = -(1+e)*discs[d1].mass*as_scalar(discs[d1].velocity.t()*k);
        if(theta > phi)
        {
            B = A*n;
            //std::cout << "1B : " << B << std::endl;        
        }
        else
        {
            //     B = ((b+1)*discs[d1].mass*(norm(discs[d1].velocity)*sin(theta) - discs[d1].radius*norm(discs[d1].angular_velocity))*(double(1)/3));   
            double factor =  (1/discs[d1].mass + ((discs[d1].radius*discs[d1].radius)/discs[d1].I));
            B = (-1)*(b+1)*as_scalar(G.t()*j)/factor;
            //std::cout << "2B : " << B << std::endl;
        }        
        //std::cout << "A : " << A << std::endl;        
        //std::cout << "iu8\n";
        //discs[d1].velocity = -1*e*discs[d1].velocity*cos(theta)*k + (discs[d1].velocity*sin(theta) - B/discs[d1].mass)*j;
        //        std::cout << "iu9\n";
        J = A*k + B*j;
        //std::cout << "J : {" << J[0] << ", "<< J[1] << ", "<< J[2] << "}\n";
        //discs[d1].angular_velocity = discs[d1].angular_velocity + (discs[d1].radius/discs[d1].I)*cross(k,J);
        if(impulses.find(d1) == impulses.end())
        {
            impulses[d1] = {{0,0,0},{0,0,0}};
            //        std::cout << "iu12\n";
        }
        impulses[d1][0] += J;
        impulses[d1][1] += cross(k,J);
        std::cout << "boundary impact\n";
        return;
    }
    ////    std::cout << "iu10\n";
    vec diff_v = discs[d1].velocity - discs[d2].velocity;
    //    std::cout << "iu11\n";
    //    std::cout << diff_v.n_rows << diff_v.n_cols << " " << k.n_rows << k.n_cols <<std::endl;
    //  std::cout << as_scalar((diff_v.t())*k) << std::endl;
    // std::cout << diff_v.n_rows << diff_v.n_cols << " " << k.n_rows << k.n_cols <<std::endl;
    //std::cout << "iu12\n";
    A = as_scalar((diff_v.t())*k)*(-1)*(1+e)*((discs[d1].mass*discs[d2].mass)/(discs[d1].mass + discs[d2].mass));
    //    std::cout << "iu12\n";
    //finding angle between G and k.
    if(theta > phi) //slipping collision
    {
        B = (-1)*n*A;
        //        std::cout << "iu12\n";
    }
    else            //sticking collision
    {   B = -(b+1)*(as_scalar((G.t())*j))/((1/discs[d1].mass) + (1/discs[d2].mass) + (pow(discs[d1].radius,2)/discs[d1].I) + (pow(discs[d2].radius,2)/discs[d2].I)); 
        //        std::cout << "iu12\n";
    }
    J = A*k + B*j;
    //std::cout << "J : {" << J[0] << ", "<< J[1] << ", "<< J[2] << "}\n";
    //    std:://cout << "iu12\n";
    if(impulses.find(d1) == impulses.end())
    {
        impulses[d1] = {{0,0,0},{0,0,0}};
        //        std::cout << "iu12\n";
    }
    if(impulses.find(d2) == impulses.end())
    {
        impulses[d2] = {{0,0,0},{0,0,0}};
        //        std::cout << "iu12\n";
    }
    //    std::cout << "iu12\n";
    impulses[d1][0] += J;
    //    std::cout << "iu13\n";
    impulses[d2][0] -= J;
    //    std::cout << "iu14\n";
    impulses[d1][1] += cross(k,J);
    //    std::cout << "iu15\n";
    impulses[d2][1] += cross(k,J);
    //    std::cout << "end impact update\n";
}

void State::velocityUpdate(int d1)
{
    //    std::cout << "3111\n";
    if(d1 != -1)
    {
        auto it = impulses.find(d1);
        if(it != impulses.end())
        {
            std::cout << "updating impulses of : " << it->first << std::endl;
            //std::cout << "impulses  " << it->second[0] << std::endl;
//            std::cout << "impulses  " << it->second[0]/discs[it->first].mass << std::endl;
//            std::cout << " impulses : " << it->second[1] << std::endl;
//            std::cout << "initial velocity : " << discs[it->first].velocity<< std::endl;
//            std::cout << "initial angular_velocity : " << discs[it->first].angular_velocity<< std::endl;
            discs[it->first].velocity += it->second[0]/discs[it->first].mass;
//            std::cout << "final velocity : " << discs[it->first].velocity<< std::endl;
//            std::cout << " angular impulses : " << (it->second[1])*(discs[it->first].radius/discs[it->first].I)<< std::endl;
            discs[it->first].angular_velocity += it->second[1]*(discs[it->first].radius/discs[it->first].I);
//            std::cout << "final angular_velocity : " << discs[it->first].angular_velocity<< std::endl;
            impulses.erase(it);
            if(norm(discs[it->first].velocity) < 1e-6)
                discs[it->first].velocity = {0,0,0};
            if(norm(discs[it->first].angular_velocity) < 1e-6)
                discs[it->first].angular_velocity = {0,0,0};
        }
        //        std::cout << "3112\n";
        return;
    }
    for(auto it = impulses.begin();it != impulses.end();it++)
    {
        std::cout << "updating impulses of : " << it->first << std::endl;
//        std::cout << "impulses  " << it->second[0] << std::endl;
//        std::cout << "impulses  " << it->second[0]/discs[it->first].mass << std::endl;
//        std::cout << " impulses : " << it->second[1] << std::endl;
//        std::cout << "initial velocity : " << discs[it->first].velocity<< std::endl;
//        std::cout << "initial angular_velocity : " << discs[it->first].angular_velocity<< std::endl;
        discs[it->first].velocity += (it->second[0]/discs[it->first].mass);
//        std::cout << "final velocity : " << discs[it->first].velocity<< std::endl;
 //       std::cout << " angular impulses : " << (it->second[1])*(discs[it->first].radius/discs[it->first].I)<< std::endl;
        discs[it->first].angular_velocity += ((it->second[1])*(discs[it->first].radius/discs[it->first].I));

//        std::cout << "final angular_velocity : " << discs[it->first].angular_velocity<< std::endl;
        if(norm(discs[it->first].velocity) < 1e-6)
            discs[it->first].velocity = {0,0,0};
        if(norm(discs[it->first].angular_velocity) < 1e-6)
            discs[it->first].angular_velocity = {0,0,0};
    }
    impulses.clear();
    //    std::cout << "3112\n";
}

void State::timeUpdate()
{
    std::cout << "before time update : \n";
    printDebug();
    for(auto&& d:discs)
    {
        std::cout << "time : " << predicted_collision.collision_instance << std::endl;
        d.position = d.position + predicted_collision.collision_instance*d.velocity - 0.5*n*g*predicted_collision.collision_instance*predicted_collision.collision_instance*normalise(d.velocity);
        d.velocity = d.velocity - n*g*predicted_collision.collision_instance*normalise(d.velocity);
    }
    std::cout << "after time update : \n";
    printDebug();
}

void State::collisionUpdate()
{
    //    std::cout << "311\n";
    timeUpdate();
    while(predicted_collision.disc_indices.size() > 0)
    {   
        std::cout << "contact pairs : \n";
        for(int m = 0;m < NO_DISCS;m++)
        {
            std::cout << "\n" << m << "=> ";
            for(auto n:discs[m].contact_pairs)
                std::cout << n << ", ";
        }
        for(auto it = predicted_collision.disc_indices.begin();it!=predicted_collision.disc_indices.end();++it)
    {
        std::cout << "DISC INDICES\n";
        std::cout << it->first << "->";
        for(int u:it->second)
            std::cout << u << ", ";
        std::cout <<"\n";
        std::cout << "381\n";
        std::vector<std::vector<int>> stack = stackGenerator(std::pair<int,vector<int>>(*it));
        std::cout << "STACK  :\n";
        for(auto i:stack)
        {
            std::cout << i[0] << ", " << i[1] << std::endl;
            //            std::cout << "391\n";
            impactUpdate(i[0],i[1]);
            //            std::cout << "3101\n";
            removeContact(i[0],i[1]);
            velocityUpdate(i[0]);
        }
    }
    //    std::cout << "3111\n";
    predicted_collision.disc_indices.clear();
    //    std::cout << "3121\n";
    //predicted_collision.collision_instance = 0;
    //    std::cout << "3151\n";
    velocityUpdate();
    for(int i = 0;i < discs.size();i++)
    {
        if(norm(discs[i].velocity) > 1e-6)
        {
            for(int j = 0;j< discs.size();j++)
            {
                if(j != i && (std::find(discs[i].contact_pairs.begin(),discs[i].contact_pairs.end(),j) == discs[i].contact_pairs.end()))
                {
                    if(inContact(i,j))
                        std::cout << "CONTACTING " << i << ", " << j << std::endl;
                }
            }
            for(int k = 0;k < discs[i].contact_pairs.size();k++)
            {
                if(isColliding(i,discs[i].contact_pairs[k]))
                {
                    std::cout << "IMMEDIATE COLLIDING : " << i << ", " << discs[i].contact_pairs[k] << std::endl;
                    predicted_collision.disc_indices[i].push_back(discs[i].contact_pairs[k]);
                    removeContact(i,discs[i].contact_pairs[k]);
                    k--;
                }
            }
            }

        //        std::cout << "3131\n";
        /*
        if(discs[i->first].contact_pairs.size() > 0)
        {
            //            std::cout << "3141\n";
            predicted_collision.disc_indices[i->first] = discs[i->first].contact_pairs;
            std::cout << "i->first : " << i->first << std::endl;
            for(int b : predicted_collision.disc_indices[i->first])
            {
                std::cout << "b : " << b << std::endl;
                removeContact(i->first,b);
            }
        }
        for(int j = 0;j< discs.size();j++)
        {
            if(j != i->first && (std::find(predicted_collision.disc_indices[i->first].begin(),predicted_collision.disc_indices[i->first].end(),j) == predicted_collision.disc_indices[i->first].end()))
            {
                std::cout << "CONTACTING " << i->first << ", " << j << std::endl;
                inContact(i->first,j);
            }
        }*/
    }
    //    std::cout << "312\n";
    }
}

bool State::inContact(const int& d1,const int& d2)
{
    //std::cout << "351\n";
    std::cout << discs[d1].position;
    std::cout << discs[d2].position;
    vec c1c2 = discs[d1].position - discs[d2].position;    
    //std::cout << "352\n";
    std::cout << "norm(c1c2) : "  << norm(c1c2) << "\n";
    std::cout << "disc radii : "  << discs[d1].radius << " " << discs[d2].radius << std::endl;
    std::cout << d1 << " " << d2 << " : " << std::abs(norm(c1c2) - (discs[d2].radius + discs[d1].radius));
    if(std::abs(norm(c1c2) - (discs[d2].radius + discs[d1].radius)) <= 1e-3)
    {
        std::cout << "IN CONTACT : " << d1 << " & " << d2 << std::endl; 
        if(find(discs[d1].contact_pairs.begin(),discs[d1].contact_pairs.end(),d2) == discs[d1].contact_pairs.end())
            discs[d1].contact_pairs.push_back(d2);   
        if(find(discs[d2].contact_pairs.begin(),discs[d2].contact_pairs.end(),d1) == discs[d2].contact_pairs.end())
            discs[d2].contact_pairs.push_back(d1);   
        return true;
    }
    return false;
    //std::cout << "353\n";
}

void State::updateState()
{
    //    std::cout << "31\n";
    if(predicted_collision.collision_instance != INIT_STATE)
    {
        //collision updating
        std::cout << "collision update!!!!!!!!!!!!!!!!!!\n"; 
        collisionUpdate();
    }
    //std::cout << "32\n"

    //create list of indices of active discs
    std::vector<int> activeIndices;
    std::vector<int> nonActiveIndices;
    for(int i = 0;i < discs.size() ;i++)
    {
        discs[i].contact_pairs.clear();
        //std::cout << "33\n";
        if(norm(discs[i].velocity) >  1e-3)
        {
            //std::cout << "34\n";
            activeIndices.push_back(i);
        }
        else
        {
           // std::cout << "35\n";
            for(auto j: nonActiveIndices)
                inContact(i,j);
            nonActiveIndices.push_back(i);
       //     std::cout << "36\n";
        }
    }
    std::cout << "ACTIVE INDICES : "; 
    for(auto i:activeIndices)
        std::cout << i << ", ";
    std::cout << std::endl;
    std::cout << "NON ACTIVE INDICES : "; 
    for(auto i:nonActiveIndices)
        std::cout << i << ", ";
    std::cout << std::endl;
    //    std::cout << "33\n";

    //choose immediate collision

    PIC minPIC(-1,-1,-1);
    for(int i = 0;i < activeIndices.size();i++)
    {
        //        std::cout << "331\n";
        PIC minPICi(PICgenerator(activeIndices[i]));
        //        std::cout << "331\n";
        for(int k = i+1;k < activeIndices.size();k++)
        {
            std::cout << "ActiveIndices[" << k << "] = "<< activeIndices[k]<< std::endl;
            //            std::cout << "332\n";
            double t = PICgenerator(activeIndices[i],activeIndices[k]);
            //            std::cout << "333\n";
                        std::cout << "i,k" << i << "," << k << " : " << t << std::endl; 
            if(t != -1)
            {
                if(fabs(t - minPICi.collision_instance) < 1e-4)
                {
                std::cout << "pushing to : " << activeIndices[i] << "->" << activeIndices[k] << std::endl;
                minPICi.disc_indices[activeIndices[i]].push_back(activeIndices[k]);
                }
                else if(t < minPICi.collision_instance || minPICi.collision_instance == -1)
                {
                    std::cout << "replacing PIC : " << activeIndices[i] << "->" << activeIndices[k] << std::endl;
                    minPICi = PIC(t,activeIndices[i],activeIndices[k]);
                }
            }
        }
        for(int k = 0; k < nonActiveIndices.size(); k++)
        {
            std::cout << "nonActiveIndices[" << k << "] = "<< nonActiveIndices[k]<< std::endl;
            double t = PICgenerator(activeIndices[i],nonActiveIndices[k]);
                        std::cout << "na,i,k" << i << activeIndices[i] << "," << k << nonActiveIndices[k] << " : " << t << std::endl;
            if(t != -1)
            {
                if(fabs(t - minPICi.collision_instance) < 1e-4)
                {
                    std::cout << "pushing to : " << activeIndices[i] << "->" << nonActiveIndices[k] << std::endl;
                    minPICi.disc_indices[activeIndices[i]].push_back(nonActiveIndices[k]);
                }
                else if(t < minPICi.collision_instance || minPICi.collision_instance == -1)
                {
                    std::cout << "replacing PIC : " << activeIndices[i] << "->" << nonActiveIndices[k] << std::endl;
                    minPICi = PIC(t,activeIndices[i],nonActiveIndices[k]);
                }
            }
            /*if(inContact(activeIndices[i],nonActiveIndices[k]) && isColliding(activeIndices[i],nonActiveIndices[k]) && predicted_collision.collision_instance > 1e-13)
            {
                std::cout << "yello man!!!!!!" << std::endl;
                if(minPICi.collision_instance == 0.0)
                    minPICi.disc_indices[activeIndices[i]].push_back(nonActiveIndices[k]);
                else
                    minPICi = PIC(0.0,activeIndices[i],nonActiveIndices[k]);
            }*/
        }
        if(minPICi.collision_instance != -1 && (minPICi.collision_instance < minPIC.collision_instance || minPIC.collision_instance == -1))
        {
            minPIC = minPICi;
        }
    }
    //    std::cout << "34\n";
    if(minPIC.collision_instance == -1)
    {
        predicted_collision.collision_instance = -1;
        computeStoppingTime(activeIndices);
        //push to state queue and stop loop
        return;
    }
    //    std::cout << "35\n";
    predicted_collision = minPIC;
    for (std::map<int,vector<int>>::iterator it=minPIC.disc_indices.begin(); it!=minPIC.disc_indices.end(); ++it)
    {
        std::cout << "collision disc : " << it->first << "-> ";
        for(int m:it->second)
            std::cout << m << ", ";
        cout <<"\n";
    }
}

void State::computeStoppingTime(std::vector<int> activeIndices_)
{
    if(activeIndices_.size() > 0)
    {
        std::cout << "1\n";
        int pos = activeIndices_[0];
        std::cout << "size of compst : " << activeIndices_.size() << std::endl;
        for(int a = 1;a < activeIndices_.size();a++)
        {
            if(norm(discs[activeIndices_[a]].velocity) > norm(discs[activeIndices_[a]].velocity)) 
            {
                pos = a;
            }
        }
        std::cout << "1\n";
        stopping_time = norm(discs[pos].velocity)/(n*g);
    }
    else
        stopping_time = 0;
    std::cout << "STOPPING TIME  : " << stopping_time << std::endl;
}

void State::removeContact(const int& d1,const int& d2)
{
    auto it = find(discs[d1].contact_pairs.begin(),discs[d1].contact_pairs.end(),d2);
    if(it != discs[d1].contact_pairs.end())
    {
        std::cout << "removing " << d2 << "from " << d1 << std::endl;
        discs[d1].contact_pairs.erase(it);
    }
    it = find(discs[d2].contact_pairs.begin(),discs[d2].contact_pairs.end(),d1);
    if(it != discs[d2].contact_pairs.end())
    {
        std::cout << "removing " << d1 << "from " << d2 << std::endl;
        discs[d2].contact_pairs.erase(it);
    }
}

vector<vector<int>> State::stackGenerator(std::pair<int,vector<int>> p)
{
    //    std::cout << "3811\n";
    std::vector<std::vector<int>> stack;
    srand(time(NULL));
    switch(p.second.size())
    {
        case 1: stack.push_back({p.first,p.second[0]});
                break;
        case 2: 
                {
                    std::cout << "3811\n";
                    int i = rand()%2;
                    stack.push_back({p.first,p.second[i]});
                    stack.push_back({p.first,p.second[(i+1)%2]});
                    removeContact(p.second[i],p.second[(i+1)%2]);
                    break;
                }
        case 3:
                {
                    std::cout << "3811\n";
                    int i = rand()%8;
                    std::cout << "3812\n";
                    orderVector(p.first,p.second);
                    std::cout << "3813\n";
                    if(0 <= i && i < 2)
                    {
                        std::cout << "3814\n";
                        stack.push_back({p.first,p.second[1]});
                        std::cout << "3815\n";
                        vector<int> t = {p.second[0],p.second[2]};
                        std::cout << "3816\n";
                        pair<int,vector<int>> o = std::make_pair(p.first,t);
                        std::cout << "3817\n";
                        vector<vector<int>> s = stackGenerator(o);
                        std::cout << "3818\n";
                        for(auto it: s)
                        {
                            stack.push_back(it);
                        }
                    }
                    else if(i < 5)
                    {
                        std::cout << "3819\n";
                        stack.push_back({p.first,p.second[0]});
                        std::cout << "3820\n";
                        stack.push_back({p.first,p.second[1]});
                        std::cout << "3821\n";
                        stack.push_back({p.first,p.second[2]});
                    }
                    else
                    {
                        std::cout << "3822\n";
                        stack.push_back({p.first,p.second[2]});
                        std::cout << "3823\n";
                        stack.push_back({p.first,p.second[1]});
                        std::cout << "3824\n";
                        stack.push_back({p.first,p.second[0]});
                    }
                    std::cout << "3825\n";
                    removeContact(p.second[0],p.second[1]);
                    std::cout << "3826\n";
                    removeContact(p.second[0],p.second[2]);
                    std::cout << "3827\n";
                    removeContact(p.second[2],p.second[1]);
                    std::cout << "3828\n";
                    break;
                }
    }
    return stack;
}

/** to generate predicted instance of collision for disc of given index and
 * the wall
 */

PIC State::PICgenerator(const int& d1)
{
    double constr1 = discs[d1].radius,constr2 = board_length - discs[d1].radius;
    std::vector<vec> a = {{constr1,-1,0},{-1,constr2,0},{constr2,-1,0},{-1,constr1,0}};
    double tMin = -1;
    int index;
    int wall_no = 5;
    for(int i = 0; i < a.size(); i++)
    {
        if(a[i][0] == -1)
            index = 1;
        else 
            index = 0;
        double coeff0 = a[i][index] - discs[d1].position[index];
        double coeff1 = discs[d1].velocity[index]*(-1);
        vec nd = normalise(discs[d1].velocity);
        double coeff2 = nd[index]*n*g*0.5;
        if(nd[index] < 0)
        {
            coeff0 *= -1;
            coeff1 *= -1;
            coeff2 *= -1;
        }
        std::cout << "coeff0 : " << coeff0 <<std::endl;
        std::cout << "coeff1 : " << coeff1 <<std::endl;
        std::cout << "coeff2 : " << coeff2 <<std::endl;
        double root = 0;
        if(solveQuadratic(coeff2,coeff1,coeff0,root))
        {
            std::cout << "root : " << root << "i : " << i << std::endl;
            std::cout << "nd : " << nd ;
            index++;
            index%=2;
            std::cout << "nd : " << nd ;
            double other_coord = n*g*(-0.5)*nd[index];
            if(n*g != 0)
            {
                std::vector<int> l = {d1};
                computeStoppingTime(l);

            }
            std::cout << "yo\n";
            other_coord =other_coord + root*discs[d1].velocity[index] + discs[d1].position[index];
            if(constr1 <= other_coord <= constr2)
            {
                if((root < tMin || tMin == -1) && root > 1e-3 && coeff0 > 1e-3)
                {    
                    if(n*g == 0  || root < stopping_time)
                    { 
                        tMin = root;
                        wall_no = i+1;
                    }
                }
            }
            std::cout << "yo\n";
        }
    }
    return PIC(tMin,d1,-1*wall_no);
}


/** to generate predicted instance of collision for 2 discs of given indices
*/
double State::PICgenerator(const int& d1,const int& d2)
{
    vec c1c2 = discs[d2].position - discs[d1].position;
    double coeff0 = dot(c1c2,c1c2) - (discs[d1].radius + discs[d2].radius)*(discs[d1].radius + discs[d2].radius);
    vec v2v1 = discs[d2].velocity - discs[d1].velocity;
    double coeff1 = 2*dot(v2v1,c1c2);
    //    std::cout << "here\n";
    vec a2a1 = (-1*n*g)*(normalise(discs[d2].velocity) - normalise(discs[d1].velocity));
    //    std::cout << "here  ok\n";
    double coeff3 = dot(a2a1,v2v1);
    double coeff2 = dot(c1c2,a2a1) + dot(v2v1,v2v1);
    double coeff4 = dot(a2a1,a2a1)/4;
        std::cout << "coeff4 : " << coeff4 << std::endl;
        std::cout << "coeff3 : " << coeff3 << std::endl;
        std::cout << "coeff2 : " << coeff2 << std::endl;
        std::cout << "coeff1 : " << coeff1 << std::endl;
        std::cout << "coeff0 : " << coeff0 << std::endl;
    double root = 0;
    bool tExists = solveQuartic(coeff4,coeff3,coeff2,coeff1,coeff0,root);
    std::cout << "tExists : " << tExists << std::endl;
    std::cout << "root found d1,d2 : " << d1 << ", " << d2 << " :" << root << std::endl;
    if(n*g != 0)
    {
        std::vector<int> l = {d1};
        computeStoppingTime(l);
        if(root >= stopping_time)
        {
            std::cout << "stopping b4 collision\n";
            return -1;
        }
    }
    if(tExists && isColliding(d1,d2))
    {
        return root;
    }
    else 
        return -1;
}

void StateQ::updateQueue()
{
    State a = q.front();
    a.printDebug();
    char ch;
    //cin >> ch;
    do
    {
        std::cout << "update loop iteration!!!!!!!!!!!!!!!!!!!!!!\n";
        //        std::cout << "3\n";
        a.updateState();
        q.push(a);
        a.printDebug();
        //cin >> ch;
        //        std::cout << "4\n";

    }while(a.predicted_collision.collision_instance >= 0 );
}

void StateQ::initialiseQueue(vec& initVel_,int noDiscs_,std::vector<vec> positionMap_)
{
    //    std::cout<< "1\n" ;
    State a(noDiscs_);
    a.discs[0].velocity = initVel_;
    a.predicted_collision.collision_instance = -2;
    a.stopping_time = -1;
    for(int o = 0;o < noDiscs_;o++)
    {
        a.discs[o].position = positionMap_[o];
    }
    q.push(a);
    /*
       for(auto d : q.front().discs)
       {
       std::cout << "velocity : " << d.velocity[0] <<","<< d.velocity[1] << "," << d.velocity[2] << std::endl; 
       std::cout << "position : " << d.position[0] <<","<< d.position[1] << "," << d.position[2] << std::endl; 
       }
       char ch;
       std::cin >> ch; 
       */
    computationThread = std::thread(&StateQ::updateQueue,this);
    std::cout<< "2\n" ;
}

void State::printDebug()
{
    for(int i = 0;i < discs.size();i++)
    {
        std::cout << "disc index :  " << i << std::endl;
        std::cout << "velocity : {" << discs[i].velocity[0] << ", " << discs[i].velocity[1] << ", " << discs[i].velocity[2] << "\n" ; 
        std::cout << "angular_velocity : {" << discs[i].angular_velocity[0] << ", " << discs[i].angular_velocity[1] << ", " << discs[i].angular_velocity[2] << "\n" ; 
        std::cout << "position : {" << discs[i].position[0] << ", " << discs[i].position[1] << ", " << discs[i].position[2] << "\n" ; 
    }
    std::cout << "pic   t: " << predicted_collision.collision_instance <<std::endl;
    for (std::map<int,vector<int>>::iterator it=predicted_collision.disc_indices.begin(); it!=predicted_collision.disc_indices.end(); ++it)
    {   std::cout << it->first << " => ";
        for(int a: it->second)
        {
            std::cout << a << ", ";
        }
        std::cout << std::endl;
    }
    std::cout << "contact pairs : \n";
    for(int m = 0;m < NO_DISCS;m++)
    {
        std::cout << "\n" << m << "=> ";
        for(auto n:discs[m].contact_pairs)
            std::cout << n << ", ";
    }
    std::cout << "STOPPING TIME : " << stopping_time << std::endl;
}

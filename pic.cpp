#include "disc.hpp"

double n = 0,g = 10,e = 1,b = 1,board_length;

map<int,vec> boundaries = {{-1,{1,0,0}},{-2,{0,-1,0}},{-3, {-1,0,0}},{-4, {0,1,0}}};

vec State::findPOC(const int& d1,const int& d2)
{
    discs[d1].position = discs[d1].position + predicted_collision.collision_instance*discs[d1].velocity - 0.5*n*g*predicted_collision.collision_instance*predicted_collision.collision_instance*normalise(discs[d1].velocity);
    if(d2 < 0)
    {
        return boundaries[d2];
    }
    discs[d2].position = discs[d2].position + predicted_collision.collision_instance*discs[d2].velocity - 0.5*n*g*predicted_collision.collision_instance*predicted_collision.collision_instance*normalise(discs[d2].velocity);
    vec c1c2 = discs[d2].position - discs[d1].position;
    return normalise(c1c2);
}

void State::updateVelocities(const int& d1,const int& d2)
{
    discs[d1].velocity -= n*g*predicted_collision.collision_instance*normalise(discs[d1].velocity);
    discs[d2].velocity -= n*g*predicted_collision.collision_instance*normalise(discs[d2].velocity);
}

void State::impactUpdate(const int& d1,const int& d2)
{
    //std::cout << "start impact update\n";
    //position vector of 2 wrt 
    //std::cout << "iu1\n";
    vec k,j;
    k = findPOC(d1,d2);
    updateVelocities(d1,d2);
    std::cout << "k : {" << k[0] << ", "<< k[1] << ", "<< k[2] << "}\n";
//    std::cout << "iu2\n";
    double phi = (3*n*(1+e))/(b+1);
//    std::cout << "iu3\n";
    vec G = discs[d1].radius*cross(discs[d1].angular_velocity,k) + discs[d1].velocity - discs[d2].velocity - discs[d2].radius*cross(discs[d2].angular_velocity,k);
    j = cross(cross(G,k),k);
//    std::cout << "iu4\n";
    j = normalise(j);
//    std::cout << "iu5\n";
    double theta = acos(as_scalar((k.t())*G));
    double A;
    double B;
    vec J;
//    std::cout << "iu6\n";
    if(d2 < 0)
    {
//        std::cout << "iu7\n";
        A = -(1+e)*discs[d1].mass*as_scalar(discs[d1].velocity*k.t());
        if(theta > phi)
            B = A*n;
        else
            B = fabs((b+1)*discs[d1].mass*(norm(discs[d1].velocity)*sin(theta) - discs[d1].radius*norm(discs[d1].angular_velocity))*(1/3));    
//        std::cout << "iu8\n";
        discs[d1].velocity = -1*e*discs[d1].velocity*cos(theta)*k + (discs[d1].velocity*sin(theta) - B/discs[d1].mass)*j;
//        std::cout << "iu9\n";
        J = A*k + B*j;
    std::cout << "J : {" << J[0] << ", "<< J[1] << ", "<< J[2] << "}\n";
        discs[d1].angular_velocity = discs[d1].angular_velocity + (discs[d1].radius/discs[d1].I)*cross(k,J);
    }
////    std::cout << "iu10\n";
    vec diff_v = discs[d1].velocity - discs[d2].velocity;
//    std::cout << "iu11\n";
    //    std::cout << diff_v.n_rows << diff_v.n_cols << " " << k.n_rows << k.n_cols <<std::endl;
    //  std::cout << as_scalar((diff_v.t())*k) << std::endl;
    // std::cout << diff_v.n_rows << diff_v.n_cols << " " << k.n_rows << k.n_cols <<std::endl;
    std::cout << "iu12\n";
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
    std::cout << "J : {" << J[0] << ", "<< J[1] << ", "<< J[2] << "}\n";
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
            discs[d1].velocity += impulses[d1][0]/discs[d1].mass;
            discs[d1].angular_velocity += impulses[d1][1]*(discs[d1].radius/discs[d1].I);
            impulses.erase(it);
        }
//        std::cout << "3112\n";
        return;
    }
    for(auto it = impulses.begin();it != impulses.end();it++)
    {
        std::cout << "updating impulses of : " << it->first << std::endl;
        std::cout << "impulses  " << it->second[0] << std::endl;
        std::cout << "impulses  " << it->second[0]/discs[it->first].mass << std::endl;
        std::cout << " impulses : " << it->second[1] << std::endl;
        std::cout << "initial velocity : " << discs[it->first].velocity<< std::endl;
        std::cout << "initial angular_velocity : " << discs[it->first].angular_velocity<< std::endl;
        discs[it->first].velocity += (it->second[0]/discs[it->first].mass);
        std::cout << "final velocity : " << discs[it->first].velocity<< std::endl;
        std::cout << " angular impulses : " << (it->second[1])*(discs[it->first].radius/discs[it->first].I)<< std::endl;
        discs[it->first].angular_velocity += ((it->second[1])*(discs[it->first].radius/discs[it->first].I));
       
        std::cout << "final angular_velocity : " << discs[it->first].angular_velocity<< std::endl;
        if(norm(discs[it->first].velocity) < 1e-6)
            discs[it->first].velocity = {0,0,0};
        if(norm(discs[it->first].angular_velocity) < 1e-6)
            discs[it->first].angular_velocity = {0,0,0};
    }
//    std::cout << "3112\n";
}

void State::collisionUpdate()
{
//    std::cout << "311\n";
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
            velocityUpdate(i[0]);
        }
    }
//    std::cout << "3111\n";
    predicted_collision.disc_indices.clear();
//    std::cout << "3121\n";
    predicted_collision.collision_instance = 0;
    for(auto i = impulses.begin();i != impulses.end();++i)
    {
//        std::cout << "3131\n";
        if(discs[i->first].contact_pairs.size() != 0)
        {
//            std::cout << "3141\n";
            predicted_collision.disc_indices[i->first] = discs[i->first].contact_pairs; 
        }
    }
//    std::cout << "3151\n";
    velocityUpdate();
//    std::cout << "312\n";
}

void State::inContact(const int& d1,const int& d2)
{
    vec c1c2 = discs[d1].position - discs[d2].position;    
    if(norm(c1c2) <= discs[d2].radius + discs[d1].radius)
    {
        discs[d1].contact_pairs.push_back(d2);   
        discs[d2].contact_pairs.push_back(d1);   
    }
}

void State::updateState()
{
//    std::cout << "31\n";
    if(predicted_collision.collision_instance != INIT_STATE)
    {
        //collision updating
//        std::cout << "collision update!!!!!!!!!!!!!!!!!!\n"; 
        collisionUpdate();
    }
//    std::cout << "32\n";

    //create list of indices of active discs
    std::vector<int> activeIndices;
    std::vector<int> nonActiveIndices;
    for(int i = 0;i < discs.size() ;i++)
    {
        if(find(activeIndices.begin(),activeIndices.end(),i) == activeIndices.end())
        {
            if(norm(discs[i].velocity)!=0)
                activeIndices.push_back(i);
            else
            {
                nonActiveIndices.push_back(i);
                for(int k = i+1;k < discs.size();k++)
                {
                    if(find(activeIndices.begin(),activeIndices.end(),i) == activeIndices.end())
                    {
                        if(norm(discs[k].velocity)!=0)
                            activeIndices.push_back(i);
                        else
                        {
                            nonActiveIndices.push_back(k);
                            inContact(i,k);
                        }
                    }
                }
            }
        }
    }
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
//            std::cout << "332\n";
            double t = PICgenerator(activeIndices[i],activeIndices[k]);
//            std::cout << "333\n";
//            std::cout << "i,k" << i << "," << k << " : " << t << std::endl; 
            if(t != -1 && t == minPICi.collision_instance || minPICi.collision_instance == -1)
            {
                minPICi.disc_indices[activeIndices[i]].push_back(activeIndices[k]);
            }
            if(t != -1 && (t < minPICi.collision_instance || minPICi.collision_instance == -1))
            {
                minPICi = PIC(t,activeIndices[i],activeIndices[k]);
            }
        }
        for(int k = 0; k < nonActiveIndices.size(); k++)
        {
            double t = PICgenerator(activeIndices[i],nonActiveIndices[k]);
//            std::cout << "na,i,k" << i << activeIndices[i] << "," << k << nonActiveIndices[k] << " : " << t << std::endl; 
            if(t != -1 && t == minPICi.collision_instance || minPICi.collision_instance == -1)
            {
                minPICi.disc_indices[activeIndices[i]].push_back(nonActiveIndices[k]);
            }
            if(t != -1 && (t < minPICi.collision_instance || minPICi.collision_instance == -1))
            {
                minPICi = PIC(t,activeIndices[i],nonActiveIndices[k]);
            }
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
    //push to state queue
    if(predicted_collision.collision_instance == -1)
        computeStoppingTime(activeIndices);
//    std::cout << "36\n";
}

void State::computeStoppingTime(std::vector<int> activeIndices_)
{
    if(activeIndices_.size() > 0)
    {
        int pos = activeIndices_[0];
        for(int a = 1;a < activeIndices_.size();a++)
        {
            if(norm(discs[activeIndices_[a]].velocity) > norm(discs[activeIndices_[a]].velocity)) 
            {
                pos = a;
            }
        }
        stopping_time = norm(discs[activeIndices_[pos]].velocity)/g;
    }
    else
        stopping_time = 0;
}

void State::removeContact(const int& d1,const int& d2)
{
    auto it = find(discs[d1].contact_pairs.begin(),discs[d1].contact_pairs.end(),d2);
    if(it != discs[d1].contact_pairs.end())
        discs[d1].contact_pairs.erase(it);
    it = find(discs[d2].contact_pairs.begin(),discs[d2].contact_pairs.end(),d1);
    if(it != discs[d2].contact_pairs.end())
        discs[d2].contact_pairs.erase(it);
}

void orderVector(vector<int> a){}

vector<vector<int>> State::stackGenerator(std::pair<int,vector<int>> p)
{
    std::vector<std::vector<int>> stack;
    srand(time(NULL));
    switch(predicted_collision.disc_indices[p.first].size())
    {
        case 1: stack.push_back({p.first,p.second[0]});
                break;
        case 2: 
                {
                    int i = rand()%2;
                    stack.push_back({p.first,p.second[i]});  
                    stack.push_back({p.first,p.second[(i+1)%2]});
                    removeContact(p.second[i],p.second[(i+1)%2]);
                    break;
                }
        case 3:
                {
                    int i = rand()%8;
                    orderVector(p.second);
                    if(0 <= i < 2)
                    {
                        stack.push_back({p.first,p.second[1]});
                        vector<int> t = {p.second[0],p.second[2]};
                        pair<int,vector<int>> o = std::make_pair(p.first,t);
                        vector<vector<int>> s = stackGenerator(o);
                        for(auto it: s)
                        {
                            stack.push_back(it);
                        }
                    }
                    else if(i < 5)
                    {
                        stack.push_back({p.first,p.second[0]});
                        stack.push_back({p.first,p.second[1]});
                        stack.push_back({p.first,p.second[2]});
                    }
                    else
                    {
                        stack.push_back({p.first,p.second[0]});
                        stack.push_back({p.first,p.second[1]});
                        stack.push_back({p.first,p.second[2]});
                    }
                    removeContact(p.second[0],p.second[1]);
                    removeContact(p.second[0],p.second[2]);
                    removeContact(p.second[2],p.second[1]);
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
    std::vector<vec> a = {{constr1,-1,0},{constr2,-1,0},{-1,constr1,0},{-1,constr2,0}};
    double tMin = -1;
    int index;
    int wall_no = 5;
    for(int i = 0; i < a.size(); i++)
    {
        if(a[i](0,0) == -1)
            index = 0;
        else 
            index = 1;
        double coeff0 = discs[d1].position(index,0) - a[i](index,0);
        double coeff1 = discs[d1].velocity(index,0);
        vec nd = normalise(discs[d1].velocity);
        double coeff2 = nd(index,0)*n*g*0.5;
        double root;
        double zero = 0;
        if(solveCubic(zero,coeff2,coeff1,coeff0,root))
        {
            index++;
            index%=2;
            double other_coord = n*g*0.5*nd(index,0);
            other_coord =other_coord + root*discs[d1].velocity(index,0) + discs[d1].position(index,0);
            if(constr1 <= other_coord <= constr2)
            {
                if(root < tMin || tMin == -1)
                {    
                    tMin = root;
                    wall_no = i+1;
                }
            }
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
    std::cout << "coeff0 : " << coeff0 << std::endl;
    vec v2v1 = discs[d2].velocity - discs[d1].velocity;
    double coeff1 = 2*dot(v2v1,c1c2);
    std::cout << "coeff1 : " << coeff1 << std::endl;
//    std::cout << "here\n";
    vec a2a1 = n*g*(normalise(discs[d2].velocity) - normalise(discs[d1].velocity));
//    std::cout << "here  ok\n";
    double coeff3 = dot(a2a1,v2v1);
    std::cout << "coeff3 : " << coeff3 << std::endl;
    double coeff2 = dot(c1c2,a2a1) + dot(v2v1,v2v1);
    std::cout << "coeff2 : " << coeff2 << std::endl;
    double coeff4 = dot(a2a1,a2a1)/2;
    std::cout << "coeff4 : " << coeff4 << std::endl;
    double root;
    bool tExists = solveQuartic(coeff4,coeff3,coeff2,coeff1,coeff0,root);
    std::cout << "tExists : " << tExists << std::endl;
    std::cout << "root found d1,d2 : " << d1 << ", " << d2 << " :" << root << std::endl;
    if(tExists && root > 0)
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
    for(auto d:discs)
    {
        std::cout << "velocity : {" << d.velocity[0] << ", " << d.velocity[1] << ", " << d.velocity[2] << "\n" ; 
        std::cout << "angular_velocity : {" << d.angular_velocity[0] << ", " << d.angular_velocity[1] << ", " << d.angular_velocity[2] << "\n" ; 
        std::cout << "position : {" << d.position[0] << ", " << d.position[1] << ", " << d.position[2] << "\n" ; 
    }
    std::cout << "pic   t: " << predicted_collision.collision_instance <<std::endl;
    for (std::map<int,vector<int>>::iterator it=predicted_collision.disc_indices.begin(); it!=predicted_collision.disc_indices.end(); ++it)
        std::cout << it->first << " => " << it->second[0] << '\n';
}

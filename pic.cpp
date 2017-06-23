#include "disc.hpp"

double n,g,e,b,board_length;

map<int,vec> boundaries = {{-1,{1,0,0}},{-2,{0,-1,0}},{-3, {-1,0,0}},{-4, {0,1,0}}};

vec State::findPOC(const int& d1,const int& d2)
{
   if(d2 < 0)
   {
       return boundaries[d2];
   }
   vec p2 = discs[d2].position + predicted_collision.collision_instance*discs[d2].velocity + 0.5*n*g*predicted_collision.collision_instance*predicted_collision.collision_instance*normalise(discs[d2].velocity);
   vec p1 = discs[d1].position + predicted_collision.collision_instance*discs[d1].velocity + 0.5*n*g*predicted_collision.collision_instance*predicted_collision.collision_instance*normalise(discs[d1].velocity);
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
    double theta = acos(as_scalar((k.t())*G));
    double A;
    double B;
    vec J;
    if(d2 < 0)
    {
        A = -(1+e)*discs[d1].mass*as_scalar(discs[d1].velocity*k.t());
        if(theta > phi)
            B = A*n;
        else
            B = fabs((b+1)*discs[d1].mass*(norm(discs[d1].velocity)*sin(theta) - discs[d1].radius*norm(discs[d1].angular_velocity))*(1/3));    
        discs[d1].velocity = -1*e*discs[d1].velocity*cos(theta)*k + (discs[d1].velocity*sin(theta) - B/discs[d1].mass)*j;
        J = A*k + B*j;
        discs[d1].angular_velocity = discs[d1].angular_velocity + (discs[d1].radius/discs[d1].I)*cross(k,J);
    }
    vec diff_v = discs[d1].velocity - discs[d2].velocity;
    std::cout << diff_v.n_rows << diff_v.n_cols << " " << k.n_rows << k.n_cols <<std::endl;
    std::cout << as_scalar((diff_v.t())*k) << std::endl;
    std::cout << diff_v.n_rows << diff_v.n_cols << " " << k.n_rows << k.n_cols <<std::endl;
    A = as_scalar((diff_v.t())*k);//*(d1.mass*d2.mass/(d1.mass + d2.mass))*(-1)*(1+e);
    //finding angle between G and k.
    if(theta > phi) //slipping collision
        B = (-1)*n*A;
    else            //sticking collision
        B = -(b+1)*(as_scalar((G.t())*j))/((1/discs[d1].mass) + (1/discs[d2].mass) + (pow(discs[d1].radius,2)/discs[d1].I) + (pow(discs[d2].radius,2)/discs[d2].I)); 
    J = A*k + B*j;
    if(impulses.find(d1) == impulses.end())
        impulses[d1] = {0,0};
    if(impulses.find(d2) == impulses.end())
        impulses[d2] = {0,0};
    impulses[d1][0] += J;
    impulses[d2][0] -= J;
    impulses[d1][1] += cross(k,J);
    impulses[d2][1] += cross(k,J);
}

void State::velocityUpdate(int d1 = -1)
{
    if(d1 != -1)
    {
        auto it = impulses.find(b);
        if(it != impulses.end())
        {
            discs[d1].velocity += impulses[d1][0]/discs[d1].mass;
            discs[d1].angular_velocity += impulses[d1][1]*(discs[d1].radius/discs[d1].I);
            impulses.erase(it);
        }
        return;
    }

    if(impulses.find(d1) != impulses.end())
    {
        discs[d1].velocity += impulses[d1][0]/discs[d1].mass;
        discs[d1].angular_velocity += impulses[d1][1]*(discs[d1].radius/discs[d1].I);
    }
}

void State::collisionUpdate()
{
    for(auto it = predicted_collision.disc_indices.begin();it!=predicted_collision.disc_indices.end();++it)
    {
        std::vector<std::vector<int>> stack = stackGenerator(std::pair<int,vector<int>>(*it));
        for(auto i:stack)
        {
           impactUpdate(i[0],i[1]);
           velocityUpdate(i[0]);
        }
    }
    predicted_collision.disc_indices.clear();
    predicted_collision.collision_instance = 0;
    for(auto i = impulses.begin();i != impulses.end();++i)
    {
        if(discs[i->first].contact_pairs.size() != 0)
        {
              predicted_collision.disc_indices[i->first] = discs[i->first].contact_pairs; 
        }
    }
    velocityUpdate();
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
    //collision updating 
    collisionUpdate();
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

    //choose immediate collision

    PIC minPIC(-1,-1,-1);
    for(int i = 0;i < activeIndices.size();i++)
    {
        PIC minPICi(PICgenerator(i));
        for(int k = i+1;k < activeIndices.size();k++)
        {
            double t = PICgenerator(i,k);
            if(t != -1 && (t < minPICi.collision_instance || minPICi.collision_instance == -1))
            {
                minPICi = PIC(t,i,k);
            }
            if(t != -1 && t == minPICi.collision_instance || minPICi.collision_instance == -1)
            {
                minPICi.disc_indices[i].push_back(k);
            }
        }
        for(int k = 0; k < nonActiveIndices.size(); k++)
        {
            double t = PICgenerator(i,k);
            if(t != -1 && (t < minPICi.collision_instance || minPICi.collision_instance == -1))
            {
                minPICi = PIC(t,i,k);
            }
            if(t != -1 && t == minPICi.collision_instance || minPICi.collision_instance == -1)
            {
                minPICi.disc_indices[i].push_back(k);
            }
        }
        if(minPICi.collision_instance != -1 && (minPICi.collision_instance < minPIC.collision_instance || minPIC.collision_instance == -1))
        {
            minPIC = minPICi;
        }
    }
    
    if(minPIC.collision_instance == -1)
    {
        predicted_collision.collision_instance = -1;
        //push to state queue and stop loop
        return;
    }
    predicted_collision = minPIC;
    //push to state queue
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
    std::vector<vec> a = {{constr1,-1},{constr2,-1},{-1,constr1},{-1,constr2}};
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
    double coeff0 = as_scalar(c1c2.t()*c1c2);
    vec v2v1 = discs[d2].velocity - discs[d1].velocity;
    double coeff1 = 2*as_scalar(v2v1.t()*c1c2);
    mat acc = {{n*g,0,0},{0,n*g,0},{0,0,n*g}};
    vec a2a1 = acc%(normalise(discs[d2].velocity) - normalise(discs[d1].velocity));
    double coeff3 = dot(a2a1,v2v1);
    double coeff2 = dot(c1c2,a2a1);
    double coeff4 = dot(a2a1,a2a1)/2;
    double root;
    bool tExists = solveQuartic(coeff4,coeff3,coeff2,coeff1,coeff0,root);
    if(tExists)
    {
        return root;
    }
    else 
        return -1;
}

void StateQ::updateQueue()
{
    State a;
    do
    {
        a.updateState();
        q.push(a);
    }while(a.predicted_collision.collision_instance >= 0 );
}

int main()
{}

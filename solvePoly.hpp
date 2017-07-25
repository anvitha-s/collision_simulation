#include <iostream>
#include <cmath>
#ifndef _SOLVE_POLY_HPP_
#define _SOLVE_POLY_HPP_

/**
 * functions solveCubic and solveQuartic to solve 4th order polynomials 
 * from :
 * https://stackoverflow.com/questions/11837078/initialize-a-constant-sized-array-in-an-initializer-list
 */

using namespace std;

const double PI = 3.14159265358979323846;
/*
//----------------------------------------------------------------------------
bool solveQuadratic(double a, double b, double c, double &root)
{
if(a == 0.0 || abs(a/b) < 1.0e-10)
{
std::cout << "ignoring coeff!!\n";
if(abs(b) < 1.0e-4)
{
std::cout << "ignoring coeff1!!\n";
return false;
}
else
{
root = -c/b;
std::cout << "a,b,c,root1 : " << a << ", " << b << ", " << c << ", " << root << std::endl;
return true;
}
}

double discriminant = b * b - 4.0 * a * c;
if(discriminant >= 0.0)
{
discriminant = sqrt(discriminant);
root = (b + discriminant) * -0.5 / a;
std::cout << "a,b,c,root1 : " << a << ", " << b << ", " << c << ", " << root << std::endl;
return true;
}

return false;
}
bool solveQuadraticOther(double a, double b, double c, double &root)
{
if(a == 0.0 || abs(a/b) < 1.0e-10)
{
std::cout << "ignoring coeff!!\n";
if(abs(b) < 1.0e-4) 
{
std::cout << "ignoring coeff1!!\n";
return false;
}
else
{
root = -c/b;
std::cout << "a,b,c,root1 : " << a << ", " << b << ", " << c << ", " << root << std::endl;
return true;
}
}

double discriminant = b * b - 4.0 * a * c;
if(discriminant >= 0.0)
{
discriminant = sqrt(discriminant);
root = (b - discriminant) * -0.5 / a;
std::cout << "a,b,c,root1 : " << a << ", " << b << ", " << c << ", " << root << std::endl;
return true;
}

return false;
}

bool solveCubic(double a, double b, double c, double d, double &root)
{
if(a == 0.0 || abs(a/b) < 1.0e-10)
{
std::cout << "ignoring coeff!!\n";
return solveQuadratic(b, c, d, root);
}

double B = b/a, C = c/a, D = d/a;

double Q = (B*B - C*3.0)/9.0, QQQ = Q*Q*Q;
double R = (2.0*B*B*B - 9.0*B*C + 27.0*D)/54.0, RR = R*R;

// 3 real roots
if(RR<QQQ)
{
    // This sqrt and division is safe, since RR >= 0, so QQQ > RR,    
    // so QQQ > 0.  The acos is also safe, since RR/QQQ < 1, and      
    // thus R/sqrt(QQQ) < 1.                                     
    double theta = acos(R/sqrt(QQQ));
    // This sqrt is safe, since QQQ >= 0, and thus Q >= 0             
    double r1, r2, r3;
    r1 = r2 = r3 = -2.0*sqrt(Q);
    r1 *= cos(theta/3.0);
    r2 *= cos((theta+2*PI)/3.0);
    r3 *= cos((theta-2*PI)/3.0);

    r1 -= B/3.0;
    r2 -= B/3.0;
    r3 -= B/3.0; 

    root = 1000000.0;

    if(r1 >= 0.0) root = r1;
    if(r2 >= 0.0 && r2 < root) root = r2;
    if(r3 >= 0.0 && r3 < root) root = r3;

    std::cout << "a,b,c,d,root1 : " << a << ", " << b << ", " << c << ", " << d << ", " << root << std::endl;
    return true;
}
// 1 real root
else
{
    double A2 = -pow(fabs(R)+sqrt(RR-QQQ),1.0/3.0);
    if (A2!=0.0) 
    {
        if (R<0.0) A2 = -A2; 
        root = A2 + Q/A2; 
    }
    root -= B/3.0;
    std::cout << "a,b,c,d,root1 : " << a << ", " << b << ", " << c << ", " << d << ", " << root << std::endl;
    return true;
}
}

bool solveQuartic(double a, double b, double c, double d, double e, double &root)
{
    // I switched to this method, and it seems to be more numerically stable.
    // http://www.gamedev.net/community/forums/topic.asp?topic_id=451048 

    // When a or (a and b) are magnitudes of order smaller than C,D,E
    // just ignore them entirely. This seems to happen because of numerical
    // inaccuracies of the line-circle algorithm. I wanted a robust solver,
    // so I put the fix here instead of there.
    if(a == 0.0 || abs(a/b) < 1.0e-10 || abs(a/c) < 1.0e-10 || abs(a/d) < 1.0e-10)
    {
        std::cout << "ignoring coeff!!\n";
        return solveCubic(b, c, d, e, root);
    }

    double B = b/a, C = c/a, D = d/a, E = e/a;
    double BB = B*B;
    double I = -3.0*BB*0.125 + C;
    double J = BB*B*0.125 - B*C*0.5 + D;
    double K = -3*BB*BB/256.0 + C*BB/16.0 - B*D*0.25 + E;

    double z;
    bool foundRoot2 = false, foundRoot3 = false, foundRoot4 = false, foundRoot5 = false;
    if(solveCubic(1.0, I+I, I*I - 4*K, -(J*J), z))
    {
        double value = z*z*z + z*z*(I+I) + z*(I*I - 4*K) - J*J;

        double p = sqrt(z);
        double r = -p;
        double q = (I + z - J/p)*0.5;
        double s = (I + z + J/p)*0.5;

        bool foundRoot = false, foundARoot;
        double aRoot;
        foundRoot = solveQuadratic(1.0, p, q, root);
        root -= B/4.0;

        foundARoot = solveQuadratic(1.0, r, s, aRoot);
        aRoot -= B/4.0;
        if((foundRoot && foundARoot && ((aRoot < root && aRoot >= 0.0) 
                        || root < 0.0)) || (!foundRoot && foundARoot)) 
        {
            root = aRoot;
            std::cout << "a,b,c,d,e,root1 : " << a << ", " << b << ", " << c << ", " << d << ", " << e << ", " << root << std::endl;
            foundRoot = true;
        }

        foundARoot = solveQuadraticOther(1.0, p, q, aRoot);
        aRoot -= B/4.0;
        if((foundRoot && foundARoot && ((aRoot < root && aRoot >= 0.0) 
                        || root < 0.0)) || (!foundRoot && foundARoot)) 
        {
            root = aRoot;
            std::cout << "a,b,c,d,e,root1 : " << a << ", " << b << ", " << c << ", " << d << ", " << e << ", " << root << std::endl;
            foundRoot = true;
        }

        foundARoot = solveQuadraticOther(1.0, r, s, aRoot);
        aRoot -= B/4.0;
        if((foundRoot && foundARoot && ((aRoot < root && aRoot >= 0.0) 
                        || root < 0.0)) || (!foundRoot && foundARoot)) 
        {
            root = aRoot;
            std::cout << "a,b,c,d,e,root1 : " << a << ", " << b << ", " << c << ", " << d << ", " << e << ", " << root << std::endl;
            foundRoot = true;
        }
        return foundRoot;
    }
    return false;
}
    */
bool solveQuadratic(double a, double b, double c, double &root)
{
    if(a == 0.0 || abs(a/b) < 1.0e-6)
    {
        std::cout << "ignoring coeff1!!\n";
        if(abs(b) < 1.0e-4) 
        {
            std::cout << "ignoring coeff2!!\n";
            return false;
        }
        else
        {
            root = -c/b;
            return true;
        }
    }

    double discriminant = b * b - 4.0 * a * c;
    if(discriminant >= 0.0)
    {
        discriminant = sqrt(discriminant);
        root = (b + discriminant) * -0.5 / a;
        return true;
    }

    return false;
}
//----------------------------------------------------------------------------
bool solveCubic(double a, double b, double c, double d, double &root)
{
    std::cout << "a : " << a << "b : " << b << "c : " << c << "d : " << d <<std::endl; 
    if(a == 0.0 || abs(a/b) < 1.0e-6)
    {
        std::cout << "ignoring coeff3!!\n";
        return solveQuadratic(b, c, d, root);
    }
    double B = b/a, C = c/a, D = d/a;

    std::cout << "B : " << B << "C : " << C << "D : " << D << std::endl;
    double Q = (B*B - C*3.0)/9.0, QQQ = Q*Q*Q;
    double R = (2.0*B*B*B - 9.0*B*C + 27.0*D)/54.0, RR = R*R;
    std::cout << "B : " << B << "Q : " << Q << "R : " << R << std::endl;
    // 3 real roots
    if(RR<QQQ)
    {
        /* This sqrt and division is safe, since RR >= 0, so QQQ > RR,    */
        /* so QQQ > 0.  The acos is also safe, since RR/QQQ < 1, and      */
        /* thus R/sqrt(QQQ) < 1.                                     */
        double theta = acos(R/sqrt(QQQ));
        /* This sqrt is safe, since QQQ >= 0, and thus Q >= 0             */
        double r1, r2, r3;
        r1 = r2 = r3 = -2.0*sqrt(Q);
        r1 *= cos(theta/3.0);
        r2 *= cos((theta+2*PI)/3.0);
        r3 *= cos((theta-2*PI)/3.0);

        r1 -= B/3.0;
        r2 -= B/3.0;
        r3 -= B/3.0; 

        root = 1000000.0;

        if(r1 >= 0.0) root = r1;
        if(r2 >= 0.0 && r2 < root) root = r2;
        if(r3 >= 0.0 && r3 < root) root = r3;

        return true;
    }
    // 1 real root
    else
    {
        double A2 = -pow(fabs(R)+sqrt(RR-QQQ),1.0/3.0);
        if (A2!=0.0) {
            if (R<0.0) A2 = -A2; 
            root = A2 + Q/A2; 
        }
        root -= B/3.0;
        std::cout << "3root : " << root << std::endl;
        return true;
    }
}
//----------------------------------------------------------------------------

bool solveQuadraticOther(double a, double b, double c, double &root)
{
    if(a == 0.0 || abs(a/b) < 1.0e-4)
    {
        std::cout << "ignoring coeff4!!\n";
        if(abs(b) < 1.0e-4) 
        {std::cout << "ignoring coeff5!!\n";
            return false;
        }
        else
        {
            root = -c/b;
            return true;
        }
    }

    double discriminant = b * b - 4.0 * a * c;
    if(discriminant >= 0.0)
    {
        discriminant = sqrt(discriminant);
        root = (b - discriminant) * -0.5 / a;
        return true;
    }

    return false;
}
//----------------------------------------------------------------------------
/*
bool solveQuartic(double a, double b, double c, double d, double e, double &root)
{
    // I switched to this method, and it seems to be more numerically stable.
    // http://www.gamedev.net/community/forums/topic.asp?topic_id=451048 

    // When a or (a and b) are magnitudes of order smaller than C,D,E
    // just ignore them entirely. This seems to happen because of numerical
    // inaccuracies of the line-circle algorithm. I wanted a robust solver,
    // so I put the fix here instead of there.
    if(a == 0.0 || abs(a/b) < 1.0e-5 || abs(a/c) < 1.0e-5 || abs(a/d) < 1e-9)
    {

        std::cout << "ignoring coeff6!!\n";
        return solveCubic(b, c, d, e, root);
    }
    double B = b/a, C = c/a, D = d/a, E = e/a;
    double BB = B*B;
    double I = -3.0*BB*0.125 + C;
    double J = BB*B*0.125 - B*C*0.5 + D;
    double K = -3*BB*BB/256.0 + C*BB/16.0 - B*D*0.25 + E;

    double z;
    bool foundRoot2 = false, foundRoot3 = false, foundRoot4 = false, foundRoot5 = false;
    if(solveCubic(1.0, I+I, I*I - 4*K, -(J*J), z))
    {
        double value = z*z*z + z*z*(I+I) + z*(I*I - 4*K) - J*J;

        double p = sqrt(z);
        double r = -p;
        double q = (I + z - J/p)*0.5;
        double s = (I + z + J/p)*0.5;

        bool foundRoot = false, foundARoot;
        double aRoot;
        foundRoot = solveQuadratic(1.0, p, q, root);
        root -= B/4.0;

        foundARoot = solveQuadratic(1.0, r, s, aRoot);
        aRoot -= B/4.0;
        if((foundRoot && foundARoot && ((aRoot < root && aRoot >= 0.0) 
                        || root < 0.0)) || (!foundRoot && foundARoot)) 
        {
            root = aRoot;
            foundRoot = true;
        }

        foundARoot = solveQuadraticOther(1.0, p, q, aRoot);
        aRoot -= B/4.0;
        if((foundRoot && foundARoot && ((aRoot < root && aRoot >= 0.0) 
                        || root < 0.0)) || (!foundRoot && foundARoot)) 
        {
            root = aRoot;
            foundRoot = true;
        }

        foundARoot = solveQuadraticOther(1.0, r, s, aRoot);
        aRoot -= B/4.0;
        if((foundRoot && foundARoot && ((aRoot < root && aRoot >= 0.0) 
                        || root < 0.0)) || (!foundRoot && foundARoot)) 
        {
            root = aRoot;
            foundRoot = true;
        }
        return foundRoot;
    }
    return false;
}
*/

bool solveQuartic(double a,double b,double c,double d,double e,double& root)
{

    cout.precision(10);
    if(abs(a) < 1e-12)
        return solveCubic(b,c,d,e,root);
    std::cout << "yo" << std::endl;
    b /= a;
    c /= a;
    d /= a;
    e /= a;
    std::cout << "a1 : " << b << "a2 : " << c << "a3 : " << d << "a4 : " << e << std::endl; 
    double coeff2 = -c;
    double coeff1 = b*d - 4*e;
    double coeff0 = 4*c*e - d*d - b*b*e;
    double root1 = 0;
    double finalroot = -1;
    std::cout << "quart to cubic : 1, " << coeff2 << ", " << coeff1 << ", " << coeff0 << std::endl;
    bool solving = solveCubic(1,coeff2,coeff1,coeff0,root1);
    std::cout << "quart to cubic : 1, " << coeff2 << ", " << coeff1 << ", " << coeff0 << std::endl; 
    std::cout << "tExists : " << solving << "root1 : " << root1 << std::endl;
    if(solving)
    {
        double root2 = 0;
        double S1 = (b*b - 4*c + 4*root1);
        double S2 = root1*root1 - 4*e;
        std::cout << "S1 : " << S1 << std::endl;
        std::cout << "S2 : " << S2 << std::endl;
        if(abs(S1) < 1e-6)
            S1 = 0;
        double qcf1 = 0.5*(b + sqrt(S1));
        double qcf0 = 0.5*(root1 + sqrt(S2));
        solving = solveQuadratic(1,qcf1,qcf0,root2);
        std::cout << "qcf1 :"<< qcf1 << ", " << "qcf0 : "<< qcf0 << std::endl; 
    std::cout << "root2 : " << root2 << std::endl;
        if(solving && root2 > 1e-13)
            finalroot = root2;
        qcf1 = 0.5*(b - sqrt(S1));
        qcf0 = 0.5*(root1 - sqrt(S2));
        std::cout << "qcf1 :"<< qcf1 << ", " << "qcf0 : "<< qcf0 << std::endl; 
        solving = solveQuadratic(1,qcf1,qcf0,root2);
    std::cout << "root2 : " << root2 << std::endl;
        if(solving)
        {
            if(root2  < finalroot && root2 > 1e-13)
                finalroot = root2;
        }
        if(finalroot >1e-13)
        {
            root = finalroot;
            std::cout << "finalroot : " << finalroot << std::endl;
            return true;
        }
    }
    return false;
}

#endif //_SOLVE_POLY_HPP_

#include <iostream>
#include <cmath>

/**
 * functions solveCubic and solveQuartic to solve 4th order polynomials 
 * from :
 * https://stackoverflow.com/questions/11837078/initialize-a-constant-sized-array-in-an-initializer-list
*/

const double PI = 3.14159265358979323846;
using namespace std;

//----------------------------------------------------------------------------
bool solveQuadratic(double a, double b, double c, double& root)
{
    if(a == 0.0 || abs(a/b) < 1.0e-6)
    {
        if(abs(b) < 1.0e-6) 
            return false;
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

bool solveQuadraticOther(double a, double b, double c, double &root)
{
    if(a == 0.0 || abs(a/b) < 1.0e-6)
    {
        if(abs(b) < 1.0e-4) 
            return false;
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

bool solveCubic(double a, double b, double c, double d, double &root)
{
    if(a == 0.0 || abs(a/b) < 1.0e-6)
        return solveQuadratic(b, c, d, root);

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
        return true;
    }
}


bool solveQuartic(double a, double b, double c, double d, double e, double &root)
{
    // I switched to this method, and it seems to be more numerically stable.
    // http://www.gamedev.n...topic_id=451048 

    // When a or (a and b) are magnitudes of order smaller than C,D,E
    // just ignore them entirely. This seems to happen because of numerical
    // inaccuracies of the line-circle algorithm. I wanted a robust solver,
    // so I put the fix here instead of there.
    if(a == 0.0 || abs(a/b) < 1.0e-6 || abs(a/c) < 1.0e-6 || abs(a/d) < 1.0e-6)
        return solveCubic(b, c, d, e, root);

    double B = b/a, C = c/a, D = d/a, E = e/a;
    double BB = B*B;
    double I = -3.0*BB*0.125 + C;
    double J = BB*B*0.125 - B*C*0.5 + D;
    double K = -3*BB*BB/256.0 + C*BB/16.0 - B*D*0.25 + E;

    double z;
    bool foundRoot2 = false, foundRoot3 = false, foundRoot4 = false, foundRoot5 = false;
    double one = 1;
    if(solveCubic(one, I+I, I*I - 4*K, -(J*J), z))
    {
        double value = z*z*z + z*z*(I+I) + z*(I*I - 4*K) - J*J;

        double p = sqrt(z);
        double r = -p;
        double q = (I + z - J/p)*0.5;
        double s = (I + z + J/p)*0.5;

        bool foundRoot = false, foundARoot;
        double aRoot;
        foundRoot = solveQuadratic(one, p, q, root);
        root -= B/4.0;

        foundARoot = solveQuadratic(one, r, s, aRoot);
        aRoot -= B/4.0;
        if((foundRoot && foundARoot && ((aRoot < root && aRoot >= 0.0) 
                        || root < 0.0)) || (!foundRoot && foundARoot)) 
        {
            root = aRoot;
            foundRoot = true;
        }

        foundARoot = solveQuadraticOther(one, p, q, aRoot);
        aRoot -= B/4.0;
        if((foundRoot && foundARoot && ((aRoot < root && aRoot >= 0.0) 
                        || root < 0.0)) || (!foundRoot && foundARoot)) 
        {
            root = aRoot;
            foundRoot = true;
        }

        foundARoot = solveQuadraticOther(one, r, s, aRoot);
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

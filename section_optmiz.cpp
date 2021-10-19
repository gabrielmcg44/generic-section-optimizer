
// Header Files
#include <iostream>
#include<conio.h>
#include<vector>
#include<cmath>
#include<limits>
#include<math.h>
using namespace std;


// Default units are in: 
// Force: kN
// Distance: cm

// Strain is in ‰



struct point {
    long double x;
    long double y;
    long double steelArea;
    long double linearSteelThickness;
    struct point *nextPoint;
};

// class point {
//     private:
//     public:
//         long double x;
//         long double y;
//         long double steelArea;
//         long double linearSteelThickness;
//         point nextPoint;
// }


struct stresses {
    long double x = 0;
    long double y = 0;
    long double normalStress;
    long double xBendingMoment;
    long double yBendingMoment;
    bool safe;
};


struct strainState {
    long double x = 0.0;
    long double y = 0.0;
    long double e0; // strain in the origin
    long double kx; // variation of strain per distance from x axis
    long double ky; // variation of strain per distance from y axis
    bool safe;
};


class polygon {
private:
public:
    point points[200];
    int numPoints = 0;

    void pushPoint(point newPoint){
        points[numPoints] = newPoint;
        numPoints++;
    }

    long double minX() {
        long double minimum = __FLT_MAX__;
        for (int i = 0; i < numPoints; i++) {
            if (points[i].x < minimum) {
                minimum = points[i].x;
            }
        }
        return minimum;
    }

    long double maxX() {
        long double maximum = - __FLT_MAX__;
        for (int i = 0; i < numPoints; i++) {
            if (points[i].x > maximum) {
                maximum = points[i].x;
            }
        }
        return maximum;
    }

    long double minY() {
        long double minimum = __FLT_MAX__;
        for (int i = 0; i < numPoints; i++) {
            if (points[i].y < minimum) {
                minimum = points[i].y;
            }
        }
        return minimum;
    }

    long double maxY() {
        long double maximum = - __FLT_MAX__;
        for (int i = 0; i < numPoints; i++) {
            if (points[i].y > maximum) {
                maximum = points[i].y;
            }
        }
        return maximum;
    }
};


long double determinant(long double** matrix, int n, bool replaceColumn = false, long double* newColumn = new long double[3] , int newColumnNumber = 0) {
    long double det = 0;

    long double** submatrix;
    long double** newmatrix;
    submatrix = new long double*[3];
    newmatrix = new long double*[3];
    for (int i = 0; i < 3; i++) {
        submatrix[i] = new long double[3];
        newmatrix[i] = new long double[3];
    }

    for (int x = 0; x < n; x++) {
        for (int y = 0; y < n; y++) {
            if(y==newColumnNumber && replaceColumn) {
                newmatrix[x][y] = newColumn[x];
            } else {
                newmatrix[x][y] = matrix[x][y];
            }
        }
    }

    if (n == 2)
    return ((newmatrix[0][0] * newmatrix[1][1]) - (newmatrix[1][0] * newmatrix[0][1]));
    else {
        for (int x = 0; x < n; x++) {
            int subi = 0;
            for (int i = 1; i < n; i++) {
                int subj = 0;
                for (int j = 0; j < n; j++) {
                    if (j == x) {
                        continue;
                    }
                    submatrix[subi][subj] = newmatrix[i][j];
                    subj++;
                }
                subi++;
            }
            det = det + (pow(-1, x) * newmatrix[0][x] * determinant( submatrix, n - 1 ));
        }
    }
    return det;
}


long double featureYieldDesign(long double fyk) {
    // featureYieldDesign 'fyd'
    return fyk / 1.15;
}


long double strainYieldDesign(long double fyk, long double E) {
    // strainYieldDesign 'epsilon yd'
    return 1000.0 * featureYieldDesign(fyk) / E;
}


long double stressCompressDesign(long double fck) {
    // stressCompressDesign 'sigma cd'
    return (0.85 / 1.4) * fck;
}


long double ec2(long double fck) {
    // 
    if (fck >= 0.0 && fck <= 5.0) {
        return 2.0;
    } else if (fck > 5.0 && fck <= 9.0) {
        return 2.0 + 0.085 * pow(fck - 5.0, 0.53);
    } else {
        return -1.0;
    }   
}


long double ecu(long double fck) {
    // 
    if (fck >= 0 && fck <= 5) {
        return 3.5;
    } else if (fck > 5.0 && fck <= 9.0) {
        return 2.6 + 35.0 * pow((9.0 - fck) / 10.0, 4);
    } else {
        return -1.0;
    }   
}


long double n(long double fck) {
    //  
    if (fck >= 0.0 && fck <= 5.0) {
        return 2.0;
    } else if (fck > 5.0 && fck <= 9.0) {
        return 1.4 + 23.4 * pow((9.0 - fck) / 10.0, 4.0);
    } else {
        return -1.0;
    }   
}


long double sigmac(long double eps, long double sigmacd, long double n, long double ec2) {
    // Calculate tension of concrete given strain
    if (eps <= 0.0) {
        return 0.0;
    } else if (eps > 0.0 && eps <= 2.0) {
        return sigmacd * (1.0 - pow(1.0 -(eps / ec2), n));
    } else {
        return sigmacd;
    }   
}


long double sigmas(long double eps, long double fyk, long double E) {
    // Calculate tension of concrete given strain
    long double eyd = strainYieldDesign(fyk, E);
    long double fyd = featureYieldDesign(fyk);
    if (eps <= - eyd) {
        return - fyd;
    } else if (eps > -eyd && eps <= eyd) {
        return fyd * (eps / eyd);
    } else {
        return fyd;
    }   
}


long double tk(long double k, long double fck, long double eps) {
    long double c = ec2(fck);
    long double n1 = n(fck);
    // cout << n1 << 'n' << '\n';
    // cout << n1 << ' ' << k << ' ' << eps << ' ' << c << ' ' << '\n';
    
    if (eps >= 0.0 && eps < c) {
        // cout << 1 - eps/c << 't' << '\n';
        // cout << 'c' << '1' << '\n';
        return (1.0/(n1 + k)) * (1.0 - pow(1 - eps/c, n1 + k));
    } else {
        // cout << 'c' << '2' << '\n';
        return 1.0/(n1 + k);
    }
}


long double I1(long double eps, long double fck) {
    long double c = ec2(fck);
    long double sigmacd = stressCompressDesign(fck);
    // cout << tk(1, fck, eps) << 't' << '\n';
    if (eps < 0.0) {
        return 0.0;
    } else {
        return (eps - c * tk(1.0, fck, eps)) * sigmacd;
    }
}


long double I2(long double eps, long double fck) {
    long double c = ec2(fck);
    long double n1 = n(fck);
    long double sigmacd = stressCompressDesign(fck);
    if (eps < 0.0) {
        return 0.0;
    } else {
        return (pow(eps, 2)/2.0 - (c/(n1 + 1.0)) * (eps - c*tk(2.0, fck, eps))) * sigmacd;
    }
}


long double I3(long double eps, long double fck) {
    long double c = ec2(fck);
    long double n1 = n(fck);
    long double sigmacd = stressCompressDesign(fck);
    if (eps < 0) {
        return 0;
    } else {
        return ((pow(eps,3))/6.0-(c/(n1+1.0))*((pow(eps,2))/2.0-(c/(n1+2.0))*(eps-c*tk(3.0,fck,eps))))*sigmacd;
    }
}


long double J1(long double eps, long double fck) {
    long double c = ec2(fck);
    long double sigmacd = stressCompressDesign(fck);
    if (eps < 0) {
        return 0;
    } else {
        return ((pow(eps,2)/2.0)-(pow(c,2))*(tk(1.0,fck,eps)-tk(2.0,fck,eps)))*sigmacd;
    }
}


long double J2(long double eps, long double fck) {
    long double c = ec2(fck);
    long double sigmacd = stressCompressDesign(fck);
    if (eps < 0.0) {
        return 0.0;
    } else {
        //cout.precision(25);
        //cout << tk(1,fck,eps) << ' ' << tk(2,fck,eps) << ' ' << tk(3,fck,eps) << ' ' << eps << ' ' << c << ' ' << sigmacd << 't' << '\n';
        //cout << ((pow(eps,3)/3)-(pow(c,3))*(tk(1,fck,eps)-2*tk(2,fck,eps)+tk(3,fck,eps)))*sigmacd << 't' << 'r' << '\n';//-(pow(c,3))*(tk(1,fck,eps)-2*tk(2,fck,eps)+tk(3,fck,eps)))*sigmacd << 't' << 'r' << '\n';
        return ((pow(eps,3)/3.0)-(pow(c,3))*(tk(1.0,fck,eps)-2.0*tk(2.0,fck,eps)+tk(3.0,fck,eps)))*sigmacd;
    }
}


long double K1(long double eps, long double fck) {
    long double c = ec2(fck);
    long double n1 = n(fck);
    long double sigmacd = stressCompressDesign(fck);
    if (eps < 0) {
        return 0;
    } else {
        return ((pow(eps,3)/3.0)-(c/(n1+1.0))*(pow(eps,2)/2.0-(pow(c,2))*(tk(2.0,fck,eps)-tk(3.0,fck,eps))))*sigmacd;
    }
}


long double f1i(long double epsIplus1, long double epsI, long double fck, long double r) {
    if (abs(epsIplus1 - epsI) < r) {
        return I1(epsI,fck);
    } else {
        return (I2(epsIplus1,fck)-I2(epsI,fck)) / (epsIplus1 - epsI);
    }
}


long double f2i(long double epsIplus1, long double epsI, long double fck, long double yI, long double yIplus1, long double r) {
    long double gI = yI*epsIplus1 - yIplus1*epsI;
    if (abs(epsIplus1 - epsI) < r) {
        return I1(epsI,fck)*((yI+yIplus1)/2);
    } else {
        return (gI*(I2(epsIplus1,fck)-I2(epsI,fck))+(yIplus1-yI)*(K1(epsIplus1,fck)-K1(epsI,fck)))/(pow(epsIplus1-epsI,2));
    }
}


long double f3i(long double epsIplus1, long double epsI, long double fck, long double yI, long double yIplus1, long double kx, long double r) {
    if (abs(epsIplus1 - epsI) < r) {
        return f2i(epsIplus1,epsI,fck,yI,yIplus1,r)+(I2(epsI,fck)/kx);
    } else {
        return f2i(epsIplus1,epsI,fck,yI,yIplus1,r)+(((I3(epsIplus1,fck)-I3(epsI,fck))/(epsIplus1-epsI))/kx);
    }
}


long double f4i(long double epsIplus1, long double epsI, long double fck, long double xI, long double xIplus1, long double r) {
    long double hI=xI*epsIplus1-xIplus1*epsI;
    if (abs(epsIplus1-epsI)<r) {
        return I1(epsI,fck)*((xI+xIplus1)/2.0);
    } else {
        return (hI*(I2(epsIplus1,fck)-I2(epsI,fck))+(xIplus1-xI)*(K1(epsIplus1,fck)-K1(epsI,fck)))/(pow(epsIplus1-epsI, 2));
    }
}


long double f5i(long double epsIplus1, long double epsI, long double fck, long double xI, long double xIplus1, long double ky, long double r) {
    if (abs(epsIplus1-epsI)<r) {
        return f4i(epsIplus1,epsI,fck,xI,xIplus1,r)-(I2(epsI,fck)/ky);
    } else {
        return f4i(epsIplus1,epsI,fck,xI,xIplus1,r)+(((I3(epsIplus1,fck)-I3(epsI,fck))/(epsIplus1-epsI))/ky);
    }
}


long double f6i(long double epsIplus1, long double epsI, long double fck, long double r) {
    long double sigmacd = stressCompressDesign(fck);    
    if (abs(epsIplus1-epsI)<r) {
        // cout << 1 << '\n';
        return sigmac(epsI,sigmacd,n(fck),ec2(fck));
    } else {
        // cout << 2 << '\n';
        // cout << I1(epsIplus1,fck) << 'I' << '1' << '\n';
        return (I1(epsIplus1,fck)-I1(epsI,fck))/(epsIplus1-epsI);
    }
}


long double f7i(long double epsIplus1, long double epsI, long double fck, long double yI, long double yIplus1, long double r) {
    long double sigmacd = stressCompressDesign(fck);   
    long double gI = yI*epsIplus1 - yIplus1*epsI; 
    if (abs(epsIplus1-epsI)<r) {
        return sigmac(epsI,sigmacd,n(fck),ec2(fck))*((yI+yIplus1)/2.0);
    } else {
        return (gI*(I1(epsIplus1,fck)-I1(epsI,fck))+(yIplus1-yI)*(J1(epsIplus1,fck)-J1(epsI,fck)))/(pow(epsIplus1-epsI,2));
    }
}


long double f8i(long double epsIplus1, long double epsI, long double fck, long double xI, long double xIplus1, long double r) {
    long double sigmacd = stressCompressDesign(fck);   
    long double hI = xI*epsIplus1 - xIplus1*epsI; 
    if (abs(epsIplus1-epsI)<r) {
        return sigmac(epsI,sigmacd,n(fck),ec2(fck))*((xI+xIplus1)/2.0);
    } else {
        return (hI*(I1(epsIplus1,fck)-I1(epsI,fck))+(xIplus1-xI)*(J1(epsIplus1,fck)-J1(epsI,fck)))/(pow(epsIplus1-epsI,2));
    }
}


long double f9i(long double epsIplus1, long double epsI, long double fck, long double yI, long double yIplus1, long double r) {
    long double sigmacd = stressCompressDesign(fck);   
    long double gI = yI*epsIplus1 - yIplus1*epsI; 
    if (abs(epsIplus1-epsI)<r) {
        return sigmac(epsI,sigmacd,n(fck),ec2(fck))*((pow(yI,2)+yI*yIplus1+pow(yIplus1,2))/3.0);
    } else {
        return ((pow(gI,2))*(I1(epsIplus1,fck)-I1(epsI,fck))+2.0*gI*(yIplus1-yI)*(J1(epsIplus1,fck)-J1(epsI,fck))+(pow(yIplus1-yI,2))*(J2(epsIplus1,fck)-J2(epsI,fck)))/(pow(epsIplus1-epsI,3));
    }
}


long double f10i(long double epsIplus1, long double epsI, long double fck, long double yI, long double yIplus1, long double xI, long double xIplus1, long double r) {
    long double gI = yI*epsIplus1 - yIplus1*epsI;
    long double hI = xI*epsIplus1 - xIplus1*epsI;
    long double sigmacd = stressCompressDesign(fck);  
    long double sigmai = sigmac(epsI, sigmacd, n(fck), ec2(fck));
    long double deltaI1 = (I1(epsIplus1,fck)-I1(epsI,fck));
    long double deltaJ1 = (J1(epsIplus1,fck)-J1(epsI,fck));
    long double deltaJ2 = (J2(epsIplus1,fck)-J2(epsI,fck));
    long double deltaY = yIplus1 - yI;
    long double deltaX = xIplus1 - xI;
    long double deltaEps = (epsIplus1 - epsI);
    if (abs(epsIplus1-epsI)<r) {
        return sigmai*((xI*yIplus1+2*(xI*yI+xIplus1*yIplus1)+xIplus1*yI)/6.0);
    } else {
        return (hI*gI*deltaI1+deltaJ1*(gI*deltaX+hI*deltaY)+deltaX*deltaY*deltaJ2)/(pow(deltaEps,3));
    }
}


long double f11i(long double epsIplus1, long double epsI, long double fck, long double xI, long double xIplus1, long double r) {
    long double hI=xI*epsIplus1-xIplus1*epsI;
    long double sigmacd = stressCompressDesign(fck);  
    if (abs(epsIplus1-epsI)<r) {
        //cout << sigmac(epsI,sigmacd,n(fck),ec2(fck))*((pow(xI,2)+xI*xIplus1+pow(xIplus1,2))/3) << 'c' << '1' << '\n';
        return sigmac(epsI,sigmacd,n(fck),ec2(fck))*((pow(xI,2)+xI*xIplus1+pow(xIplus1,2))/3.0);
    } else {
        //cout << hI << ' ' << I1(epsIplus1,fck) << ' ' << I1(epsI,fck) << ' ' << epsIplus1  << ' ' << fck  << ' ' << epsI  << ' ' <<  xIplus1  << ' ' << xI  << ' ' << J1(epsIplus1,fck) << ' ' <<  J1(epsI,fck) << ' ' << J2(epsIplus1,fck) << ' ' <<  J2(epsI,fck) << 'C' << '2' << '\n';
        //cout << ((pow(hI,2))*(I1(epsIplus1,fck)-I1(epsI,fck))+2*hI*(xIplus1-xI)*(J1(epsIplus1,fck)-J1(epsI,fck))+(pow(xIplus1-xI,2))*(J2(epsIplus1,fck)-J2(epsI,fck)))/(pow(epsIplus1-epsI,3)) << 'R' << '2' << '\n';
        return ((pow(hI,2))*(I1(epsIplus1,fck)-I1(epsI,fck))+2*hI*(xIplus1-xI)*(J1(epsIplus1,fck)-J1(epsI,fck))+(pow(xIplus1-xI,2))*(J2(epsIplus1,fck)-J2(epsI,fck)))/(pow(epsIplus1-epsI,3));
    }
}


long double S1(long double eps, long double eyd, long double fyd) {
    if (abs(eps) < eyd) {
        return (1.0/2.0)*(fyd/eyd)*pow(eps,2);
    } else {
        long double result = abs(eps*fyd) - eyd*fyd/2.0;
        return abs(eps*fyd) - eyd*fyd/2.0;
    }
}


long double KS1(long double eps, long double eyd, long double fyd) {
    if (eps < -eyd) {
        return -(1.0/2.0)*(fyd)*pow(eps,2)+(fyd*pow(eyd,2))/6.0;
    } else if (eps >= -eyd && eyd > eps){
        return (pow(eps,3)/3.0)*(fyd/eyd);
    } else {
        return (1.0/2.0)*(fyd)*pow(eps,2)-(fyd*pow(eyd,2))/6.0;
    }
}


long double Dc(long double eps, long double fck) {
    long double c = ec2(fck);
    long double nf = n(fck);
    long double sigmacd = stressCompressDesign(fck); 
    if (eps<c && eps >= 0.0) {
        return sigmacd * (nf / c) * pow(1 - eps / c, nf - 1.0);
    } else {
        return 0.0;
    }
}


long double Ds(long double eps, long double fyk, long double E) {
    long double eyd = strainYieldDesign(fyk, E); 
    if (eps < eyd && eps >= -eyd) {
        return E / 1000.0;
    } else {
        return 0.0;
    }
}


long double NsEps0(long double lambda, long double DsI, long double dist, long double deltaEps, long double deltaSigmas, long double r) {
    if (abs(deltaEps) < r) {
        return lambda * dist * DsI;
    } else {
        return lambda * (dist / deltaEps) * deltaSigmas;
    } 
}


long double NsKx(long double lambda, long double DsI, long double dist, long double deltaEps, long double deltaSigmas, long double gi, long double deltaY, long double deltaSigmaEps, long double deltaS1, long double Yi, long double r) {
    if (abs(deltaEps) < r) {
        return - lambda * dist * DsI * (Yi + deltaY / 2.0);
    } else {
        return - lambda * (dist / pow(deltaEps, 2)) * (gi * deltaSigmas + deltaY * deltaSigmaEps - deltaY * deltaS1);
    }
}


long double NsKy(long double lambda, long double DsI, long double dist, long double deltaEps, long double deltaSigmas, long double hi, long double deltaX, long double deltaSigmaEps, long double deltaS1, long double Xi, long double r) {
    if (abs(deltaEps) < r) {
        return lambda * dist * DsI * (Xi + deltaX / 2.0);
    } else {
        return lambda * (dist / pow(deltaEps, 2)) * (hi * deltaSigmas + deltaX * deltaSigmaEps - deltaX * deltaS1);
    }
}


long double MsxKx (long double lambda, long double DsI, long double dist, long double deltaEps, long double deltaSigmas, long double gi, long double deltaY, long double deltaSigmaEps, long double deltaSigmaEps2, long double deltaS1, long double deltaKS1, long double Yi, long double r) {
    if (abs(deltaEps) < r) {
        return lambda * dist * DsI * (pow(Yi, 2) + Yi * deltaY + pow(deltaY, 2)/3.0);
    } else {
        return lambda * (dist/pow(deltaEps,3)) * (pow(gi, 2) * deltaSigmas + 2.0 * gi * deltaY * (deltaSigmaEps - deltaS1) + pow(deltaY, 2) * (deltaSigmaEps2 - 2.0 * deltaKS1));
    }
}


long double MsyKy (long double lambda, long double DsI, long double dist, long double deltaEps, long double deltaSigmas, long double hi, long double deltaX, long double deltaSigmaEps, long double deltaSigmaEps2, long double deltaS1, long double deltaKS1, long double Xi, long double r) {
    if (abs(deltaEps) < r) {
        return lambda * dist * DsI * (pow(Xi, 2) + Xi * deltaX + pow(deltaX, 2)/3.0);
    } else {
        return lambda * (dist/pow(deltaEps,3)) * (pow(hi, 2) * deltaSigmas + 2.0 * hi * deltaX * (deltaSigmaEps - deltaS1) + pow(deltaX, 2) * (deltaSigmaEps2 - 2.0 * deltaKS1));
    }
}


long double MsxKy (long double lambda, long double DsI, long double dist, long double deltaEps, long double deltaSigmas, long double gi, long double hi, long double deltaX, long double deltaY, long double deltaSigmaEps, long double deltaSigmaEps2, long double deltaS1, long double deltaKS1, long double Xi, long double Yi, long double r) {
    if (abs(deltaEps) < r) {
        return -lambda * dist * DsI * (Xi*Yi + (Yi*deltaX + Xi*deltaY)/2.0 + deltaX*deltaY/3.0);
    } else {
        return -lambda * (dist/pow(deltaEps, 3)) * (gi*hi*deltaSigmas + (gi*deltaX + hi*deltaY) * (deltaSigmaEps - deltaS1) + deltaX*deltaY*(deltaSigmaEps2 - 2.0*deltaKS1));
    }
}


class crossSection {
private:
    polygon concreteSection;
    polygon steelLinearDist;

    //Access - Specifier
public:
    long double fyk = 50.0; //featureYieldKnown
    long double fck = 4.0; //featureCompressionKnown
    long double E = 21000.0; //youngModulus
    long double r = 0.00001; //tolerance
    long double p = 0.00000001;
    long double t = 0.000000000001;
    //Variable Declaration
    void defineConcreteSection(polygon section){
        concreteSection = section;
    }
    
    void defineSteelLinearDist(polygon linearDist) {
        steelLinearDist = linearDist;
    }

    void showCrossSection() {
        for (int i = 0; i < concreteSection.numPoints; i++) {
            cout << concreteSection.points[i].x << ',' << concreteSection.points[i].y << '\n';
        }
    } 

    long double Ac() {
        long double soma = 0;
        for (int i = 0; i < concreteSection.numPoints; i++) {
            int Iplus1;
            if (i != concreteSection.numPoints - 1) {
                Iplus1 = i + 1;
            } else {
                Iplus1 = 0;
            }
            long double xI = concreteSection.points[i].x;
            long double xIplus1 = concreteSection.points[Iplus1].x;
            long double yI = concreteSection.points[i].y;
            long double yIplus1 = concreteSection.points[Iplus1].y;          
            long double aI = xI*yIplus1 - xIplus1*yI;
            soma += aI;
        }
        return soma / 2.0;
    }
    
    long double Sy() {
        long double soma = 0.0;
        for (int i = 0; i < concreteSection.numPoints; i++) {
            int Iplus1;
            if (i != concreteSection.numPoints - 1) {
                Iplus1 = i + 1;
            } else {
                Iplus1 = 0;
            }
            long double xI = concreteSection.points[i].x;
            long double xIplus1 = concreteSection.points[Iplus1].x;
            long double yI = concreteSection.points[i].y;
            long double yIplus1 = concreteSection.points[Iplus1].y;          
            long double aI = xI*yIplus1 - xIplus1*yI;
            soma += aI * (xI + xIplus1);
        }
        return soma / 6.0;
    }
    
    long double Sx() {
        long double soma = 0.0;
        for (int i = 0; i < concreteSection.numPoints; i++) {
            int Iplus1;
            if (i != concreteSection.numPoints - 1) {
                Iplus1 = i + 1;
            } else {
                Iplus1 = 0;
            }
            long double xI = concreteSection.points[i].x;
            long double xIplus1 = concreteSection.points[Iplus1].x;
            long double yI = concreteSection.points[i].y;
            long double yIplus1 = concreteSection.points[Iplus1].y;          
            long double aI = xI*yIplus1 - xIplus1*yI;
            soma += aI * (yI + yIplus1);
        }
        return soma / 6.0;
    }
    
    long double Iyy() {
        long double soma = 0.0;
        for (int i = 0; i < concreteSection.numPoints; i++) {
            int Iplus1;
            if (i != concreteSection.numPoints - 1) {
                Iplus1 = i + 1;
            } else {
                Iplus1 = 0;
            }
            long double xI = concreteSection.points[i].x;
            long double xIplus1 = concreteSection.points[Iplus1].x;
            long double yI = concreteSection.points[i].y;
            long double yIplus1 = concreteSection.points[Iplus1].y;          
            long double aI = xI*yIplus1 - xIplus1*yI;
            soma += aI * (pow(xI,2) +  xI*xIplus1 + pow(xIplus1,2));
        }
        return soma / 12.0;
    } 
    
    long double Ixy() {
        long double soma = 0.0;
        for (int i = 0; i < concreteSection.numPoints; i++) {
            int Iplus1;
            if (i != concreteSection.numPoints - 1) {
                Iplus1 = i + 1;
            } else {
                Iplus1 = 0;
            }
            long double xI = concreteSection.points[i].x;
            long double xIplus1 = concreteSection.points[Iplus1].x;
            long double yI = concreteSection.points[i].y;
            long double yIplus1 = concreteSection.points[Iplus1].y;          
            long double aI = xI*yIplus1 - xIplus1*yI;
            soma += aI * (xI*yIplus1+2.0*(xI*yI+xIplus1*yIplus1)+xIplus1*yI);
        }
        return soma / 24.0;
    }

    long double Ixx() {
        long double soma = 0.0;
        for (int i = 0; i < concreteSection.numPoints; i++) {
            int Iplus1;
            if (i != concreteSection.numPoints - 1) {
                Iplus1 = i + 1;
            } else {
                Iplus1 = 0;
            }
            long double xI = concreteSection.points[i].x;
            long double xIplus1 = concreteSection.points[Iplus1].x;
            long double yI = concreteSection.points[i].y;
            long double yIplus1 = concreteSection.points[Iplus1].y;          
            long double aI = xI*yIplus1 - xIplus1*yI;
            soma += aI * (pow(yI,2) +  yI*yIplus1 + pow(yIplus1,2));
        }
        return soma / 12.0;
    }

    long double NcFOC(strainState e0kxky) {
        long double sigmacd = stressCompressDesign(fck); 
        long double e0 = e0kxky.e0;
        long double kx = e0kxky.kx;
        long double ky = e0kxky.ky;
        long double soma = 0.0;
        if (abs(kx)<r && abs(ky)<r) {
            return sigmac(e0,sigmacd,n(fck),ec2(fck))*Ac();
        } else if (abs(kx) >= r) {
            for (int i = 0; i < concreteSection.numPoints; i++) {
                int Iplus1;
                if (i != concreteSection.numPoints - 1) {
                    Iplus1 = i + 1;
                } else {
                    Iplus1 = 0;
                }
                long double xI = concreteSection.points[i].x;
                long double xIplus1 = concreteSection.points[Iplus1].x;
                long double yI = concreteSection.points[i].y;
                long double yIplus1 = concreteSection.points[Iplus1].y;
                long double epsI = e0 + ky*xI - kx*yI;
                long double epsIplus1 = e0 + ky*xIplus1 - kx*yIplus1;    
                soma += (xIplus1 - xI) * f1i(epsIplus1, epsI, fck, r);

            }
            return (1.0 / kx) * soma;
        } else {
            for (int i = 0; i < concreteSection.numPoints; i++) {
                int Iplus1;
                if (i != concreteSection.numPoints - 1) {
                    Iplus1 = i + 1;
                } else {
                    Iplus1 = 0;
                }
                long double xI = concreteSection.points[i].x;
                long double xIplus1 = concreteSection.points[Iplus1].x;
                long double yI = concreteSection.points[i].y;
                long double yIplus1 = concreteSection.points[Iplus1].y;
                long double epsI = e0 + ky*xI - kx*yI;
                long double epsIplus1 = e0 + ky*xIplus1 - kx*yIplus1;                 
                soma += (yIplus1 - yI) * f1i(epsIplus1, epsI, fck, r); 

            }
            return (1.0 / ky) * soma;            
        }
    }
    
    long double McxFOC(strainState e0kxky) {
        long double sigmacd = stressCompressDesign(fck); 
        long double e0 = e0kxky.e0;
        long double kx = e0kxky.kx;
        long double ky = e0kxky.ky;
        long double soma = 0.0;
        if (abs(kx)<r && abs(ky)<r) {
            return sigmac(e0,sigmacd,n(fck),ec2(fck))*Sx();
        } else if (abs(kx) >= r) {
            for (int i = 0; i < concreteSection.numPoints; i++) {
                int Iplus1;
                if (i != concreteSection.numPoints - 1) {
                    Iplus1 = i + 1;
                } else {
                    Iplus1 = 0;
                }
                long double xI = concreteSection.points[i].x;
                long double xIplus1 = concreteSection.points[Iplus1].x;
                long double yI = concreteSection.points[i].y;
                long double yIplus1 = concreteSection.points[Iplus1].y;
                long double epsI = e0 + ky*xI - kx*yI;
                long double epsIplus1 = e0 + ky*xIplus1 - kx*yIplus1;
                soma += (xIplus1 - xI) * f3i(epsIplus1, epsI, fck, yI, yIplus1, kx, r);
            }
            return -(1.0 / kx) * soma;
        } else {
            for (int i = 0; i < concreteSection.numPoints; i++) {
                int Iplus1;
                if (i != concreteSection.numPoints - 1) {
                    Iplus1 = i + 1;
                } else {
                    Iplus1 = 0;
                }
                long double xI = concreteSection.points[i].x;
                long double xIplus1 = concreteSection.points[Iplus1].x;
                long double yI = concreteSection.points[i].y;
                long double yIplus1 = concreteSection.points[Iplus1].y;
                long double epsI = e0 + ky*xI - kx*yI;
                long double epsIplus1 = e0 + ky*xIplus1 - kx*yIplus1;
                soma += (yIplus1 - yI) * f2i(epsIplus1, epsI, fck, yI, yIplus1, r);
            }
            return -(1.0 / ky) * soma;            
        }
    }

    long double McyFOC(strainState e0kxky) {
        long double sigmacd = stressCompressDesign(fck); 
        long double e0 = e0kxky.e0;
        long double kx = e0kxky.kx;
        long double ky = e0kxky.ky;
        long double soma = 0.0;
        if (abs(kx)<r && abs(ky)<r) {
            return sigmac(e0,sigmacd,n(fck),ec2(fck))*Sy();
        } else if (abs(kx) >= r) {
            for (int i = 0; i < concreteSection.numPoints; i++) {
                int Iplus1;
                if (i != concreteSection.numPoints - 1) {
                    Iplus1 = i + 1;
                } else {
                    Iplus1 = 0;
                }
                long double xI = concreteSection.points[i].x;
                long double xIplus1 = concreteSection.points[Iplus1].x;
                long double yI = concreteSection.points[i].y;
                long double yIplus1 = concreteSection.points[Iplus1].y;
                long double epsI = e0 + ky*xI - kx*yI;
                long double epsIplus1 = e0 + ky*xIplus1 - kx*yIplus1;
                soma += (xIplus1 - xI) * f4i(epsIplus1, epsI, fck, xI, xIplus1, r);
            }
            return (1.0 / kx) * soma;
        } else {
            for (int i = 0; i < concreteSection.numPoints; i++) {
                int Iplus1;
                if (i != concreteSection.numPoints - 1) {
                    Iplus1 = i + 1;
                } else {
                    Iplus1 = 0;
                }
                long double xI = concreteSection.points[i].x;
                long double xIplus1 = concreteSection.points[Iplus1].x;
                long double yI = concreteSection.points[i].y;
                long double yIplus1 = concreteSection.points[Iplus1].y;
                long double epsI = e0 + ky * xI - kx * yI;
                long double epsIplus1 = e0 + ky*xIplus1 - kx*yIplus1;
                soma += (yIplus1 - yI) * f5i(epsIplus1, epsI, fck, yI, yIplus1, ky, r);
            }
            return (1.0 / ky) * soma;            
        }
    }   

    long double NsFOC(strainState e0kxky) {
        long double eyd = strainYieldDesign(fyk, E);
        long double fyd = featureYieldDesign(fyk);
        long double e0 = e0kxky.e0;
        long double kx = e0kxky.kx;
        long double ky = e0kxky.ky;   
        long double total = 0.0;
        for (int i = 0; i < steelLinearDist.numPoints; i++) {
            long double xI = steelLinearDist.points[i].x;
            long double yI = steelLinearDist.points[i].y;            
            long double epsI = e0 + ky * xI - kx * yI;
            point *nextPoint = steelLinearDist.points[i].nextPoint;
            if (nextPoint != NULL) {
                long double lambda = steelLinearDist.points[i].linearSteelThickness;
                long double deltaX = nextPoint->x - xI;
                long double deltaY = nextPoint->y - yI;
                long double dist = sqrt(pow(deltaX, 2) + pow(deltaY, 2));
                long double epsIplus1 = e0 + ky * nextPoint->x - kx * nextPoint->y;
                long double deltaEps = epsIplus1 - epsI;
                long double S1I = S1(epsI, eyd, fyd);
                long double S1Iplus1 = S1(epsIplus1, eyd, fyd);
                long double deltaS1 = S1Iplus1 - S1I;
                if (abs(deltaEps) < r) {
                    total += lambda * dist * sigmas(epsI, fyk, E);
                } else {
                    total += lambda * (dist / deltaEps) * deltaS1;
                }
            }
        }  
        return total;   
    }

    long double MsxFOC(strainState e0kxky) {
        long double eyd = strainYieldDesign(fyk, E);
        long double fyd = featureYieldDesign(fyk);
        long double e0 = e0kxky.e0;
        long double kx = e0kxky.kx;
        long double ky = e0kxky.ky;   
        long double total = 0.0;
        for (int i = 0; i < steelLinearDist.numPoints; i++) {
            long double xI = steelLinearDist.points[i].x;
            long double yI = steelLinearDist.points[i].y;            
            long double epsI = e0 + ky * xI - kx * yI;
            point *nextPoint = steelLinearDist.points[i].nextPoint;
            if (nextPoint != NULL) {
                long double lambda = steelLinearDist.points[i].linearSteelThickness;
                long double deltaX = nextPoint->x - xI;
                long double deltaY = nextPoint->y - yI;
                long double dist = sqrt(pow(deltaX, 2) + pow(deltaY, 2));
                long double epsIplus1 = e0 + ky * nextPoint->x - kx * nextPoint->y;
                long double deltaEps = epsIplus1 - epsI;
                long double S1I = S1(epsI, eyd, fyd);
                long double S1Iplus1 = S1(epsIplus1, eyd, fyd);
                long double deltaS1 = S1Iplus1 - S1I;
                long double gi = yI * epsIplus1 - nextPoint->y * epsI;              
                long double KS1I = KS1(epsI, eyd, fyd);
                long double KS1Iplus1 = KS1(epsIplus1, eyd, fyd);
                long double deltaKS1 = KS1Iplus1 - KS1I;
                if (abs(deltaEps) < r) {
                    total += lambda * dist * sigmas(epsI, fyk, E) * (yI + nextPoint->y) / 2.0;
                } else {
                    total += lambda * (dist / pow(deltaEps, 2)) * (gi * deltaS1 + deltaY * deltaKS1);
                }
            }          
        }
        return -total;  
    }

    long double MsyFOC(strainState e0kxky) {
        long double eyd = strainYieldDesign(fyk, E);
        long double fyd = featureYieldDesign(fyk);
        long double e0 = e0kxky.e0;
        long double kx = e0kxky.kx;
        long double ky = e0kxky.ky;   
        long double total = 0.0;
        for (int i = 0; i < steelLinearDist.numPoints; i++) {
            long double xI = steelLinearDist.points[i].x;
            long double yI = steelLinearDist.points[i].y;            
            long double epsI = e0 + ky * xI - kx * yI;
            point *nextPoint = steelLinearDist.points[i].nextPoint;
            if (nextPoint != NULL) {
                long double lambda = steelLinearDist.points[i].linearSteelThickness;
                long double deltaX = nextPoint->x - xI;
                long double deltaY = nextPoint->y - yI;
                long double dist = sqrt(pow(deltaX, 2) + pow(deltaY, 2));
                long double epsIplus1 = e0 + ky * nextPoint->x - kx * nextPoint->y;
                long double deltaEps = epsIplus1 - epsI;
                long double S1I = S1(epsI, eyd, fyd);
                long double S1Iplus1 = S1(epsIplus1, eyd, fyd);
                long double deltaS1 = S1Iplus1 - S1I;
                long double hi = xI * epsIplus1 - nextPoint->x * epsI;
                long double KS1I = KS1(epsI, eyd, fyd);
                long double KS1Iplus1 = KS1(epsIplus1, eyd, fyd);
                long double deltaKS1 = KS1Iplus1 - KS1I;
                if (abs(deltaEps) < r) {
                    total += lambda * dist * sigmas(epsI, fyk, E) * (xI + nextPoint->x) / 2.0;
                } else {
                    total += lambda * (dist / pow(deltaEps, 2)) * (hi * deltaS1 + deltaX * deltaKS1);
                }
            }           
        }
        return total; 
    }

    long double NtotFOC(strainState e0kxky) {
        long double Nc = NcFOC(e0kxky);
        long double Ns = NsFOC(e0kxky);
        return Nc + Ns;
    }

    long double MxtotFOC(strainState e0kxky) {
        long double Mcx = McxFOC(e0kxky);
        long double Msx = MsxFOC(e0kxky);
        return Mcx + Msx;
    }   

    long double MytotFOC(strainState e0kxky) {
        long double Mcy = McyFOC(e0kxky);
        long double Msy = MsyFOC(e0kxky);
        return Mcy + Msy;
    }      

    long double** Jc(strainState e0kxky) {
        long double** sum2D;
        sum2D = new long double*[3];
        for (int i = 0; i < 3; i++) {
            sum2D[i] = new long double[3];
            for (int j = 0; j < 3; j++) {
                sum2D[i][j] = 0;
            }
        }

        long double e0 = e0kxky.e0;
        long double kx = e0kxky.kx;
        long double ky = e0kxky.ky;
        // cout << kx << 'k' << 'x' << '\n';
        // cout << ky << 'k' << 'y' << '\n';
        if (abs(kx) > r) {
            cout << '1' << '\n';
            for (int i = 0; i < concreteSection.numPoints; i++) {
                long double xI = concreteSection.points[i].x;
                long double yI = concreteSection.points[i].y;
                long double epsI = e0 + ky * xI - kx * yI;
                
                int Iplus1;
                if (i != concreteSection.numPoints - 1) {
                    Iplus1 = i + 1;
                } else {
                    Iplus1 = 0;
                }
                long double xIplus1 = concreteSection.points[Iplus1].x;
                long double yIplus1 = concreteSection.points[Iplus1].y;                
                long double epsIplus1 = e0 + ky * xIplus1 - kx * yIplus1;

                long double dx = xIplus1 - xI;
                // cout << dx << 'd' << 'x' << '\n';

                long double f6if = f6i(epsIplus1, epsI, fck, r);
                long double f7if = f7i(epsIplus1, epsI, fck, yI, yIplus1, r);
                long double f8if = f8i(epsIplus1, epsI, fck, xI, xIplus1, r);
                long double f9if = f9i(epsIplus1, epsI, fck, yI, yIplus1, r);
                long double f10if = f10i(epsIplus1, epsI, fck, yI, yIplus1, xI, xIplus1, r);
                long double f11if = f11i(epsIplus1, epsI, fck, xI, xIplus1, r);

                // cout << epsIplus1 << ' ' <<  epsI << ' ' <<  fck << ' ' <<  yI << ' ' <<  yIplus1 << ' ' <<  xI << ' ' <<  xIplus1 << ' ' <<  r << '\n';

                long double sum[3][3] = {{dx*f6if, -dx*f7if, dx*f8if}, 
                                   {-dx*f7if, dx*f9if, -dx*f10if}, 
                                   {dx*f8if, -dx*f10if, dx*f11if}};

                for (int i = 0; i < 3; i++) {
                    for (int j = 0; j < 3; j++) {
                        sum2D[i][j] += sum[i][j];
                    }
                }
                   
            }
            long double Ncf = NcFOC(e0kxky);
            long double Mcxf = McxFOC(e0kxky);
            long double Mcyf = McyFOC(e0kxky);

            long double sub[3][3] = {{0, Ncf, 0}, {Ncf, 2*Mcxf, Mcyf}, {0, Mcyf, 0}};

            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    sum2D[i][j] = (1.0/kx) * (sum2D[i][j] - sub[i][j]);
                }
            }            
            return sum2D;

        } else if (abs(ky) > r) {
            cout << '2' << '\n';
            for (int i = 0; i < concreteSection.numPoints; i++) {
                long double xI = concreteSection.points[i].x;
                long double yI = concreteSection.points[i].y;
                long double epsI = e0 + ky * xI - kx * yI;
                int Iplus1;
                if (i != concreteSection.numPoints - 1) {
                    Iplus1 = i + 1;
                } else {
                    Iplus1 = 0;
                }
                long double xIplus1 = concreteSection.points[Iplus1].x;
                long double yIplus1 = concreteSection.points[Iplus1].y;                
                long double epsIplus1 = e0 + ky * xIplus1 - kx * yIplus1;

                long double dy = yIplus1 - yI;

                long double f6if = f6i(epsIplus1, epsI, fck, r);
                long double f7if = f7i(epsIplus1, epsI, fck, yI, yIplus1, r);
                long double f8if = f8i(epsIplus1, epsI, fck, xI, xIplus1, r);
                long double f9if = f9i(epsIplus1, epsI, fck, yI, yIplus1, r);
                long double f10if = f10i(epsIplus1, epsI, fck, yI, yIplus1, xI, xIplus1, r);
                long double f11if = f11i(epsIplus1, epsI, fck, xI, xIplus1, r);
                //cout << epsIplus1 << ' ' << epsI << ' ' << fck << ' ' << xI << ' ' << xIplus1 << '\n';
                //cout << f6if << ' ' << f7if << ' ' << f8if << ' ' << f9if << ' ' << f10if << ' ' << f11if << '\n';

                long double sum[3][3] = {{dy*f6if, -dy*f7if, dy*f8if}, 
                                   {-dy*f7if, dy*f9if, -dy*f10if}, 
                                   {dy*f8if, -dy*f10if, dy*f11if}};    

                for (int i = 0; i < 3; i++) {
                    for (int j = 0; j < 3; j++) {
                        sum2D[i][j] += sum[i][j];
                    }
                }           
            
            }
            long double Ncf = NcFOC(e0kxky);
            long double Mcxf = McxFOC(e0kxky);
            long double Mcyf = McyFOC(e0kxky);

            long double sub[3][3] = {{0, 0, Ncf}, {0, 0, Mcxf}, {Ncf, Mcxf, 2*Mcyf}};

            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    sum2D[i][j] = (1.0/ky) * (sum2D[i][j] - sub[i][j]);
                }
            }                  
            return sum2D;
        } else {
            long double dc = Dc(e0, fck);
            long double result[3][3] = {{dc*Ac(), dc*Sx(), dc*Sy()}, {dc*Sx(), dc*Ixx(), -dc*Ixy()}, {dc*Sy(), -dc*Ixy(), dc*Iyy()}};
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    sum2D[i][j] = result[i][j];
                }
            }               
            return sum2D;
        }
    }

    long double** Js(strainState e0kxky) {
        long double** sum2D;
        long double eyd = strainYieldDesign(fyk, E);
        long double fyd = featureYieldDesign(fyk);
        sum2D = new long double*[3];
        for (int i = 0; i < 3; i++) {
            sum2D[i] = new long double[3];
            for (int j = 0; j < 3; j++) {
                sum2D[i][j] = 0.0;
            }
        }

        long double e0 = e0kxky.e0;
        long double kx = e0kxky.kx;
        long double ky = e0kxky.ky;

        for (int i = 0; i < steelLinearDist.numPoints; i++) {
            long double xI = steelLinearDist.points[i].x;
            long double yI = steelLinearDist.points[i].y;
            long double lambda = steelLinearDist.points[i].linearSteelThickness;
            long double epsI = e0 + ky * xI - kx * yI;
            int Iplus1;
            if (i != concreteSection.numPoints - 1) {
                Iplus1 = i + 1;
            } else {
                Iplus1 = 0;
            }
            long double xIplus1 = steelLinearDist.points[Iplus1].x;
            long double yIplus1 = steelLinearDist.points[Iplus1].y;  
            long double epsIplus1 = e0 + ky * xIplus1 - kx * yIplus1; 
            long double deltaX = xIplus1 - xI;
            long double deltaY = yIplus1 - yI;            
            long double dist = sqrt(pow(deltaX, 2) + pow(deltaY, 2));           
            long double deltaEps = epsIplus1 - epsI;            
            long double sigmaI = sigmas(epsI, fyk, E);
            long double sigmaIplus1 = sigmas(epsIplus1, fyk, E);
            long double deltaSigmas = sigmaIplus1 - sigmaI;        
            long double DsI = Ds(epsI, fyk, E);
            long double S1I = S1(epsI, eyd, fyd);
            long double S1Iplus1 = S1(epsIplus1, eyd, fyd);
            long double deltaS1 = S1Iplus1 - S1I;
            long double sigmaEpsI = sigmaI * epsI;
            long double sigmaEpsIplus1 = sigmaIplus1 * epsIplus1;
            long double deltaSigmaEps = sigmaEpsIplus1 - sigmaEpsI;
            long double sigmaEps2I = sigmaI * pow(epsI, 2);
            long double sigmaEps2Iplus1 = sigmaIplus1 * pow(epsIplus1, 2);
            long double deltaSigmaEps2 = sigmaEps2Iplus1 - sigmaEps2I;
            long double gi = yI * epsIplus1 - yIplus1 * epsI;
            long double hi = xI * epsIplus1 - xIplus1 * epsI;
            long double KS1I = KS1(epsI, eyd, fyd);
            long double KS1Iplus1 = KS1(epsIplus1, eyd, fyd);
            long double deltaKS1 = KS1Iplus1 - KS1I;
            long double NsEps0I = NsEps0(lambda,DsI,dist,deltaEps,deltaSigmas,r);
            long double NsKxI = NsKx(lambda,DsI,dist,deltaEps,deltaSigmas,gi,deltaY,deltaSigmaEps,deltaS1,yI,r); 
            long double NsKyI = NsKy(lambda,DsI,dist,deltaEps,deltaSigmas,hi,deltaX,deltaSigmaEps,deltaS1,xI,r);
            long double MsxKxI = MsxKx(lambda,DsI,dist,deltaEps,deltaSigmas,gi,deltaY,deltaSigmaEps,deltaSigmaEps2,deltaS1,deltaKS1,yI,r);
            long double MsyKyI = MsyKy(lambda,DsI,dist,deltaEps,deltaSigmas,hi,deltaX,deltaSigmaEps,deltaSigmaEps2,deltaS1,deltaKS1,xI,r);
            long double MsxKyI = MsxKy(lambda,DsI,dist,deltaEps,deltaSigmas,gi,hi,deltaX,deltaY,deltaSigmaEps,deltaSigmaEps2,deltaS1,deltaKS1,xI,yI,r);
            // cout << NsKyI << 'N' << '\n';
            long double sub[3][3] = {{NsEps0I, NsKxI, NsKyI}, {NsKxI, MsxKxI, MsxKyI}, {NsKyI, MsxKyI, MsyKyI}};

            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    sum2D[i][j] += sub[i][j];
                }
            } 
        }

        return sum2D;
    }

    strainState e0kxky (stresses nMxMy) {
        long double Nd = nMxMy.normalStress;
        long double Mxd = nMxMy.xBendingMoment;
        long double Myd = nMxMy.yBendingMoment;
        long double sigmacd = stressCompressDesign(fck);

        strainState e0kxky = {0.0, 0.0, 0.0, 0.0, 0.0};

        bool over = false;
        int counter = 0;
        long double f_checkpoint = __FLT_MAX__;

        while (!over) {
            
            long double Nc = NcFOC(e0kxky);
            long double Mcx = McxFOC(e0kxky);
            long double Mcy = McyFOC(e0kxky);
        
            long double Ns = NsFOC(e0kxky);
            long double Msx = MsxFOC(e0kxky);
            long double Msy = MsyFOC(e0kxky);


            long double Nt = Ns + Nc;
            long double Mxt = Mcx + Msx;
            long double Myt = Mcy + Msy;

            cout << 'N' << 'd' << Nd << '\n';
            cout << 'M' << 'x' << 'd' << Mxd << '\n';
            cout << 'M' << 'y' << 'd'  << Myd << '\n';

            cout << 'N' << Nt << '\n';
            cout << 'M' << 'x' << Mxt << '\n';
            cout << 'M' << 'y' << Myt << '\n';

            long double Ac_ = this->Ac();

            long double maxY = - __FLT_MAX__;
            long double minY = __FLT_MAX__;

            for (int i = 0; i < concreteSection.numPoints; i++) {
                long double yI = concreteSection.points[i].y;
                if (yI > maxY) {
                    maxY = yI;
                }
                if (yI < minY) {
                    minY = yI;
                }
            }

            long double h = maxY - minY;
            long double f = sqrt(pow(((Nd-Nt)/(sigmacd*Ac_)),2)+pow(((Mxd-Mxt)/(sigmacd*Ac_*h)),2)+pow(((Myd-Myt)/(sigmacd*Ac_*h)),2));
            cout << 'F' << f << ' ' << Ac_ << ' ' << sigmacd << ' ' << '\n';
            if (f <= p) {
                over = true;
                bool broke = brokeSection(e0kxky);
                if (broke) {
                    e0kxky.safe = false;
                    cout << 'k' << '1' << '\n';
                    return e0kxky; 
                } else {
                    e0kxky.safe = true;
                    cout <<'k' << '2' << '\n';
                    return e0kxky;
                }
            } else {
                cout.precision(35);

                cout << e0kxky.e0 << '\n';
                cout << e0kxky.kx << '\n';
                cout << e0kxky.ky << '\n';
                cout << '\n';

                long double** Jcf = Jc(e0kxky);
                long double** Jsf = Js(e0kxky);

                long double** Jt = new long double*[3];
                for (int i = 0; i < 3; i++) {
                    Jt[i] = new long double[3];
                    for (int j = 0; j < 3; j++) {
                        Jt[i][j] = Jsf[i][j] + Jcf[i][j];
                        cout << Jsf[i][j];
                        cout << ' ';
                    }
                    cout << '\n';
                }
                cout << '\n';

                long double G = determinant(Jt, 3.0);
                long double* dStresses = new long double[3];
                dStresses[0] = Nd - Nt;
                dStresses[1] = Mxd - Mxt;
                dStresses[2] = Myd - Myt;

                if (abs(G) > t) {

                    long double Dx = determinant(Jt, 3.0, true, dStresses, 0.0);
                    long double delE0 = Dx / G;
                    long double Dy = determinant(Jt, 3.0, true, dStresses, 1.0);
                    long double delKx = Dy / G;
                    long double Dz = determinant(Jt, 3.0, true, dStresses, 2.0);
                    long double delKy = Dz / G;

                    e0kxky.e0 += delE0;
                    e0kxky.kx += delKx;
                    e0kxky.ky += delKy;
                    cout << '\n' << G << '\n' << '\n';

                    cout << delE0 << '\n';
                    cout << delKx << '\n';
                    cout << delKy << '\n';
                    cout << '\n';

                } else {
                    over = true;
                    e0kxky.safe = false;
                    cout << 'k' << '3' << '\n';
                    return e0kxky; 
                }
            }
            counter++;
            cout << 'C' << 'T' << counter << '\n';
            if (counter % 1000 == 0) {
                if (f > 0.9*f_checkpoint) {
                    over = true;
                    e0kxky.safe = false;
                    cout << counter << 'k' <<  '4' << '\n';
                    return e0kxky;                
                }
                f_checkpoint = f;

                if (counter > 100000) {
                    over = true;
                    e0kxky.safe = false;
                    cout << 'k' << '5' << '\n';
                    return e0kxky;                       
                }
            }
        }
    }  
    
    long double* maxminStrainsConcreteSteel(strainState e0kxky, long double *extremes) {
        extremes[0] = - __FLT_MAX__, extremes[2] = extremes[0];
        extremes[1] = __FLT_MAX__, extremes[3] = extremes[1];
        
        for (int i = 0; i < concreteSection.numPoints; i++) {
            long double x = concreteSection.points[i].x;
            long double y = concreteSection.points[i].y;
            long double strain = e0kxky.e0 + e0kxky.ky * x - e0kxky.kx * y;
            if (strain > extremes[0]) {
                extremes[0] = strain;
            }
            if (strain < extremes[1]) {
                extremes[1] = strain;
            }
        }

        for (int i = 0; i < steelLinearDist.numPoints; i++) {
            long double x = steelLinearDist.points[i].x;
            long double y = steelLinearDist.points[i].y;
            long double strain = e0kxky.e0 + e0kxky.ky * x - e0kxky.kx * y;
            if (strain > extremes[2]) {
                extremes[2] = strain;
            }
            if (strain < extremes[3]) {
                extremes[3] = strain;
            }
        }

        return extremes;
    }

    bool brokeSection(strainState e0kxky) {
        bool broke = false;
        long double arr[4];
        long double* extremes = maxminStrainsConcreteSteel(e0kxky, arr);
        long double u = ecu(fck);
        long double c = ec2(fck);
        if (extremes[0] > u){
            return true;
        }
        if ((u - c) * extremes[0] + c * extremes[1] > u * c){
            return true;
        }
        if (extremes[3] < -10.0){
            return true;
        }
        return false;
    }
};



//Main Function

int main() {
    // Object Creation For Class

    point point1;
    point1.x = -10;
    point1.y = -50;

    point point2;
    point2.x = 10;
    point2.y = -50;

    point point3;
    point3.x = 10;
    point3.y = 50;

    point point4;
    point4.x = -10;
    point4.y = 50;

    polygon PointPolygonConcrete;
    
    PointPolygonConcrete.pushPoint(point1);
    PointPolygonConcrete.pushPoint(point2);
    PointPolygonConcrete.pushPoint(point3);
    PointPolygonConcrete.pushPoint(point4);


    point point5;
    point5.x = -7;
    point5.y = -47;
    

    point point6;
    point6.x = 7;
    point6.y = -47;
    point5.nextPoint = &point6;
    point5.linearSteelThickness = 1;

    point point7;
    point7.x = 7;
    point7.y = 47;
    point6.nextPoint = &point7;
    point6.linearSteelThickness = 1;

    point point8;
    point8.x = -7;
    point8.y = 47;
    point7.nextPoint = &point8;
    point7.linearSteelThickness = 1;
    point8.nextPoint = &point5;
    point8.linearSteelThickness = 1;

    polygon PointPolygonSteel;

    PointPolygonSteel.pushPoint(point5);
    PointPolygonSteel.pushPoint(point6);
    PointPolygonSteel.pushPoint(point7);
    PointPolygonSteel.pushPoint(point8);


    crossSection mySection;
    mySection.defineConcreteSection(PointPolygonConcrete);
    mySection.defineSteelLinearDist(PointPolygonSteel);
    mySection.showCrossSection();

    stresses solicitingStresses = {0, 0, 4, 42, 2.03};
    cout << solicitingStresses.xBendingMoment << 'X';

    strainState e0kxky = {0, 0, 0.000106461, 2.4232e-07, 8.22006e-06};

    cout << '\n';

    long double** Jsf = mySection.Js(e0kxky);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            cout << Jsf[i][j] << ' ';
        }
        cout << '\n';
    }

    cout << '\n';

    strainState deform = mySection.e0kxky(solicitingStresses);

    cout << deform.e0 << '\n';
    cout << deform.kx << '\n';
    cout << deform.ky << '\n';
    cout << deform.safe << '\n';

    // cout << mySection.Jc(e0kxky)

    //getchar();
    return 0;
}


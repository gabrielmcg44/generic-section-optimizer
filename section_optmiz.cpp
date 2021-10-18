
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
    double x;
    double y;
    double steelArea;
    double linearSteelThickness;
    struct point *nextPoint;
};

// class point {
//     private:
//     public:
//         double x;
//         double y;
//         double steelArea;
//         double linearSteelThickness;
//         point nextPoint;
// }


struct stresses {
    double x = 0;
    double y = 0;
    double normalStress;
    double xBendingMoment;
    double yBendingMoment;
    bool safe;
};


struct strainState {
    double x = 0;
    double y = 0;
    double e0; // strain in the origin
    double kx; // variation of strain per distance from x axis
    double ky; // variation of strain per distance from y axis
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

    double minX() {
        double minimum = __FLT_MAX__;
        for (int i = 0; i < numPoints; i++) {
            if (points[i].x < minimum) {
                minimum = points[i].x;
            }
        }
        return minimum;
    }

    double maxX() {
        double maximum = - __FLT_MAX__;
        for (int i = 0; i < numPoints; i++) {
            if (points[i].x > maximum) {
                maximum = points[i].x;
            }
        }
        return maximum;
    }

    double minY() {
        double minimum = __FLT_MAX__;
        for (int i = 0; i < numPoints; i++) {
            if (points[i].y < minimum) {
                minimum = points[i].y;
            }
        }
        return minimum;
    }

    double maxY() {
        double maximum = - __FLT_MAX__;
        for (int i = 0; i < numPoints; i++) {
            if (points[i].y > maximum) {
                maximum = points[i].y;
            }
        }
        return maximum;
    }
};


double determinant(double** matrix, int n, bool replaceColumn = false, double* newColumn = new double[3] , int newColumnNumber = 0) {
    double det = 0;

    double** submatrix;
    double** newmatrix;
    submatrix = new double*[3];
    newmatrix = new double*[3];
    for (int i = 0; i < 3; i++) {
        submatrix[i] = new double[3];
        newmatrix[i] = new double[3];
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


double featureYieldDesign(double fyk) {
    // featureYieldDesign 'fyd'
    return fyk / 1.15;
}


double strainYieldDesign(double fyk, double E) {
    // strainYieldDesign 'epsilon yd'
    return 1000 * featureYieldDesign(fyk) / E;
}


double stressCompressDesign(double fck) {
    // stressCompressDesign 'sigma cd'
    return (0.85 / 1.4) * fck;
}


double ec2(double fck) {
    // 
    if (fck >= 0.0 && fck <= 5.0) {
        return 2.0;
    } else if (fck > 5.0 && fck <= 9.0) {
        return 2.0 + 0.085 * pow(fck - 5.0, 0.53);
    } else {
        return -1;
    }   
}


double ecu(double fck) {
    // 
    if (fck >= 0 && fck <= 5) {
        return 3.5;
    } else if (fck > 5 && fck <= 9) {
        return 2.6 + 35 * pow((9 - fck) / 10, 4);
    } else {
        return -1;
    }   
}


double n(double fck) {
    //  
    if (fck >= 0.0 && fck <= 5.0) {
        return 2.0;
    } else if (fck > 5.0 && fck <= 9.0) {
        return 1.4 + 23.4 * pow((9.0 - fck) / 10.0, 4.0);
    } else {
        return -1;
    }   
}


double sigmac(double eps, double sigmacd, double n, double ec2) {
    // Calculate tension of concrete given strain
    if (eps <= 0) {
        return 0;
    } else if (eps > 0 && eps <= 2) {
        return sigmacd * (1 - pow(1 -(eps / ec2), n));
    } else {
        return sigmacd;
    }   
}


double sigmas(double eps, double fyk, double E) {
    // Calculate tension of concrete given strain
    double eyd = strainYieldDesign(fyk, E);
    double fyd = featureYieldDesign(fyk);
    if (eps <= - eyd) {
        return - fyd;
    } else if (eps > -eyd && eps <= eyd) {
        return fyd * (eps / eyd);
    } else {
        return fyd;
    }   
}


double tk(double k, double fck, double eps) {
    double c = ec2(fck);
    double n1 = n(fck);
    // cout << n1 << 'n' << '\n';
    // cout << n1 << ' ' << k << ' ' << eps << ' ' << c << ' ' << '\n';
    
    if (eps >= 0 && eps < c) {
        // cout << 1 - eps/c << 't' << '\n';
        // cout << 'c' << '1' << '\n';
        return (1/(n1 + k)) * (1 - pow(1 - eps/c, n1 + k));
    } else {
        // cout << 'c' << '2' << '\n';
        return 1/(n1 + k);
    }
}


double I1(double eps, double fck) {
    double c = ec2(fck);
    double sigmacd = stressCompressDesign(fck);
    // cout << tk(1, fck, eps) << 't' << '\n';
    if (eps < 0) {
        return 0;
    } else {
        return (eps - c * tk(1, fck, eps)) * sigmacd;
    }
}


double I2(double eps, double fck) {
    double c = ec2(fck);
    double n1 = n(fck);
    double sigmacd = stressCompressDesign(fck);
    if (eps < 0) {
        return 0;
    } else {
        return (pow(eps, 2)/2 - (c/(n1 + 1)) * (eps - c*tk(2, fck, eps))) * sigmacd;
    }
}


double I3(double eps, double fck) {
    double c = ec2(fck);
    double n1 = n(fck);
    double sigmacd = stressCompressDesign(fck);
    if (eps < 0) {
        return 0;
    } else {
        return ((pow(eps,3))/6-(c/(n1+1))*((pow(eps,2))/2-(c/(n1+2))*(eps-c*tk(3,fck,eps))))*sigmacd;
    }
}


double J1(double eps, double fck) {
    double c = ec2(fck);
    double sigmacd = stressCompressDesign(fck);
    if (eps < 0) {
        return 0;
    } else {
        return ((pow(eps,2)/2)-(pow(c,2))*(tk(1,fck,eps)-tk(2,fck,eps)))*sigmacd;
    }
}


double J2(double eps, double fck) {
    double c = ec2(fck);
    double sigmacd = stressCompressDesign(fck);
    if (eps < 0) {
        return 0;
    } else {
        return ((pow(eps,3)/3)-(pow(c,3))*(tk(1,fck,eps)-2*tk(2,fck,eps)+tk(3,fck,eps)))*sigmacd;
    }
}


double K1(double eps, double fck) {
    double c = ec2(fck);
    double n1 = n(fck);
    double sigmacd = stressCompressDesign(fck);
    if (eps < 0) {
        return 0;
    } else {
        return ((pow(eps,3)/3)-(c/(n1+1))*(pow(eps,2)/2-(pow(c,2))*(tk(2,fck,eps)-tk(3,fck,eps))))*sigmacd;
    }
}


double f1i(double epsIplus1, double epsI, double fck, double r) {
    if (abs(epsIplus1 - epsI) < r) {
        return I1(epsI,fck);
    } else {
        return (I2(epsIplus1,fck)-I2(epsI,fck)) / (epsIplus1 - epsI);
    }
}


double f2i(double epsIplus1, double epsI, double fck, double yI, double yIplus1, double r) {
    double gI = yI*epsIplus1 - yIplus1*epsI;
    if (abs(epsIplus1 - epsI) < r) {
        return I1(epsI,fck)*((yI+yIplus1)/2);
    } else {
        return (gI*(I2(epsIplus1,fck)-I2(epsI,fck))+(yIplus1-yI)*(K1(epsIplus1,fck)-K1(epsI,fck)))/(pow(epsIplus1-epsI,2));
    }
}


double f3i(double epsIplus1, double epsI, double fck, double yI, double yIplus1, double kx, double r) {
    if (abs(epsIplus1 - epsI) < r) {
        return f2i(epsIplus1,epsI,fck,yI,yIplus1,r)+(I2(epsI,fck)/kx);
    } else {
        return f2i(epsIplus1,epsI,fck,yI,yIplus1,r)+(((I3(epsIplus1,fck)-I3(epsI,fck))/(epsIplus1-epsI))/kx);
    }
}


double f4i(double epsIplus1, double epsI, double fck, double xI, double xIplus1, double r) {
    double hI=xI*epsIplus1-xIplus1*epsI;
    if (abs(epsIplus1-epsI)<r) {
        return I1(epsI,fck)*((xI+xIplus1)/2);
    } else {
        return (hI*(I2(epsIplus1,fck)-I2(epsI,fck))+(xIplus1-xI)*(K1(epsIplus1,fck)-K1(epsI,fck)))/(pow(epsIplus1-epsI, 2));
    }
}


double f5i(double epsIplus1, double epsI, double fck, double xI, double xIplus1, double ky, double r) {
    if (abs(epsIplus1-epsI)<r) {
        return f4i(epsIplus1,epsI,fck,xI,xIplus1,r)-(I2(epsI,fck)/ky);
    } else {
        return f4i(epsIplus1,epsI,fck,xI,xIplus1,r)+(((I3(epsIplus1,fck)-I3(epsI,fck))/(epsIplus1-epsI))/ky);
    }
}


double f6i(double epsIplus1, double epsI, double fck, double r) {
    double sigmacd = stressCompressDesign(fck);    
    if (abs(epsIplus1-epsI)<r) {
        // cout << 1 << '\n';
        return sigmac(epsI,sigmacd,n(fck),ec2(fck));
    } else {
        // cout << 2 << '\n';
        // cout << I1(epsIplus1,fck) << 'I' << '1' << '\n';
        return (I1(epsIplus1,fck)-I1(epsI,fck))/(epsIplus1-epsI);
    }
}


double f7i(double epsIplus1, double epsI, double fck, double yI, double yIplus1, double r) {
    double sigmacd = stressCompressDesign(fck);   
    double gI = yI*epsIplus1 - yIplus1*epsI; 
    if (abs(epsIplus1-epsI)<r) {
        return sigmac(epsI,sigmacd,n(fck),ec2(fck))*((yI+yIplus1)/2);
    } else {
        return (gI*(I1(epsIplus1,fck)-I1(epsI,fck))+(yIplus1-yI)*(J1(epsIplus1,fck)-J1(epsI,fck)))/(pow(epsIplus1-epsI,2));
    }
}


double f8i(double epsIplus1, double epsI, double fck, double xI, double xIplus1, double r) {
    double sigmacd = stressCompressDesign(fck);   
    double hI = xI*epsIplus1 - xIplus1*epsI; 
    if (abs(epsIplus1-epsI)<r) {
        cout << 'c' << '1' << '\n';
        return sigmac(epsI,sigmacd,n(fck),ec2(fck))*((xI+xIplus1)/2);
    } else {
        cout << 'c' << '2' << '\n';
        cout << (J1(epsIplus1,fck)-J1(epsI,fck)) << '\n';
        return (hI*(I1(epsIplus1,fck)-I1(epsI,fck))+(xIplus1-xI)*(J1(epsIplus1,fck)-J1(epsI,fck)))/(pow(epsIplus1-epsI,2));
    }
}


double f9i(double epsIplus1, double epsI, double fck, double yI, double yIplus1, double r) {
    double sigmacd = stressCompressDesign(fck);   
    double gI = yI*epsIplus1 - yIplus1*epsI; 
    if (abs(epsIplus1-epsI)<r) {
        return sigmac(epsI,sigmacd,n(fck),ec2(fck))*((pow(yI,2)+yI*yIplus1+pow(yIplus1,2))/3);
    } else {
        return ((pow(gI,2))*(I1(epsIplus1,fck)-I1(epsI,fck))+2*gI*(yIplus1-yI)*(J1(epsIplus1,fck)-J1(epsI,fck))+(pow(yIplus1-yI,2))*(J2(epsIplus1,fck)-J2(epsI,fck)))/(pow(epsIplus1-epsI,3));
    }
}


double f10i(double epsIplus1, double epsI, double fck, double yI, double yIplus1, double xI, double xIplus1, double r) {
    double gI = yI*epsIplus1 - yIplus1*epsI;
    double hI = xI*epsIplus1 - xIplus1*epsI;
    double sigmacd = stressCompressDesign(fck);  
    double sigmai = sigmac(epsI, sigmacd, n(fck), ec2(fck));
    double deltaI1 = (I1(epsIplus1,fck)-I1(epsI,fck));
    double deltaJ1 = (J1(epsIplus1,fck)-J1(epsI,fck));
    double deltaJ2 = (J2(epsIplus1,fck)-J2(epsI,fck));
    double deltaY = yIplus1 - yI;
    double deltaX = xIplus1 - xI;
    double deltaEps = (epsIplus1 - epsI);
    if (abs(epsIplus1-epsI)<r) {
        return sigmai*((xI*yIplus1+2*(xI*yI+xIplus1*yIplus1)+xIplus1*yI)/6);
    } else {
        return (hI*gI*deltaI1+deltaJ1*(gI*deltaX+hI*deltaY)+deltaX*deltaY*deltaJ2)/(pow(deltaEps,3));
    }
}


double f11i(double epsIplus1, double epsI, double fck, double xI, double xIplus1, double r) {
    double hI=xI*epsIplus1-xIplus1*epsI;
    double sigmacd = stressCompressDesign(fck);  
    if (abs(epsIplus1-epsI)<r) {
        return sigmac(epsI,sigmacd,n(fck),ec2(fck))*((pow(xI,2)+xI*xIplus1+pow(xIplus1,2))/3);
    } else {
        return ((pow(hI,2))*(I1(epsIplus1,fck)-I1(epsI,fck))+2*hI*(xIplus1-xI)*(J1(epsIplus1,fck)-J1(epsI,fck))+(pow(xIplus1-xI,2))*(J2(epsIplus1,fck)-J2(epsI,fck)))/(pow(epsIplus1-epsI,3));
    }
}


double S1(double eps, double eyd, double fyd) {
    if (abs(eps) < eyd) {
        return (1.0/2.0)*(fyd/eyd)*pow(eps,2);
    } else {
        double result = abs(eps*fyd) - eyd*fyd/2;
        return abs(eps*fyd) - eyd*fyd/2;
    }
}


double KS1(double eps, double eyd, double fyd) {
    if (eps < -eyd) {
        return -(1.0/2.0)*(fyd)*pow(eps,2)+(fyd*pow(eyd,2))/6;
    } else if (eps >= -eyd && eyd > eps){
        return (pow(eps,3)/3)*(fyd/eyd);
    } else {
        return (1.0/2.0)*(fyd)*pow(eps,2)-(fyd*pow(eyd,2))/6;
    }
}


double Dc(double eps, double fck) {
    double c = ec2(fck);
    double nf = n(fck);
    double sigmacd = stressCompressDesign(fck); 
    if (eps<c && eps >= 0) {
        return sigmacd * (nf / c) * pow(1 - eps / c, nf - 1);
    } else {
        return 0;
    }
}


double Ds(double eps, double fyk, double E) {
    double eyd = strainYieldDesign(fyk, E); 
    if (eps < eyd && eps >= -eyd) {
        return E / 1000;
    } else {
        return 0;
    }
}


double NsEps0(double lambda, double DsI, double dist, double deltaEps, double deltaSigmas, double r) {
    if (abs(deltaEps) < r) {
        return lambda * dist * DsI;
    } else {
        return lambda * (dist / deltaEps) * deltaSigmas;
    } 
}


double NsKx(double lambda, double DsI, double dist, double deltaEps, double deltaSigmas, double gi, double deltaY, double deltaSigmaEps, double deltaS1, double Yi, double r) {
    if (abs(deltaEps) < r) {
        return - lambda * dist * DsI * (Yi + deltaY / 2);
    } else {
        return - lambda * (dist / pow(deltaEps, 2)) * (gi * deltaSigmas + deltaY * deltaSigmaEps - deltaY * deltaS1);
    }
}


double NsKy(double lambda, double DsI, double dist, double deltaEps, double deltaSigmas, double hi, double deltaX, double deltaSigmaEps, double deltaS1, double Xi, double r) {
    if (abs(deltaEps) < r) {
        return - lambda * dist * DsI * (Xi + deltaX / 2);
    } else {
        return - lambda * (dist / pow(deltaEps, 2)) * (hi * deltaSigmas + deltaX * deltaSigmaEps - deltaX * deltaS1);
    }
}


double MsxKx (double lambda, double DsI, double dist, double deltaEps, double deltaSigmas, double gi, double deltaY, double deltaSigmaEps, double deltaSigmaEps2, double deltaS1, double deltaKS1, double Yi, double r) {
    if (abs(deltaEps) < r) {
        return lambda * dist * DsI * (pow(Yi, 2) + Yi * deltaY + pow(deltaY, 2)/3);
    } else {
        return lambda * (dist/pow(deltaEps,3)) * (pow(gi, 2) * deltaSigmas + 2 * gi * deltaY * (deltaSigmaEps - deltaS1) + pow(deltaY, 2) * (deltaSigmaEps2 - 2 * deltaKS1));
    }
}


double MsyKy (double lambda, double DsI, double dist, double deltaEps, double deltaSigmas, double hi, double deltaX, double deltaSigmaEps, double deltaSigmaEps2, double deltaS1, double deltaKS1, double Xi, double r) {
    if (abs(deltaEps) < r) {
        return lambda * dist * DsI * (pow(Xi, 2) + Xi * deltaX + pow(deltaX, 2)/3);
    } else {
        return lambda * (dist/pow(deltaEps,3)) * (pow(hi, 2) * deltaSigmas + 2 * hi * deltaX * (deltaSigmaEps - deltaS1) + pow(deltaX, 2) * (deltaSigmaEps2 - 2 * deltaKS1));
    }
}


double MsxKy (double lambda, double DsI, double dist, double deltaEps, double deltaSigmas, double gi, double hi, double deltaX, double deltaY, double deltaSigmaEps, double deltaSigmaEps2, double deltaS1, double deltaKS1, double Xi, double Yi, double r) {
    if (abs(deltaEps) < r) {
        return -lambda * dist * DsI * (Xi*Yi + (Yi*deltaX + Xi*deltaY)/2 + deltaX*deltaY/3);
    } else {
        return -lambda * (dist/pow(deltaEps, 3)) * (gi*hi*deltaSigmas + (gi*deltaX + hi*deltaY) * (deltaSigmaEps - deltaS1) + deltaX*deltaY*(deltaSigmaEps2 - 2*deltaKS1));
    }
}


class crossSection {
private:
    polygon concreteSection;
    polygon steelLinearDist;

    //Access - Specifier
public:
    double fyk = 50; //featureYieldKnown
    double fck = 4; //featureCompressionKnown
    double E = 21000; //youngModulus
    double r = 0.00001; //tolerance
    double p = 0.000001;
    double t = 0.000000000001;
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

    double Ac() {
        double soma = 0;
        for (int i = 0; i < concreteSection.numPoints; i++) {
            int Iplus1;
            if (i != concreteSection.numPoints - 1) {
                Iplus1 = i + 1;
            } else {
                Iplus1 = 0;
            }
            double xI = concreteSection.points[i].x;
            double xIplus1 = concreteSection.points[Iplus1].x;
            double yI = concreteSection.points[i].y;
            double yIplus1 = concreteSection.points[Iplus1].y;          
            double aI = xI*yIplus1 - xIplus1*yI;
            soma += aI;
        }
        return soma / 2;
    }
    
    double Sy() {
        double soma = 0;
        for (int i = 0; i < concreteSection.numPoints; i++) {
            int Iplus1;
            if (i != concreteSection.numPoints - 1) {
                Iplus1 = i + 1;
            } else {
                Iplus1 = 0;
            }
            double xI = concreteSection.points[i].x;
            double xIplus1 = concreteSection.points[Iplus1].x;
            double yI = concreteSection.points[i].y;
            double yIplus1 = concreteSection.points[Iplus1].y;          
            double aI = xI*yIplus1 - xIplus1*yI;
            soma += aI * (xI + xIplus1);
        }
        return soma / 6;
    }
    
    double Sx() {
        double soma = 0;
        for (int i = 0; i < concreteSection.numPoints; i++) {
            int Iplus1;
            if (i != concreteSection.numPoints - 1) {
                Iplus1 = i + 1;
            } else {
                Iplus1 = 0;
            }
            double xI = concreteSection.points[i].x;
            double xIplus1 = concreteSection.points[Iplus1].x;
            double yI = concreteSection.points[i].y;
            double yIplus1 = concreteSection.points[Iplus1].y;          
            double aI = xI*yIplus1 - xIplus1*yI;
            soma += aI * (yI + yIplus1);
        }
        return soma / 6;
    }
    
    double Iyy() {
        double soma = 0;
        for (int i = 0; i < concreteSection.numPoints; i++) {
            int Iplus1;
            if (i != concreteSection.numPoints - 1) {
                Iplus1 = i + 1;
            } else {
                Iplus1 = 0;
            }
            double xI = concreteSection.points[i].x;
            double xIplus1 = concreteSection.points[Iplus1].x;
            double yI = concreteSection.points[i].y;
            double yIplus1 = concreteSection.points[Iplus1].y;          
            double aI = xI*yIplus1 - xIplus1*yI;
            soma += aI * (pow(xI,2) +  xI*xIplus1 + pow(xIplus1,2));
        }
        return soma / 12;
    } 
    
    double Ixy() {
        double soma = 0;
        for (int i = 0; i < concreteSection.numPoints; i++) {
            int Iplus1;
            if (i != concreteSection.numPoints - 1) {
                Iplus1 = i + 1;
            } else {
                Iplus1 = 0;
            }
            double xI = concreteSection.points[i].x;
            double xIplus1 = concreteSection.points[Iplus1].x;
            double yI = concreteSection.points[i].y;
            double yIplus1 = concreteSection.points[Iplus1].y;          
            double aI = xI*yIplus1 - xIplus1*yI;
            soma += aI * (xI*yIplus1+2*(xI*yI+xIplus1*yIplus1)+xIplus1*yI);
        }
        return soma / 24;
    }

    double Ixx() {
        double soma = 0;
        for (int i = 0; i < concreteSection.numPoints; i++) {
            int Iplus1;
            if (i != concreteSection.numPoints - 1) {
                Iplus1 = i + 1;
            } else {
                Iplus1 = 0;
            }
            double xI = concreteSection.points[i].x;
            double xIplus1 = concreteSection.points[Iplus1].x;
            double yI = concreteSection.points[i].y;
            double yIplus1 = concreteSection.points[Iplus1].y;          
            double aI = xI*yIplus1 - xIplus1*yI;
            soma += aI * (pow(yI,2) +  yI*yIplus1 + pow(yIplus1,2));
        }
        return soma / 12;
    }

    double NcFOC(strainState e0kxky) {
        double sigmacd = stressCompressDesign(fck); 
        double e0 = e0kxky.e0;
        double kx = e0kxky.kx;
        double ky = e0kxky.ky;
        double soma = 0;
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
                double xI = concreteSection.points[i].x;
                double xIplus1 = concreteSection.points[Iplus1].x;
                double yI = concreteSection.points[i].y;
                double yIplus1 = concreteSection.points[Iplus1].y;
                double epsI = e0 + ky*xI - kx*yI;
                double epsIplus1 = e0 + ky*xIplus1 - kx*yIplus1;    
                soma += (xIplus1 - xI) * f1i(epsIplus1, epsI, fck, r);

            }
            return (1 / kx) * soma;
        } else {
            for (int i = 0; i < concreteSection.numPoints; i++) {
                int Iplus1;
                if (i != concreteSection.numPoints - 1) {
                    Iplus1 = i + 1;
                } else {
                    Iplus1 = 0;
                }
                double xI = concreteSection.points[i].x;
                double xIplus1 = concreteSection.points[Iplus1].x;
                double yI = concreteSection.points[i].y;
                double yIplus1 = concreteSection.points[Iplus1].y;
                double epsI = e0 + ky*xI - kx*yI;
                double epsIplus1 = e0 + ky*xIplus1 - kx*yIplus1;                 
                soma += (yIplus1 - yI) * f1i(epsIplus1, epsI, fck, r); 

            }
            return (1 / ky) * soma;            
        }
    }
    
    double McxFOC(strainState e0kxky) {
        double sigmacd = stressCompressDesign(fck); 
        double e0 = e0kxky.e0;
        double kx = e0kxky.kx;
        double ky = e0kxky.ky;
        double soma = 0;
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
                double xI = concreteSection.points[i].x;
                double xIplus1 = concreteSection.points[Iplus1].x;
                double yI = concreteSection.points[i].y;
                double yIplus1 = concreteSection.points[Iplus1].y;
                double epsI = e0 + ky*xI - kx*yI;
                double epsIplus1 = e0 + ky*xIplus1 - kx*yIplus1;
                soma += (xIplus1 - xI) * f3i(epsIplus1, epsI, fck, yI, yIplus1, kx, r);
            }
            return -(1 / kx) * soma;
        } else {
            for (int i = 0; i < concreteSection.numPoints; i++) {
                int Iplus1;
                if (i != concreteSection.numPoints - 1) {
                    Iplus1 = i + 1;
                } else {
                    Iplus1 = 0;
                }
                double xI = concreteSection.points[i].x;
                double xIplus1 = concreteSection.points[Iplus1].x;
                double yI = concreteSection.points[i].y;
                double yIplus1 = concreteSection.points[Iplus1].y;
                double epsI = e0 + ky*xI - kx*yI;
                double epsIplus1 = e0 + ky*xIplus1 - kx*yIplus1;
                soma += (yIplus1 - yI) * f2i(epsIplus1, epsI, fck, yI, yIplus1, r);
            }
            return -(1 / ky) * soma;            
        }
    }

    double McyFOC(strainState e0kxky) {
        double sigmacd = stressCompressDesign(fck); 
        double e0 = e0kxky.e0;
        double kx = e0kxky.kx;
        double ky = e0kxky.ky;
        double soma = 0;
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
                double xI = concreteSection.points[i].x;
                double xIplus1 = concreteSection.points[Iplus1].x;
                double yI = concreteSection.points[i].y;
                double yIplus1 = concreteSection.points[Iplus1].y;
                double epsI = e0 + ky*xI - kx*yI;
                double epsIplus1 = e0 + ky*xIplus1 - kx*yIplus1;
                soma += (xIplus1 - xI) * f4i(epsIplus1, epsI, fck, xI, xIplus1, r);
            }
            return (1 / kx) * soma;
        } else {
            for (int i = 0; i < concreteSection.numPoints; i++) {
                int Iplus1;
                if (i != concreteSection.numPoints - 1) {
                    Iplus1 = i + 1;
                } else {
                    Iplus1 = 0;
                }
                double xI = concreteSection.points[i].x;
                double xIplus1 = concreteSection.points[Iplus1].x;
                double yI = concreteSection.points[i].y;
                double yIplus1 = concreteSection.points[Iplus1].y;
                double epsI = e0 + ky * xI - kx * yI;
                double epsIplus1 = e0 + ky*xIplus1 - kx*yIplus1;
                soma += (yIplus1 - yI) * f5i(epsIplus1, epsI, fck, yI, yIplus1, ky, r);
            }
            return (1 / ky) * soma;            
        }
    }   

    double NsFOC(strainState e0kxky) {
        double eyd = strainYieldDesign(fyk, E);
        double fyd = featureYieldDesign(fyk);
        double e0 = e0kxky.e0;
        double kx = e0kxky.kx;
        double ky = e0kxky.ky;   
        double total = 0;
        for (int i = 0; i < steelLinearDist.numPoints; i++) {
            double xI = steelLinearDist.points[i].x;
            double yI = steelLinearDist.points[i].y;            
            double epsI = e0 + ky * xI - kx * yI;
            point *nextPoint = steelLinearDist.points[i].nextPoint;
            if (nextPoint != NULL) {
                double lambda = steelLinearDist.points[i].linearSteelThickness;
                double deltaX = nextPoint->x - xI;
                double deltaY = nextPoint->y - yI;
                double dist = sqrt(pow(deltaX, 2) + pow(deltaY, 2));
                double epsIplus1 = e0 + ky * nextPoint->x - kx * nextPoint->y;
                double deltaEps = epsIplus1 - epsI;
                double S1I = S1(epsI, eyd, fyd);
                double S1Iplus1 = S1(epsIplus1, eyd, fyd);
                double deltaS1 = S1Iplus1 - S1I;
                if (abs(deltaEps) < r) {
                    total += lambda * dist * sigmas(epsI, fyk, E);
                } else {
                    total += lambda * (dist / deltaEps) * deltaS1;
                }
            }
        }  
        return total;   
    }

    double MsxFOC(strainState e0kxky) {
        double eyd = strainYieldDesign(fyk, E);
        double fyd = featureYieldDesign(fyk);
        double e0 = e0kxky.e0;
        double kx = e0kxky.kx;
        double ky = e0kxky.ky;   
        double total = 0;
        for (int i = 0; i < steelLinearDist.numPoints; i++) {
            double xI = steelLinearDist.points[i].x;
            double yI = steelLinearDist.points[i].y;            
            double epsI = e0 + ky * xI - kx * yI;
            point *nextPoint = steelLinearDist.points[i].nextPoint;
            if (nextPoint != NULL) {
                double lambda = steelLinearDist.points[i].linearSteelThickness;
                double deltaX = nextPoint->x - xI;
                double deltaY = nextPoint->y - yI;
                double dist = sqrt(pow(deltaX, 2) + pow(deltaY, 2));
                double epsIplus1 = e0 + ky * nextPoint->x - kx * nextPoint->y;
                double deltaEps = epsIplus1 - epsI;
                double S1I = S1(epsI, eyd, fyd);
                double S1Iplus1 = S1(epsIplus1, eyd, fyd);
                double deltaS1 = S1Iplus1 - S1I;
                double gi = yI * epsIplus1 - nextPoint->y * epsI;              
                double KS1I = KS1(epsI, eyd, fyd);
                double KS1Iplus1 = KS1(epsIplus1, eyd, fyd);
                double deltaKS1 = KS1Iplus1 - KS1I;
                if (abs(deltaEps) < r) {
                    total += lambda * dist * sigmas(epsI, fyk, E) * (yI + nextPoint->y) / 2.0;
                } else {
                    total += lambda * (dist / pow(deltaEps, 2)) * (gi * deltaS1 + deltaY * deltaKS1);
                }
            }          
        }
        return -total;  
    }

    double MsyFOC(strainState e0kxky) {
        double eyd = strainYieldDesign(fyk, E);
        double fyd = featureYieldDesign(fyk);
        double e0 = e0kxky.e0;
        double kx = e0kxky.kx;
        double ky = e0kxky.ky;   
        double total = 0;
        for (int i = 0; i < steelLinearDist.numPoints; i++) {
            double xI = steelLinearDist.points[i].x;
            double yI = steelLinearDist.points[i].y;            
            double epsI = e0 + ky * xI - kx * yI;
            point *nextPoint = steelLinearDist.points[i].nextPoint;
            if (nextPoint != NULL) {
                double lambda = steelLinearDist.points[i].linearSteelThickness;
                double deltaX = nextPoint->x - xI;
                double deltaY = nextPoint->y - yI;
                double dist = sqrt(pow(deltaX, 2) + pow(deltaY, 2));
                double epsIplus1 = e0 + ky * nextPoint->x - kx * nextPoint->y;
                double deltaEps = epsIplus1 - epsI;
                double S1I = S1(epsI, eyd, fyd);
                double S1Iplus1 = S1(epsIplus1, eyd, fyd);
                double deltaS1 = S1Iplus1 - S1I;
                double hi = xI * epsIplus1 - nextPoint->x * epsI;
                double KS1I = KS1(epsI, eyd, fyd);
                double KS1Iplus1 = KS1(epsIplus1, eyd, fyd);
                double deltaKS1 = KS1Iplus1 - KS1I;
                if (abs(deltaEps) < r) {
                    total += lambda * dist * sigmas(epsI, fyk, E) * (xI + nextPoint->x) / 2;
                } else {
                    total += lambda * (dist / pow(deltaEps, 2)) * (hi * deltaS1 + deltaX * deltaKS1);
                }
            }           
        }
        return total; 
    }

    double NtotFOC(strainState e0kxky) {
        double Nc = NcFOC(e0kxky);
        double Ns = NsFOC(e0kxky);
        return Nc + Ns;
    }

    double MxtotFOC(strainState e0kxky) {
        double Mcx = McxFOC(e0kxky);
        double Msx = MsxFOC(e0kxky);
        return Mcx + Msx;
    }   

    double MytotFOC(strainState e0kxky) {
        double Mcy = McyFOC(e0kxky);
        double Msy = MsyFOC(e0kxky);
        return Mcy + Msy;
    }      

    double** Jc(strainState e0kxky) {
        double** sum2D;
        sum2D = new double*[3];
        for (int i = 0; i < 3; i++) {
            sum2D[i] = new double[3];
            for (int j = 0; j < 3; j++) {
                sum2D[i][j] = 0;
            }
        }

        double e0 = e0kxky.e0;
        double kx = e0kxky.kx;
        double ky = e0kxky.ky;
        // cout << kx << 'k' << 'x' << '\n';
        // cout << ky << 'k' << 'y' << '\n';
        if (abs(kx) > r) {
            for (int i = 0; i < concreteSection.numPoints; i++) {
                double xI = concreteSection.points[i].x;
                double yI = concreteSection.points[i].y;
                double epsI = e0 + ky * xI - kx * yI;
                
                int Iplus1;
                if (i != concreteSection.numPoints - 1) {
                    Iplus1 = i + 1;
                } else {
                    Iplus1 = 0;
                }
                double xIplus1 = concreteSection.points[Iplus1].x;
                double yIplus1 = concreteSection.points[Iplus1].y;                
                double epsIplus1 = e0 + ky * xIplus1 - kx * yIplus1;

                double dx = xIplus1 - xI;
                // cout << dx << 'd' << 'x' << '\n';

                double f6if = f6i(epsIplus1, epsI, fck, r);
                double f7if = f7i(epsIplus1, epsI, fck, yI, yIplus1, r);
                double f8if = f8i(epsIplus1, epsI, fck, xI, xIplus1, r);
                double f9if = f9i(epsIplus1, epsI, fck, yI, yIplus1, r);
                double f10if = f10i(epsIplus1, epsI, fck, yI, yIplus1, xI, xIplus1, r);
                double f11if = f11i(epsIplus1, epsI, fck, xI, xIplus1, r);

                // cout << epsIplus1 << ' ' <<  epsI << ' ' <<  fck << ' ' <<  yI << ' ' <<  yIplus1 << ' ' <<  xI << ' ' <<  xIplus1 << ' ' <<  r << '\n';
                cout.precision(17);
                cout << f6if << ' ' <<  f7if << ' ' <<  f8if << ' ' <<  f9if << ' ' <<  f10if << ' ' <<  f11if << ' ' << '\n';

                double sum[3][3] = {{dx*f6if, -dx*f7if, dx*f8if}, 
                                   {-dx*f7if, dx*f9if, -dx*f10if}, 
                                   {dx*f8if, -dx*f10if, dx*f11if}};

                for (int i = 0; i < 3; i++) {
                    for (int j = 0; j < 3; j++) {
                        sum2D[i][j] += sum[i][j];
                    }
                }
                   
            }
            double Ncf = NcFOC(e0kxky);
            double Mcxf = McxFOC(e0kxky);
            double Mcyf = McyFOC(e0kxky);

            double sub[3][3] = {{0, Ncf, 0}, {Ncf, 2*Mcxf, Mcyf}, {0, Mcyf, 0}};

            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    sum2D[i][j] = (1/kx) * (sum2D[i][j] - sub[i][j]);
                }
            }            
            return sum2D;

        } else if (abs(ky) > r) {
            for (int i = 0; i < concreteSection.numPoints; i++) {
                double xI = concreteSection.points[i].x;
                double yI = concreteSection.points[i].y;
                double epsI = e0 + ky * xI - kx * yI;
                int Iplus1;
                if (i != concreteSection.numPoints - 1) {
                    Iplus1 = i + 1;
                } else {
                    Iplus1 = 0;
                }
                double xIplus1 = concreteSection.points[Iplus1].x;
                double yIplus1 = concreteSection.points[Iplus1].y;                
                double epsIplus1 = e0 + ky * xIplus1 - kx * yIplus1;

                double dy = yIplus1 - yI;

                double f6if = f6i(epsIplus1, epsI, fck, r);
                double f7if = f7i(epsIplus1, epsI, fck, yI, yIplus1, r);
                double f8if = f8i(epsIplus1, epsI, fck, xI, xIplus1, r);
                double f9if = f9i(epsIplus1, epsI, fck, yI, yIplus1, r);
                double f10if = f10i(epsIplus1, epsI, fck, yI, yIplus1, xI, xIplus1, r);
                double f11if = f11i(epsIplus1, epsI, fck, xI, xIplus1, r);

                double sum[3][3] = {{dy*f6if, -dy*f7if, dy*f8if}, 
                                   {-dy*f7if, dy*f9if, -dy*f10if}, 
                                   {dy*f8if, -dy*f10if, dy*f11if}};    

                for (int i = 0; i < 3; i++) {
                    for (int j = 0; j < 3; j++) {
                        sum2D[i][j] += sum[i][j];
                    }
                }           
            
            }
            double Ncf = NcFOC(e0kxky);
            double Mcxf = McxFOC(e0kxky);
            double Mcyf = McyFOC(e0kxky);

            double sub[3][3] = {{0, 0, Ncf}, {0, 0, Mcxf}, {Ncf, Mcxf, 2*Mcyf}};

            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    sum2D[i][j] = (1/ky) * (sum2D[i][j] - sub[i][j]);
                }
            }                  
            return sum2D;
        } else {
            double dc = Dc(e0, fck);
            double result[3][3] = {{dc*Ac(), dc*Sx(), dc*Sy()}, {dc*Sx(), dc*Ixx(), -dc*Ixy()}, {dc*Sy(), -dc*Ixy(), dc*Iyy()}};
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    sum2D[i][j] = result[i][j];
                }
            }               
            return sum2D;
        }
    }

    double** Js(strainState e0kxky) {
        double** sum2D;
        double eyd = strainYieldDesign(fyk, E);
        double fyd = featureYieldDesign(fyk);
        sum2D = new double*[3];
        for (int i = 0; i < 3; i++) {
            sum2D[i] = new double[3];
            for (int j = 0; j < 3; j++) {
                sum2D[i][j] = 0;
            }
        }

        double e0 = e0kxky.e0;
        double kx = e0kxky.kx;
        double ky = e0kxky.ky;

        for (int i = 0; i < steelLinearDist.numPoints; i++) {
            double xI = steelLinearDist.points[i].x;
            double yI = steelLinearDist.points[i].y;
            double lambda = steelLinearDist.points[i].linearSteelThickness;
            double epsI = e0 + ky * xI - kx * yI;
            int Iplus1;
            if (i != concreteSection.numPoints - 1) {
                Iplus1 = i + 1;
            } else {
                Iplus1 = 0;
            }
            double xIplus1 = steelLinearDist.points[Iplus1].x;
            double yIplus1 = steelLinearDist.points[Iplus1].y;  
            double epsIplus1 = e0 + ky * xIplus1 - kx * yIplus1; 
            double deltaX = xIplus1 - xI;
            double deltaY = yIplus1 - yI;            
            double dist = sqrt(pow(deltaX, 2) + pow(deltaY, 2));           
            double deltaEps = epsIplus1 - epsI;            
            double sigmaI = sigmas(epsI, fyk, E);
            double sigmaIplus1 = sigmas(epsIplus1, fyk, E);
            double deltaSigmas = sigmaIplus1 - sigmaI;        
            double DsI = Ds(epsI, fyk, E);
            double S1I = S1(epsI, eyd, fyd);
            double S1Iplus1 = S1(epsIplus1, eyd, fyd);
            double deltaS1 = S1Iplus1 - S1I;
            double sigmaEpsI = sigmaI * epsI;
            double sigmaEpsIplus1 = sigmaIplus1 * epsIplus1;
            double deltaSigmaEps = sigmaEpsIplus1 - sigmaEpsI;
            double sigmaEps2I = sigmaI * pow(epsI, 2);
            double sigmaEps2Iplus1 = sigmaIplus1 * pow(epsIplus1, 2);
            double deltaSigmaEps2 = sigmaEps2Iplus1 - sigmaEps2I;
            double gi = yI * epsIplus1 - yIplus1 * epsI;
            double hi = xI * epsIplus1 - xIplus1 * epsI;
            double KS1I = KS1(epsI, eyd, fyd);
            double KS1Iplus1 = KS1(epsIplus1, eyd, fyd);
            double deltaKS1 = KS1Iplus1 - KS1I;
            double NsEps0I = NsEps0(lambda,DsI,dist,deltaEps,deltaSigmas,r);
            double NsKxI = NsKx(lambda,DsI,dist,deltaEps,deltaSigmas,gi,deltaY,deltaSigmaEps,deltaS1,yI,r); 
            double NsKyI = NsKy(lambda,DsI,dist,deltaEps,deltaSigmas,hi,deltaX,deltaSigmaEps,deltaS1,xI,r);
            double MsxKxI = MsxKx(lambda,DsI,dist,deltaEps,deltaSigmas,gi,deltaY,deltaSigmaEps,deltaSigmaEps2,deltaS1,deltaKS1,yI,r);
            double MsyKyI = MsyKy(lambda,DsI,dist,deltaEps,deltaSigmas,hi,deltaX,deltaSigmaEps,deltaSigmaEps2,deltaS1,deltaKS1,xI,r);
            double MsxKyI = MsxKy(lambda,DsI,dist,deltaEps,deltaSigmas,gi,hi,deltaX,deltaY,deltaSigmaEps,deltaSigmaEps2,deltaS1,deltaKS1,xI,yI,r);

            double sub[3][3] = {{NsEps0I, NsKxI, NsKyI}, {NsKxI, MsxKxI, MsxKyI}, {NsKyI, MsxKyI, MsyKyI}};

            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    sum2D[i][j] += sub[i][j];
                }
            } 
        }

        return sum2D;
    }

    strainState e0kxky (stresses nMxMy) {
        double Nd = nMxMy.normalStress;
        double Mxd = nMxMy.xBendingMoment;
        double Myd = nMxMy.yBendingMoment;
        double sigmacd = stressCompressDesign(fck);

        strainState e0kxky = {0, 0, 0, 0, 0};

        bool over = false;
        int counter = 0;
        double f_checkpoint = __FLT_MAX__;

        while (!over) {
            
            double Nc = NcFOC(e0kxky);
            double Mcx = McxFOC(e0kxky);
            double Mcy = McyFOC(e0kxky);
        
            double Ns = NsFOC(e0kxky);
            double Msx = MsxFOC(e0kxky);
            double Msy = MsyFOC(e0kxky);


            double Nt = Ns + Nc;
            double Mxt = Mcx + Msx;
            double Myt = Mcy + Msy;

            double Ac_ = this->Ac();

            double maxY = - __FLT_MAX__;
            double minY = __FLT_MAX__;

            for (int i = 0; i < concreteSection.numPoints; i++) {
                double yI = concreteSection.points[i].y;
                if (yI > maxY) {
                    maxY = yI;
                }
                if (yI < minY) {
                    minY = yI;
                }
            }

            double h = maxY - minY;
            double f = sqrt(pow(((Nd-Nt)/(sigmacd*Ac_)),2)+pow(((Mxd-Mxt)/(sigmacd*Ac_*h)),2)+pow(((Myd-Myt)/(sigmacd*Ac_*h)),2));

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
                double** Jcf = Jc(e0kxky);
                double** Jsf = Js(e0kxky);

                double** Jt = new double*[3];
                for (int i = 0; i < 3; i++) {
                    Jt[i] = new double[3];
                    for (int j = 0; j < 3; j++) {
                        Jt[i][j] = Jsf[i][j] + Jcf[i][j];
                    }
                }

                double G = determinant(Jt, 3);
                double* dStresses = new double[3];
                dStresses[0] = Nd - Nt;
                dStresses[1] = Mxd - Mxt;
                dStresses[2] = Myd - Myt;

                if (abs(G) > t) {

                    double Dx = determinant(Jt, 3, true, dStresses, 0);
                    double delE0 = Dx / G;
                    double Dy = determinant(Jt, 3, true, dStresses, 1);
                    double delKx = Dy / G;
                    double Dz = determinant(Jt, 3, true, dStresses, 2);
                    double delKy = Dz / G;

                    e0kxky.e0 += delE0;
                    e0kxky.kx += delKx;
                    e0kxky.ky += delKy;

                } else {
                    over = true;
                    e0kxky.safe = false;
                    cout << 'k' << '3' << '\n';
                    return e0kxky; 
                }
            }
            counter++;
            if (counter % 1000 == 0) {
                if (f > 0.9*f_checkpoint) {
                    over = true;
                    e0kxky.safe = false;
                    cout << 'k' <<  '4' << '\n';
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
    
    double* maxminStrainsConcreteSteel(strainState e0kxky, double *extremes) {
        extremes[0] = - __FLT_MAX__, extremes[2] = extremes[0];
        extremes[1] = __FLT_MAX__, extremes[3] = extremes[1];
        
        for (int i = 0; i < concreteSection.numPoints; i++) {
            double x = concreteSection.points[i].x;
            double y = concreteSection.points[i].y;
            double strain = e0kxky.e0 + e0kxky.ky * x - e0kxky.kx * y;
            if (strain > extremes[0]) {
                extremes[0] = strain;
            }
            if (strain < extremes[1]) {
                extremes[1] = strain;
            }
        }

        for (int i = 0; i < steelLinearDist.numPoints; i++) {
            double x = steelLinearDist.points[i].x;
            double y = steelLinearDist.points[i].y;
            double strain = e0kxky.e0 + e0kxky.ky * x - e0kxky.kx * y;
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
        double arr[4];
        double* extremes = maxminStrainsConcreteSteel(e0kxky, arr);
        double u = ecu(fck);
        double c = ec2(fck);
        if (extremes[0] > u){
            return true;
        }
        if ((u - c) * extremes[0] + c * extremes[1] > u * c){
            return true;
        }
        if (extremes[3] < -10){
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

    stresses solicitingStresses = {0, 0, 1, 1, 1};

    strainState e0kxky = {0, 0, 0.00000000, 0.00000000, 0.00000000};

    double** Jcf = mySection.Jc(e0kxky);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            cout << Jcf[i][j] << ' ';
        }
        cout << '\n';
    }

    cout << '\n';

    strainState deform = mySection.e0kxky(solicitingStresses);

    cout << deform.e0 << '\n';
    cout << deform.kx << '\n';
    cout << deform.ky << '\n';

    // cout << mySection.Jc(e0kxky)

    //getchar();
    return 0;
}


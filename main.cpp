//
// Created by Vitto Resnick on 2/7/24.
//
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <fstream>
#include <armadillo>

using namespace std;
using namespace arma;

// Gaussian G(X) centered at X
double G(double x, double X, double alpha, double l){
    return pow(x-X,l)*exp(-alpha*pow(x-X,2));
}

// 1D Overlap Integral Integrand dSxAB
double dSxAB_1D(double x,double XA,double XB,double alpha,double beta,double lA,double lB){
    return G(x,XA,alpha,lA)*G(x,XB,beta,lB);
}

// Numerical Integration: Extended Trapezoidal Rule - Attempt 1 (Works better)
double trap_SxAB_1D(double a, double b, int n,double XA,double XB,double alpha,double beta,double lA,double lB)
{
    double h = (b-a)/(n-1);
// Evaluate endpoints
    double value = 0.5*( dSxAB_1D(a,XA,XB,alpha,beta,lA,lB)+ dSxAB_1D(b,XA,XB,alpha,beta,lA,lB));
// Now the midpoints
    for(int k=2; k < n; k++){
        value+=dSxAB_1D( a + h*(k-1) ,XA,XB,alpha,beta,lA,lB);
    }
    value*=h;
    return value;
}

// Numerical Integration: Extended Trapezoidal Rule - Attempt 2 (Works worse), based on Numerical Recipes
double trapzd_SxAB(double a, double b, double n,double XA,double XB,double alpha,double beta,double lA,double lB){
    double x,tnm,sum,del;
    static float s;
    int it,j;
    if (n == 1) {
        return (s = 0.5*(b-a)*(dSxAB_1D(a,XA,XB,alpha,beta,lA,lB) + dSxAB_1D(b,XA,XB,alpha,beta,lA,lB)));
    } else {
        for (it = 1, j = 1; j < n-1; j++) {
            it <<= 1;
            tnm=it;
            del=(b-a)/tnm; //This is the spacing of the points to be added.
            x = a+0.5*del;
            for (sum = 0.0, j = 1; j <= it; j++, x+=del) {
                sum += dSxAB_1D(x,XA,XB,alpha,beta,lA,lB);
            }
            s = 0.5*(s+(b-a)*sum/tnm); //This replaces s by its refined value.
        }

        return s;
    }
}

// Definite 1D Overlap Integral Evaluation with Bounds of Integration Calculation
double OverlapIntegral1D(double XA,double XB,double alpha,double beta,double lA,double lB){
    // Calculate a & b bounds of integration, exploit gaussian
    double a, b, stdA, stdB;
    // Gaussian: f(x) = a*exp(-(x-b)^2/(2c^2))
    stdA = sqrt(1/(2*alpha));
    stdB = sqrt(1/(2*beta ));

    if ((XA-4*stdA) <= (XB-4*stdB)){
        a = XA-4*stdA;
    }else {
        a = XB-4*stdB;
    }
    if ((XA+4*stdA) >= (XB+4*stdB)){
        b = XA+4*stdA;
    }else {
        b = XB+4*stdB;
    }
    double n = 10000; // Found from checking different n values
    return trap_SxAB_1D(a,b,n,XA,XB,alpha,beta,lA,lB);
}

// Count words in input file to know which question algorithm to run
int count_words(string file_name){
    ifstream inFile;
    inFile.open(file_name);

    string line;
    int numWords = 0;

    while(getline(inFile, line)){
        stringstream lineStream(line);
        while(getline(lineStream, line, ' ')){
            numWords++;
        }
    }

    inFile.close();
    return numWords;
}

// For given shell, find subshells/functions as all combinations of angular momentum components
vector<vector<double>> findTriplets(int L){
    vec bank = linspace<vec>(L, 0,L+1);
    vector<vector<double>> triplets;
    int n = bank.n_elem;
    for         (int i = 0; i < n; i++){
        for     (int j = 0; j < n; j++){
            for (int k = 0; k < n; k++){
                if (bank[i]+bank[j]+bank[k] == L){
                    vector<double> triplet = {bank(i),bank(j),bank(k)};
                    triplets.push_back(triplet);
                }
            }
        }
    }
    return triplets;
}

// Factorial function
int factorial(int n){
    int res = 1,i;
    for (i=1;i<=n;i++){
        res *= i;
    }
    return res;
}

// Double factorial
int dfact(int n){
    int i;double res=1.0;
    for(i=n;i>=1;i-=2){
        res *=i;
    }
    return res;
}

// Bionomial Coefficient Example: m choose n
int binomCoeff(int m, int n){
    return factorial(m)/(factorial(n)*factorial(m-n));
}

// Calculation one of three directional components for SAB
double SxABComp(double XA, double XB, double alpha, double beta, double lA, double lB){

    double P  = exp(-alpha*beta*pow(XA-XB,2)/(alpha + beta)); // Calculate prefactor
    double XP = (alpha*XA + beta*XB)/(alpha + beta);
    double doubleSum = 0; // Initialize double sum

    // Compute  double sound
    for     (int i = 0; i < lA+1; i++){
        double innerSum = 0;
        for (int j = 0; j < lB+1; j++){
            if ((i+j)% 2 == 0){ // Only do even i+j terms
                double summand = binomCoeff(lA,i)*binomCoeff(lB,j)*dfact(i+j-1)*pow(XP-XA,lA-i)*pow(XP-XB,lB-j)/pow(2*(alpha+beta),(i+j)/2);
                innerSum += summand;
            }
        }
        doubleSum += innerSum;
    }
    return P*sqrt(M_PI/(alpha + beta))*doubleSum;
}

int main() {
    string file_name = "/Users/vittor/Documents/CLASSES/SPRING 2024/CHEM_179_HW2/6.txt";

    // Count number of words in file to determine which question to do.
    if (count_words(file_name) < 7){
        // Do question 1

        // Read in input file
        ifstream inputFile(file_name);

        // Throw error if file was not opened correctly
        if (!inputFile) {
            cerr << "Error opening file." << endl;
        }

        // Read in variables
        double XA, XB, alpha, beta, lA, lB; // Initialize variables
        inputFile >> XA >> alpha >> lA;
        inputFile >> XB >> beta  >> lB;
        inputFile.close();       // Close the txt file

        // Calculate overlap integral
        double OI1 = OverlapIntegral1D(XA,XB,alpha,beta,lA,lB);
        cout << "1d numerical overlap integral between Gaussian functions is " << OI1 << endl;
    } else {
        // Do question 2

        // Read in input file
        ifstream inputFile(file_name);

        // Throw error if file was not opened correctly
        if (!inputFile) {
            cerr << "Error opening file." << endl;
        }

        // Read in variables
        double XA, YA, ZA, XB, YB, ZB, alpha, beta, LA, LB; // Initialize variables
        inputFile >> XA >> YA >> ZA >> alpha >> LA;
        inputFile >> XB >> YB >> ZB >> beta  >> LB;
        inputFile.close();       // Close the txt file

        // Calculate number of subshells/functions for the shell
        double NumOrbFuncA1 = (LA + 1)*(LA + 2)/2;
        double NumOrbFuncB2 = (LB + 1)*(LB + 2)/2;

        // Echo input
        cout << "Shell 1 has " << NumOrbFuncA1 << " functions." << endl;
        cout << "This shell info: R(" << XA << ", "<< YA << ", "<<  ZA << "), with angular momentum: " << LA << ", coefficient: " << alpha << endl;
        cout << "Shell 2 has " << NumOrbFuncB2 << " functions." << endl;
        cout << "This shell info: R(" << XB << ", "<< YB << ", "<<  ZB << "), with angular momentum: " << LB << ", coefficient: " << beta  << endl;

        // Find subshells/functions' angular momentum components for each shell
        vector<vector<double>> tripletsA = findTriplets(LA);
        vector<vector<double>> tripletsB = findTriplets(LB);

        // Calculate Overlap Integral Matrix
        mat OverlapIntegrals(NumOrbFuncA1,NumOrbFuncB2, fill::zeros); // Initialize matrix
        for     (int i = 0; i < NumOrbFuncA1; i++){ // Iterate through shell 1's subshells/functions
            for (int j = 0; j < NumOrbFuncB2; j++){ // Iterate through shell 2's subshells/functions

                // Extra angular momentum components for this shell 1 subshell/function
                double lA = tripletsA[i][0], mA = tripletsA[i][1], nA = tripletsA[i][2];
                // Extra angular momentum components for this shell 2 subshell/function
                double lB = tripletsB[j][0], mB = tripletsB[j][1], nB = tripletsB[j][2];

                double SAB = 1; // Initialize SAB product
                SAB *= SxABComp(XA,XB,alpha,beta,lA,lB); // Compute SxAB
                SAB *= SxABComp(YA,YB,alpha,beta,mA,mB); // Compute SyAB
                SAB *= SxABComp(ZA,ZB,alpha,beta,nA,nB); // Compute SzAB
                OverlapIntegrals(i,j) = SAB; // Add 3D overlap integral value to matrix
            }
        }
        // Print Analytical 3D Overlap Integral Matrix
        OverlapIntegrals.print("Overlap integral between Shell 1 and Shell 2");
        // Echo Shell 1 and Shell 2's subshells/functions
        cout << "The components of angular momentum (l, m, n) for the matrix column, from top to bottom, are listed sequentially as: (";
        for (int i = 0; i < tripletsA.size(); i++){
            cout << tripletsA[i][0] << ", " << tripletsA[i][1]  << ", " << tripletsA[i][2];
            if (i != tripletsA.size()-1){
                cout << "), (";
            } else {
                cout << ")." << endl;
            }
        }
        cout << "The components of angular momentum (l, m, n) for the matrix row, from left to right, are listed sequentially as: (";
        for (int i = 0; i < tripletsB.size(); i++){
            cout << tripletsB[i][0] << ", " << tripletsB[i][1]  << ", " << tripletsB[i][2];
            if (i != tripletsB.size()-1){
                cout << "), (";
            } else {
                cout << ")." << endl;
            }
        }
    }
    return 0;
}



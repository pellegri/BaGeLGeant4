
using namespace std;

////////////////////////////////////
////        BiRelKin        ////
////////////////////////////////////
//
//  A calculator for Binary Relativistic Kinematics
//  Kevin C.W. Li, kcwli@sun.ac.za
//  10/02/15


////////////////////////
////    Constants   ////
////////////////////////

////  Speed of Light
double c2 = 931.494;     // MeV/u, c^2
double c4 = c2*c2;  // (MeV/u)^2, c^4

////    Lab Variables
double Q, Etotal;
double theta_lab, phi_lab;

////    Quadratic equation variables
double a = 0.0;
double b = 0.0;
double c = 0.0;


void BiRelKin(double *m, double *T, double *E, double *p, double ThetaSCAT,  double &PhiSCAT, double Ex) {
    
    ////    Q-value calculation
    Q = (m[2] + m[3])*c2 - (m[0] + m[1])*c2; // MeV
    
    ////    Conversion of ejectile scattering angle from degrees to radians
    theta_lab = ThetaSCAT*0.017453292; // radians
    phi_lab = PhiSCAT*0.017453292; // radians

    ////    Initial Total Energy Calculation
    E[0] = T[0] + (m[0]*c2);
    E[1] = T[1] + (m[1]*c2);
    Etotal = E[0] + E[1] + Q - Ex;

    ////    Initial Momentum Calculation
    p[0] = (1/sqrt(c2))*sqrt((E[0]*E[0]) - (m[0]*m[0]*c4));
    p[1] = (1/sqrt(c2))*sqrt((E[1]*E[1]) - (m[1]*m[1]*c4));

    
    ////////////////////////////////////////////////////////////////////////////////
    ////    E[2] is now solved via the solution of a quadratic equation

    a = 4*p[0]*p[0]*c2*cos(theta_lab)*cos(theta_lab) - 4*Etotal*Etotal;
    b = (4*Etotal*Etotal*Etotal) - (4*p[0]*p[0]*c2*Etotal) + (4*m[2]*m[2]*c4*Etotal) - (4*m[3]*m[3]*c4*Etotal);
    c = (2*p[0]*p[0]*c2*Etotal*Etotal) - (2*m[2]*m[2]*c4*Etotal*Etotal) + (2*m[2]*m[2]*c4*p[0]*p[0]*c2) + (2*m[3]*m[3]*c4*Etotal*Etotal) - (2*m[3]*m[3]*c4*p[0]*p[0]*c2) + (2*m[3]*m[3]*c4*m[2]*m[2]*c4) - (Etotal*Etotal*Etotal*Etotal) - (p[0]*p[0]*p[0]*p[0]*c4) - (m[2]*m[2]*m[2]*m[2]*c4*c4) - (m[3]*m[3]*m[3]*m[3]*c4*c4) - (4*m[2]*m[2]*c4*p[0]*p[0]*c2*cos(theta_lab)*cos(theta_lab));

    ////    Total Energy, Kinetic Energy and Momentum of Ejectile
    E[2] = (- b - sqrt((b*b) - (4*a*c)))/(2*a);
    T[2] = E[2] - m[2]*c2;
    p[2] = (1/sqrt(c2))*sqrt((E[2]*E[2]) - (m[2]*m[2]*c4));

    ////    Total Energy, Kinetic Energy and Momentum of Recoil
    E[3] = Etotal - E[2];
    T[3] = E[3] - m[3]*c2;
    p[3] = (1/sqrt(c2))*sqrt((E[3]*E[3]) - (m[3]*m[3]*c4));

    ////    Angle between beam axis and Recoil velocity
    phi_lab = asin((p[2]/p[3])*sin(theta_lab)); // radians
    PhiSCAT = phi_lab/0.0174532925; // deg
    
    bool printResults = false;
    
    if (printResults)
    {
        cout << "///////////////////////////////" << endl;
        cout << "////   Particle Massses    ////" << endl;
        
        cout << "m[0]:    " << m[0] << " u"  << endl;
        cout << "m[1]:    " << m[1] << " u"  << endl;
        cout << "m[2]:    " << m[2] << " u"  << endl;
        cout << "m[3]:    " << m[3] << " u"  << endl;
        
        cout << "///////////////////////" << endl;
        cout << "////   Ejectile    ////" << endl;

        cout << "Total Energy, E[2]:    " << E[2] << " MeV"  << endl;
        cout << "Kinetic Energy, T[2]:    " << T[2] << " MeV"  << endl;
        cout << "Momentum, p[2]:    " << p[2] << " MeV/c"  << endl;

        
        cout << "///////////////////////" << endl;
        cout << "////   Recoil    ////" << endl;
        cout << "Total Energy, E[3]:    " << E[3] << " MeV"  << endl;
        cout << "Kinetic Energy, T[3]:    " << T[3] << " MeV"  << endl;
        cout << "Momentum, p[3]:    " << p[3] << " MeV/c"  << endl;
        
        
        cout << "   " << endl;
        cout << "ThetaSCAT (angle of ejectile w.r.t. beam axis):    " << ThetaSCAT << " deg"  << endl;
        cout << "PhiSCAT (angle of recoil w.r.t. beam axis):    " << PhiSCAT << " deg"  << endl;
    }
    
    
    
    
}

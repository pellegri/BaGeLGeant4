#include "Randomize.hh"
#include "G4ThreeVector.hh"

#include "BiRelKin.hh"
#include "DCS_PDR_1minus.h"
#include "DCS_PDR_2plus.h"

std::vector<std::vector<double>> differentialCrossSectionVec;

//============================================================================================================

bool    setDCSDistribution;
bool    setSingleDCSAngle;

double  range_DCS_AcceptanceRejection;
double  ejectileThetaMin, ejectileThetaMax;
double  singleEjectileTheta;

//============================================================================================================

bool    setGammaDecayDistribution;
bool    setSingleGammaDecayAngle;

int     elecMultipolarity;
double  range_AngularDistribution_GammaDecay;
double  singleGammaDecayTheta;

//============================================================================================================

double Evaluate_CrossSection(double thetaValue);
double Sample_CrossSection_Ejectile();
double Evaluate_DifferentialCrossSection(double thetaValue);

//============================================================================================================

void InitialiseVariables_DCS() {
    
    setDCSDistribution = false;
    setSingleDCSAngle = false;
    
    range_DCS_AcceptanceRejection = 0.0;
    ejectileThetaMin = 0.0;
    ejectileThetaMax = 0.0;
    singleEjectileTheta = 0.0;
    
    DefineDCS_PDF_1minus();
    DefineDCS_PDF_2plus();
}

void InitialiseVariables_GammaDecayDistribution() {
    
    setGammaDecayDistribution = false;
    setSingleGammaDecayAngle = false;
    
    elecMultipolarity = 0;
    range_AngularDistribution_GammaDecay = 0.0;
}

//============================================================================================================

void SetupEventGenerator_DifferentialCrossSection(int distN_DCS, double ejectileTheta_Min, double ejectileTheta_Max) {
    
    InitialiseVariables_DCS();
    
    //--------------------------------------------------------------------
    if(distN_DCS==1)
    {
        setDCSDistribution = true;
        differentialCrossSectionVec = DCS_PDR_2plus;
    }
    else if(distN_DCS==2)
    {
        setDCSDistribution = true;
        differentialCrossSectionVec = DCS_PDR_2plus;
    }
    else
    {
        std::cout << "WARNING: Undefined differential cross section number: " << distN_DCS << std::endl;
    }
    
    ejectileThetaMin = ejectileTheta_Min;
    ejectileThetaMax = ejectileTheta_Max;
    
    //--------------------------------------------------------------------
    int nTestPoints = 1000;
    double interval = (ejectileTheta_Max-ejectileTheta_Min)/(nTestPoints-1);

    for(int i=0; i<nTestPoints; i++)
    {
        double thetaValue = (i*interval) + ejectileTheta_Min;
        //double crossSectionValue = Evaluate_CrossSection((i*interval) + ejectileTheta_Min);
        double crossSectionValue = Evaluate_DifferentialCrossSection(thetaValue);

        if(crossSectionValue>range_DCS_AcceptanceRejection)
        {
            range_DCS_AcceptanceRejection = crossSectionValue;
        }
    }
    
    //--------------------------------------------------------------------
    //      A safety of 5% to ensure the full
    range_DCS_AcceptanceRejection *= 1.05;
}

//============================================================================================================

void SetupEventGenerator_SingleEjectileAngle(double ejectileTheta) {
    
    InitialiseVariables_DCS();
    
    //--------------------------------------------------------------------
    setSingleDCSAngle = true;
    singleEjectileTheta = ejectileTheta;
}

//============================================================================================================

double Evaluate_DifferentialCrossSection(double thetaValue) {
    
    double functionValue = 0.0;
    bool validRange = false;
    
    //--------------------------------------------------------------------
    for(int i=0; i<(int) differentialCrossSectionVec.size()-1; i++)
    {
        if(thetaValue>=differentialCrossSectionVec[i][0] && thetaValue<differentialCrossSectionVec[i+1][0])
        {
            validRange = true;
            
            double x1, x2, y1, y2;
            double m, c;

            x1 = differentialCrossSectionVec[i][0];
            x2 = differentialCrossSectionVec[i+1][0];
            y1 = differentialCrossSectionVec[i][1];
            y2 = differentialCrossSectionVec[i+1][1];

            m = (y2-y1)/(x2-x1);
            c = y2 - (m*x2);
            
            functionValue = (m*thetaValue) + c;
        }
    }
    
    if(!validRange)
    {
        functionValue = 0.0;
    }
    
    return functionValue;
}

//============================================================================================================

double Evaluate_CrossSection(double thetaValue) {
    
    double functionValue = sin(thetaValue*deg)*Evaluate_DifferentialCrossSection(thetaValue);
    
    return functionValue;
}

//============================================================================================================

double Sample_CrossSection_Ejectile() {
    
    double thetaValue = 0.0;
    
    if(setSingleDCSAngle)
    {
        thetaValue = singleEjectileTheta;
    }
    else if(setDCSDistribution)
    {
        //----------------------------
        bool found = false;
        
        while(!found)
        {
            thetaValue = G4RandFlat::shoot(ejectileThetaMin, ejectileThetaMax);
            
            double rangeDist = G4RandFlat::shoot(0.0, range_DCS_AcceptanceRejection);
            double crossSectionValue = Evaluate_CrossSection(thetaValue);
            
            //std::cout << "rangeDist: " << rangeDist << std::endl;
            //std::cout << "crossSectionValue: " << rangeDist << std::endl;

            if(rangeDist<=crossSectionValue)
            {
                found = true;
            }
        }
    }
    else
    {
        thetaValue = 0.0;
    }
    
    return thetaValue;
}

//============================================================================================================

void SetupEventGenerator_AngularDistribution_GammaDecay(int eMultipolarity) {

    InitialiseVariables_GammaDecayDistribution();
    
    elecMultipolarity = eMultipolarity;

    if(elecMultipolarity==-1)
    {
        setGammaDecayDistribution = true;
    }
    else if(elecMultipolarity==1)
    {
        setGammaDecayDistribution = true;
        range_AngularDistribution_GammaDecay = 3.0;
    }
    else if(elecMultipolarity==2)
    {
        setGammaDecayDistribution = true;
        range_AngularDistribution_GammaDecay = (15.0/4.0);
    }
    else
    {
        std::cout << "WARNING: Undefined Electric transition multipolarity: " << eMultipolarity << std::endl;
    }
}

//============================================================================================================

void SetupEventGenerator_SingleAngle_GammaDecay(double angle) {

    InitialiseVariables_GammaDecayDistribution();

    setSingleGammaDecayAngle = true;
    singleGammaDecayTheta = angle;
}

//============================================================================================================

double EvalAngularCorrelation_GammaDecay(double theta) {
    
    double functionValue = 0.0;
    
    if(elecMultipolarity==1)
    {
        functionValue = 3*(1-pow(cos(theta*deg), 2.0)); // m = 0;
        //functionValue += (3.0/2.0)*(1+pow(cos(theta*deg), 2.0)); // m = -1;
        //functionValue += (3.0/2.0)*(1+pow(cos(theta*deg), 2.0)); // m = +1;
    }
    else if(elecMultipolarity==2)
    {
        functionValue = (5.0/2.0)*(6.0*pow(cos(theta*deg), 2.0) - 6.0*pow(cos(theta*deg), 4.0)); // m = 0
        //functionValue += (5.0/2.0)*(1.0 - 3.0*pow(cos(theta*deg), 2.0) + 4.0*pow(cos(theta*deg), 4.0)); // m = -1
        //functionValue += (5.0/2.0)*(1.0 - 3.0*pow(cos(theta*deg), 2.0) + 4.0*pow(cos(theta*deg), 4.0)); // m = +1
        //functionValue += (5.0/2.0)*(1.0 - 4.0*pow(cos(theta*deg), 4.0)); // m = -2
        //functionValue += (5.0/2.0)*(1.0 - 4.0*pow(cos(theta*deg), 4.0)); // m = +2
    }
    
    return functionValue;
}

//============================================================================================================

double EvalAngularDistribution_GammaDecay(double theta) {
    
    double functionValue = sin(theta*deg)*EvalAngularCorrelation_GammaDecay(theta);
    
    return functionValue;
}

//============================================================================================================

double Sample_AngularDistribution_GammaDecay() {
    
    double thetaValue = 0.0;
    
    if(setSingleGammaDecayAngle)
    {
        thetaValue = singleGammaDecayTheta;
    }
    else if(setGammaDecayDistribution)
    {
        //----------------------------
        if(elecMultipolarity==-1)
        {
            thetaValue = acos(1 - (2.0*G4UniformRand()))/deg; // 0.0->180.0
        }
        else
        {
            bool found = false;
            while(!found)
            {
                thetaValue = G4RandFlat::shoot(0.0, 180.0);
                
                double rangeDist = G4RandFlat::shoot(0.0, range_AngularDistribution_GammaDecay);
                
                if(rangeDist<=EvalAngularDistribution_GammaDecay(thetaValue))
                {
                    found = true;
                }
            }
        }
    }
    else
    {
        thetaValue = 0.0;
    }
    
    return thetaValue;
}

//============================================================================================================

void Sample_AngularDistribution_GammaDecay_LAB(double A0, double A1, double A2, double A3, double beamEnergy, double excitationEnergy, double &thetaGamma_LAB, double &thetaGamma_COM, double &gammaEnergy)
{
    //------------------------------------------------
    double m[4], T[4], E[4], p[4];
    double theta_ejectile = 0.0;
    double phi_ejectile = 0.0;
    double theta_recoil = 0.0;
    double phi_recoil = 0.0;
    
    for(G4int i=0; i<4; i++)
    {
        T[i] = 0.0;
        E[i] = 0.0;
        p[i] = 0.0;
    }
    
    //----------------------------------------------
    //      The masses of the binary reaction
    m[0] = A0; // u
    m[1] = A1; // u
    m[2] = A2; // u
    m[3] = A3; // u
    
    //----------------------------
    //      Projectile Energy
    T[0] = beamEnergy; // MeV
    
    //----------------------------
    //      Ejectile Energy
    T[1] = 0.0; // MeV
    
    //----------------------------------------
    //      Choosing the ejectile angles
    theta_ejectile = Sample_CrossSection_Ejectile(); // deg
    phi_ejectile = G4RandFlat::shoot(0.0, 360.0); // deg
    
    //----------------------------------------------------------
    //      Calculating the relativistic binary kinematics
    BiRelKin(m, T, E, p, theta_ejectile, theta_recoil, excitationEnergy);
    
    //------------------------------------------------
    //      Choosing the gamma-ray decay angles
    thetaGamma_COM = Sample_AngularDistribution_GammaDecay(); // deg

    //------------------------------------------------------------
    //      Calculating the kinematics of the recoil nucleus
    double v_recoil = (p[3]/E[3]);
    double beta_recoil = v_recoil/sqrt(c2);
    double gamma_recoil = 1.0/sqrt(1-pow(beta_recoil, 2.0));
    
    double theta_GammaDecay_LAB_relativeToEjectile = acos((cos(thetaGamma_COM*deg) + beta_recoil)/(1 + beta_recoil*cos(thetaGamma_COM*deg)))/deg; // deg
    double phi_GammaDecay_LAB_relativeToEjectile = G4RandFlat::shoot(0.0, 360.0); // deg

    G4double mx, my, mz;
    mx = sin(theta_GammaDecay_LAB_relativeToEjectile*deg)*cos(phi_GammaDecay_LAB_relativeToEjectile*deg);
    my = sin(theta_GammaDecay_LAB_relativeToEjectile*deg)*sin(phi_GammaDecay_LAB_relativeToEjectile*deg);
    mz = cos(theta_GammaDecay_LAB_relativeToEjectile*deg);
    
    G4ThreeVector gamma_Direction_relativeToEjectile(mx, my, mz);
    gamma_Direction_relativeToEjectile.unit();

    //----------------------------------------------------
    //      Determing the angle of the recoil nucleus
    phi_recoil = G4RandFlat::shoot(0.0, 360.0);
    
    G4ThreeVector requiredFinalXaxis, requiredFinalYaxis, requiredFinalZaxis;
    
    if(theta_recoil==0.0)
    {
        requiredFinalZaxis = G4ThreeVector(0.0, 0.0, 1.0);
        requiredFinalXaxis = G4ThreeVector(1.0, 0.0, 0.0);
        requiredFinalYaxis = G4ThreeVector(0.0, 1.0, 0.0);
    }
    else if(theta_recoil==180.0)
    {
        requiredFinalZaxis = -G4ThreeVector(0.0, 0.0, 1.0);
        requiredFinalXaxis = G4ThreeVector(1.0, 0.0, 0.0);
        requiredFinalYaxis = -G4ThreeVector(0.0, 1.0, 0.0);
    }
    else
    {
        requiredFinalZaxis = G4ThreeVector(sin(theta_recoil*deg)*cos(phi_recoil*deg), sin(theta_recoil*deg)*sin(phi_recoil*deg), -sin(theta_recoil*deg)).unit();
        requiredFinalXaxis = requiredFinalZaxis.orthogonal().unit();
        requiredFinalYaxis = requiredFinalZaxis.cross(requiredFinalXaxis).unit();
    }
    
    G4RotationMatrix rotMatrixToLabFrame;
    rotMatrixToLabFrame.rotateAxes(requiredFinalXaxis, requiredFinalYaxis, requiredFinalZaxis);
    
    //----------------------------
    G4ThreeVector gammaDirection_lab = rotMatrixToLabFrame*gamma_Direction_relativeToEjectile;
    
    //----------------------------
    thetaGamma_LAB = acos(gammaDirection_lab.z()/gammaDirection_lab.mag())/deg;
    gammaEnergy = gammaEnergy/(gamma_recoil*(1 + beta_recoil*cos(theta_GammaDecay_LAB_relativeToEjectile)));
}

//============================================================================================================

void Sample_AngularDistribution_GammaDecay_LAB(double beta, double thetaRecoil_LAB, double &thetaGamma_LAB, double &thetaGamma_COM, double &gammaEnergy)
{
    //------------------------------------------------
    double theta_ejectile = 0.0;
    double phi_ejectile = 0.0;
    double theta_recoil = thetaRecoil_LAB;
    double phi_recoil = 0.0;
    
    //----------------------------------------
    //      Choosing the ejectile angles
    theta_ejectile = Sample_CrossSection_Ejectile(); // deg
    phi_ejectile = G4RandFlat::shoot(0.0, 360.0); // deg
    
    //------------------------------------------------
    //      Choosing the gamma-ray decay angles
    thetaGamma_COM = Sample_AngularDistribution_GammaDecay(); // deg
    
    //------------------------------------------------------------
    //      Calculating the kinematics of the recoil nucleus
    double beta_recoil = beta;
    double gamma_recoil = 1.0/sqrt(1-pow(beta_recoil, 2.0));
    double theta_GammaDecay_LAB_relativeToEjectile = acos((cos(thetaGamma_COM*deg) + beta_recoil)/(1 + beta_recoil*cos(thetaGamma_COM*deg)))/deg; // deg
    double phi_GammaDecay_LAB_relativeToEjectile = G4RandFlat::shoot(0.0, 360.0); // deg

    G4double mx, my, mz;
    mx = sin(theta_GammaDecay_LAB_relativeToEjectile*deg)*cos(phi_GammaDecay_LAB_relativeToEjectile*deg);
    my = sin(theta_GammaDecay_LAB_relativeToEjectile*deg)*sin(phi_GammaDecay_LAB_relativeToEjectile*deg);
    mz = cos(theta_GammaDecay_LAB_relativeToEjectile*deg);
    
    G4ThreeVector gamma_Direction_relativeToEjectile(mx, my, mz);
    gamma_Direction_relativeToEjectile.unit();
    
    //----------------------------------------------------
    //      Determing the angle of the recoil nucleus
    phi_recoil = G4RandFlat::shoot(0.0, 360.0);
    
    G4ThreeVector requiredFinalXaxis, requiredFinalYaxis, requiredFinalZaxis;
    
    if(theta_recoil==0.0)
    {
        requiredFinalZaxis = G4ThreeVector(0.0, 0.0, 1.0);
        requiredFinalXaxis = G4ThreeVector(1.0, 0.0, 0.0);
        requiredFinalYaxis = G4ThreeVector(0.0, 1.0, 0.0);
    }
    else if(theta_recoil==180.0)
    {
        requiredFinalZaxis = -G4ThreeVector(0.0, 0.0, 1.0);
        requiredFinalXaxis = G4ThreeVector(1.0, 0.0, 0.0);
        requiredFinalYaxis = -G4ThreeVector(0.0, 1.0, 0.0);
    }
    else
    {
        requiredFinalZaxis = G4ThreeVector(sin(theta_recoil*deg)*cos(phi_recoil*deg), sin(theta_recoil*deg)*sin(phi_recoil*deg), -sin(theta_recoil*deg)).unit();
        requiredFinalXaxis = requiredFinalZaxis.orthogonal().unit();
        requiredFinalYaxis = requiredFinalZaxis.cross(requiredFinalXaxis).unit();
    }
    
    G4RotationMatrix rotMatrixToLabFrame;
    rotMatrixToLabFrame.rotateAxes(requiredFinalXaxis, requiredFinalYaxis, requiredFinalZaxis);
    
    //----------------------------
    G4ThreeVector gammaDirection_lab = rotMatrixToLabFrame*gamma_Direction_relativeToEjectile;

    //----------------------------
    thetaGamma_LAB = acos(gammaDirection_lab.z()/gammaDirection_lab.mag())/deg;
    gammaEnergy = gammaEnergy/(gamma_recoil*(1 + beta_recoil*cos(theta_GammaDecay_LAB_relativeToEjectile)));
}

//============================================================================================================










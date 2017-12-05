int main()
{
    
    gROOT->Reset();
  
    TFile f("Analyse.root","recreate");

    double nrofThreads = 24;

    //char rootFileDirectory[] = "K600-buildMT";
    char rootFileDirectory[] = "K600-buildMT2";
    //char rootFileDirectory[] = "K600-buildMT-1332";
    //char rootFileDirectory[] = "/Users/kevinchingweili/Physics/GEANT4/K600/master/K600-build";
    
    char name[512];
    char condition[512];
    char command[512];

    
    TChain *DataTreeChain = new TChain("DataTreeSim","DataTreeSim");
    TChain *GeometryAnalysisTree = new TChain("GeometryAnalysisTree","GeometryAnalysisTree");
    TChain *InputVariableTree = new TChain("InputVariableTree","InputVariableTree");

    for(Int_t i=0; i<nrofThreads; i++)
    {
        sprintf(name,"K600Output_t%d.root", i);
        //sprintf(name,"%s/K600Output.root", rootFileDirectory, i);
        
        DataTreeChain->Add(name);
        GeometryAnalysisTree->Add(name);
        InputVariableTree->Add(name);
    }
    

    
    TCanvas * c1 = new TCanvas("c1", "c1", 1000, 700);
    
    
    ////================================================////
    
    std::vector<std::shared_ptr<TH1F>> hEnergy_energy;
    std::vector<double> integrals;
    std::vector<double> efficiencies;
    
    for(int i=0; i<9; i++)
    {
        double initialParticleKineticEnergy = 0.0;
        
        if(i==0)
        {
            initialParticleKineticEnergy = 0.25;
        }
        else if(i==1)
        {
            initialParticleKineticEnergy = 0.5;
        }
        else if(i==2)
        {
            initialParticleKineticEnergy = 1;
        }
        else if(i==3)
        {
            initialParticleKineticEnergy = 1.332;
        }
        else if(i==4)
        {
            initialParticleKineticEnergy = 2;
        }
        else if(i==5)
        {
            initialParticleKineticEnergy = 5;
        }
        else if(i==6)
        {
            initialParticleKineticEnergy = 10;
        }
        else if(i==7)
        {
            initialParticleKineticEnergy = 15;
        }
        else if(i==8)
        {
            initialParticleKineticEnergy = 20;
        }

        //initialParticleKineticEnergy *= 1.0e+03;
        
        sprintf(name,"hEnergy_energy%d", i);
        hEnergy_energy.push_back(std::shared_ptr<TH1F>(new TH1F(name, "", 22000, 0.0, 22000.0)));
        
        double lowerE, higherE;
        lowerE = initialParticleKineticEnergy-0.010;
        higherE = initialParticleKineticEnergy+0.010;
        
        int lowerEBin, higherEBin;
        lowerEBin = (int) 1000*lowerE;
        higherEBin = (int) 1000*higherE;
        
        sprintf(condition,"(%f)<InitialParticleKineticEnergy && InitialParticleKineticEnergy<(%f)", lowerE, higherE);
        sprintf(command,"CLOVER_Energy>>%s", name);
        
        DataTreeChain->Draw(command, condition, "");
        //std::cout << "condition: " << condition << std::endl;
        //std::cout << "lowerEBin: " << lowerEBin << std::endl;
        //std::cout << "higherEBin: " << higherEBin << std::endl;

        DataTreeChain->Draw(command, condition, "");

        integrals.push_back(hEnergy_energy.back()->Integral(lowerEBin, higherEBin));
        efficiencies.push_back(100.0*integrals.back()/1000000.0);
        
        std::cout << "Energy: " << initialParticleKineticEnergy << " [MeV], Efficiency: " << efficiencies.back() << std::endl;
    }
    

    f.Write();
    c1->Close();
    
    return 0;
    
}

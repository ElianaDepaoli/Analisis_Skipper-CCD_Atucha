/*
Autora: Eliana Depaoli

Con este código se leen recursivamente los archivos *.root almacenados en el directorio indicado en la ejecución. 

Dentro de cada archivo leído se guarda un nuevo TREE "calc" que contiene los siguientes branches:
SEE (Single Electron Event per pixel per day), 
DEE (Double Electron Event per pixel per day), 
OCCUPANCY (total amount of occupied pixels per pixel per day), 
"oproba" (probabilidad de cada evento del tree hitSumm de ser un neutrino).
SER (Single Electron Event per pixel) y RN (Readout Noise) de los cuadrantes 1 y 2 obtenidos por ajuste de los histogramas de carga en los intervalos de 0 y 1 e-. Errores en SER_1,2 y RN_1,2.
Tiempo de exposición en días reducido por un factor debido a la electura entre imágenes sucesivas.
Cantidad de columnas brillantes extraídas..

La proba se obtiene de la f(Varx,Vary), pdf de la forma de los eventos de neutrinos en la CCD, 
calculada partir de una simulación de eventos de energía 16 < E < 10000 eV uniformemente distribuida, 
distribuidos uniformemente en volumen y luego transportados a la superficie y clusterizados 
con el procesamiento de SENSEI [Mariano Cababié (12/09/2022)].

Compilación: g++ -o addtreetoimage.exe addtreetoimage.cpp `root-config --cflags --glibs`

Ejecución: ./addtreetoimage [ruta/al/directorio/con/los*.root]

Nota importante: en la variable filename_sim hay que poner la ruta completa al archivo 'atucha1.root', que tiene los resultados de las simulaciones necesario para los cálculos que hace la función 'ProcessImage'

Ultima modificacion: 01/06/2025 Branches were added

*/
#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
#include <fstream>
#include <chrono>
#include <thread>
#include <cstdlib>
#include <string>
#include <algorithm> 
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TSystem.h>
#include <TSystemDirectory.h>
#include "TMath.h"
#include "TFitResult.h"
#include "TF1.h"
#include "TH2D.h"

using std::cout; using std::endl; using std::vector; using std::string; using std :: copy; using std::setprecision;
//Function definition ········································
void Enable_and_Set_Branches(TTree* & tree);
void Enable_and_Set_Branches_hS(TTree* & tree);
void Enable_and_Set_Branches_he(TTree* & tree);
void Enable_and_Set_Branches_cp(TTree* & tree);
double imagedurationindays (Char_t DATESTART[0], Char_t DATEEND[0]);
double imagedurationindays (Char_t DATESTART[0], Char_t DATEEND[0], float factor);
void ProcessImage(const char* filename);
void ProcessDirectory(const char* dirname);
void extracthc(int vmin[], int vmax[],int arraysize, bool &isnthotcol, Int_t xm , Int_t xM);

double Nmax{};
//To calulate SER % Readout Noise ·························
Double_t poisson_gauss_conv(Double_t *x, Double_t *par){//lambda [e-/pxl] , mu, RN, normalización

        double xx = x[0];
        double y = 0.0;

        for(int i{0}; i < Nmax+1; ++i){
            double pois = pow(par[0],i)/TMath::Factorial(i);
            double gauss = TMath::Exp(-0.5*pow((xx-i-par[1])/par[3],2)); 
            y += pois*gauss;
        }

        return TMath::Exp(-par[0])*y/(par[3] * TMath::Sqrt(2 * TMath::Pi()))*par[2];//mejora el p-valor si está par[2] (factor de normalización)
}

////////////////////////////////////////////////////////////////////////////////
string filename_sim{"/home/dario/Documentos/Atucha/Scripts/cutcatalog/May_25/atucha1.root"};
double factorVarx{pow(100,0.5)};//este factor arregla que las mediciones tienen un bineado x10 en las columnas mientras que las simulaciones no lo tienen 
//Variances histogram
int nxbin{50};double xl{0.0}; double xu{1.6}; int nybin{50};double yl{0.0}; double yu{1.6};
double pasox{xu/nxbin}; double pasoy{yu/nxbin};
//Probability Histogram
int npbin{25};double pl{0.0}; double pu{0.017};

TH2D * var2D_sim_h = new TH2D("Simulated events", "Simulated events with uniform distribution energy", nxbin, xl, xu, nybin,yl,yu);

//Branches in simulation tree ······
double oVarX; double oVarY;

////Branches in Atucha's trees ······
//hitSumm ·····
Float_t xVar; Float_t yVar; Int_t runID; float e; float n; 

//headerTree_0
Char_t DATESTART [256]; 
Char_t DATEEND [256];
Char_t RUNID_head;

// calPixTree 
int x_cp; int y_cp; int RUNID_cp; int ohdu_cp;
double ePix_cp;//tiene los eventos de 0 e-

//Quality cuts ······
float factor = 0.5;//corrección por lectura entre dos imágenes sucesivas
// Hot Columns definition: Superan el cuantil 0.95 de la distribución de columnas con eventos de 1 e y 1 píxel, y se repiten en varias carpetas de mediciones (+ de 15 días de estadística cada carpeta)
int hotcol_ohdu_1_min[] = {4, 99, 174, 185, 201, 206, 231, 276, 286, 307};
int hotcol_ohdu_1_max[] = {6, 100, 178, 187, 202, 207, 233, 277, 288, 307};
int hotcol_ohdu_2_min[] = {6, 53,  58, 106, 110, 128, 138, 199, 307};
int hotcol_ohdu_2_max[] = {9, 54,  60, 109, 116, 129, 139, 199, 307};
bool isnthotcol;//to remove hot columns
int sizehotcol_1_min{sizeof(hotcol_ohdu_1_min)/sizeof(hotcol_ohdu_1_min[0])};
int sizehotcol_2_min{sizeof(hotcol_ohdu_2_min)/sizeof(hotcol_ohdu_2_min[0])};

//Borders of the active area
int xBaryMin=3; // 
int xBaryMax=305; // 5 on the left from the overscan
int yBaryMin=3; 
int yBaryMax=545;
unsigned int hot_col_amount = 26 + 26;//run <= 31: 20 + 10;
double amount_cutted_pix = (xBaryMax - xBaryMin - hot_col_amount)*20*(yBaryMax - yBaryMin);//amount of pixels ohdu1 + ohdu2

// usefull quadrants
int ohdu_1=1;
int ohdu_2=2;

// SER & Readout noise calculation ········
TH1D * epix_ohdu1_actar_hist = new TH1D("epix active area ohdu1", "epix active area", 50, -0.5, 1.5);
TH1D * epix_ohdu2_actar_hist = new TH1D("epix active area ohdu2", "epix active area", 150, -0.5, 1.5);


////////////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[]) {//argc (argument count) represents the number of command-line arguments
  //How to proceed if only ./addtreeimage.cpp were provided with no address 
  if (argc < 2) {
    cout << "Uso: ./programa directorio_raiz" << endl;
    return 1;
  }
  //

  const char* rootDir = argv[1];
  
  // Procesa el directorio raíz
  ProcessDirectory(rootDir);
  //ProcessImage(rootDir);
  return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ProcessImage(const char* filename) {
    // Abre el archivo ROOT :::::: DO WE NEED THIS STEP? ::::::
    //::::::
    TFile* file = TFile::Open(filename);
    if (!file || file->IsZombie()) {
        std::cout << "Error opening file: " << filename << std::endl;
        return;
    }
    //::::::

    //Openning simulation file
    TFile * f_sim = TFile::Open(filename_sim.c_str(),"READ");//c_str(): string -> array of characters & terminates it with a null character. 
    
    if (!f_sim->IsOpen()) {std::cerr << "Failed to load file" << filename_sim << std::endl;}
    
    //Enabling simulation tree
    TTree * tsim = (TTree*) f_sim->Get("hitSumm");//retrieves the tree from the file
    if(!tsim){
        cout << "Tree hitSumm not found in simulation file" << endl;
        f_sim->Close();
    }

    //Openning Atucha file ···
    TFile * f_atu = TFile::Open(filename,"update");//
    //Enabling Atucha trees
    if (!f_atu->IsOpen()) {std::cerr << "Failed to load file: " << filename << std::endl;}
    if ( f_atu->IsOpen()) {std::cerr << "Loading file:        " << filename << std::endl;}

    //Enabling Atucha trees ····

    TTree * hStree = (TTree*) f_atu->Get("ThitSumm");//retrieves the tree from the file ThitSumm
    if(!hStree){
        cout << "Tree ThitSumm not found in Data file" << endl;
        f_atu->Close();
    }

    //TTree * htree = (TTree*) f_atu->Get("headerTree_0");
    TTree * htree = (TTree*) f_atu->Get("Theader");
    TTree * cptree = (TTree*) f_atu->Get("TcalPixTree");

    //··························
    //NEW TREE & BRANCHES ······
    //··························
    TTree* new_tree = new TTree("calc","calc");
    Double_t oproba{};//probability of being a neutrino
    Double_t SEE{};//[e/time in days/amount of pixels]
    Double_t DEE{};//[e/time in days/amount of pixels]
    Double_t OCCUPANCY{};//[e/time in days/amount of pixels]
    double RN_1{};//[e] Ajuste
    double err_RN_1{};//[e] Ajuste
    double RN_2{};//[e] Ajuste
    double err_RN_2{};//[e] Ajuste
    double SER_1{};//[e/pxl] Ajuste
    double err_SER_1{};//[e/pxl] Ajuste
    double SER_2{};//[e/pxl] Ajuste
    double err_SER_2{};//[e/pxl] Ajuste
    double HOTCOL_amount{};//[#columns]
    double time_expo{};//[days]

    TBranch *proba_new_branch = new_tree->Branch("oproba",&oproba,"oproba/D");
    TBranch *SEE_new_branch = new_tree->Branch("SEE",&SEE,"SEE/D");
    TBranch *DEE_new_branch = new_tree->Branch("DEE",&DEE,"DEE/D");
    TBranch *OCCUPANCY_new_branch = new_tree->Branch("OCCUPANCY",&OCCUPANCY,"OCCUPANCY/D");
    TBranch *RN_1_new_branch = new_tree->Branch("RN_1",&RN_1,"RN_1/D");
    TBranch *RN_2_new_branch = new_tree->Branch("RN_2",&RN_2,"RN_2/D");
    TBranch *err_RN_1_new_branch = new_tree->Branch("err_RN_1",&err_RN_1,"err_RN_1/D");
    TBranch *err_RN_2_new_branch = new_tree->Branch("err_RN_2",&err_RN_2,"err_RN_2/D");
    TBranch *SER_1_new_branch = new_tree->Branch("SER_1",&SER_1,"SER_1/D");
    TBranch *SER_2_new_branch = new_tree->Branch("SER_2",&SER_2,"SER_2/D");
    TBranch *err_SER_1_new_branch = new_tree->Branch("err_SER_1",&err_SER_1,"err_SER_1/D");
    TBranch *err_SER_2_new_branch = new_tree->Branch("err_SER_2",&err_SER_2,"err_SER_2/D");
    TBranch *HOTCOL_amount_new_branch = new_tree->Branch("HOTCOL_amount",&HOTCOL_amount,"HOTCOL_amount/D");
    TBranch *time_expo_new_branch = new_tree->Branch("time_expo",&time_expo,"time_expo/D");

    Enable_and_Set_Branches(tsim);//enable branches in simulation
    Enable_and_Set_Branches_hS(hStree);//hitSumm
    Enable_and_Set_Branches_he(htree);//Header
    Enable_and_Set_Branches_cp(cptree);//calPixTree

    //Filling 2D variances histogram simulation////////////////////////////////////////////////////////////////////////////
    for (Long64_t entry{0}; entry < tsim ->GetEntries();++entry){
        tsim -> GetEntry(entry);
        var2D_sim_h -> Fill(factorVarx*oVarX,oVarY);
    }

    var2D_sim_h->Scale(1/var2D_sim_h->Integral("width")); //normalizing 

    //Calculation in Atucha data using HitSumm ···········
    unsigned int che {0};
    double clusters_1e {};//SEE
    int occupancy {0};
    double clusters_2e {0};//DEE
    Double_t see_aux{};
    Double_t dee_aux{};
    Double_t occ_aux{};
    cout << "hitSumm -> GetEntries() = " << hStree ->GetEntries()  << '\n';

    for (Long64_t entry{0}; entry < hStree ->GetEntries() ;++entry)//100
    {
        hStree -> GetEntry(entry);
        htree ->GetEntry(che);
        if (n == 1.0 && e == 1.0){++clusters_1e;}//n: pixels occupied, e: electrons in cluster
        if (n == 1.0 && e == 2.0){++clusters_2e;}//
        occupancy+=n;
        
        if(entry == (hStree ->GetEntries()-1)){
            cout << "clusters_1e =" << clusters_1e << '\n';
            cout << "clusters_2e =" << clusters_2e << '\n';
            time_expo = imagedurationindays(&DATESTART[0],&DATEEND[0], factor);
            cout << "imageduration in days = " << setprecision(4)<< time_expo << '\n';
            see_aux=clusters_1e/time_expo/amount_cutted_pix;
            dee_aux=clusters_2e/time_expo/amount_cutted_pix;
            occ_aux=occupancy/time_expo/amount_cutted_pix;
            clusters_1e = 0;
            clusters_2e = 0;
            occupancy = 0;
        }

    }

    // Calculation in Atucha data using caPixTre ···········
    for(int j{0}; j < cptree -> GetEntries(); ++j){
        cptree->GetEntry(j);
        //SER & Readout noise ····
        //epix en área activa

        isnthotcol=true;
        if(x_cp < 307 && ePix_cp < 10 && ePix_cp > -1){
            if (ohdu_cp == ohdu_1 || ohdu_cp == ohdu_2){
                if (ohdu_cp==ohdu_1){//events outside hot columns
                    extracthc(hotcol_ohdu_1_min, hotcol_ohdu_1_max,sizehotcol_1_min, isnthotcol, x_cp, x_cp);
                    if (isnthotcol==true){//si el pixel no está en una columna brillante
                        epix_ohdu1_actar_hist->Fill(ePix_cp);
                    }
                } else if (ohdu_cp==ohdu_2){
                    extracthc(hotcol_ohdu_2_min, hotcol_ohdu_2_max,sizehotcol_2_min, isnthotcol, x_cp, x_cp);
                    if (isnthotcol==true){//si el pixel no está en una columna brillante
                        epix_ohdu2_actar_hist->Fill(ePix_cp);
                    }
                }
            }
        }
    }

    //Histograms normalization ·····
    epix_ohdu1_actar_hist->Scale( 1./epix_ohdu1_actar_hist->Integral(),"WIDTH");//histogram describe a probability density function. 
    epix_ohdu2_actar_hist->Scale( 1./epix_ohdu2_actar_hist->Integral(),"WIDTH");
    //··················
    //Fittings      ····
    //··················
    //SER: Poisson convolucionada con gaussiana
    TF1 *pgc = new TF1("poisson_gauss_conv", poisson_gauss_conv, -1, 2,4);//(,, x_lim_inf, x_lim_sup, #parameters)
    Nmax= 1;
    //ohdu 1 ····
    cout << '\n' << endl;
    cout << " --------- OHDU 1 ACTIVE AREA --------- " << endl;
    pgc->SetParameters(0.2,0.03, 1.0, 0.21); //lambda, mu, normalización, sigma. Inicializo.
    auto result02 = epix_ohdu1_actar_hist->Fit("poisson_gauss_conv", "R+S");
    //double pvalue02 = ROOT::Math::chisquared_cdf_c(result02->Chi2(), result02->Ndf());
    SER_1 = result02 -> Parameter(0);
    err_SER_1 = result02 -> ParError(0);
    RN_1 = result02 -> Parameter(3);
    err_RN_1 = result02 -> ParError(3);
    //ohdu 2 ····
    cout << '\n' << endl;
    cout << " --------- OHDU 2 ACTIVE AREA -------- " << endl;
    result02 = epix_ohdu2_actar_hist->Fit("poisson_gauss_conv", "R+S");
    //pvalue02 = ROOT::Math::chisquared_cdf_c(result02->Chi2(), result02->Ndf());
    SER_2 = result02 -> Parameter(0);
    err_SER_2 = result02 -> ParError(0);
    RN_2 = result02 -> Parameter(3);
    err_RN_2 = result02 -> ParError(3);

    //Saving new branches into new tree ······
    for (Long64_t entry{0}; entry < hStree ->GetEntries() ;++entry)//100
    {
        hStree -> GetEntry(entry);
        oproba=(var2D_sim_h->GetBinContent(trunc(xVar/pasox)+1,trunc(yVar/pasoy)+1))*pasox*pasoy;
        SEE = see_aux;
        DEE = dee_aux;
        OCCUPANCY = occ_aux;
        HOTCOL_amount = hot_col_amount;
        new_tree->Fill();
    }

    cout << "SEE = " << SEE << '\n';
    cout << "DEE = " << DEE << '\n';
    cout << "OCCUPANCY = " << OCCUPANCY << '\n';
    cout << HOTCOL_amount << '\n';
    cout << SER_1 << '\n';
    cout << SER_2 << '\n';
    cout << time_expo << 'n';

    new_tree->Write();
    hStree->AddFriend(new_tree);
    hStree->Write();

    //Closing files
    f_sim -> Close();
    f_atu ->Close();
    file ->Close();
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
// Simulations ···
void Enable_and_Set_Branches(TTree* & tree){
  tree->SetBranchStatus("*",0); //disable all branches
  tree->SetBranchStatus("oVarX",1);
  tree->SetBranchStatus("oVarY",1);
  tree->SetBranchAddress ("oVarY",&oVarY);
  tree->SetBranchAddress ("oVarX",&oVarX);
  
  /*tree->SetBranchStatus("eventID",1); 
  tree->SetBranchStatus("hdu",1);
  tree->SetBranchStatus("oEnergy",1);
  tree->SetBranchStatus("oMuX",1);
  tree->SetBranchStatus("oMuY",1);
  tree->SetBranchStatus("oMuXAfterBin",1);
  tree->SetBranchStatus("oMuYAfterBin",1);
  tree->SetBranchStatus("oQ",1);
  tree->SetBranchStatus("oCharge",1);
  tree->SetBranchStatus("oSigma",1);
  tree->SetBranchStatus("flag",1);
  tree->SetBranchAddress ("eventID",&eventID);//links de branch 'eventID' (1st) in tree with the variable 'eventID' where data will be store
  tree->SetBranchAddress ("hdu",&hdu);
  tree->SetBranchAddress ("oEnergy",&oEnergy);
  tree->SetBranchAddress ("oMuX",&oMuX);
  tree->SetBranchAddress ("oMuY",&oMuY);
  tree->SetBranchAddress ("oMuXAfterBin",&oMuXAfterBin);
  tree->SetBranchAddress ("oMuYAfterBin",&oMuYAfterBin);
  tree->SetBranchAddress ("oCharge",&oCharge);
  tree->SetBranchAddress ("oSigma",&oSigma);
  tree->SetBranchAddress ("flag",&oflag);*/
}

//hitSumm ···
void Enable_and_Set_Branches_hS(TTree* & tree){
    tree->SetBranchStatus("*",0); //disable all branches
    //tree->SetBranchStatus("flag",1);
    tree->SetBranchStatus("runID",1);
    tree->SetBranchStatus("xVar",1);
    tree->SetBranchStatus("yVar",1);
    tree->SetBranchStatus("e",1);
    tree->SetBranchStatus("n",1);

    //tree->SetBranchAddress ("flag",&flag);
    tree->SetBranchAddress ("runID",&runID);
    tree->SetBranchAddress ("e",&e);
    tree->SetBranchAddress ("n",&n);
    tree->SetBranchAddress ("xVar",&xVar);
    tree->SetBranchAddress ("yVar",&yVar);
}

//Tree Header
void Enable_and_Set_Branches_he(TTree* & tree)
{
    tree->SetBranchStatus("*",0); //disable all branches
    tree->SetBranchStatus("DATESTART",1); 
    tree->SetBranchStatus("DATEEND",1); 
    tree->SetBranchStatus("RUNID",1);
    tree->SetBranchAddress ("DATESTART",&DATESTART);
    tree->SetBranchAddress ("DATEEND",&DATEEND);
    tree->SetBranchAddress ("RUNID",&RUNID_head);
}

//calPixTree
void Enable_and_Set_Branches_cp(TTree* & tree){
    tree->SetBranchStatus("*",0); //disable all branches
    tree->SetBranchStatus("ohdu",1);
    tree->SetBranchStatus("RUNID",1);
    tree->SetBranchStatus("x",1);
    tree->SetBranchStatus("y",1);
    tree->SetBranchStatus("ePix",1); 

    tree->SetBranchAddress ("RUNID",&RUNID_cp);
    tree->SetBranchAddress ("ohdu",&ohdu_cp);
    tree->SetBranchAddress ("x",&x_cp);
    tree->SetBranchAddress ("y",&y_cp);
    tree->SetBranchAddress ("ePix",&ePix_cp);
}

//Time calculation//////////////////////////////////////////////////////////////////////////////
double imagedurationindays (Char_t DATESTART[0], Char_t DATEEND[0])
{
    double difference;
    struct tm itm {};
    struct tm ftm {};
    strptime(DATEEND,"%Y-%m-%dT%H:%M:%S",&ftm);//converting string to date/time
    strptime(DATESTART,"%Y-%m-%dT%H:%M:%S",&itm);
    time_t idate = mktime(&itm);//converting tm structure to time_t format 
    time_t fdate = mktime(&ftm);//the number of seconds since 00:00, Jan 1 1970 UTC, not counting leap seconds
    difference = difftime(fdate,idate)/60/60/24;//end - start in seconds and then to days
    return difference;
}

//Timming calculation ························································································
double imagedurationindays (Char_t DATESTART[0], Char_t DATEEND[0],float factor)
{
    double difference;
    struct tm itm {};
    struct tm ftm {};
    strptime(DATEEND,"%Y-%m-%dT%H:%M:%S",&ftm);//converting string to date/time
    strptime(DATESTART,"%Y-%m-%dT%H:%M:%S",&itm);
    time_t idate = mktime(&itm);//converting tm structure to time_t format 
    time_t fdate = mktime(&ftm);//the number of seconds since 00:00, Jan 1 1970 UTC, not counting leap seconds
    difference = difftime(fdate,idate)/60/60/24;//end - start in seconds and then to days
    return difference*factor;
}

////Función para calcular SER ·························
/*Double_t poisson_gauss_conv(Double_t *x, Double_t *par){//lambda [e-/pxl] , mu, RN, normalización

        double xx = x[0];
        double y = 0.0;

        for(int i{0}; i < Nmax+1; ++i){
            double pois = pow(par[0],i)/TMath::Factorial(i);
            double gauss = TMath::Exp(-0.5*pow((xx-i-par[1])/par[3],2)); 
            y += pois*gauss;
        }

        return TMath::Exp(-par[0])*y/(par[3] * TMath::Sqrt(2 * TMath::Pi()))*par[2];//mejora el p-valor si está par[2] (factor de normalización)
}
*/
// Extracting hot columns ·················································································
void extracthc(int vmin[], int vmax[],int arraysize, bool &isnthotcol, Int_t xm , Int_t xM){
    int j = 0;
    while (isnthotcol == true && j < arraysize){
        if (xM >= vmin[j] && xm <= vmax[j]) isnthotcol = 0;
        ++j;
    }
}

//
void ProcessDirectory(const char* dirname) {
  // Abre el directorio
  TSystemDirectory dir(dirname, dirname);
  
  // Recorre los elementos del directorio
  TList* files = dir.GetListOfFiles();
  if (files) {
    TSystemFile* file;
    TString filename;
    TIter next(files);
    while ((file = dynamic_cast<TSystemFile*>(next()))) {
      filename = file->GetName();
      
      if (!file->IsDirectory() && filename.EndsWith(".root") && filename != "." && filename != "..") {
        ProcessImage(Form("%s/%s", dirname, filename.Data()));
        cout<<"Agregando branches al archivo: "<<Form("%s/%s", dirname, filename.Data())<<endl<<endl;
      }
    }
  }
}
///////////////////////////////////////////////////////////////////////////////
//Explicaciones
//g++ -o addtreetoimage.exe addtreetoimage.cpp `root-config --cflags --glibs`
//-o <file>                Place the output into <file>.
//Without -o addtreetoimage.cpp 'd be compiled to "a.out" file 

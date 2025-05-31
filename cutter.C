/* 
Script to cut root files produced by skExtracct and create a new file reduced according to a set of quality cuts

-> clusters registered in quadrants with single electron counting (ohdu 1 & ohdu 2)
-> clusters outside the edges of CCD
-> clusters outside hot columns
-> clusters with variance lower than a selected limit

The files are stored in a folder created by the script named out_dirname = "/Cutted_sin_extra_branches_new_hot_col".

Authors: Dario Rodrigues & Eliana Depaoli
Last modification: 
-> August 22th, 2024: Data type of variables bleedX & bleedY switched to Float
-> May 27th, 2025: New hot columns added. The script must be into the folder containing the files to cut. 
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
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TObject.h"
#include <TBranch.h>
#include <TSystem.h>
#include <TSystemDirectory.h>
#include "TH2D.h"
//using std::cout; using std::vector; using std::string; using std :: copy;
using namespace std;

// Input Directory name //////////////////////////////////////////////
const char* script_path = gSystem->Which(".", __FILE__);  // Obtiene la ruta del script
const char* main_dirname = gSystem->DirName(script_path); // Extrae el directorio
const char* sub_dirname = "/hits";

// Output filename ///////////////////////////////////////////////////
const char* out_dirname = "/Cutted_calPixTree_HOT_COL_new_EDGES_xVarMin_0_xVarMax_4_yVarMin_0yVarMax_4"; //"/Cutted_sin_extra_branches"

// Trees
string treename_cp = "calPixTree";
string treename_he = "headerTree_0";
string treename_hs = "hitSumm";
//variables to store data from the leaves
//hitSumm ···
Int_t runID; Int_t ohdu; Int_t flag; 
Int_t xMin; Int_t xMax;
Int_t yMin; Int_t yMax; 
Float_t bleedX; Float_t bleedY;
Float_t ddistance; 
Float_t clusterDist; float e; float n; 
Float_t xBary;Float_t yBary;
Float_t xVar; Float_t yVar;
Float_t eRaw; Float_t xy1Var; Float_t xy2Var;
Float_t xy1Width; Float_t xy2Width;
Int_t nSavedPix;
Int_t expoStart;
const Int_t kMaxTrack = 2e8;//nSavedPix;//
Float_t ePix[kMaxTrack];//[]-size-array all with null value//epix upper value according to the TBrowser b plot of this leaf 
Int_t xPix[kMaxTrack];//[]-size-array all with null value//xPix upper value according to the "TBrowser b" plot of this leaf 
Int_t yPix[kMaxTrack];//[]-size-array all with null value//yPix upper value = 550; it doesn't work

//headerTree_0 ···
Char_t DATESTART [256]; 
Char_t DATEEND [256];
Char_t RUNID_head;

// calPixTree ···
int x_cp; int y_cp; int RUNID_cp; int ohdu_cp;
double ePix_cp;//tiene los eventos de 0 e-
//Quality cuts ///////////////////////////////////////////////////

// Hot Columns definition ///////////
// Superan el cuantil 0.95 de la distribución de columnas con eventos de 1 e y 1 píxel, y se repiten en varias carpetas de mediciones (+ de 15 días de estadística cada carpeta) 
// 1 electron events

int hotcol_ohdu_1_min[] = {4, 99, 174, 185, 201, 206, 231, 276, 286, 307};
int hotcol_ohdu_1_max[] = {6, 100, 178, 187, 202, 207, 233, 277, 288, 307};
int hotcol_ohdu_2_min[] = {6, 53,  58, 106, 110, 128, 138, 199, 307};
int hotcol_ohdu_2_max[] = {9, 54,  60, 109, 116, 129, 139, 199, 307};
/*
int hotcol_ohdu_1_min[] = {3, 99, 174, 201, 206, 231, 276, 286, 307};
int hotcol_ohdu_1_max[] = {6, 99, 177, 202, 206, 232, 277, 288, 307};
int hotcol_ohdu_2_min[] = {6, 53,  58, 199, 307};
int hotcol_ohdu_2_max[] = {9, 54,  60, 199, 307};
*/
bool isnthotcol;//to remove hot columns
int sizehotcol_1_min{sizeof(hotcol_ohdu_1_min)/sizeof(hotcol_ohdu_1_min[0])};
int sizehotcol_2_min{sizeof(hotcol_ohdu_2_min)/sizeof(hotcol_ohdu_2_min[0])};

// remove events partially inside the active area
int xBaryMin=3; // 
int xBaryMax=305; // 5 on the left from the overscan
int yBaryMin=3; 
int yBaryMax=545;

// usefull quadrants
int ohdu_1=1;
int ohdu_2=2;

// Event selection ////////////////////////////////////////////////
double xVarMin=pow(0.0,2);
double xVarMax=pow(2.0,2);

//selection of bulk events
double yVarMin=pow(0.0,2);
double yVarMax=pow(2.0,2);

//////////////////////////////////////////////////////////////////////////////////////////
void Enable_and_Set_Branches(TTree* & tree);
void Enable_and_Set_Branches_2(TTree* & tree);
void Enable_and_Set_Branches_3(TTree* & tree);

//////////////////////////////////////////////////////////////////////////////////////////
void extracthc(int vmin[], int vmax[],int arraysize, bool &isnthotcol, Int_t xm , Int_t xM){
    int j = 0;
    while (isnthotcol == true && j < arraysize){
        if (xM >= vmin[j] && xm <= vmax[j]) isnthotcol = 0;
        ++j;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
void ProcessImage(const char* filename,const char* output_filename){

  TFile * f_exp = TFile::Open(filename);//open root file

  if (!f_exp->IsOpen()) {std::cerr << "ERROR: cannot open the root file with experimental data" << std::endl;}

  //Trees to open ······················································
  TTree * hStree = (TTree*) f_exp->Get("hitSumm");
  TTree * htree = (TTree*) f_exp->Get("headerTree_0");
  TTree * cptree = (TTree*) f_exp->Get(treename_cp.c_str());

  TFile * f_exp_cutted = new TFile(output_filename, "RECREATE");//creates root file to save new trees
  
  //Creates trees to fill
  TTree * t_hitSumm = new TTree("ThitSumm", "T_hitSumm");
  TTree * t_Head = new TTree("Theader","T_header");
  TTree * t_calPixTree = new TTree("TcalPixTree", treename_cp.c_str());

  Int_t Entries_hStree = hStree -> GetEntries();
  Int_t Entries_htree = htree -> GetEntries(); 
  Int_t Entries_cptree = cptree -> GetEntries(); 

  Enable_and_Set_Branches(hStree);
  Enable_and_Set_Branches_2(htree);
  Enable_and_Set_Branches_3(cptree);

  //Creating branches in Tree
  //hitSumm ---
  t_hitSumm->Branch ("runID",&runID,"runID/I");
  t_hitSumm->Branch("ohdu",&ohdu,"ohdu/I");
  t_hitSumm->Branch("expoStart",&expoStart,"expoStart/I");
  t_hitSumm->Branch ("flag",&flag,"flag/I");
  t_hitSumm->Branch ("xMin",&xMin,"xMin/I");
  t_hitSumm->Branch ("xMax",&xMax,"xMax/I");
  t_hitSumm->Branch ("yMin",&yMin,"yMin/I");
  t_hitSumm->Branch ("yMax",&yMax,"yMax/I");
  t_hitSumm->Branch ("eRaw",&eRaw,"eRaw/F"); 
  t_hitSumm->Branch ("distance",&ddistance,"distance/F");
  t_hitSumm->Branch ("bleedX",&bleedX,"bleedX/I");
  t_hitSumm->Branch ("bleedY",&bleedY,"bleedY/I");
  t_hitSumm->Branch ("clusterDist",&clusterDist,"clusterDist/F");
  t_hitSumm->Branch ("e",&e,"e/F");
  t_hitSumm->Branch ("n",&n,"n/F");
  t_hitSumm->Branch ("xBary",&xBary,"xBary/F");
  t_hitSumm->Branch ("yBary",&yBary,"yBary/F");
  t_hitSumm->Branch ("xVar",&xVar,"xVar/F");
  t_hitSumm->Branch ("yVar",&yVar,"yVar/F");
  t_hitSumm->Branch ("xy1Var",&xy1Var,"xy1Var/F");
  t_hitSumm->Branch ("xy2Var",&xy2Var,"xy2Var/F");
  t_hitSumm->Branch ("xy1Width",&xy1Width,"xy1Width/F");
  t_hitSumm->Branch ("xy2Width",&xy2Width,"xy2Width/F");
  t_hitSumm->Branch ("nSavedPix",&nSavedPix,"nSavedPix/I");
  t_hitSumm->Branch ("yPix",&yPix,"yPix[nSavedPix]/I");
  t_hitSumm->Branch ("xPix",&xPix,"xPix[nSavedPix]/I");
  t_hitSumm->Branch ("ePix",&ePix,"ePix[nSavedPix]/F");
  
  //header ---
  t_Head->Branch("DATESTART",&DATESTART,"string/C");
  t_Head->Branch("DATEEND",&DATEEND,"string/C");
  t_Head->Branch("RUNID",&RUNID_head,"string/C");
  
  // calPixTree --
  t_calPixTree->Branch ("runID",&RUNID_cp,"runID/I");
  t_calPixTree->Branch("ohdu",&ohdu_cp,"ohdu/I");
  t_calPixTree->Branch ("x",&x_cp,"x/I");
  t_calPixTree->Branch ("y",&y_cp,"y/I");
  t_calPixTree->Branch ("ePix",&ePix_cp,"ePix/D");
  
  ///////////////////////////////////////////////////////////////////////////
  //Saving in HitSumm Tree ···
  for(Int_t i_event=0;i_event<Entries_hStree; i_event++){
    hStree->GetEntry(i_event);//fill all objects linked to branches with status 1
    isnthotcol=true;
      
    if (ohdu==ohdu_1 || ohdu==ohdu_2){//events in quadrants with single electron counting

      //cout << " acá estoy " << endl;
      if (xBary<xBaryMax && xBary>xBaryMin && yBary<yBaryMax && yBary>yBaryMin){//events inside edges of CCD

        if (ohdu==ohdu_1){//events outside hot columns
            extracthc(hotcol_ohdu_1_min, hotcol_ohdu_1_max,sizehotcol_1_min, isnthotcol, xMin, xMax);
        } else if (ohdu==ohdu_2){extracthc(hotcol_ohdu_2_min, hotcol_ohdu_2_max,sizehotcol_2_min, isnthotcol, xMin, xMax);}
        
        if (isnthotcol==true){
          if(xVarMin<=xVar && xVar<=xVarMax && yVarMin<=yVar && yVar<=yVarMax) {
            t_hitSumm->Fill();
          }
        }
      }
    }
  }

  //Saving in Header Tree ···
  for(int j=0;j<Entries_htree;j++)  {
    htree->GetEntry(j);
    t_Head->Fill();
  }

  //Saving in  calPixTree ···
  for(int j=0;j<Entries_cptree; j++){
    cptree->GetEntry(j);//fill all objects linked to branches with status 1
    isnthotcol=true;
      
    if (ohdu_cp == ohdu_1 || ohdu_cp == ohdu_2){
        if (ohdu_cp==ohdu_1){//events outside hot columns
            extracthc(hotcol_ohdu_1_min, hotcol_ohdu_1_max,sizehotcol_1_min, isnthotcol, x_cp, x_cp);
        } else if (ohdu_cp==ohdu_2){extracthc(hotcol_ohdu_2_min, hotcol_ohdu_2_max,sizehotcol_2_min, isnthotcol, x_cp, x_cp);}

        if (isnthotcol==true){t_calPixTree->Fill();}
    }
  
  }


  f_exp_cutted->Write();
  f_exp_cutted->Close();
  f_exp->Close();


}

//////////////////////////////////////////////////////////////////////////////////////////
void cutter() {

  //Carpetas salida
  std::string fulloutdirname = std::string(main_dirname) + out_dirname;
  const char* cfulloutdirname = fulloutdirname.c_str();
  cout << "cfulloutdirname " << cfulloutdirname << endl;
  gSystem->mkdir(cfulloutdirname);//gsystem allows to execute system command Creo la carpeta donde se guardan los archivos nuevos

  //Carpetas entrada
  std::string fulldirname = std::string(main_dirname) + sub_dirname;
  const char* cfulldirname = fulldirname.c_str();
  TSystemDirectory dir(cfulldirname, cfulldirname);
  
  //TDirectory *dir1 = parentdir->Get(out_dirname); if (!dir1) dir1 = parentdir->mkdir(out_dirname,out_dirname);

  // Recorre los elementos de la carpeta donde están los archivos a recortar
  TList* files = dir.GetListOfFiles();
  if (files) {
    TSystemFile* file;
    TString filename;
    TIter next(files);
    while ((file = dynamic_cast<TSystemFile*>(next()))) {
      filename = file->GetName();
      cout << "filename " << filename << endl;

      if (!file->IsDirectory() && filename.EndsWith(".root") && filename.Contains("hits_corr_proc_run")  && filename != "." && filename != "..") {
        
        ProcessImage(Form("%s/%s", cfulldirname, filename.Data()),Form("%s/Cutted_%s", cfulloutdirname, filename.Data()));
        cout<<"Creating file: "<<Form("%s/Cutted_%s", cfulloutdirname, filename.Data())<<endl;
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////////////////
//Turn on a leaf and link it with a named-after-it variable

//////////////////////////////////////////////////////////////////////////////////////////
//  hitsumm
//////////////////////////////////////////////////////////////////////////////////

void Enable_and_Set_Branches(TTree* & tree){

  //tree->SetBranchStatus("*",0); //disable all branches

  //tree->SetBranchStatus("expoStart",1);
  //tree->SetBranchAddress ("nSat",&nSat);
  tree->SetBranchAddress ("runID",&runID);
  tree->SetBranchAddress ("ohdu",&ohdu);
  tree->SetBranchAddress ("flag",&flag);
  tree->SetBranchAddress ("xMin",&xMin);
  tree->SetBranchAddress ("xMax",&xMax);
  tree->SetBranchAddress ("yMin",&yMin);
  tree->SetBranchAddress ("yMax",&yMax); 
  tree->SetBranchAddress ("distance",&ddistance);
  tree->SetBranchAddress ("eRaw",&eRaw);
  tree->SetBranchAddress ("xy1Var",&xy1Var); 
  tree->SetBranchAddress ("xy2Var",&xy2Var);
  tree->SetBranchAddress ("xy1Width",&xy1Width); 
  tree->SetBranchAddress ("xy2Width",&xy2Width);
  tree->SetBranchAddress ("bleedX",&bleedX);
  tree->SetBranchAddress ("bleedY",&bleedY);
  tree->SetBranchAddress ("clusterDist",&clusterDist);
  tree->SetBranchAddress ("e",&e);
  tree->SetBranchAddress ("n",&n);
  tree->SetBranchAddress ("xBary",&xBary);
  tree->SetBranchAddress ("yBary",&yBary);
  tree->SetBranchAddress ("xVar",&xVar);
  tree->SetBranchAddress ("yVar",&yVar);
  tree->SetBranchAddress ("nSavedPix",&nSavedPix);
  tree->SetBranchAddress ("xPix",&xPix);
  tree->SetBranchAddress ("yPix",&yPix);
  tree->SetBranchAddress ("ePix",&ePix);
  tree->SetBranchAddress ("expoStart",&expoStart);
}

//////////////////////////////////////////////////////////////////////////////////////////
//  Header
//////////////////////////////////////////////////////////////////////////////////
void Enable_and_Set_Branches_2(TTree* & tree){
  tree->SetBranchStatus("*",0); //disable all branches
  tree->SetBranchStatus("DATESTART",1); 
  tree->SetBranchAddress ("DATESTART",&DATESTART);
  tree->SetBranchStatus("DATEEND",1); 
  tree->SetBranchAddress ("DATEEND",&DATEEND);
  tree->SetBranchStatus("RUNID",1); 
  tree->SetBranchAddress ("RUNID",&RUNID_head);
}

//////////////////////////////////////////////////////////////////////////////////////////
//  calPixTree
//////////////////////////////////////////////////////////////////////////////////
void Enable_and_Set_Branches_3(TTree* & tree){
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
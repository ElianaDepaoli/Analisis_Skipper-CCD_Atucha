/*
Nota importante: Tener cuidado con que en esta version se harcodeo un factor 2 en la exposure de reactor OFF 2023
para corregir el hecho de que en 2023 hubo limpieza entre imagenes.

Este archivo debe estar en la carpeta que contiene a las subcarpetas con los root procesados por cutter.C y addtreetoimage.cpp

Autores: Eliana Depaoli y Dario Rodrigues
Con este codigo se grafican histogramas de energia. 
Los eventos fueron seleccionados con cutter.C
Los branches SEE, DEE, ocupancia de la imagen y probabilidad de cada evento de ser neutrino fueron agregados con addtreetoimage.cpp

Se eliminaron eventos con:
1) Algun pixel dentro de las columnas brillantes 
2) (x,y) Baricentro dentro de los bordes de la CCD
3) 0 <= xVar <= 2 ; 0 <= yVar <= 2 

Se comparan los espectros producidos con imagenes en distintos intervalos de SEE.
Se comparan los espectros producidos con eventos con distinta probabilidad de provenir de una dstribución 
de energia "tipo nu". La distribucion (VarX, Vary) se obtuvo mediante simulaciones.

Se grafican histogramas antes y despues de los cortes.
Además se grafica SEE vs tiempo. 
Se etiqueta el eje horizontal con la fecha de la imagen de la cual proviene la corriente SEE. 
Además se genera un histograma de SEE.

Ultima modificacion: 04/12/2024 Factor de exposición adentro de la función que genera los histogramas
20 de diciembre de 2023 ---------> ON2022 vs OFF2023
01 de junio de 2025 Adaptación a nuevos branches en tree calc
Tiempo de exposición en días reducido por un factor de exposición debido a la electura entre imágenes sucesivas.

Para generar una lista de los archivos ordenada por número de run y luego por número de imagen
ls | awk -F'[_]' '{split($0,a,"run_"); split($0,b,"img"); print a[2], b[2], $0}' | sort -n -k1,1 -k2,2 | cut -d' ' -f3- > OFF_2023.txt

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
#include <sstream>
#include <TROOT.h>
#include "TFile.h"
#include "TTree.h"
#include "TObject.h"
#include "TString.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TLine.h"
#include "TLegend.h"
#include "TSystem.h"
#include "TGaxis.h"
#include "TObjString.h"
#include "TH2D.h"
#include "TGraphErrors.h"
#include "TChain.h"
#include <time.h>
using std::cout; using std::vector; using std::string; using std :: copy; using std :: cerr;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void concatenateRootTrees(const string & filelist, const string & path, TChain & chain, const string & treename);
double imagedurationindays (Char_t DATESTART[0], Char_t DATEEND[0]);
double imagedurationindays (Char_t DATESTART[0], Char_t DATEEND[0],float factor);
void Enable_and_Set_Branches_hs(TChain & tree);
void Enable_and_Set_Branches_he(TChain & tree);
void Enable_and_Set_Branches_calc(TChain & tree);

void generate_energy_histo_with_quality_cuts(TChain & texp,TChain & header, TChain &calc, TH1D*h_e_bulk, 
	double calibration, double & total_day_time_bc, double & total_day_time_ac, 
    vector<string> & startdate_vec, vector <double> &SEE_v, float see_min, float see_max, double exposure_factor);

void SEE_evolution_plot(vector<string> vecx, vector<Double_t> vecy,const int paso, const char *canv_name,const char *canv_title);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////Branches in Atucha tree////////////////////////////////////////////////////////////////
//hitSumm
Int_t runID; Int_t ohdu; Int_t flag; 
Int_t xMin; Int_t xMax;
Int_t yMin; Int_t yMax;

Float_t ddistance; 
//Float_t bleedX; Float_t bleedY;
Int_t bleedX; Int_t bleedY;

Float_t clusterDist; float e; float n; 
Float_t xBary;Float_t yBary;
Float_t xVar; Float_t yVar;
Float_t eRaw; Float_t xy1Var; Float_t xy2Var;
Float_t xy1Width; Float_t xy2Width;

Int_t nSavedPix;
Int_t expoStart;
const Int_t kMaxTrack = 2e7;//nSavedPix;//
Float_t ePix[kMaxTrack]; //[]-size-array all with null value//epix upper value according to the TBrowser b plot of this leaf 
Int_t xPix[kMaxTrack];   //[]-size-array all with null value//xPix upper value according to the "TBrowser b" plot of this leaf 
Int_t yPix[kMaxTrack];   //[]-size-array all with null value//yPix upper value = 550; it doesn't work*/

//header
Char_t DATESTART [512]; 
Char_t DATEEND [512];
Char_t RUNID_head;

//calc Tree
Double_t SEE{};//see/time in days/amount of pixels
Double_t DEE{};
Double_t OCCUPANCY{};
Double_t oproba{};//probability of being a neutrino
Double_t SER_1{};
Double_t SER_2{};
Double_t RN_1{};
Double_t RN_2{};
Double_t HOTCOL_amount{};
Double_t time_expo{};

// Quality cuts
// Hot Columns definition: Superan el cuantil 0.95 de la distribución de columnas con eventos de 1 e y 1 píxel, y se repiten en varias carpetas de mediciones (+ de 15 días de estadística cada carpeta)

////Edges of the active area
int xBaryMin=3; // 
int xBaryMax=305; // 5 on the left from the overscan
int yBaryMin=3; 
int yBaryMax=509;

// usefull quadrants
int ohdu_1=1;
int ohdu_2=2;

// Event selection ////////////////////////////////////////////////
double xVarMin_GM = pow(0.0,2);
double xVarMax_GM = pow(0.3,2); //1.6

//selection of bulk events
//0.26<yVar<1.1 se quda con el 80% de los eventos de acuerdo a la simulacion en atucha1.root
//double yVarMin_GM = 0.26; 
//double yVarMax_GM = 1.10; 

double yVarMin_GM = pow(0.3,2); // 0.3 removes serial-register events
double yVarMax_GM = pow(1.0,2);

// An extra quality cuts on compatibility with neutrino geometry
double proba_min = 0.0;

// Expousure factor between images with cleanning in between
double exposure_factor_A = 2.0;
double exposure_factor_B = 2.0;

//Calculation////////////////////////////////////////////////////////////////////////////////////////////
double calibration = 1.0163134;//correction with 8.048 keV peak from Cu fluorescence
unsigned int hot_col_amount = 26 + 26; //esto está en los archivos root, sería mejor obtenerlo de allí
double amount_cutted_pix = 2*(xBaryMax - xBaryMin - hot_col_amount)*(yBaryMax - yBaryMin)*10; //amount of pixels ohdu1 + ohdu2

//Timming
vector<string> startdate_A; //to save DATESTART from images for plotting purposes
vector<string> startdate_B;
double total_day_time_bseec{};
double total_day_time_aseec{};

//mass
double rho_Si = 0.00233; // kg/cm3
double height = 0.06750; // [cm]
double row_num = 1024/2;
double col_num = 6144;
double CCDNPRES = 7;
// [cm2] 1024*15e-4*6144*15e-4 = 14.16 cm2
// decidimos no restar hot_col_num porque esto deberia arreglarse con la eficiencia
double area = (xBaryMax - xBaryMin)*20*(yBaryMax - yBaryMin)*15*15*1e-8; 
double mass_in_kg = area*height*rho_Si; //0.70/1000;

/////////////////////////////////////////////////////////////////////////////
//ENERGY
//float eminbin = 0.015; //[keV]
//float emaxbin = 0.215; //[keV]

float eminbin = 0.00;  //[keV]
float emaxbin = 1.00;//0.340; //[keV]

int ebines = 17; 

double x_error_binsize = (emaxbin-eminbin)/ebines/2; // in keV
double bin_size_in_keV = x_error_binsize*2;

TH1D * e_bulk_A_hist = new TH1D("Data SET A", "Reactor A ", ebines, eminbin, emaxbin);
TH1D * e_bulk_B_hist = new TH1D("Data SET B", "Reactor B ", ebines, eminbin, emaxbin);

vector<double> SEE_A; vector<double> SEE_B;
vector<double> SER_A; vector<double> SER_B;
vector<double> RN_A; vector<double> RN_B;

////////////////////////////////////////////////////////////////////////////////////////////////////
void fill_histo_with_asymetric_errors(int numbin, double to_kdru, double efficiency, TH1D *h, Double_t (x_p)[],Double_t (y_p)[],Double_t (e_x)[],Double_t (e_y_u)[],Double_t (e_y_d)[])
{
    for (int j{0}; j<numbin;++j) {

        x_p[j]=h->GetXaxis()->GetBinCenter(j+1); //GetBinLowEdge(j+1);
        y_p[j]=h->GetBinContent(j+1);
        e_x[j]=x_error_binsize;
       
        if(y_p[j]>0 && y_p[j] < 1e2) {
	        e_y_d[j]=h->GetBinContent(j+1) - ROOT::MathMore::chisquared_quantile(0.16,2*(h->GetBinContent(j+1)))/2;
	        e_y_u[j]=ROOT::MathMore::chisquared_quantile(1-0.16,2*(h->GetBinContent(j+1))+2)/2 - h->GetBinContent(j+1);
        }else if(y_p[j]>=1e2){
            e_y_u[j]=pow(y_p[j],0.5);
            e_y_d[j]=pow(y_p[j],0.5);
        }else if(y_p[j]==0){
            e_y_d[j]=0;
            e_y_u[j]=ROOT::MathMore::chisquared_quantile(1-0.16,2*(h->GetBinContent(j+1))+2)/2 - h->GetBinContent(j+1);
        }
        
        //cout<< x_p[j]*1000<<"   "<<y_p[j]<<"    "<<e_y_u[j]<<endl;

        /*
        std::cout << std::fixed << std::setprecision(0) << j << " " 
        << std::fixed << std::setprecision(3) << x_p[j] - std::floor(x_p[j]) << " " 
        << std::fixed << std::setprecision(0) << y_p[j] << " " 
        << std::fixed << std::setprecision(3) << e_y_d[j] - std::floor(e_y_d[j]) << " " 
        << std::fixed << std::setprecision(3) << e_y_u[j] - std::floor(e_y_u[j]) << std::endl;
        */

        e_y_d[j]/=to_kdru*efficiency;
        e_y_u[j]/=to_kdru*efficiency;
        y_p[j]/=to_kdru*efficiency;
    }
        cout<<"To correct by exposure and get the results in dru, counts and its error must be divided by: "<<to_kdru*efficiency<<endl<<endl;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void spectra_AvsB(){
    /////////////////PATH A LOS ARCHIVOS //////////////////////////////////////////////////
	const char* script_path = gSystem->Which(".", __FILE__);  // Obtiene la ruta del script
    const char* path = gSystem->DirName(script_path); // Extrae el directorio //"/media/dario/0e0c1913-ca43-4252-b2bb-783703f57ea5/proc/Cutted_con_extra_branches_con_varianzas/";//"/mnt/Data/Atucha_Data/";
	
    cout << "script_path " << script_path << endl;
    cout << "path: " << path << endl;

	string folder_A = "/Cutted_con_extra_branches_HOT_COL_new_EDGES_xVarMin_0_xVarMax_4_yVarMin_0yVarMax_4/";//"OFF_2023/";//"ON_2023/";//"OFF_2024/";
	string folder_B = "/Cutted_con_extra_branches_EDGES_xVarMin_0_xVarMax_4_yVarMin_0yVarMax_4/";//"OFF_2024/"; 

	string filelist_A = "OFF_R33_full_shield_2024__HC_Edges.txt";//"OFF_2023.txt";//"ON_2023.txt";//"OFF_2023.txt";
    string filelist_B = "OFF_R33_full_shield_2024_Edges.txt";//"OFF_2024.txt";

    cout << "filelist_A: "<< filelist_A << endl;
    TString title_A = filelist_A.substr(0, filelist_A.find('.')).c_str();
    TString title_B = filelist_B.substr(0, filelist_B.find('.')).c_str();

    TString filename_spect = "spectra" + title_A + title_B + ".png";

    //////////////////////////////////////////////////////////////////////////////////////
    double efficiency_A=0.4910; // Eficiencia despues de limpiar, 2023
    //double efficiency_B=0.4155; // Eficiencia antes de limpiar, 2022
    double efficiency_B=0.4910;

	string treename_hS = "ThitSumm"; //Antes se llamaba hitSumm
	string treename_he = "Theader"; // Antes se llamaba header
	string treename_c = "calc";

	float SEE_min_A {0.01};
	float SEE_max_A {0.90};

    float SEE_min_B {0.01}; 
    float SEE_max_B {0.90}; 
	
	cout << endl;
    cout << "mass [g] = " << mass_in_kg*1000 << endl;
    cout << "bin size [eV] = " << bin_size_in_keV*1000 << endl;

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

	TChain chain_hS_A(treename_hS.c_str());//c_str() converts string to char*;
	TChain chain_he_A(treename_he.c_str());
	TChain chain_c_A(treename_c.c_str());

    cout << endl;
    cout << "---- Data SET A -----" << endl;
	concatenateRootTrees(filelist_A,path+folder_A,chain_hS_A, treename_hS);
	concatenateRootTrees(filelist_A,path+folder_A,chain_he_A, treename_he);
	concatenateRootTrees(filelist_A,path+folder_A,chain_c_A, treename_c);

    generate_energy_histo_with_quality_cuts(chain_hS_A, chain_he_A, chain_c_A, e_bulk_A_hist, 
    	calibration, total_day_time_bseec,total_day_time_aseec, startdate_A, SEE_A, SEE_min_A, SEE_max_A, exposure_factor_A);

    // Expousure calculation for Data Set A ///////////////////////////

    double expousure_time_in_days_A = total_day_time_aseec;
    total_day_time_aseec = 0.0;
    double to_kdru_A = mass_in_kg*bin_size_in_keV*expousure_time_in_days_A*1000;
    
    cout << "Expousure time [days] after  SEE selection = " << expousure_time_in_days_A << endl;
    cout << "to_kdru_A = " << to_kdru_A << endl;
    
    double integralA = e_bulk_A_hist->Integral(5,50)/to_kdru_A/efficiency_A;
    cout<<"Integral Data SET A: "<<integralA<<endl;

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

    TChain chain_hS_B(treename_hS.c_str());//c_str() converts string to char*;
    TChain chain_he_B(treename_he.c_str());
    TChain chain_c_B(treename_c.c_str());

    cout << endl;
    cout << "------ Data SET B -------" << endl;
    concatenateRootTrees(filelist_B,path+folder_B,chain_hS_B, treename_hS);
    concatenateRootTrees(filelist_B,path+folder_B,chain_he_B, treename_he);
    concatenateRootTrees(filelist_B,path+folder_B,chain_c_B, treename_c);
                                                                    
    generate_energy_histo_with_quality_cuts(chain_hS_B, chain_he_B, chain_c_B, e_bulk_B_hist, 
        calibration, total_day_time_bseec,total_day_time_aseec, startdate_B, SEE_B, SEE_min_B, SEE_max_B,exposure_factor_B);

    // Exposure_factor aparece para corregir expousure por limpieza entre lecturas
    //for (int i = 0; i < SEE_B.size(); ++i) SEE_B[i]=exposure_factor*SEE_B[i]; 
    //metí esta cuenta dentro de la función generate_energy_histo_with_quality_cuts

    // Expousure calculation for Data SET B ///////////////////////////
    
    //double expousure_time_in_days_B=total_day_time_aseec/exposure_factor; 
    // el exposure_factor fue agregado debido a la limpieza. metí esta cuenta dentro de la función generate_energy_histo_with_quality_cuts
    double expousure_time_in_days_B=total_day_time_aseec;
    total_day_time_aseec=0.0;
    double to_kdru_B=mass_in_kg*bin_size_in_keV*expousure_time_in_days_B*1000;
    
    cout << "Expousure time [days] after  SEE selection = " << expousure_time_in_days_B << endl;
    cout << "to_kdru_B = " << to_kdru_B << endl;

    double integralB = e_bulk_B_hist->Integral(5,50)/to_kdru_B/efficiency_B;
    cout<<"Integral Data SET B: "<<integralB<<endl<<endl<<endl;

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Histograms content to arrays & asymetric errors calculation with Neyman belt 

    Double_t x_plot[ebines]; //energy
    Double_t e_x[ebines];    //energy error will be set equal to 0

    Double_t e_y_A_d[ebines];
    Double_t e_y_A_u[ebines];
    Double_t e_plot_A[ebines];//OFF AQC ASEEC Soft [counts/rates]
    
    Double_t e_y_B_d[ebines];
    Double_t e_y_B_u[ebines];
    Double_t e_plot_B[ebines];//OFF AQC BSEEC [counts/rates]

    cout<<"For Data SET A: "<<endl;
    fill_histo_with_asymetric_errors(ebines, to_kdru_A, efficiency_A, e_bulk_A_hist, x_plot,e_plot_A,e_x,e_y_A_u,e_y_A_d);
    
    cout<<"For Data SET B: "<<endl;
    fill_histo_with_asymetric_errors(ebines, to_kdru_B, efficiency_B, e_bulk_B_hist, x_plot,e_plot_B,e_x,e_y_B_u,e_y_B_d);
    
    TCanvas* c1 = new TCanvas("canvas", "Spectra", 1200, 600);
    c1->SetGrid();
    //c1->SetLogy();

    TGraphAsymmErrors * grA = new TGraphAsymmErrors(ebines,x_plot,e_plot_A,e_x,e_x,e_y_A_d,e_y_A_u);
    TGraphAsymmErrors * grB = new TGraphAsymmErrors(ebines,x_plot,e_plot_B,e_x,e_x,e_y_B_d,e_y_B_u);
    
    grB->SetMarkerStyle(20);
    grB->SetLineColor(kOrange+1);
    grB->SetMarkerColor(kOrange+1);
    grB->SetMarkerSize(1.3);

    grB->GetXaxis()->SetTitle("Energy (keV)");
    grB->GetYaxis()->SetTitle("Events [kDRU]");
    grB->GetXaxis()->SetTitleSize(0.05);
    grB->GetYaxis()->SetTitleSize(0.05);
    grB->GetXaxis()->SetRangeUser(0.0, 10.5);
    grB->GetYaxis()->SetRangeUser(0.0, 200);
    grB->SetTitle("Atucha II");
    grB->SetTitle("");
    grB->Draw("AP");
    
    grA->SetMarkerStyle(23);
    grA->SetLineColor(kAzure);//kTeal+4
    grA->SetMarkerColor(kAzure);
    grA->SetMarkerSize(1.3);
    grA->Draw("P same");

    auto legend2 = new TLegend(0.75,0.75,0.99,0.95); 
    legend2->AddEntry(grA,title_A,"p"); // .c_str()
    legend2->AddEntry(grB,title_B,"p"); // .c_str()
    legend2->Draw();

    /////////////////////////////////////////////////////////////////////////

    // Crea el pad para el inset
    double inset_x1 = 0.10; // Coordenada x inferior izquierda del inset
    double inset_y1 = 0.40; // Coordenada y inferior izquierda del inset
    double inset_x2 = 0.60; // Coordenada x superior derecha del inset
    double inset_y2 = 0.93; // Coordenada y superior derecha del inset

    TPad *pad_inset = new TPad("pad_inset", "Inset", inset_x1, inset_y1, inset_x2, inset_y2);
    pad_inset->SetGrid();
    //pad_inset->SetLogy();

    pad_inset->SetFillStyle(4000); // Transparente
    pad_inset->Draw();
    pad_inset->cd();

    // Clonar y ajustar los TGraphAsymmErrors en el inset
    TGraphAsymmErrors *cloneGraph1 = dynamic_cast<TGraphAsymmErrors *>(grA->Clone());
    TGraphAsymmErrors *cloneGraph2 = dynamic_cast<TGraphAsymmErrors *>(grB->Clone());

    // Cambia el tamaño del texto en los ejes
    const double label_size =0.05;
    const double title_size =0.05;
    TAxis* xAxis = cloneGraph1->GetXaxis();
    TAxis* yAxis = cloneGraph1->GetYaxis();
    xAxis->SetLabelSize(label_size);
    yAxis->SetLabelSize(label_size);
    xAxis->SetTitleSize(title_size);
    yAxis->SetTitleSize(title_size);

    // Dibujar los gráficos clonados en el inset
    cloneGraph1->Draw("AP");
    cloneGraph2->Draw("P same");
    cloneGraph1->SetTitle("");

    c1->Update();

    //c1->SaveAs("Energy_spectrum_0p26yVar1p10.png");
    c1->SaveAs(filename_spect);//.c_str()

	SEE_evolution_plot(startdate_A, SEE_A, 50, title_A,title_A);
    SEE_evolution_plot(startdate_B, SEE_B, 50, title_B, title_B);

}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// Open files /////////////////////////////////////////////////////////////////////////////////////////
void concatenateRootTrees(const string & filelist, const string & path, TChain & chain, const string & treename)
{
	std::ifstream inputFileList(filelist.c_str());
	
	if(!inputFileList.is_open()) cerr << "File list missing " << filelist << endl ;
	
	string filename;
	while(std::getline(inputFileList,filename)) //leer una linea de inputFileList y almacenarla en filename
	{
		string fullpath = path + filename;
		//cout << fullpath << endl;
		chain.Add(fullpath.c_str());
	}

	inputFileList.close();

	cout << "Entries in " << treename << " tree " << chain.GetEntries() << endl;

}


//////////////////////////////////////////////////////////////////////////////////////////////////// 
void generate_energy_histo_with_quality_cuts(TChain & texp,TChain & header, TChain &calc, TH1D*h_e_bulk, 
	double calibration, double & total_day_time_bc, double & total_day_time_ac, 
    vector<string> & startdate_vec, vector <double> &SEE_v, float see_min, float see_max, double exposure_factor){

	double count=0;
	cout<<"Mass in kg: "<<mass_in_kg<<endl;

	TH2F* h_variance = new TH2F("hist2D", "Variance histogram", 100, 0, 0.1, 100, 0, 1.0);
	TH1F* h_xVar = new TH1F("hist1D", "xVar histogram", 200, -0.05, 0.35);
	TH1F* h_yVar = new TH1F("hist1D", "yVar histogram", 200, 0.25, 1.05);

    int Entries_exp = texp.GetEntries();//int Entries_exp = 1e6;// 1 Entries_exp = 1 cluster
    cout << "Entries in experimental data file: " << Entries_exp << endl;
    
    Enable_and_Set_Branches_hs(texp); //activates branches of interest
    Enable_and_Set_Branches_he(header);
    Enable_and_Set_Branches_calc(calc);

    unsigned int che {0};  //to move along the tree header
    double electrones {0}; //to fill the histogram
    double energia {0};    //to fill the histogram
    int aux_images{};      //amount of images leftover after SEE cut   
    Int_t runID_OLD {0};   //image the event i_event-1 belong to

    for (int i_event = 0; i_event < Entries_exp; ++i_event){
    //if (i_event==1e4) break;
    	
        texp.GetEntry(i_event); //reads activated branches in tree hitSumm in i_event
        calc.GetEntry(i_event); //reads activated branches in tree calc in i_event

        // If this is a new image, count the previous one
        if ((runID != runID_OLD)) {
        	header.GetEntry(che); //reads activated branches in tree header in i_event
        	total_day_time_bc += imagedurationindays(&DATESTART[0],&DATEEND[0]);//expousure before SEE cut
        	//cout << "runID ="<< runID << endl;
        	
        	if(see_min<SEE && SEE<see_max){ 
        		total_day_time_ac += imagedurationindays(&DATESTART[0],&DATEEND[0]); //expousure after SEE cut
        		++aux_images;		
        	}
        	            
            startdate_vec.push_back(DATESTART);
            SEE_v.push_back(SEE*exposure_factor);
            runID_OLD=runID;
            ++che;           
        }

        //if (e >= 16.0/3.75/calibration && e<=1.0e4/3.75/calibration){
        if (yMax<670){
        	if (xVar<=1.6 && yVar<=1.6 && xVar>=0.0 && yVar>=0.0) {//selection of bulk events according to simulation
	            electrones=0;
	            //It calculates total number of electrons in the cluster
	            for (int i = 0; i < nSavedPix; ++i) electrones+=ePix[i];
	             
	            energia=electrones*3.75*calibration*0.001; //keV 
	            
	            //Events selection - Quality cuts
	        	if(see_min<SEE && SEE<see_max)	{

					// Event selection ////////////////////////////////////////////////

	                if(xVar>=xVarMin_GM && xVar<xVarMax_GM && yVar>yVarMin_GM && yVar<yVarMax_GM){
		                if (oproba>=proba_min){
		                	h_e_bulk->Fill(energia);
		                	//h_variance->Fill(-log(xVar*10),-log(yVar));
		                	h_variance->Fill(xVar,yVar);
		                	h_xVar->Fill(sqrt(xVar));
		                	if (xVar==0) count+=1;
		                	h_yVar->Fill(sqrt(yVar));
		                } 
	                	
	                }
	        	}
            }
        }
        //}
    }

    total_day_time_ac = total_day_time_ac/exposure_factor;
    cout << endl;
    cout << "# images before SEE cut = " << che << endl;
    cout << "# images after SEE cut = " << aux_images << endl<< endl;
    che=0;aux_images=0;runID_OLD=0;
/*
    TCanvas* canvas = new TCanvas("canvas", "Histogramas de varianzas", 1200, 600);
    gStyle->SetOptStat(0);
    //canvas->SetLogx();
    //canvas->SetLogy();
    //canvas->SetLogz();
    canvas->Divide(2);
    //canvas->cd(1);
    //gStyle->SetPalette(kViridis);
    //h_variance->Draw("colz");
    canvas->cd(1);
    canvas->cd(1)->SetLogy();
    canvas->cd(1)->SetGridx();
    canvas->cd(1)->SetGridy();
    h_xVar->Draw();
    h_xVar->SetTitle("");
	h_xVar->GetXaxis()->SetTitle("#sigma_{X}");
	h_xVar->GetYaxis()->SetTitle("events");

    canvas->cd(2);
    canvas->cd(2)->SetLogy();
    canvas->cd(2)->SetGridx();
    canvas->cd(2)->SetGridy();
    h_yVar->Draw();
    h_yVar->SetTitle("");
	h_yVar->GetXaxis()->SetTitle("#sigma_{Y}");
	h_yVar->GetYaxis()->SetTitle("events");
    cout<<"Counter: "<<count<<endl<<endl;
    canvas->SaveAs("Sigma_x_Sigma_y.png");
*/
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
double imagedurationindays (Char_t DATESTART[0], Char_t DATEEND[0]){
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

/////////////////////////////////////////////////////////////////////////////////////////////////////
void SEE_evolution_plot(vector<string> vecx, vector<Double_t> vecy,const int paso, const char *canv_name,const char *canv_title){
    TCanvas * c = new TCanvas(canv_name,canv_title, 1000, 750);
    c->cd();
	c->SetGrid(2);
    //fill the TGraph
   	int n = TMath::Min(vecy.size(), vecx.size());

	TGraph* tg = new TGraph(n);
	for (Int_t i = 0; i <n; i++) tg->SetPoint(i, i + 1., vecy[i]); //Set x and y values for point number i.
	//Set x axis label
	auto h = new TH1F("h","h",n,0,n);
   	tg->SetHistogram(h);
	
	int k{0};
	while(k<n){
		h->GetXaxis()->SetBinLabel(k + 1, vecx[k].c_str());
		k +=paso;
	}

	h->GetXaxis()-> LabelsOption("v");
	h->SetStats(0);
	h->Draw("AXIS");
	
	//tg->GetXaxis()->SetTimeDisplay(1); not working 'cause must put labels in seconds. Later on.
	tg->SetMarkerColor(kAzure-1);//(kAzure-1);//kRed);(kOrange+1)
	tg->SetMarkerStyle(33);
	//tg->SetLineColor(kRed);
	tg->SetTitle("Atucha");//not working don't know why
	tg->GetYaxis()->SetTitle("#SEE [1e events/day/pix]");
	//tg->GetYaxis()->SetRangeUser(0.03, 0.2);
	tg->Draw("P same");

	auto legend = new TLegend(0.50,0.75,0.85,0.95); 
    legend->AddEntry(tg,canv_title,"p"); 
    legend->Draw();
	
	gPad->RedrawAxis("g");
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
void Enable_and_Set_Branches_hs(TChain & tree){
    tree.SetBranchStatus("*",1); //enable all branches
    tree.SetBranchAddress ("runID",&runID);
    tree.SetBranchAddress ("ohdu",&ohdu);
    tree.SetBranchAddress ("flag",&flag);
    tree.SetBranchAddress ("xMin",&xMin);
    tree.SetBranchAddress ("xMax",&xMax);
    tree.SetBranchAddress ("yMin",&yMin);
    tree.SetBranchAddress ("yMax",&yMax); 
    tree.SetBranchAddress ("distance",&ddistance);
    tree.SetBranchAddress ("eRaw",&eRaw);
    tree.SetBranchAddress ("xy1Var",&xy1Var); 
    tree.SetBranchAddress ("xy2Var",&xy2Var);
    tree.SetBranchAddress ("xy1Width",&xy1Width); 
    tree.SetBranchAddress ("xy2Width",&xy2Width);
    tree.SetBranchAddress ("bleedX",&bleedX);
    tree.SetBranchAddress ("bleedY",&bleedY);
    //tree.SetBranchAddress ("clusterDist",&clusterDist);
    tree.SetBranchAddress ("e",&e);
    tree.SetBranchAddress ("n",&n);
    tree.SetBranchAddress ("xBary",&xBary);
    tree.SetBranchAddress ("yBary",&yBary);
    tree.SetBranchAddress ("xVar",&xVar);
    tree.SetBranchAddress ("yVar",&yVar);
    tree.SetBranchAddress ("nSavedPix",&nSavedPix);
    tree.SetBranchAddress ("xPix",&xPix);
    tree.SetBranchAddress ("yPix",&yPix);
    tree.SetBranchAddress ("ePix",&ePix);
    tree.SetBranchAddress ("expoStart",&expoStart);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
void Enable_and_Set_Branches_he(TChain & tree){
    tree.SetBranchStatus("*",1); //enable all branches
    /*tree.SetBranchStatus("DATESTART",1); 
    tree.SetBranchStatus("DATEEND",1); 
    tree.SetBranchStatus("RUNID",1);*/
    tree.SetBranchAddress ("DATESTART",&DATESTART);
    tree.SetBranchAddress ("DATEEND",&DATEEND);
    tree.SetBranchAddress ("RUNID",&RUNID_head);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
void Enable_and_Set_Branches_calc(TChain & tree){
    tree.SetBranchStatus("*",1); //enable all branches
    tree.SetBranchAddress ("oproba",&oproba);
    tree.SetBranchAddress ("SEE",&SEE);
    tree.SetBranchAddress ("DEE",&DEE);
    tree.SetBranchAddress ("OCCUPANCY",&OCCUPANCY);
    tree.SetBranchAddress ("RN_1",&RN_1);
    tree.SetBranchAddress ("RN_2",&RN_2);
    tree.SetBranchAddress ("SER_1",&SER_1);
    tree.SetBranchAddress ("SER_2",&SER_2);
    tree.SetBranchAddress ("HOTCOL_amount",&HOTCOL_amount);
    tree.SetBranchAddress ("time_expo",&time_expo);
}
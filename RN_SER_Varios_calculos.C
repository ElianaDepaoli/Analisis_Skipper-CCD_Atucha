/* Autora: E.Depaoli Fecha inicial: 24/04/2025
Código para calcular el ruido de lectura y el SEE de las imágenes de Atucha 
utilizando los catálogos hits_corr_proc_run*.root. Calcula el p-valor del ajuste.
-> Histograma de carga sin clusterizar (ePix en calPixTree) restringida a los dos overscans, X e Y, para cada extensión. 
No se contabiliza carga en píxeles de columnas brillantes.
Ajusta cada histograma con la convolución entre Poisson y Gauss. Parámetros que se obtienen de este ajuste: el ruido de lectura, SER*tiempo.
-> Calcula cluster_1e como la suma del número de eventos que ocupan 1 pxl y tienen 1 electrón (n == 1.0 && e == 1.0) 
en todas las extensiones. Calcula cluster_almenos_1e +=n 
SEE = cluster_1e/N_tot_pxls/t_exposición
Occupancy = cluster_almenos_1e/N_tot_pxls/t_exposición_imagen
-> Calcula umbrales a partir de una relación que se obtiene al igualar las densidades de probabilidad de los picos de
0 y de 1 electrón. Esta cuenta surge de igualar los errores de tipo 1 y luego derivar aplicando el teorema fundamental
del cálculo. No está bien este resultado, como lo demuestra el hecho de que no se cumple que los errores de tipo 1 
efectivamente den iguales. La diferencia entre ellos es del 30 al 60 %. HAY QUE IMPLEMENTAR LA TABLA (xc, lambda, RN)
Observaciones:
-> El pico de 0 e- del OHDU 1 tiene menor altura que el del OHDU 2 porque su sigma 
es mayor

*/

#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
#include <fstream>
#include <sys/stat.h> // Para verificar si el archivo existe
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
using std::cout; using std::vector; using std::string; using std :: copy;

// Function declaration -------------------------------------------------------------
void Enable_and_Set_Branches_A(TTree* & tree);
void Enable_and_Set_Branches_B(TTree* & tree);
void Enable_and_Set_Branches_C(TTree* & tree);
double imagedurationindays (Char_t DATESTART[0], Char_t DATEEND[0], float factor);
void extracthc(int vmin[], int vmax[],int arraysize, bool &isnthotcol, Int_t xm , Int_t xM);
bool fileExists(const std::string& filename);
int hot_col_amount(int vmin[], int vmax[], int arraysize);
// ----------------------------------------------------------------------------------

//Configuracion global para los gráficos --------------------------
void SetGlobalStyle(){
	TStyle *st1 = new TStyle("st1","my style");
	st1->SetTitleSize(0.05,"XYZ");
	st1->SetLabelSize(0.05, "XYZ");
	st1->SetPadBorderMode(0);
	st1->SetLegendTextSize(0.035);
	st1->SetLegendBorderSize(-1);
	st1->SetGridColor(0);
	st1->SetLegendFillColor(0);
	gStyle->SetLegendFillColor(0);
	gROOT->SetStyle("st1");//
	gROOT->ForceStyle();//esto no funciona: el estilo se aplica a todos los objetos creados después
	st1->cd ();//This is now the current style gStyle
}

void SetOhdu1Style(TH1 *h1){
	h1->GetYaxis()->SetTitleOffset(0.8);
	h1->SetMarkerColor(kBlue+1);
	h1->SetMarkerStyle(20);
	h1->SetMarkerSize(1);
	h1->SetMarkerStyle(20);
	h1->SetMarkerSize(1);
	h1->SetLineColor(kBlue+1);
	h1->SetLineWidth(2);
	h1->GetXaxis()->SetTitle("ePix [e-]");
	h1->GetYaxis()->SetTitle("density");
}

void SetOhdu2Style(TH1 *h1){
	h1->SetMarkerColor(kCyan+1);
	h1->SetMarkerStyle(20);
	h1->SetMarkerSize(1);
	h1->SetMarkerStyle(20);
	h1->SetMarkerSize(1);
	h1->SetLineColor(kCyan+1);
	h1->SetLineWidth(2);
	h1->GetXaxis()->SetTitle("ePix [e-]");
	h1->GetYaxis()->SetTitle("density");
}
	
// Global Variables ---------------------------------------------------------------------
//string filename{"/media/dario/0e0c1913-ca43-4252-b2bb-783703f57ea5/proc/run_31_19Dic2023/hits_corr_proc_run_31_19Dic2023__EXP1_NSAMP300_VSUB70_img1005.root"};
//string filename{"/media/dario/0e0c1913-ca43-4252-b2bb-783703f57ea5/proc/run_40/hits/hits_corr_proc_run_40_01Feb2025__EXP1_NSAMP300_VSUB70_img99.root"};
//string filename{"/home/eliana/Documentos/Atucha/images/run_40/hits/hits_corr_proc_run_40_01Feb2025__EXP1_NSAMP300_VSUB70_img99.root"};
string filename{"/home/eliana/Documentos/Atucha/images/hits_corr_proc_run_43_31Mar2025__EXP1_NSAMP300_VSUB70_img17.root"};

double xc_oh1{0.655};
double xc_oh2{0.610};
double xc_ohdu1_ovx{0.0};
double xc_ohdu2_ovx{0.0};

double Nmax{};
//Estimadores ·············································································
//Active area ····
double RN_1{}; double err_RN_1{};
double RN_2{}; double err_RN_2{};
double lamb_1{}; double err_lamb_1{};
double lamb_2{}; double err_lamb_2{};
double N1_N0_actar_oh1{}; double N1_N0_actar_oh2{}; double N1_N0_actar_oh3{}; double N1_N0_actar_oh4{};
//Overscan ·····
double RN_1_ovx{}; double err_RN_1_ovx{};
double RN_2_ovx{}; double err_RN_2_ovx{};
double lamb_1_ov{}; double err_lamb_1_ov{};
double lamb_2_ov{}; double err_lamb_2_ov{};
double N1_N0_ovx_oh1{}; double N1_N0_ovx_oh2{}; double N1_N0_ovx_oh3{}; double N1_N0_ovx_oh4{};
double N1_N0_ovy_oh1{}; double N1_N0_ovy_oh2{}; double N1_N0_ovy_oh3{}; double N1_N0_ovy_oh4{};

// Timming ····
double time_expo{0};
//Tree Variables ···
// Trees
string treename_cp = "calPixTree";
string treename_he = "headerTree_0";
string treename_hs = "hitSumm";

// calPixTree --
int x; int y; int RUNID; int ohdu;
double ePix;//tiene los eventos de 0 e-

//header ---
Char_t DATESTART [256]; //[256];
Char_t DATEEND [256];
Char_t RUNID_head[256];

//hitSumm ---	
int NSAMP = 300; int SR = 307;
Int_t runID; float e; float n; 
Int_t ohdu_hs; Float_t xBary; Float_t yBary;
Int_t xMin; Int_t xMax; Int_t yMin; Int_t yMax;
const Int_t kMaxTrack = 2e7;
int xPix[kMaxTrack];
int yPix[kMaxTrack];

//n: occupied pixels, e: electrons in cluster, discretized
//ePix: no contiene los eventos de 0 e-
//Quality cuts ///////////////////////////////////////////////////
// Hot Columns definition ///////////
// median +- sigma "hits_corr_proc_Reactor_ON_xVarMax_1_yVarMax_2.root"
// 1 electron events

int hotcol_ohdu_1_min[] = {3, 99, 174, 201, 206, 231, 276, 286, 307};
int hotcol_ohdu_1_max[] = {6, 99, 177, 202, 206, 232, 277, 288, 307};
int hotcol_ohdu_2_min[] = {6, 53,  58, 199, 307};
int hotcol_ohdu_2_max[] = {9, 54,  60, 199, 307};

bool isnthotcol;//to remove hot columns
int sizehotcol_1_min{sizeof(hotcol_ohdu_1_min)/sizeof(hotcol_ohdu_1_min[0])};
int sizehotcol_2_min{sizeof(hotcol_ohdu_2_min)/sizeof(hotcol_ohdu_2_min[0])};

// remove events partially inside the active area
//Borders of the active area
int xBaryMin=3; // 
int xBaryMax=305; // 5 on the left from the overscan
int yBaryMin=3; 
int yBaryMax=545;

// usefull quadrants
int ohdu_1=1;
int ohdu_2=2;

float factor = 0.5;//corrección por lectura entre dos imágenes sucesivas

//Archivo de salida
string ofilename = "ruido_de_lectura_Atucha.txt";
string onoise_DCfilename = "ruido_de_lectura_SEE_Atucha.txt";
string oumbralextract_filename = "th_for_skextract";

//Función para calcular SER ·························

Double_t poisson_gauss_conv(Double_t *x, Double_t *par){//lambda [e-/pxl] , mu, normalización, sigma

		double xx = x[0];
		double y = 0.0;

		for(int i{0}; i < Nmax+1; ++i){
			double pois = pow(par[0],i)/TMath::Factorial(i);
			double gauss = TMath::Exp(-0.5*pow((xx-i-par[1])/par[3],2)); 
			y += pois*gauss;
		}

		return TMath::Exp(-par[0])*y/(par[3] * TMath::Sqrt(2 * TMath::Pi()))*par[2];//mejora el p-valor si está par[2] (factor de normalización)
}

// Rutina  ----------------------------------------------------------------------------------
void RN_SER_Varios_calculos(){
	SetGlobalStyle();

	//int hc_amount=hot_col_amount(hotcol_ohdu_1_min, hotcol_ohdu_1_max, sizehotcol_1_min) + hot_col_amount(hotcol_ohdu_2_min, hotcol_ohdu_2_max, sizehotcol_2_min);// = 21 + 106;
	//cout << "hc_amount = " << hc_amount << endl;
	double amount_useful_pix_oh1 = (xBaryMax - xBaryMin - hot_col_amount(hotcol_ohdu_1_min, hotcol_ohdu_1_max, sizehotcol_1_min))*(yBaryMax - yBaryMin);//amount of pixels ohdu1 + ohdu2
	double amount_useful_pix_oh2 = (xBaryMax - xBaryMin - hot_col_amount(hotcol_ohdu_2_min, hotcol_ohdu_2_max, sizehotcol_2_min))*(yBaryMax - yBaryMin);
	double amount_useful_pix = (xBaryMax - xBaryMin)*(yBaryMax - yBaryMin);
	// Genero histogramas 	········
	int nbins_epix_overscan = 35;
	float lower_bin = -1.0;//0.66;//
	float upper_bin = 1.8;//0.66;//
	TH1D * epix_ohdu1_ovx_hist = new TH1D("epix ohdu1", "Overscan X", nbins_epix_overscan, lower_bin, upper_bin);// ohdu1
	TH1D * epix_ohdu2_ovx_hist = new TH1D("epix ohdu2", "Overscan X", nbins_epix_overscan, lower_bin, upper_bin);
	TH1D * epix_ohdu1_ovy_hist = new TH1D("epix ohdu1", "Overscan Y", 50, -0.5, 1.5);
	TH1D * epix_ohdu2_ovy_hist = new TH1D("epix ohdu2", "Overscan Y", 150, -0.5, 1.5);

	//Open root catalogue
	TFile * file = TFile::Open(filename.c_str(),"READ");
	//c_str() converts a string to an array of characters & terminates it with a null character. No parameters are allowed, a pointer to the array is returned.
	if (!file->IsOpen()) {std::cerr << "Failed to load file" << filename << std::endl;}
	//Retrieve the tree from the file
	TTree * cptree = (TTree*) file->Get(treename_cp.c_str());
	TTree * htree = (TTree*) file->Get(treename_he.c_str());
	TTree * hStree = (TTree*) file->Get(treename_hs.c_str());

	int Entries_cptree =  cptree -> GetEntries();//cantidad de eventos en el árbol
	int Entries_htree =  htree -> GetEntries();
	int Entries_hStree =  hStree -> GetEntries();

	Enable_and_Set_Branches_A(cptree);
	Enable_and_Set_Branches_B(htree);
	Enable_and_Set_Branches_C(hStree);

	htree->GetEntry(0);
	cout << "Starting date" << &DATESTART[0] << endl;
	cout << "Ending date" << &DATEEND[0] << endl;
		
	//SEE. Viejo. ····
    double clusters_1e_oh1 {0};double clusters_1e_oh2 {0};double clusters_1e_oh3 {0};double clusters_1e_oh4 {0};//SEE
    double clusters_2e_oh1 {0};double clusters_2e_oh2 {0};double clusters_2e_oh3 {0};double clusters_2e_oh4 {0};
    Double_t see_oh1{0};Double_t see_oh2{0}; Double_t see_oh3{0}; Double_t see_oh4{0};
    Double_t occ_ohd1{0}; Double_t occ_ohd2{0}; Double_t occ_ohd3{0}; Double_t occ_ohd4{0};
	// hitSumm Tree ----------------------------------------------------------------
    cout << "Cantidad de entradas en hitSumm = " << Entries_hStree << endl;

    for (Long64_t entry{0}; entry < Entries_hStree ; ++entry)
    {
        hStree -> GetEntry(entry);
    	isnthotcol=true;
      
    	if (xBary<xBaryMax && xBary>xBaryMin && yBary<yBaryMax && yBary>yBaryMin){//events inside edges of CCD
    		if (ohdu_hs == ohdu_1 || ohdu_hs == ohdu_2){//events in quadrants with single electron counting
    			//Extraer columnas brillantes no reduce el valor final de clusters_1e porque estos no contienen eventos con n > 1. Pero debería reducir la ocupancia ¡y lo hace!
    			if (ohdu_hs==ohdu_1){//events outside hot columns
        		    extracthc(hotcol_ohdu_1_min, hotcol_ohdu_1_max,sizehotcol_1_min, isnthotcol, xMin, xMax);
		        } else if (ohdu_hs==ohdu_2){extracthc(hotcol_ohdu_2_min, hotcol_ohdu_2_max,sizehotcol_2_min, isnthotcol, xMin, xMax);}

		 	}       
    			if (isnthotcol==true){//si el evento no contiene ni toca una columna brillante
			        if(yMax < 512 && xMax < 307.9 && xMin >= 1){
			        	if (n == 1.0 && e == 1.0 ){//n: pixels occupied, e: electrons in cluster
				        	switch(ohdu_hs){
								case 1:
									++clusters_1e_oh1;
								case 2:
									++clusters_1e_oh2;
								case 3:
									++clusters_1e_oh3;

								case 4:
									++clusters_1e_oh4;
								}
						}

						switch(ohdu_hs){
							case 1:
								occ_ohd1+=n;
								if (n == 1.0 && e == 2.0 ) ++clusters_2e_oh1;
							case 2:
								occ_ohd2+=n;
								if (n == 1.0 && e == 2.0 ) ++clusters_2e_oh2;
							case 3:
								occ_ohd3+=n;
								if (n == 1.0 && e == 2.0 ) ++clusters_2e_oh3;
							case 4:
								occ_ohd4+=n;
								if (n == 1.0 && e == 2.0 ) ++clusters_2e_oh4;
						}

					}
			    }
	    }

    //··················
	//Cálculos 		····
	//··················

    //Luego de pasar por el último evento, ocupancia y DC 
        if(entry > (Entries_hStree-2)){
        	//cout << "final entry = " << entry << endl;
            time_expo = imagedurationindays(&DATESTART[0],&DATEEND[0], factor);
            //cout << " clusters_1e_oh2 " << clusters_1e_oh2/imagedurationindays(&DATESTART[0],&DATEEND[0],factor)/amount_useful_pix << endl;
        }
    }

    cout << " ················· Active area  - Clusters ··········· " << endl;
    cout << "imagedurationindays = " << setprecision(4)<< time_expo << '\n';
    cout << "amount of useful pixels in ohdu 1 = " << amount_useful_pix_oh1 << endl;
    cout << "amount of useful pixels in ohdu 2 = " << amount_useful_pix_oh2 << endl;
    cout << "amount of useful pixels in ohdu 3 = " << amount_useful_pix << endl;
            
    double readout_time_per_pxl_per_sample = 2*time_expo*86400/(1024*(6144/20+12.8)*NSAMP/2);//
    double t_media_ov_px = NSAMP*SR*readout_time_per_pxl_per_sample/(24*60*60);//67e-6/86400
    //double t_expo_overscan_y = NSAMP*readout_time_per_pxl_per_sample*(SR+1)/2 + time_expo*24*60*60;
    double t_expo_overscan_y = NSAMP*readout_time_per_pxl_per_sample*(SR+1)/(2*24*60*60) + time_expo;
    cout << "readout_time_per_pxl_per_sample = " << 1e6*readout_time_per_pxl_per_sample << " microsec" << endl;
    cout << "tiempo de exposición promedio de un pixel x-overscan  = " << t_media_ov_px << "day " << endl; //" sec"  
    cout << "tiempo de exposición promedio de un pixel y-overscan  = " << t_expo_overscan_y <<  "day " << endl; //" sec"

    occ_ohd1 = static_cast<float>(occ_ohd1)/ static_cast<float>(amount_useful_pix_oh1) / time_expo; // / (time_expo*24*60*60); ///;//imagedurationindays(&DATESTART[0],&DATEEND[0],factor)/
    occ_ohd2 = static_cast<float>(occ_ohd2)/ static_cast<float>(amount_useful_pix_oh2)  / time_expo; // / (time_expo*24*60*60);//occ_ohd2;///amount_useful_pix_oh2;
    occ_ohd3 = static_cast<float>(occ_ohd3)/ static_cast<float>(amount_useful_pix)  / time_expo; // / (time_expo*24*60*60);//occ_ohd3;///amount_useful_pix;
    occ_ohd4 = static_cast<float>(occ_ohd4)/ static_cast<float>(amount_useful_pix)  / time_expo; // / (time_expo*24*60*60); //occ_ohd4;///amount_useful_pix;
    see_oh1 = clusters_1e_oh1/amount_useful_pix_oh1  /time_expo; // / (time_expo*24*60*60);
    see_oh2 = clusters_1e_oh2/amount_useful_pix_oh2  /time_expo; // / (time_expo*24*60*60);
    see_oh3 = clusters_1e_oh3/amount_useful_pix  /time_expo; // / (time_expo*24*60*60);
    see_oh4 = clusters_1e_oh4/amount_useful_pix  /time_expo; // / (time_expo*24*60*60);


    cout << "see_oh1 = " << see_oh1 << " ---- " << "occ_ohd1 = " << occ_ohd1 <<  " e/pix/day " << endl;
    cout << "see_oh2 = " << see_oh2 << " ---- " << "occ_ohd2 = " << occ_ohd2 <<  " e/pix/day " << endl;
    cout << "see_oh3 = " << see_oh3 << " ---- " << "occ_ohd3 = " << occ_ohd3 <<  " e/pix/day " << endl;
    cout << "see_oh4 = " << see_oh4 << " ---- " << "occ_ohd4 = " << occ_ohd4 <<  " e/pix/day " << endl;
    cout << "lamb_oh1 = 2*N_2/N_1 = " << 2*clusters_2e_oh1/clusters_1e_oh1  <<  " e/pix" << endl;
    cout << "lamb_oh2 = 2*N_2/N_1 = " << 2*clusters_2e_oh2/clusters_1e_oh2  <<  " e/pix"  << endl;
    cout << "lamb_oh3 = 2*N_2/N_1 = " << 2*clusters_2e_oh3/clusters_1e_oh3  <<  " e/pix" << endl;
    cout << "lamb_oh4 = 2*N_2/N_1 = " << 2*clusters_2e_oh4/clusters_1e_oh4  <<  " e/pix" << endl;
    
    int count_0e_oh1{0}; int count_1e_oh1{0}; int count_masde2e_oh1{0};
    int count_0e_oh1_ovy{0}; int count_1e_oh1_ovy{0};
    int count_0e_oh2{0}; int count_1e_oh2{0};
    int count_0e_oh2_ovy{0}; int count_1e_oh2_ovy{0};

	// CalPixTree ---------------------------------------------------------------
	for(int j{0}; j < Entries_cptree; ++j){
		cptree->GetEntry(j);
		//Overscan x ·····
		if(x > 307.9 && ePix < 2 && ePix > -1){//
			switch(ohdu){
			case 1:
				epix_ohdu1_ovx_hist->Fill(ePix);
				if(ePix < 0.7) count_0e_oh1+=1;
				if((ePix >= 0.7) && (ePix < 1.5)) count_1e_oh1+=1;
				if(ePix >= 1.5) count_masde2e_oh1+=1;
			case 2:
				epix_ohdu2_ovx_hist->Fill(ePix);
				if(ePix < 0.7) count_0e_oh2+=1;
				if((ePix >= 0.7) && (ePix < 1.5)) count_1e_oh2+=1;
				
			}
		}

		//SER ····
		isnthotcol = true;
		//Overscan y ····
		if(x <= 307.9  && x > 0.8 && y > 512 && ePix < 2 && ePix > -1){
			//Transforma isnthotcol = false si el pixel está en una columna dañada
			switch(ohdu){
			case 1:
				extracthc(hotcol_ohdu_1_min, hotcol_ohdu_1_max,sizehotcol_1_min, isnthotcol, x, x);
			case 2: 
				extracthc(hotcol_ohdu_2_min, hotcol_ohdu_2_max,sizehotcol_2_min, isnthotcol, x, x);
			}
	        
        	if (isnthotcol==true){//si el pixel no está en una columna brillante
				
				switch(ohdu){
				case 1: 
					epix_ohdu1_ovy_hist->Fill(ePix);
					if(ePix < 0.7) count_0e_oh1_ovy+=1;
					if((ePix >= 0.7) && (ePix < 1.5)) count_1e_oh1_ovy+=1;
				case 2:
					epix_ohdu2_ovy_hist->Fill(ePix);
					if(ePix < 0.7) count_0e_oh2_ovy+=1;
					if((ePix >= 0.7) && (ePix < 1.5)) count_1e_oh2_ovy+=1;
				
	    		}

    	    }
		}

	}

	cout << '\n' << endl;
	cout << "·············calPixTree ············" << endl;
	cout << "Overscan x ···············" << endl;
    cout << "N0  = " << count_0e_oh1 << endl;
    cout << "N1  = " << count_1e_oh1 << endl;
    cout << "N>1 = " << count_masde2e_oh1 << endl;
    cout << "Ohdu 1 N1/N0 = " << static_cast<float>(count_1e_oh1)/ static_cast<float>(count_0e_oh1) << endl;
    cout << "Ohdu 2 N1/N0 = " << static_cast<float>(count_1e_oh2)/ static_cast<float>(count_0e_oh2) << endl;
	cout << "Overscan y ···············" << endl;
	cout << "Ohdu 1 N1/N0 = " << static_cast<float>(count_1e_oh1_ovy)/ static_cast<float>(count_0e_oh1_ovy) << endl;
    cout << "Ohdu 2 N1/N0 = " << static_cast<float>(count_1e_oh2_ovy)/ static_cast<float>(count_0e_oh2_ovy) << endl;
	
	//Normalizo histogramas ····
	epix_ohdu1_ovy_hist->Scale( 1./epix_ohdu1_ovy_hist->Integral(),"WIDTH");//histogram describe a probability density function. 
	epix_ohdu1_ovx_hist->Scale(1./epix_ohdu1_ovx_hist->Integral(), "WIDTH");
	epix_ohdu2_ovy_hist->Scale( 1./epix_ohdu2_ovy_hist->Integral(),"WIDTH");
	epix_ohdu2_ovx_hist->Scale(1./epix_ohdu2_ovx_hist->Integral(), "WIDTH");
	//Clones de ohdu histograms ···
	//TH1*ohdu1_overscan_hc = (TH1*)epix_ohdu1_ovx_hist->Clone();
	//TH1*ohdu2_overscan_hc = (TH1*)epix_ohdu2_ovx_hist->Clone();
	
	//··················
	//Ajustes 		····
	//··················
	//Solo Ruido de lectura: Suma de dos gaussianas con igual varianza ·······················
	TF1*total_gauss = new TF1("G0+G1","[0]*exp(-0.5*pow((x-[1])/[2],2)) + [3]*exp(-0.5*pow((x-[4])/[2],2))", -1.0,1.5);//
	total_gauss->SetLineColor(kMagenta-4);
	//Inicializo parámetros
	total_gauss->SetParameters(1.5,0.0003, 0.22, 3.0, 1.0);
	total_gauss->SetParLimits(2,0.17, 0.25);
	
	//“Q” Quiet mode (minimum printing); “S” The result of the fit is returned in the TFitResultPtr; "R" to restrict the fit to the range specified in the TF1 constructor
	// SER + Ruido de lectura: Poisson convolucionada con gaussiana ·······················
	TF1 *pgc = new TF1("poisson_gauss_conv", poisson_gauss_conv, -1, 2,4);//(,, x_lim_inf, x_lim_sup, #parameters)
	Nmax= 1;
	//pgc->SetParameters(0.2,0.03, 1.0); //lambda, mu, normalización. Inicializo.
	pgc->SetParameters(0.2,0.03, 1.0, 0.21); //lambda, mu, normalización, sigma. Inicializo.
	pgc->SetParNames("lambda", "mu", "norm", "RN");
	pgc->FixParameter(1,0);
	pgc->SetLineColor(kMagenta-4);//kOrange-3);
	
	//OVERSCAN ··············
	cout << '\n' << endl;
	cout << " --------- OVERSCAN X --------- " << endl;
	// --------- Convolution Poisson-Gauss---------;

	cout << " OHDU 1 " << endl;
	auto resultovx = epix_ohdu1_ovx_hist->Fit("poisson_gauss_conv", "R+S");//+Q
	double pvalueovx = ROOT::Math::chisquared_cdf_c(resultovx->Chi2(), resultovx->Ndf()); //	cout << "pvalue OVERSCAN  = " << pvalueovx << endl;
	lamb_1_ov = resultovx -> Parameter(0); err_lamb_1_ov = resultovx -> ParError(0); RN_1_ovx = resultovx -> Parameter(3); err_RN_1_ovx = resultovx -> ParError(3); xc_ohdu1_ovx = 0.5 - pow(RN_1_ovx,2)*TMath::Log(lamb_1_ov);
	xc_ohdu1_ovx = 0.5 - pow(RN_1_ovx,2)*TMath::Log(lamb_1_ov);	//Ver cálculo en cuaderno
	//Error tipo 1 resultante ····· 
	double error_tipo1_N0 = ROOT::Math::normal_cdf_c(xc_ohdu1_ovx,RN_1_ovx);//, 1);// (double x, double sigma=1, double x0=0)
	double error_tipo1_N1 = lamb_1_ov*(ROOT::Math::normal_cdf(xc_ohdu1_ovx,RN_1_ovx,1));//

	cout << "pvalue = " << pvalueovx << endl;
	cout << "SER ohdu 1 = " << std::setprecision(9) << lamb_1_ov / t_media_ov_px << " +- " << err_lamb_1_ov / t_media_ov_px << endl;
	cout << "xc ohdu 1 = " << xc_ohdu1_ovx << endl;
	cout << "error_tipo1_N0 = " << error_tipo1_N0 << endl;
	cout << "error_tipo1_N1 = " << error_tipo1_N1 << endl;
	
	cout << '\n' << endl;
	cout << " OHDU 2 " << endl;
	resultovx = epix_ohdu2_ovx_hist->Fit("poisson_gauss_conv", "R+S");//+Q
	pvalueovx = ROOT::Math::chisquared_cdf_c(resultovx->Chi2(), resultovx->Ndf());
	lamb_2_ov = resultovx -> Parameter(0); err_lamb_2_ov = resultovx -> ParError(0); RN_2_ovx = resultovx -> Parameter(3);	err_RN_2_ovx = resultovx -> ParError(3); xc_ohdu2_ovx = 0.5 - pow(RN_2_ovx,2)*TMath::Log(lamb_2_ov);
	xc_ohdu2_ovx = 0.5 - pow(RN_2_ovx,2)*TMath::Log(lamb_2_ov);	//Ver cálculo en cuaderno
	//Error tipo 1 resultante ····· 
	error_tipo1_N0 = ROOT::Math::normal_cdf_c(xc_ohdu2_ovx,RN_2_ovx);//, 1);// (double x, double sigma=1, double x0=0)
	error_tipo1_N1 = lamb_2_ov*(ROOT::Math::normal_cdf(xc_ohdu2_ovx,RN_2_ovx,1));//
	
	cout << "pvalue = " << pvalueovx << endl;
	cout << "SER ohdu 2 = " << std::setprecision(9) << lamb_2_ov/ t_media_ov_px  << " +- " << err_lamb_2_ov/ t_media_ov_px  << endl;
	cout << "xc ohdu 2 = " << xc_ohdu2_ovx << endl;
	cout << "error_tipo1_N0 = " << error_tipo1_N0 << endl;
	cout << "error_tipo1_N1 = " << error_tipo1_N1 << endl;
	
	cout << '\n' << endl;

	//Grafico OVERSCAN X ···
	//Poisson convolucionada con Gauss ·····························
	//gStyle->SetOptFit(0111);
    TCanvas * c1 = new TCanvas("c1", "Charge overscan histogram", 1200, 800);
	gStyle->SetOptStat(11);
	gStyle->SetOptFit(0111);
	c1 -> SetLogy();
	c1-> Divide(2,1);
	//ohdu 1 ····
	c1->cd(1);
	gPad->SetLogy();
    gPad->SetGrid();
    //gPad->SetLeftMargin(0.0);
    //epix_ohdu1_ovx_hist->GetYaxis()->SetTitleOffset(0.8);
	SetOhdu1Style(epix_ohdu1_ovx_hist);
	epix_ohdu1_ovx_hist->Draw();
	//ohdu 2 ····
	c1->cd(2);
	gPad->SetLogy();
    gPad->SetGrid();
	SetOhdu1Style(epix_ohdu2_ovx_hist);
	epix_ohdu2_ovx_hist->Draw();

	TObjArray *tokens = TString(filename).Tokenize("_, .");
	TString OVX_out = ((TObjString*)tokens->At(3))->GetString() + "_" + ((TObjString*)tokens->At(4))->GetString()+ "_" + ((TObjString*)tokens->At(9))->GetString() + "_OVX" + ".png";
	TString OVY_out = ((TObjString*)tokens->At(3))->GetString() + "_" + ((TObjString*)tokens->At(4))->GetString()+ "_" + ((TObjString*)tokens->At(9))->GetString() + "_OVY" + ".png";
	//std::cout << OVX_out << std::endl;
	delete tokens;
	c1->SaveAs(OVX_out.Data());

	cout << " --------- OVERSCAN Y --------- " << endl;
	cout << " OHDU 1 " << endl;
	
	//Grafico y ajusto OVERSCAN Y ···
	//Poisson convolucionada con Gauss ·····························
	//gStyle->SetOptFit(0111);
    TCanvas * c2 = new TCanvas("c2", "Charge overscan histogram", 1200, 800);
	gStyle->SetOptStat(11);
	gStyle->SetOptFit(0111);
	c2-> Divide(2,1);
	c2->cd(1);
	//ohdu 1 ····
	auto resultovy = epix_ohdu1_ovy_hist->Fit("poisson_gauss_conv", "R+S");
	double pvalueovy = ROOT::Math::chisquared_cdf_c(resultovy->Chi2(), resultovy->Ndf());
	lamb_1 = resultovy -> Parameter(0); err_lamb_1 = resultovy -> ParError(0); RN_1 = resultovy -> Parameter(3); err_RN_1 = resultovy -> ParError(3);
	xc_oh1 = 0.5 - pow(RN_1,2)*TMath::Log(lamb_1);
	//Error tipo 1 ····	
	error_tipo1_N0 = ROOT::Math::normal_cdf_c(xc_oh1, RN_1,0);
	error_tipo1_N1 = lamb_1*(ROOT::Math::normal_cdf(xc_oh1, RN_1,1));

	cout << "pvalue = " << pvalueovy << endl;
	cout << "SER ohdu 1 = " << std::setprecision(9) << lamb_1 / t_expo_overscan_y << " +- " << err_lamb_1 / t_expo_overscan_y << endl;
	cout << "xc ohdu 1 = " << xc_oh1 << endl;	
	cout << "error_tipo1_N0 = " << error_tipo1_N0 << endl;
	cout << "error_tipo1_N1 = " << error_tipo1_N1 << endl;
	gPad->SetLogy();
    gPad->SetGrid();
	SetOhdu1Style(epix_ohdu1_ovy_hist);
	epix_ohdu1_ovy_hist->Draw();
		
	//ohdu 2 ····
	c2->cd(2);
	gPad->SetLogy();
    gPad->SetGrid();
	//N(0,noise) y N(1,noise) ··········
	/*TF1 * N0noise = new TF1("normal_0_noise", "gaus(0)", -0.5,1.5);//
	TF1 * N1noise = new TF1("normal_1_noise", "gaus(0)", -0.5,1.5);//gaus(0) = [0]*exp(-0.5*((x-[1])/[2])**2) and (0) means start numbering parameters at 0
	TF1 * N0noise_oh2 = new TF1("normal_0_noise", "gaus(0)", -0.5,1.5);
	TF1 * N1noise_oh2 = new TF1("normal_1_noise", "gaus(0)", -0.5,1.5);
	N0noise->SetLineColor(kGray+2);
	N1noise->SetLineColor(kOrange+2);
	N0noise_oh2->SetLineColor(kGray+2);
	N1noise_oh2->SetLineColor(kOrange+2);
	N0noise->SetParameters(pow(resultovy -> Parameter(3),-1)*pow(2*TMath::Pi(),-0.5), resultovy -> Parameter(1), resultovy -> Parameter(3));//resultovy -> Parameter(2)
	N1noise->SetParameters(resultovy->Parameter(0)*pow(resultovy -> Parameter(3),-1)*pow(2*TMath::Pi(),-0.5), 1+resultovy -> Parameter(1), resultovy -> Parameter(3));
	//double lamb_recal = ROOT::Math::normal_cdf(-xc/RN_1)/ROOT::Math::normal_cdf((xc-1)/RN_1);
	//cout << "lamb_recal = " << lamb_recal << endl;
	cout << '\n' << endl;
	*/
	cout << " OHDU 2 " << endl;
	resultovy = epix_ohdu2_ovy_hist->Fit("poisson_gauss_conv", "R+S");
	pvalueovy = ROOT::Math::chisquared_cdf_c(resultovy->Chi2(), resultovy->Ndf());
	lamb_2 = resultovy -> Parameter(0); err_lamb_2 = resultovy -> ParError(0); RN_2 = resultovy -> Parameter(3); err_RN_2 = resultovy -> ParError(3);
	//N0noise_oh2->SetParameters(pow(resultovy -> Parameter(3),-1)*pow(2*TMath::Pi(),-0.5), resultovy -> Parameter(1), resultovy -> Parameter(3));
	//N1noise_oh2->SetParameters(resultovy->Parameter(0)*pow(resultovy -> Parameter(3),-1)*pow(2*TMath::Pi(),-0.5), 1+resultovy -> Parameter(1), resultovy -> Parameter(3));
	xc_oh2 = 0.5 - pow(RN_2,2)*TMath::Log(lamb_2);
	//Error tipo 1 ····	
	error_tipo1_N0 = ROOT::Math::normal_cdf_c(xc_oh2, RN_2,0);
	error_tipo1_N1 = lamb_2*(ROOT::Math::normal_cdf(xc_oh2, RN_2,1));

	cout << "pvalue = " << pvalueovy << endl;
	cout << "SER ohdu 2 = " << std::setprecision(9) << lamb_2 / t_expo_overscan_y << " +- " << err_lamb_2 / t_expo_overscan_y << endl;
	cout << "xc ohdu 2 = " << xc_oh2 << endl;	
	cout << "error_tipo1_N0 = " << error_tipo1_N0 << endl;
	cout << "error_tipo1_N1 = " << error_tipo1_N1 << endl;
	SetOhdu1Style(epix_ohdu2_ovy_hist);
	epix_ohdu2_ovy_hist->Draw();
	c2->SaveAs(OVY_out.Data());
}

//Information in leaves to local or global variables  ·································································
// calPixTree --
void Enable_and_Set_Branches_A(TTree* & tree){
	tree->SetBranchStatus("*",0); //disable all branches
	tree->SetBranchStatus("ohdu",1);
	tree->SetBranchStatus("RUNID",1);
  	tree->SetBranchStatus("x",1);
  	tree->SetBranchStatus("y",1);
	tree->SetBranchStatus("ePix",1); 

	tree->SetBranchAddress ("RUNID",&RUNID);
	tree->SetBranchAddress ("ohdu",&ohdu);
	tree->SetBranchAddress ("x",&x);
  	tree->SetBranchAddress ("y",&y);
  	tree->SetBranchAddress ("ePix",&ePix);
}

//header ---	
void Enable_and_Set_Branches_B(TTree* & tree){
    tree->SetBranchStatus("*",0); //disable all branches
    tree->SetBranchStatus("DATESTART",1);//enable branch 
    tree->SetBranchStatus("DATEEND",1); 
    tree->SetBranchStatus("RUNID",1);

    tree->SetBranchAddress ("DATESTART",&DATESTART);//save branch content into variable
    tree->SetBranchAddress ("DATEEND",&DATEEND);
    tree->SetBranchAddress ("RUNID",&RUNID_head);
}

//hitSumm ---	
void Enable_and_Set_Branches_C(TTree* & tree){
    tree->SetBranchStatus("*",0); //disable all branches
    //tree->SetBranchStatus("flag",1);
	//tree->SetBranchStatus("xPix",1);
	//tree->SetBranchStatus("yPix",1);    
    tree->SetBranchStatus("runID",1);
    tree->SetBranchStatus("e",1);
    tree->SetBranchStatus("n",1);
    tree->SetBranchStatus("ohdu",1);
    tree->SetBranchStatus("xBary",1); 
	tree->SetBranchStatus("yBary",1); 
	tree->SetBranchStatus("xMin",1); 
	tree->SetBranchStatus("xMax",1);
	tree->SetBranchStatus("yMin",1); 
	tree->SetBranchStatus("yMax",1);
    
    //tree->SetBranchAddress ("flag",&flag);
	//tree->SetBranchAddress ("xPix",&xPix);
	//tree->SetBranchAddress ("yPix",&yPix);
    tree->SetBranchAddress ("runID",&runID);
    tree->SetBranchAddress ("e",&e);
    tree->SetBranchAddress ("n",&n);
    tree->SetBranchAddress ("ohdu",&ohdu_hs);
    tree->SetBranchAddress ("xBary",&xBary);
	tree->SetBranchAddress ("yBary",&yBary);
	tree->SetBranchAddress ("xMin",&xMin);
	tree->SetBranchAddress ("xMax",&xMax);
	tree->SetBranchAddress ("yMin",&yMin);
	tree->SetBranchAddress ("yMax",&yMax);

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

// Extracting hot columns ·················································································
void extracthc(int vmin[], int vmax[],int arraysize, bool &isnthotcol, Int_t xm , Int_t xM){
    int j = 0;
    while (isnthotcol == true && j < arraysize){
        if (xM >= vmin[j] && xm <= vmax[j]) isnthotcol = 0;
        ++j;
    }
}

// Función para verificar si un archivo existe ·······························································
bool fileExists(const std::string& filename) {
    struct stat buffer;
    return (stat(filename.c_str(), &buffer) == 0);
}

// Calcular cantidad de píxeles útiles por cuadrante ·············
int hot_col_amount(int vmin[], int vmax[], int arraysize){
	int xm{0};
	int j = 0;
	while(j < arraysize){
		xm += vmax[j] - vmin[j] + 1;
		++j;
		if (j == arraysize){cout << xm << endl;}
	}

	return xm;
}

//gaus(0) is a substitute for [0]*exp(-0.5*((x-[1])/[2])**2) and (0) means start numbering parameters at 0. expo(3) is a substitute for exp([3]+[4]*x)

//https://root.cern.ch/root/htmldoc/guides/users-guide/FittingHistograms.html

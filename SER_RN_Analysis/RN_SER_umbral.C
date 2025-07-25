/* Autora: E.Depaoli Fecha inicial: 24/04/2025
Código para calcular el ruido de lectura y el SEE de las imágenes de Atucha 
utilizando los catálogos hits_corr_proc_run*.root. Calcula el p-valor del ajuste.
-> Histograma de ePix en calPixTree restringida al overscan, para cada extensión. 
Ajusta cada histograma con una suma de 2 gaussianas. Guarda en un archivo el sigma (ruido de
lectura), su incerteza y el p-valor del ajuste.
-> Calcula cluster_1e como la suma del número de eventos que ocupan 1 pxl y tienen
1 electrón (n == 1.0 && e == 1.0) en las extensiones 1 y 2. 
SEE = cluster_1e/N_tot_pxls/t_exposición
-> Histograma de carga en el área activa de los cuadrantes 1 y 2 separados (los que tienen menor ruido de lectura) usando
ePix en calPixTree, habiendo eliminado las columnas brillantes que venimos usando. Ajusto con la convolución entre 
Poisson y Gauss. Parámetros que se obtienen de este ajuste: el ruido de lectura, SER*tiempo.

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
string filename{"/home/eliana/Documentos/Atucha/images/run_40/hits/hits_corr_proc_run_40_01Feb2025__EXP1_NSAMP300_VSUB70_img99.root"};
double xc{0.655};
double xc_oh2{0.610};
double xc_ohdu1_ov{0.0};
double xc_ohdu2_ov{0.0};

double Nmax{};
double RN_1{};
double err_RN_1{};
double RN_2{};
double err_RN_2{};
double lamb_1{};
double err_lamb_1{};
double lamb_2{};
double err_lamb_2{};

double lamb_1_ov{};
double err_lamb_1_ov{};
double RN_1_ov{};
double err_RN_1_ov{};

double lamb_2_ov{};
double err_lamb_2_ov{};
double RN_2_ov{};
double err_RN_2_ov{};

// Variables for exporting to archivo
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
Int_t runID; float e; float n; 
Int_t ohdu_hs; Float_t xBary; Float_t yBary;
Int_t xMin; Int_t xMax;
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

unsigned int hot_col_amount = 71 + 106;
double amount_cutted_pix = (xBaryMax - xBaryMin - hot_col_amount)*20*(yBaryMax - yBaryMin);//amount of pixels ohdu1 + ohdu2
float factor = 0.5;//corrección por lectura entre dos imágenes sucesivas

//Archivo de salida
string ofilename = "ruido_de_lectura_Atucha.txt";
string onoise_DCfilename = "ruido_de_lectura_SEE_Atucha.txt";
string oumbralextract_filename = "th_for_skextract";

//Función para calcular SER ·························

Double_t poisson_gauss_conv(Double_t *x, Double_t *par){//lambda [e-/pxl] , mu, normalización

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
void RN_SER_umbral(){
	SetGlobalStyle();

	// Genero histogramas 	········
	int nbins_epix_overscan = 35;
	float lower_bin = -1.0;
	float upper_bin = 1.8;
	TH1D * epix_ohdu1_overscan_hist = new TH1D("epix overscan ohdu1", "epix overscan", nbins_epix_overscan, lower_bin, upper_bin);// ohdu1
	TH1D * epix_ohdu2_overscan_hist = new TH1D("epix overscan ohdu2", "epix overscan", nbins_epix_overscan, lower_bin, upper_bin);
	TH1D * epix_ohdu1_actar_hist = new TH1D("epix active area ohdu1", "epix active area", 50, -0.5, 1.5);
	TH1D * epix_ohdu2_actar_hist = new TH1D("epix active area ohdu2", "epix active area", 150, -0.5, 1.5);

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
		
	// hitSumm Tree ----------------------------------------------------------------
	//SEE. Viejo. ····
    double clusters_1e {0};//SEE
    int occupancy {0};
    Double_t see_aux{0};
    Double_t occ_aux{0};
    cout << "Cantidad de entradas en hitSumm = " << Entries_hStree << endl;

    for (Long64_t entry{0}; entry < Entries_hStree ; ++entry)
    {
        hStree -> GetEntry(entry);
    	isnthotcol=true;
      
    	if (ohdu_hs == ohdu_1 || ohdu_hs == ohdu_2){//events in quadrants with single electron counting
    		if (xBary<xBaryMax && xBary>xBaryMin && yBary<yBaryMax && yBary>yBaryMin){//events inside edges of CCD
    			//Extraer columnas brillantes no reduce el valor final de clusters_1e porque estos no contienen eventos con n > 1. Pero debería reducir la ocupancia ¡y lo hace!
    			if (ohdu_hs==ohdu_1){//events outside hot columns
        		    extracthc(hotcol_ohdu_1_min, hotcol_ohdu_1_max,sizehotcol_1_min, isnthotcol, xMin, xMax);
		        } else if (ohdu_hs==ohdu_2){extracthc(hotcol_ohdu_2_min, hotcol_ohdu_2_max,sizehotcol_2_min, isnthotcol, xMin, xMax);}

    			if (isnthotcol==true){//si el evento no contiene ni toca una columna brillante
			        if (n == 1.0 && e == 1.0){++clusters_1e;}//n: pixels occupied, e: electrons in cluster
			        occupancy+=n;
			        //cout << n << endl;
			    }
		 	}       
	    }

    //··················
	//Cálculos 		····
	//··················

    //Luego de pasar por el último evento, ocupancia y DC 
        if(entry > (Entries_hStree-2)){
        	//cout << "final entry = " << entry << endl;
            cout << "clusters_1e =" << clusters_1e << " [pixels]" <<'\n';
            cout << "occupancy = " << occupancy  << " [pixels]" << endl;
            time_expo = imagedurationindays(&DATESTART[0],&DATEEND[0], factor);
            cout << "imagedurationindays = " << setprecision(4)<< time_expo << '\n';
            cout << "amount of useful pixels = " << amount_cutted_pix << endl;
            see_aux=clusters_1e/imagedurationindays(&DATESTART[0],&DATEEND[0],factor)/amount_cutted_pix;
            occ_aux=occupancy/amount_cutted_pix;//imagedurationindays(&DATESTART[0],&DATEEND[0],factor)/
			//variables que usaré en la próxima ronda se inicializan         
            clusters_1e = 0;
            occupancy = 0;
        }
    }

    cout << "see_aux = " << see_aux << endl;
    cout << "occ_aux = " << occ_aux << endl;

	// CalPixTree ---------------------------------------------------------------
	for(int j{0}; j < Entries_cptree; ++j){
		//Ruido de lectura ····
		cptree->GetEntry(j);
		//epix en overscan
		if(x > 307 && ePix < 10 && ePix > -1){//
			switch(ohdu){
			case 1:
				epix_ohdu1_overscan_hist->Fill(ePix);
			case 2:
				epix_ohdu2_overscan_hist->Fill(ePix);
			}
		}

		//SER ····
		//epix en área activa
		//Extraer columnas brillantes, funciona.
		isnthotcol=true;
		if(x < 307 && ePix < 10 && ePix > -1){
			if (ohdu == ohdu_1 || ohdu == ohdu_2){
				if (ohdu==ohdu_1){//events outside hot columns
				    extracthc(hotcol_ohdu_1_min, hotcol_ohdu_1_max,sizehotcol_1_min, isnthotcol, x, x);
				    if (isnthotcol==true){//si el pixel no está en una columna brillante
						epix_ohdu1_actar_hist->Fill(ePix);
			    	}
		        } else if (ohdu==ohdu_2){
		        	extracthc(hotcol_ohdu_2_min, hotcol_ohdu_2_max,sizehotcol_2_min, isnthotcol, x, x);
		        	if (isnthotcol==true){//si el pixel no está en una columna brillante
						epix_ohdu2_actar_hist->Fill(ePix);
			    	}
	    	    }
			}
		}
	}

	//Normalizo histogramas ····
	//epix_ohdu1_actar_hist->Scale( 1./epix_ohdu1_actar_hist->Integral());//Normalize to a probability distribution. Sum of histogram content is equal to 1.
	epix_ohdu1_actar_hist->Scale( 1./epix_ohdu1_actar_hist->Integral(),"WIDTH");//histogram describe a probability density function. 
	epix_ohdu1_overscan_hist->Scale(1./epix_ohdu1_overscan_hist->Integral(), "WIDTH");
	epix_ohdu2_actar_hist->Scale( 1./epix_ohdu2_actar_hist->Integral(),"WIDTH");
	epix_ohdu2_overscan_hist->Scale(1./epix_ohdu2_overscan_hist->Integral(), "WIDTH");
	//Clones de ohdu histograms ···
	TH1*ohdu1_overscan_hc = (TH1*)epix_ohdu1_overscan_hist->Clone();
	TH1*ohdu2_overscan_hc = (TH1*)epix_ohdu2_overscan_hist->Clone();
	
	//··················
	//Ajustes 		····
	//··················
	//OVERSCAN ··············· Solo Ruido de lectura: Suma de dos gaussianas con igual varianza ·······················
	TF1*total_gauss = new TF1("G0+G1","[0]*exp(-0.5*pow((x-[1])/[2],2)) + [3]*exp(-0.5*pow((x-[4])/[2],2))", -1.0,1.5);//
	total_gauss->SetLineColor(kMagenta-4);
	//Inicializo parámetros
	total_gauss->SetParameters(1.5,0.0003, 0.22, 3.0, 1.0);
	total_gauss->SetParLimits(2,0.17, 0.25);
	//total_gauss->FixParameter(1,0);
	//total_gauss->FixParameter(4,1);
	
	// ohdu 1 ···
	/*cout << '\n' << endl;
	cout << " ---------- OHDU 1 OVERSCAN G0+G1 -------------" << endl;
	//“Q” Quiet mode (minimum printing); “S” The result of the fit is returned in the TFitResultPtr; "R" to restrict the fit to the range specified in the TF1 constructor
	auto result01=epix_ohdu1_overscan_hist->Fit("G0+G1", "R+S");
	double pvalue01 = ROOT::Math::chisquared_cdf_c(result01->Chi2(), result01->Ndf());
	RN_1 = result01 -> Parameter(2);
	err_RN_1 = result01 -> ParError(2);
	cout << '\n' << endl;
	cout << "pvalue OVERSCAN = " << pvalue01 << endl;
	cout << "" << RN_1 << endl;
	cout << err_RN_1 << endl;

	// ohdu 2 ···
	cout << '\n' << endl;
	cout << " ------------- OHDU 2 OVERSCAN G0+G1 ------------" << endl;
	result01=epix_ohdu2_overscan_hist->Fit("G0+G1", "R+S");
	pvalue01 = ROOT::Math::chisquared_cdf_c(result01->Chi2(), result01->Ndf());
	RN_2 = result01 -> Parameter(2);
	err_RN_2 = result01 -> ParError(2);
	cout << '\n' << endl;
	cout << "pvalue OVERSCAN  = " << pvalue01 << endl;
	cout << RN_2 << endl;
	cout << err_RN_2 << endl;
*/
	// SER + Ruido de lectura: Poisson convolucionada con gaussiana ·······················
	TF1 *pgc = new TF1("poisson_gauss_conv", poisson_gauss_conv, -1, 2,4);//(,, x_lim_inf, x_lim_sup, #parameters)
	Nmax= 1;
	//pgc->SetParameters(0.2,0.03, 1.0); //lambda, mu, normalización. Inicializo.
	pgc->SetParameters(0.2,0.03, 1.0, 0.21); //lambda, mu, normalización, sigma. Inicializo.
	pgc->SetLineColor(kMagenta-4);//kOrange-3);
	
	//OVERSCAN ··············
	//ohdu 1 ····
	cout << '\n' << endl;
	cout << " --------- OHDU 1 OVERSCAN Convolution Poisson-Gauss--------- " << endl;
	auto result03 = ohdu1_overscan_hc->Fit("poisson_gauss_conv", "R+S");
	double pvalue03 = ROOT::Math::chisquared_cdf_c(result03->Chi2(), result03->Ndf());
	lamb_1_ov = result03 -> Parameter(0);
	err_lamb_1_ov = result03 -> ParError(0);
	RN_1_ov = result03 -> Parameter(3);
	err_RN_1_ov = result03 -> ParError(3);
	xc_ohdu1_ov = 0.5 - pow(RN_1_ov,2)*TMath::Log(lamb_1_ov);

	cout << '\n' << endl;
	cout << " OHDU 1 " << endl;
	cout << "pvalue OVERSCAN  = " << pvalue03 << endl;
	cout << "lambda ohdu 1 = " << std::setprecision(9) << lamb_1_ov << " +- " << err_lamb_1_ov << endl;
	cout << "sigma ohdu 1 = " << std::setprecision(9) << RN_1_ov << " +- " << err_RN_1_ov << endl;
	cout << "xc = " << xc_ohdu1_ov << endl;
	
	//ohdu 2 ····
	cout << '\n' << endl;
	cout << " --------- OHDU 2 OVERSCAN  Convolution Poisson-Gauss --------- " << endl;
	result03 = ohdu2_overscan_hc->Fit("poisson_gauss_conv", "R+S");
	pvalue03 = ROOT::Math::chisquared_cdf_c(result03->Chi2(), result03->Ndf());
	lamb_2_ov = result03 -> Parameter(0);
	err_lamb_2_ov = result03 -> ParError(0);
	RN_2_ov = result03 -> Parameter(3);
	err_RN_2_ov = result03 -> ParError(3);
	xc_ohdu2_ov = 0.5 - pow(RN_2_ov,2)*TMath::Log(lamb_2_ov);

	cout << '\n' << endl;
	cout << " OHDU 1 " << endl;
	cout << "pvalue OVERSCAN  = " << pvalue03 << endl;
	cout << "lambda ohdu 2 = " << std::setprecision(9) << lamb_2_ov << " +- " << err_lamb_2_ov << endl;
	cout << "sigma ohdu 2 = " << std::setprecision(9) << RN_2_ov << " +- " << err_RN_2_ov << endl;
	cout << "xc = " << xc_ohdu2_ov << endl;
	

	//N(0,noise) y N(1,noise) ··········
	TF1 * N0noise = new TF1("normal_0_noise", "gaus(0)", -0.5,1.5);//
	TF1 * N1noise = new TF1("normal_1_noise", "gaus(0)", -0.5,1.5);//gaus(0) = [0]*exp(-0.5*((x-[1])/[2])**2) and (0) means start numbering parameters at 0
	TF1 * N0noise_oh2 = new TF1("normal_0_noise", "gaus(0)", -0.5,1.5);
	TF1 * N1noise_oh2 = new TF1("normal_1_noise", "gaus(0)", -0.5,1.5);
	N0noise->SetLineColor(kGray+2);
	N1noise->SetLineColor(kOrange+2);
	N0noise_oh2->SetLineColor(kGray+2);
	N1noise_oh2->SetLineColor(kOrange+2);
	
	//AREA ACTIVA ··············
	//ohdu 1 ····
	cout << '\n' << endl;
	cout << " --------- OHDU 1 ACTIVE AREA --------- " << endl;
	auto result02 = epix_ohdu1_actar_hist->Fit("poisson_gauss_conv", "R+S");
	double pvalue02 = ROOT::Math::chisquared_cdf_c(result02->Chi2(), result02->Ndf());
	lamb_1 = result02 -> Parameter(0);
	err_lamb_1 = result02 -> ParError(0);
	RN_1 = result02 -> Parameter(3);
	err_RN_1 = result02 -> ParError(3);
	xc = 0.5 - pow(RN_1,2)*TMath::Log(lamb_1);

	cout << '\n' << endl;
	cout << " OHDU 1 " << endl;
	cout << "pvalue active area  = " << pvalue02 << endl;
	cout << "lambda ohdu 1 = " << std::setprecision(9) << lamb_1 << " +- " << err_lamb_1 << endl;
	cout << "sigma ohdu 1 = " << std::setprecision(9) << RN_1 << " +- " << err_RN_1 << endl;
	cout << "xc = " << xc << endl;
	
	//For plotting ·····
	N0noise->SetParameters(pow(result02 -> Parameter(3),-1)*pow(2*TMath::Pi(),-0.5), result02 -> Parameter(1), result02 -> Parameter(3));//result02 -> Parameter(2)
	N1noise->SetParameters(result02->Parameter(0)*pow(result02 -> Parameter(3),-1)*pow(2*TMath::Pi(),-0.5), 1+result02 -> Parameter(1), result02 -> Parameter(3));
	//Error tipo 1 resultante ····· 
	double error_tipo1_N0 = ROOT::Math::normal_cdf_c(xc,RN_1);//, 1);// (double x, double sigma=1, double x0=0)
	double error_tipo1_N1 = lamb_1*(ROOT::Math::normal_cdf(xc,RN_1,1));//
	cout << "error_tipo1_N0 = " << error_tipo1_N0 << endl;
	cout << "error_tipo1_N1 = " << error_tipo1_N1 << endl;
	cout << '\n' << endl;
	double lamb_recal = ROOT::Math::normal_cdf(-xc/RN_1)/ROOT::Math::normal_cdf((xc-1)/RN_1);
	cout << "lamb_recal = " << lamb_recal << endl;

	//ohdu 2 ····
	cout << '\n' << endl;
	cout << " --------- OHDU 2 ACTIVE AREA -------- " << endl;
	result02 = epix_ohdu2_actar_hist->Fit("poisson_gauss_conv", "R+S");
	pvalue02 = ROOT::Math::chisquared_cdf_c(result02->Chi2(), result02->Ndf());
	lamb_2 = result02 -> Parameter(0);
	err_lamb_2 = result02 -> ParError(0);
	RN_2 = result02 -> Parameter(3);
	err_RN_2 = result02 -> ParError(3);
	N0noise_oh2->SetParameters(pow(result02 -> Parameter(3),-1)*pow(2*TMath::Pi(),-0.5), result02 -> Parameter(1), result02 -> Parameter(3));
	N1noise_oh2->SetParameters(result02->Parameter(0)*pow(result02 -> Parameter(3),-1)*pow(2*TMath::Pi(),-0.5), 1+result02 -> Parameter(1), result02 -> Parameter(3));
	//Error tipo 1 ····	
	xc = 0.5 - pow(RN_2,2)*TMath::Log(lamb_2);	
	error_tipo1_N0 = ROOT::Math::normal_cdf_c(xc, result02 -> Parameter(3),0);
	error_tipo1_N1 = result02->Parameter(0)*(ROOT::Math::normal_cdf(xc, result02 -> Parameter(3),1));

	cout << '\n' << endl;
	cout << " OHDU 2 " << endl;
	cout << "pvalue active area  = " << pvalue02 << endl;
	cout << "lambda ohdu 2 = " << std::setprecision(9) << lamb_2 << " +- " << err_lamb_2 << endl;
	cout << "sigma ohdu 2 = " << std::setprecision(9) << RN_2 << " +- " << err_RN_2 << endl;
	cout << "xc = " << xc << endl;
	cout << "error_tipo1_N0 = " << error_tipo1_N0 << endl;
	cout << "error_tipo1_N1 = " << error_tipo1_N1 << endl;
	cout << '\n' << endl;
	
	//Graficar resultados
	/*TCanvas * c_pgc = new TCanvas("c_pgc", "SEE_Sim", 1200, 800);
	c_pgc->SetGrid();
	pgc->Draw();
	*/

	//Grafico OVERSCAN ···
	TCanvas * c1 = new TCanvas("c1", "Charge histogram", 1200, 800);
	//gStyle->SetOptStat(1);
	c1 ->SetLogy();
	c1->Divide(2,1);
	//ohdu 1 ····
	c1->cd(1);
	c1->SetGrid();
	SetOhdu1Style(epix_ohdu1_overscan_hist);
	epix_ohdu1_overscan_hist->SetTitle("Overscan");//calPixTree
	epix_ohdu1_overscan_hist->GetYaxis()->SetRangeUser(1e-4,1e1);
	epix_ohdu1_overscan_hist->SetStats(0);
	gPad->SetLogy();
    gPad->SetGrid();

	epix_ohdu1_overscan_hist->Draw();
	auto l1_RN = new TLegend(0.55,0.65,0.90,0.80); 
    l1_RN->AddEntry(epix_ohdu1_overscan_hist,"Measured charge","p");//ohdu 1
    l1_RN->AddEntry(total_gauss,"N(#mu_{0}, #sigma_{RN}) + N(#mu_{1}, #sigma_{RN})","l");
    l1_RN->Draw();
    //ohdu 2 ····
	c1->cd(2);
	SetOhdu2Style(epix_ohdu2_overscan_hist);
	gPad ->SetLogy();
    gPad->SetGrid();
	epix_ohdu2_overscan_hist->SetTitle("Overscan");
	epix_ohdu2_overscan_hist->GetYaxis()->SetRangeUser(1e-4,1e1);
	epix_ohdu2_overscan_hist->Draw();//"same"

	auto l2_RN = new TLegend(0.55,0.65,0.90,0.75); 
    l2_RN->AddEntry(epix_ohdu2_overscan_hist,"ohdu 2","l");
    l2_RN->AddEntry(total_gauss,"N(0,RN) + N(1,RN)");
    l2_RN->Draw();

    //Poisson convolucionada con Gauss ·····························
    TCanvas * c2 = new TCanvas("c2", "Charge histogram", 1200, 800);
	//gStyle->SetOptStat(1);
	c2 ->SetLogy();
	c2->Divide(2,1);
	//ohdu 1 ····
	c2->cd(1);
	c2->SetGrid();
	ohdu1_overscan_hc->Draw();

	//ohdu 2 ····
	c2->cd(2);
	c2->SetGrid();
	ohdu2_overscan_hc->Draw();

	//Graficos en area activa ·············	
	TCanvas * c1AA = new TCanvas("c1AA", "Charge histogram active area", 1200, 800);
    gStyle->SetOptStat(1);
	c1AA->Divide(2,1);
	//ohdu 1 ····
	c1AA->cd(1);
	//gPad->SetLogy();
    gPad->SetGrid();

    epix_ohdu1_actar_hist->SetTitle("calPixTree");
    SetOhdu1Style(epix_ohdu1_actar_hist);
    epix_ohdu1_actar_hist->GetYaxis() -> SetRangeUser(1e-3, 2);
	epix_ohdu1_actar_hist->Draw();
	N0noise->Draw("same");
	//Fill area beneath the curve
	TF1* N0Range = (TF1*)N0noise->Clone("N0Range");
	N0Range->SetRange(xc, 1.5);
 	N0Range->SetFillColorAlpha(kGray+1,0.35);
 	N0Range->SetFillStyle(1001);
 	N0Range->Draw("SAME FC");

	N1noise->Draw("same");
	//Fill area beneath the curve
	TF1* N1Range = (TF1*)N1noise->Clone("N1Range");
	N1Range->SetRange(-0.5, xc);
 	N1Range->SetFillColorAlpha(kOrange+1,0.35);
 	N1Range->SetFillStyle(1001);
 	N1Range->Draw("SAME FC");

 	TLine*xc_line = new TLine(xc,2e-3,xc,1);
 	xc_line -> SetLineColor(kRed);
 	xc_line ->SetLineWidth(2);
 	xc_line -> SetLineStyle(kDashed);

 	xc_line->Draw("same");

	auto l1_DC = new TLegend(0.25, 0.15,0.55,0.30); 
    l1_DC->AddEntry(epix_ohdu1_actar_hist,"ohdu 1","l");
    l1_DC->AddEntry(pgc, "ajuste", "l");
    l1_DC->AddEntry(N0noise, "N(0,#sigma_{RN})", "l");
    l1_DC->AddEntry(N1noise, "N(1,#sigma_{RN})", "l");
    l1_DC->AddEntry(xc_line, Form("x_{c} = %.3f",xc), "l");
    l1_DC->Draw();

    //ohdu 2 ····
	c1AA->cd(2);
	//gPad ->SetLogy();
    gPad->SetGrid();

	SetOhdu2Style(epix_ohdu2_actar_hist);
    epix_ohdu2_actar_hist->SetTitle("calPixTree");
    epix_ohdu2_actar_hist->GetYaxis() -> SetRangeUser(1e-3, 2);
	epix_ohdu2_actar_hist->Draw();	
	N0noise_oh2->Draw("same");
	//Fill area beneath the curve
	TF1* N0Range_oh2 = (TF1*)N0noise_oh2->Clone("N0Range_oh2");
	N0Range_oh2->SetRange(xc_oh2, 1.5);
 	N0Range_oh2->SetFillColorAlpha(kGray+1,0.35);
 	N0Range_oh2->SetFillStyle(1001);
 	N0Range_oh2->Draw("SAME FC");

	N1noise_oh2->Draw("same");
	//Fill area beneath the curve
	TF1* N1Range_oh2 = (TF1*)N1noise_oh2->Clone("N1Range_oh2");
	N1Range_oh2->SetRange(-0.5, xc_oh2);
 	N1Range_oh2->SetFillColorAlpha(kOrange+1,0.35);
 	N1Range_oh2->SetFillStyle(1001);
 	N1Range_oh2->Draw("SAME FC");

 	TLine*xc_line_oh2 = new TLine(xc_oh2,1e-3,xc_oh2,1);
 	xc_line_oh2 -> SetLineColor(kRed);
 	xc_line_oh2 ->SetLineWidth(2);
 	xc_line_oh2 -> SetLineStyle(kDashed);
 	xc_line_oh2->Draw("same");

	auto l2_DC = new TLegend(0.25, 0.15,0.55,0.30); 
    l2_DC->AddEntry(epix_ohdu2_actar_hist,"ohdu 2","l");
    l2_DC->AddEntry(pgc, "ajuste", "l");
    l2_DC->AddEntry(N0noise_oh2, "N(0,#sigma_{RN})", "l");
    l2_DC->AddEntry(N1noise_oh2, "N(1,#sigma_{RN})", "l");
    l2_DC->AddEntry(xc_line, Form("x_{c} = %.3f",xc_oh2), "l");
    l2_DC->Draw();

	// Resultados a archivo

	// Verifico si ya existe el archivo donde se guardarán los resultados
    bool ofile_already_exists = fileExists(ofilename);
    bool onoise_DCfile_already_exists = fileExists(onoise_DCfilename);
    //cout << ofile_already_exists << endl;
    //cout << onoise_DCfile_already_exists << endl;

	ofstream ofile(ofilename, std::ios::app);
	ofstream onoise_DCfile(onoise_DCfilename, std::ios::app);


	// Si el archivo no existía, escribir la línea de encabezado
    /*if (!ofile_already_exists) {
        ofile << "RUNID " << '\t' << "DATESTART"<< '\t' << "ohdu " << '\t' << "N [pxl]" << '\t' << "Error [pxl]" << '\t' << "Mean [e-]" << '\t' << "Error [e-]" << '\t' << "Std Dev [e-]" << '\t' << "Error [e-]" << '\t' << "p-valor" << '\t' << "ohdu " << '\t' << "N [pxl]" << '\t' << "Error [pxl]" << '\t' << "Mean [e-]" << '\t' << "Error [e-]" << '\t' << "Std Dev [e-]" << '\t' << "Error [e-]" << '\t' << "p-valor" << '\t' << "ohdu " << '\t' << "N [pxl]" << '\t' << "Error [pxl]" << '\t' << "Mean [e-]" << '\t' << "Error [e-]" << '\t' << "Std Dev [e-]" << '\t' << "Error [e-]" << '\t' << "p-valor" << '\t' << "ohdu " << '\t' << "N [pxl]" << '\t' << "Error [pxl]" << '\t' << "Mean [e-]" << '\t' << "Error [e-]" << '\t' << "Std Dev [e-]" << '\t' << "Error [e-]" << '\t' << "p-valor" << '\n';
    }
		
	ofile << RUNID_head << '\t' << &DATESTART[0] << '\t' << "1" << '\t' << result->Parameter(0) << '\t' << result->ParError(0) << '\t' << result->Parameter(1) << '\t' << result->ParError(1) << '\t' << result->Parameter(2) << '\t' << result->ParError(2) << '\t' << pvalue0 << '\t'; 
	
	if(!onoise_DCfile_already_exists) onoise_DCfile << "RUNID " << '\t' << "DATESTART"<< '\t' << "ohdu " << '\t' << "Readout Noise [e-]" << '\t' << "Error [e-]" << '\t' << "ohdu " << '\t' << "Readout Noise [e-]" << '\t' << "Error [e-]"<< '\t' << "SEE[e-/pxl/day]" << '\t' << "OCC [n/pxl/day]" << '\t';
	
	onoise_DCfile << RUNID_head << '\t' << &DATESTART[0] << '\t'  << "1" << '\t' << result->Parameter(2) << '\t' << result->ParError(2) << '\t';	
*/
	// ohdu 2 ----------------------------------------------------------------------
	/*G0->SetParameters(1400,epix_ohdu2_overscan_hist->GetMean(),epix_ohdu2_overscan_hist->GetStdDev());
	result = epix_ohdu2_overscan_hist->Fit("G0", "SQ");
	pvalue0 = ROOT::Math::chisquared_cdf_c(result->Chi2(), result->Ndf());
	// Información a archivo
	ofile << "2" << '\t' << result->Parameter(0) << '\t' << result->ParError(0) << '\t' << result->Parameter(1) << '\t' << result->ParError(1) << '\t' << result->Parameter(2) << '\t' << result->ParError(2) << '\t' << pvalue0 << '\t';
	onoise_DCfile << "2" << '\t' << result->Parameter(2) << '\t' << result->ParError(2) << '\t';
	// Información a archivo
	onoise_DCfile << see_aux << '\t' << occ_aux << '\n';
*/
	// Cierro el archivo
    /*ofile.close();
    onoise_DCfile.close();
    file->Close();
*/
}

//Information in leaves to local or global variables  ·································································
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

void Enable_and_Set_Branches_B(TTree* & tree){
    tree->SetBranchStatus("*",0); //disable all branches
    tree->SetBranchStatus("DATESTART",1);//enable branch 
    tree->SetBranchStatus("DATEEND",1); 
    tree->SetBranchStatus("RUNID",1);

    tree->SetBranchAddress ("DATESTART",&DATESTART);//save branch content into variable
    tree->SetBranchAddress ("DATEEND",&DATEEND);
    tree->SetBranchAddress ("RUNID",&RUNID_head);
}

void Enable_and_Set_Branches_C(TTree* & tree){
    tree->SetBranchStatus("*",0); //disable all branches
    //tree->SetBranchStatus("flag",1);
    tree->SetBranchStatus("runID",1);
    tree->SetBranchStatus("e",1);
    tree->SetBranchStatus("n",1);
    tree->SetBranchStatus("ohdu",1);
    tree->SetBranchStatus("xBary",1); 
	tree->SetBranchStatus("yBary",1); 
	tree->SetBranchStatus("xMin",1); 
	tree->SetBranchStatus("xMax",1);    
    
    //tree->SetBranchAddress ("flag",&flag);
    tree->SetBranchAddress ("runID",&runID);
    tree->SetBranchAddress ("e",&e);
    tree->SetBranchAddress ("n",&n);
    tree->SetBranchAddress ("ohdu",&ohdu_hs);
    tree->SetBranchAddress ("xBary",&xBary);
	tree->SetBranchAddress ("yBary",&yBary);
	tree->SetBranchAddress ("xMin",&xMin);
	tree->SetBranchAddress ("xMax",&xMax);
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

//gaus(0) is a substitute for [0]*exp(-0.5*((x-[1])/[2])**2) and (0) means start numbering parameters at 0. expo(3) is a substitute for exp([3]+[4]*x)

//https://root.cern.ch/root/htmldoc/guides/users-guide/FittingHistograms.html

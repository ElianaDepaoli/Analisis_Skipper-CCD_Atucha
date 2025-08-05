/* Autora: E.Depaoli Fecha inicial: 01/08/2025
Última modificación: 04/08/2025
Código para graficar histogramas de carga de las imágenes de Atucha utilizando los catálogos hits_corr_proc_run*.root
de una carpeta. 
Procesa los archivos que están en dirname+subdirname.

Importante: 
Correr con la siguiente linea en consola para que no imprima en pantalla cada gráfico:
root -b -q histogramas_de_carga_app.C

Si no se usa la linea anterior, se bloquea la pantalla por la emergencia de gráficos.

Para usarlo
MODO CARPETA COMPLETA
main_dirname = ruta completa a la carpeta donde están los catálogos a leer; 
sub_dirname = nombre de la carpeta que contiene los catálogos;
Comentar las lineas con comentario://to run just one file using filename
Descomentar las lineas con comentario: //to run over a full folder.
MODO ÚNICO ARCHIVO
Descomentar las lineas con comentario://to run just one file using filename
Comentar las lineas con comentario: //to run over a full folder.

Lo que produce:
-> Histogramas de carga las extensiones 1 y 2, sobre overscan x y area activa. Los guarda como *.png

Detalle:

-> Histograma de carga sin clusterizar (ePix en calPixTree). 
-> No se contabiliza carga en píxeles de columnas brillantes.

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
//Función para calcular SER ·························
double Nmax{};
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

//Funcion para normalizar histograma 2D ············
void normaliza_minmax(TH2D * h2){

	// Copiar todos los valores de los bins del TH2D a un TH1D temporal
	TH1D h1d_temp("h1d_temp", "Temp histogram", 20e6, 0, 0); // Bins dummy

	for (int i = 1; i <= h2->GetNbinsX(); ++i) {
	    for (int j = 1; j <= h2->GetNbinsY(); ++j) {
	        double content = h2->GetBinContent(i, j);
	        if (content != 0) { // Ignorar bins vacíos (opcional)
	            h1d_temp.Fill(content);
	        }
	    }
	}

	//double q_min, q_max;
	Double_t probl[1] = {0.01};
	Double_t probh[1] = {0.99};

	//h1d_temp.GetQuantiles(1, &q_min, probl); // Percentil 1 %
	//h1d_temp.GetQuantiles(1, &q_max, probh); // Percentil 99 %

	double q_min = h2->GetMinimum();
	double q_max = h2->GetMaximum();

		for (int i{1}; i <= h2->GetNbinsX(); ++i){
			for(int j{1}; j <= h2->GetNbinsY(); ++j){
				h2->AddBinContent(i,j,-q_min);
			}
		}

		h2->Scale(1./(q_max-q_min));
	}	
// ----------------------------------------------------------------------------------

//Configuracion global para los gráficos --------------------------
void SetGlobalStyle(){
	TStyle *st1 = new TStyle("st1","my style");
	st1->SetTitleSize(0.04,"XYZ");
	st1->SetTitleX(0.5);     // Posición X normalizada (0.5 = centro)
    //st1->SetTitleY(0.98);    // Posición Y normalizada (0.98 = cerca del borde superior)
	st1->SetLabelSize(0.04, "XYZ");
	st1->SetPadBorderMode(0);
	st1->SetCanvasColor(0);
	//Títulos
	st1->SetTitleStyle(0);               // Fondo transparente (0 = sin relleno)
    st1->SetTitleBorderSize(0);          // Sin borde
    st1->SetTitleFillColor(0);          // Fondo transparente (color 0 = blanco/transparente)

	st1->SetStatColor(0);
	st1->SetTitleOffset(1.3);
	gStyle->SetStatStyle(0);
	//st1->SetLegendTextSize(0.035);
	st1->SetLegendBorderSize(-1);
	st1->SetGridColor(0);
	gStyle->SetLegendFillColor(0);

	// Ajustar márgenes para mejor visualización
    st1->SetPadLeftMargin(0.12);   // Aumentar margen izquierdo (default suele ser ~0.10)
    st1->SetPadRightMargin(0.05);  // Margen derecho
    st1->SetPadTopMargin(0.06);    // Margen superior
    st1->SetPadBottomMargin(0.10); // Margen inferior
	
	gROOT->ForceStyle();//
	gROOT->SetStyle("st1");//
	st1->cd ();//This is now the current style gStyle
}

void SetOhdu1Style(TH1 *h1){
	h1->GetYaxis()->SetTitleOffset(1.5);
	h1->SetMarkerColor(kBlue+1);
	h1->SetMarkerStyle(20);
	h1->SetMarkerSize(1);
	h1->SetMarkerStyle(20);
	h1->SetMarkerSize(1);
	h1->SetLineColor(kBlue+1);
	h1->SetFillColor(kBlue+1);
	h1->SetLineWidth(2);
	h1->GetYaxis()->SetRangeUser(0.1, 5e5);
	h1->GetXaxis()->SetTitle("ePix [e-]");
	h1->GetYaxis()->SetTitle("Number");
}

void SetStyle2(TH1 *h1){
	h1->SetMarkerColor(kOrange+1);
	h1->SetMarkerStyle(20);
	h1->SetMarkerSize(1);
	h1->SetMarkerStyle(20);
	h1->SetMarkerSize(1);
	h1->SetLineColor(kOrange+1);
	h1->SetFillColor(kOrange+1);
	h1->SetLineWidth(2);
	h1->GetXaxis()->SetTitle("ePix [e-]");
	h1->GetYaxis()->SetTitle("Number");
}
	
// Global Variables ---------------------------------------------------------------------
//Archivos y carpetas de entrada ········
//string filename{"/home/eliana/Documentos/Atucha/images/run_46_29Jul2025/proc/hits/hits_corr_proc_run_46_29Jul2025__EXP1_NSAMP300_VSUB70_NROW100_img11.root"};
//string filename{"/home/eliana/Documentos/Atucha/images/run_48_01Ago2025_test_vl-0_5/proc/hits/hits_corr_proc_run_48_01Ago2025_test_vl-0_5__EXP1_NSAMP300_VSUB70_NROW100_img30.root"};
//const char* main_dirname = "/media/eliana/0e0c1913-ca43-4252-b2bb-783703f57ea5/proc/";
const char* main_dirname = "/home/eliana/Documentos/Atucha/images/";
//const char* sub_dirname = "run_46_29Jul2025/proc/hits/";//"run_39/hits/subset/";
const char* sub_dirname = "run_47_31Jul2025_test_vl-1_0/proc/hits/";
//const char* sub_dirname = "run_48_01Ago2025_test_vl-0_5/proc/hits/";
	
//Carpetas salida
//const char* cfulloutdirname = gSystem->pwd();
//Creo una carpeta con igual nombre que el run que proceso
TString run_folder = TString(sub_dirname);
run_folder = run_folder(0, run_folder.First('/'));

// Crear carpeta nueva dentro del directorio actual
TString cfulloutdirname = TString(gSystem->pwd()) + "/" + run_folder;
//gSystem->mkdir(cfulloutdirname);  // kTRUE crea también si ya existe

//Estimadores ·············································································
// Timming ····
double image_duration_day{0};
double exposure_active_area{0};
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
Char_t NCOL[256];//se llena de [0] a [2], las demás posiciones quedan vacías
Char_t NROW[256];
Char_t NSAMP[256];
Char_t CCDNCOL [256];

//hitSumm ---	
int SR = 307; int filas_aa = 512;
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
int hotcol_ohdu_2_min[] = {6, 53, 58, 199, 307};
int hotcol_ohdu_2_max[] = {9, 54, 60, 199, 307};

bool isnthotcol;//to remove hot columns
int sizehotcol_1_min{sizeof(hotcol_ohdu_1_min)/sizeof(hotcol_ohdu_1_min[0])};
int sizehotcol_2_min{sizeof(hotcol_ohdu_2_min)/sizeof(hotcol_ohdu_2_min[0])};

// remove events partially inside the active area
//Borders of the active area
int xBaryMin=3; // 
int xBaryMax=305; // 5 on the left from the overscan
int yBaryMin=3; 
int yBaryMax=545;

float factor = 1.0;//0.5;//corrección por lectura entre dos imágenes sucesivas

// Rutina  ----------------------------------------------------------------------------------
//void histogramas_de_carga_app() {//to run just one file using filename
void ProcessImage(const char* filename){//to run over a full folder
//SetGlobalStyle();

	//Open root catalogue
	//cout << "filename: " << cfulloutdirname << endl;
	TFile * file = TFile::Open(filename,"READ");//to run over a full folder
	//TFile * file = TFile::Open(filename.c_str(),"READ");// //to run just one file using filename
	//c_str() converts a string to an array of characters & terminates it with a null character. No parameters are allowed, a pointer to the array is returned.
	if (!file->IsOpen()) {std::cerr << "Failed to load file" << filename << std::endl;}

	////Archivo de salida
	TObjArray *tokens = TString(filename).Tokenize("_, .");
	TString OVX_out = ((TObjString*)tokens->At(12))->GetString() + ".png";// ((TObjString*)tokens->At(18))->GetString()
	TString h2_out = "imagen_" + ((TObjString*)tokens->At(12))->GetString() + ".png";
	delete tokens;

/*	
	std::cout << "sub_dirname: " << sub_dirname << std::endl;
	std::cout << "run_folder extraído: " << run_folder << std::endl;
	cout << "OVX_out " << OVX_out << endl; 
	cout << "h2_out " << h2_out << endl;
*/
	
	double amount_useful_pix_oh1 = (xBaryMax - xBaryMin - hot_col_amount(hotcol_ohdu_1_min, hotcol_ohdu_1_max, sizehotcol_1_min))*(yBaryMax - yBaryMin)*10;//amount of pixels ohdu1 + ohdu2
	double amount_useful_pix_oh2 = (xBaryMax - xBaryMin - hot_col_amount(hotcol_ohdu_2_min, hotcol_ohdu_2_max, sizehotcol_2_min))*(yBaryMax - yBaryMin)*10;
	double amount_useful_pix = (xBaryMax - xBaryMin)*(yBaryMax - yBaryMin)*10;
	
	// Genero histogramas 	········
	int nbins_epix_overscan = 35;
	float lower_bin = -1.0;//0.66;//
	float upper_bin = 10.0;//1.8;//0.66;//
	float ePix_h=upper_bin;//{2.0};
    float ePix_l=lower_bin;//{-1.0};
	TH1D * epix_ohdu1_ovx_hist = new TH1D("epix ohdu1 OVX", "Overscan X", nbins_epix_overscan, lower_bin, upper_bin);// ohdu1
	TH1D * epix_ohdu2_ovx_hist = new TH1D("epix ohdu2 OVX", "Overscan X", nbins_epix_overscan, lower_bin, upper_bin);
	TH1D * epix_ohdu1_actar_hist = new TH1D("epix ohdu1", " Active Area", nbins_epix_overscan, lower_bin, upper_bin);
	TH1D * epix_ohdu2_actar_hist = new TH1D("epix ohdu2", " Active Area", nbins_epix_overscan, lower_bin, upper_bin);
	TH1D * epix_ohdu3_actar_hist = new TH1D("epix ohdu3", " Active Area", nbins_epix_overscan, lower_bin, upper_bin);
	TH1D * epix_ohdu4_actar_hist = new TH1D("epix ohdu4", " Active Area", nbins_epix_overscan, lower_bin, upper_bin);
	
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
	cout << " NCOL = " << NCOL 	<< endl;
	cout << "CCDNCOL = " << CCDNCOL << endl;
	int OVX_INIT_COL = atoi(CCDNCOL)/2.0;
	int OVX_COL = atoi(NCOL) - OVX_INIT_COL;
	cout << "OVX_COL = " << OVX_COL << endl;
	cout << "Inicio overscan = " << OVX_INIT_COL +1 << endl;

	TH2D * h2D_oh1 = new TH2D( "h2D_oh1", " Ohdu 1",/* X-dim */ atoi(NCOL), 0.0, atoi(NCOL), /* Y-dim */ atoi(NROW), 0.0, atoi(NROW));
	TH2D * h2D_oh2 = new TH2D( "h2D_oh2", " Ohdu 2",/* X-dim */ atoi(NCOL), 0.0, atoi(NCOL), /* Y-dim */ atoi(NROW), 0.0, atoi(NROW));


	// CalPixTree ---------------------------------------------------------------
	for(int j{0}; j < Entries_cptree; ++j){
		cptree->GetEntry(j);
		isnthotcol = true;
		//2D histogram to inspect events ·····
		switch(ohdu){
			case 1:
				//cout << "ohdu " << ohdu << endl;
				h2D_oh1->Fill(x,y,ePix);
				break;
			case 2:
				//cout << "ohdu " << ohdu << endl;
				h2D_oh2->Fill(x,y,ePix);
				break;

		}

		//Overscan x ·····
		if(x >= OVX_INIT_COL && ePix < ePix_h && ePix > ePix_l){//

			switch(ohdu){
			case 1:
				epix_ohdu1_ovx_hist->Fill(ePix);
				break;
			case 2:
				epix_ohdu2_ovx_hist->Fill(ePix);
				break;
				
			}

		}

		//Active Area - ····
		if(x < OVX_INIT_COL  && x > 7 && ePix < ePix_h && ePix > ePix_l){
			
			//Transforma isnthotcol = false si el pixel está en una columna dañada
			switch(ohdu){
			case 1:
				extracthc(hotcol_ohdu_1_min, hotcol_ohdu_1_max,sizehotcol_1_min, isnthotcol, x, x);
				break;
			case 2: 
				extracthc(hotcol_ohdu_2_min, hotcol_ohdu_2_max,sizehotcol_2_min, isnthotcol, x, x);
				break;
			}
	        
        	if (isnthotcol==true){//si el pixel no está en una columna brillante
				
				switch(ohdu){
				case 1: 
					epix_ohdu1_actar_hist->Fill(ePix);
					break;
				case 2:
					epix_ohdu2_actar_hist->Fill(ePix);
					break;
				case 3:
					epix_ohdu3_actar_hist->Fill(ePix);
					break;
				case 4:
					epix_ohdu4_actar_hist->Fill(ePix);
					break;
	    		}

    	    }

		}

	}
	
	//Normalizo histogramas ····
	//Probando ····
	/*
	cout << "Antes " << h2D_oh1->GetBinContent(1020,13) << endl;
	cout << "Minimum " << h2D_oh1->GetMinimum() << endl;
	h2D_oh1->AddBinContent(1020,13, -h2D_oh1->GetMinimum());
	cout << "Después " << h2D_oh1->GetBinContent(1020,13) << endl;
	*/

	//Aplico normalización min-max como en ds9
	normaliza_minmax(h2D_oh1);
	normaliza_minmax(h2D_oh2);
	//epix_ohdu1_ovx_hist->Scale(1./epix_ohdu1_ovx_hist->Integral(), "WIDTH");
	//epix_ohdu2_ovx_hist->Scale(1./epix_ohdu2_ovx_hist->Integral(), "WIDTH");
	
	//··················
	//Gráficos 		····
	//··················
	gStyle->SetOptStat(0);
	//gStyle->SetOptFit(0111);
	
	TCanvas * c2 = new TCanvas("c2", "2D Charge histogram", 1200, 800 );
	//c2 -> SetLogz();
	c2-> Divide(1,2);
	//ohdu 1 ····
	c2->cd(1);
	h2D_oh1->Draw("COL2Z");//
	//ohdu 2 ····
	c2->cd(2);
	h2D_oh2->Draw("COL2Z");//
	gPad->Update();
/*
    TCanvas * c1 = new TCanvas("c1", "Charge histogram", 1200, 800);
	gStyle->SetOptStat(0);//(11);
	c1 -> SetLogy();
	c1-> Divide(2,1);
	//ohdu 1 ····
	c1->cd(1);
	gPad->SetLogy();
    gPad->SetGrid();
	SetOhdu1Style(epix_ohdu1_actar_hist);
	SetStyle2(epix_ohdu1_ovx_hist);
    epix_ohdu1_actar_hist->SetTitle("ohdu 1");//
	epix_ohdu1_actar_hist->Draw();
	epix_ohdu1_ovx_hist->Draw("same");
	auto leg = new TLegend(0.55,0.75,0.90,0.90);;//(0.75,0.10,0.90,0.30);
	leg->SetTextSize(0.05);
    leg->AddEntry(epix_ohdu1_actar_hist,"Active Area", "l");
    leg->AddEntry(epix_ohdu1_ovx_hist,"Overscan x", "l");
    leg->Draw();
	//ohdu 2 ····
	c1->cd(2);
	gPad->SetLogy();
    gPad->SetGrid();
	SetOhdu1Style(epix_ohdu2_actar_hist);
	SetStyle2(epix_ohdu2_ovx_hist);
	epix_ohdu2_actar_hist->SetTitle("ohdu 2");
	epix_ohdu2_actar_hist->Draw();
	epix_ohdu2_ovx_hist->Draw("same");
	leg->Draw();
*/
	// Resultados a archivo ····
	//c1->SaveAs(Form("%s/%s", cfulloutdirname.Data(), OVX_out.Data()));
	//c2->SaveAs(Form("%s/%s", cfulloutdirname.Data(), h2_out.Data()));

}

//··················
// Rutine		····
//··················
//to run just one file comment the following

void histogramas_de_carga_app() {

  //Carpetas entrada

	TString fulldirname = TString(main_dirname) + sub_dirname;
	TString cmd = "ls \"" + fulldirname + "\" | sort -V";
	TString output = gSystem->GetFromPipe(cmd);
	TObjArray* files = output.Tokenize("\n");

	for (int i = 0; i < files->GetEntries(); ++i) {
	  TString fname = ((TObjString*)files->At(i))->GetString();
	  if (fname.EndsWith(".root") && fname.Contains("hits_corr_proc_run")) {
	    cout << fulldirname + fname << std::endl;
	    ProcessImage((fulldirname + "/" + fname).Data());
		}
	}

delete files;//esto es importante para liberar memoria porque TObjArray //to run over a full folder.

}

//··············································
//Function definition 						····
//··············································

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
    tree->SetBranchStatus("NROW",1);
    tree->SetBranchStatus("NCOL",1);
    tree->SetBranchStatus("CCDNCOL",1);
    tree->SetBranchStatus("NSAMP",1);

    tree->SetBranchAddress ("DATESTART",&DATESTART);//save branch content into variable
    tree->SetBranchAddress ("DATEEND",&DATEEND);
    tree->SetBranchAddress ("RUNID",&RUNID_head);
    tree->SetBranchAddress ("NROW",&NROW);
    tree->SetBranchAddress ("NCOL",&NCOL);
    tree->SetBranchAddress ("CCDNCOL",&CCDNCOL);
    tree->SetBranchAddress ("NSAMP",&NSAMP);
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
    return (stat(filename.c_str(), &buffer));
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


//https://root.cern.ch/root/htmldoc/guides/users-guide/FittingHistograms.html

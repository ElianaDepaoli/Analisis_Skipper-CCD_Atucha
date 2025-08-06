/* Autora: E.Depaoli Fecha inicial: 09/07/2025
Última modificación: 
06/08/25 Corregí regiones (area activa y OVX), corregí factores de normalización. 
05/08/2025. Arreglé bug en todos los switch.
Código para calcular el ruido de lectura y el SEE de las imágenes de Atucha utilizando los catálogos hits_corr_proc_run*.root
de una carpeta. Procesa los archivos que están en dirname+subdirname en orden numérico según el campo "img" del nombre completo
del catálogo.

Importante: 
Correr con la siguiente linea en consola para que no imprima en pantalla cada gráfico:
root -b -q RN_SER_Varios_calculos_app.C

Si no se usa la linea anterior, se bloquea la pantalla por la emergencia de gráficos.

Para usarlo
ofilename = nombre de un archivo de salida *.txt;
MODO CARPETA COMPLETA
main_dirname = ruta completa a la carpeta donde están los catálogos a leer; 
sub_dirname = nombre de la carpeta que contiene los catálogos;
Comentar las lineas con comentario://to run just one file using filename
Descomentar las lineas con comentario: //to run over a full folder.
MODO ÚNICO ARCHIVO
Descomentar las lineas con comentario://to run just one file using filename
Comentar las lineas con comentario: //to run over a full folder.

Lo que produce:
-> Un archivo *.txt conteniendo la tasa de eventos de 1 electrón en [e/pxl/day] sobre overscan x, overscan y y sobre el
área activa, el ruido de lectura en las extensiones 1 y 2. Almacena donde está el script.
-> Histogramas con su ajustes en las extensiones 1 y 2, sobre overscan x y overscan y. Almacena en carpeta que crea el 
propio script con el nombre de la carpeta procesada.

Detalle:

-> Histograma de carga sin clusterizar (ePix en calPixTree) restringida a los dos overscans, X e Y, para cada extensión. 
No se contabiliza carga en píxeles de columnas brillantes.
Ajusta cada histograma con la convolución entre Poisson y Gauss. 
Parámetros que se obtienen de este ajuste: el ruido de lectura, SER*tiempo. Cuadrantes 1 y 2.
-> Calcula cluster_1e como la suma del número de eventos que ocupan 1 pxl y tienen 1 electrón (n == 1.0 && e == 1.0) 
en todas las extensiones. Calcula cluster_almenos_1e +=n 
SEE = cluster_1e/N_tot_pxls/t_exposición_imagen
Occupancy = cluster_almenos_1e/N_tot_pxls/t_exposición_imagen

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
	
// Global Variables ---------------------------------------------------------------------
//string filename{"/home/eliana/Documentos/Atucha/images/hits_corr_proc_run_40_01Feb2025__EXP1_NSAMP300_VSUB70_img99.root"};
//string filename{"/media/eliana/0e0c1913-ca43-4252-b2bb-783703f57ea5/proc/run_41/hits/hits_corr_proc_run_41_19Feb2025__EXP1_NSAMP300_VSUB70_img633.root"};
string filename{"/home/eliana/Documentos/Atucha/images/run_46_29Jul2025/proc/hits/hits_corr_proc_run_46_29Jul2025__EXP1_NSAMP300_VSUB70_NROW100_img7.root"};
//const char* main_dirname = "/media/eliana/0e0c1913-ca43-4252-b2bb-783703f57ea5/proc/";
const char* main_dirname = "/home/eliana/Documentos/Atucha/images/";
const char* sub_dirname = "run_48_01Ago2025_test_vl-0_5/proc/hits/"; //"run_47_31Jul2025_test_vl-1_0/proc/hits/";// "run_46_29Jul2025/proc/hits/";//"run_39/hits/subset/";
//Archivo de salida
string ofilename = "RN_SER_Atucha_r48.txt";//;//"PROBANDO_r41_img633.txt";//"RN_SER_Atucha_r39_subset.txt";

//Carpetas salida
//const char* cfulloutdirname = gSystem->pwd();
//Creo una carpeta con igual nombre que el run que proceso
TString run_folder = TString(sub_dirname);
run_folder = run_folder(0, run_folder.First('/'));
// Crear carpeta nueva dentro del directorio actual
TString cfulloutdirname = TString(gSystem->pwd()) + "/" + run_folder;
gSystem->mkdir(cfulloutdirname);  // kTRUE crea también si ya existe


double Nmax{};
//Estimadores ·············································································
//Overscan Y ····
double RN_1_ovy{}; double err_RN_1_ovy{};
double RN_2_ovy{}; double err_RN_2_ovy{};
double lamb_1_ovy{}; double err_lamb_1_ovy{};
double lamb_2_ovy{}; double err_lamb_2_ovy{};
//Overscan X·····
double RN_1_ovx{}; double err_RN_1_ovx{};
double RN_2_ovx{}; double err_RN_2_ovx{};
double lamb_1_ovx{}; double err_lamb_1_ovx{};
double lamb_2_ovx{}; double err_lamb_2_ovx{};

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
Char_t CCDNCOL[256];
Char_t CCDNROW[256];

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

float factor = 1.0;//0.5;//corrección por lectura entre dos imágenes sucesivas

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
//void RN_SER_Varios_calculos_app() {//to run just one file using filename
void ProcessImage(const char* filename){//to run over a full folder
	//SetGlobalStyle();
	//Open root catalogue
	cout << filename << endl;
	TFile * file = TFile::Open(filename,"READ");//to run over a full folder
	//TFile * file = TFile::Open(filename.c_str(),"READ");// //to run just one file using filename
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
	cout << "NCOL = " << NCOL 	<< endl;
	cout << "CCDNCOL = " << CCDNCOL << endl;
	cout << "CCDNROW = " << atoi(&CCDNROW [0])<< endl;
	cout << "NROW = " << atoi(&NROW[0]) << endl;
    cout << "NSAMP = " << atoi(&NSAMP[0]) << endl; //<< "···· NCOL = " << atoi(&NCOL[0]) << "···· NROW = " << atoi(&NROW[0]) << endl;

	int OVX_INIT_COL = atoi(CCDNCOL)/2.0; //cantidad de columnas del area activa
	int OVX_COL = atoi(NCOL) - OVX_INIT_COL;
	cout << "OVX_COL = " << OVX_COL << endl;//cantidad de columnas en el overscan
	cout << "Inicio overscan = " << OVX_INIT_COL +1 << endl;

	// remove events partially inside the active area
	//Borders of the active area
	int xBaryMin=3; // 
	int xBaryMax= OVX_INIT_COL-50;//3072-50;//305; // 5 on the left from the overscan
	int yBaryMin=3; 
	int yBaryMax=atoi(&NROW[0])-5;//95;//545;

	cout << "xBaryMax = " << xBaryMax << endl;
	cout << "yBaryMax = " << yBaryMax << endl;

	//double amount_useful_pix_oh1 = (xBaryMax - xBaryMin - hot_col_amount(hotcol_ohdu_1_min, hotcol_ohdu_1_max, sizehotcol_1_min))*(yBaryMax - yBaryMin)*10;//amount of pixels ohdu1 + ohdu2
	//double amount_useful_pix_oh2 = (xBaryMax - xBaryMin - hot_col_amount(hotcol_ohdu_2_min, hotcol_ohdu_2_max, sizehotcol_2_min))*(yBaryMax - yBaryMin)*10;
	//double amount_useful_pix = (xBaryMax - xBaryMin)*(yBaryMax - yBaryMin)*10;
	
	//Leo menos filas que las reales y sin binneado en x ················································
	double amount_useful_pix_oh1 = (xBaryMax - xBaryMin - hot_col_amount(hotcol_ohdu_1_min, hotcol_ohdu_1_max, sizehotcol_1_min))*(yBaryMax - yBaryMin);//atoi(&NROW[0])
	double amount_useful_pix_oh2 = (xBaryMax - xBaryMin - hot_col_amount(hotcol_ohdu_2_min, hotcol_ohdu_2_max, sizehotcol_2_min))*(yBaryMax - yBaryMin);
	double amount_useful_pix = (xBaryMax - xBaryMin)*(yBaryMax - yBaryMin);

	// Genero histogramas 	········
	int nbins_epix_overscan = 100;
	float lower_bin = -1.0;//0.66;//
	float upper_bin = 1.8;//0.66;//
	TH1D * epix_ohdu1_ovx_hist = new TH1D("epix ohdu1 OVX", "Overscan X", nbins_epix_overscan, lower_bin, upper_bin);// ohdu1
	TH1D * epix_ohdu2_ovx_hist = new TH1D("epix ohdu2 OVX", "Overscan X", nbins_epix_overscan, lower_bin, upper_bin);
	TH1D * epix_ohdu1_ovy_hist = new TH1D("epix ohdu1 OVY", "Overscan Y", 50, -0.5, 1.5);
	TH1D * epix_ohdu2_ovy_hist = new TH1D("epix ohdu2 OVY", "Overscan Y", 150, -0.5, 1.5);
		
	//SEE. Viejo. ····
    double clusters_1e_oh1 {0};double clusters_1e_oh2 {0};double clusters_1e_oh3 {0};double clusters_1e_oh4 {0};//SEE
    Double_t see_oh1{0};Double_t see_oh2{0}; Double_t see_oh3{0}; Double_t see_oh4{0};
    Double_t occ_ohd1{0}; Double_t occ_ohd2{0}; Double_t occ_ohd3{0}; Double_t occ_ohd4{0};
	// hitSumm Tree ----------------------------------------------------------------
    cout << "Cantidad de entradas en hitSumm = " << Entries_hStree << endl;

    for (Long64_t entry{0}; entry < Entries_hStree ; ++entry)
    {
        hStree -> GetEntry(entry);
    	isnthotcol=true;
      
    	if (xBary<xBaryMax && xBary>xBaryMin && yBary<yBaryMax && yBary>yBaryMin){//events inside edges of CCD
    		/*if (ohdu_hs == 1 || ohdu_hs == 2){//events in quadrants with single electron counting
    			
    			if (ohdu_hs == 1){//events outside hot columns
        		    extracthc(hotcol_ohdu_1_min, hotcol_ohdu_1_max,sizehotcol_1_min, isnthotcol, xMin, xMax);
		        } else if (ohdu_hs == 2){extracthc(hotcol_ohdu_2_min, hotcol_ohdu_2_max,sizehotcol_2_min, isnthotcol, xMin, xMax);}

		 	} */

		 	//Transforma isnthotcol = false si el pixel está en una columna dañada. //Extraer columnas brillantes no reduce el valor final de clusters_1e porque estos no contienen eventos con n > 1. Pero debería reducir la ocupancia ¡y lo hace!
			switch(ohdu_hs){

				case 1:
					extracthc(hotcol_ohdu_1_min, hotcol_ohdu_1_max,sizehotcol_1_min, isnthotcol, xMin, xMax);
					break;
				case 2: 
					extracthc(hotcol_ohdu_2_min, hotcol_ohdu_2_max,sizehotcol_2_min, isnthotcol, xMin, xMax);
					break;

			}      

		 	//Active Area ····
			if (isnthotcol==true){//si el evento no contiene ni toca una columna brillante
		        if(yMax <= atoi(NROW) && xMax <= OVX_INIT_COL && xMin >= 8){// filas activas, antes del overscan, después del prescan
		        	if (n == 1.0 && e == 1.0 ){//n: pixels occupied, e: electrons in cluster
			        	switch(ohdu_hs){
							case 1:
								++clusters_1e_oh1;
								break;
							case 2:
								++clusters_1e_oh2;
								break;
							case 3:
								++clusters_1e_oh3;
								break;
							case 4:
								++clusters_1e_oh4;
								break;
							}
					}

					switch(ohdu_hs){
						case 1:
							occ_ohd1+=n;
							break;
						case 2:
							occ_ohd2+=n;
							break;
						case 3:
							occ_ohd3+=n;
							break;
						case 4:
							occ_ohd4+=n;
							break;
					}

				}
			}
	    }

    //··················
	//Cálculos 		····
	//··················

    //Luego de pasar por el último evento, ocupancia y DC 
        if(entry > (Entries_hStree-2)){
            image_duration_day = imagedurationindays(&DATESTART[0],&DATEEND[0], factor);
        }
    }

    cout << "-----------------------------------------------------------------------------------" << endl;

    cout << "amount of useful pixels in ohdu 1 = " << amount_useful_pix_oh1 << endl;
    cout << "amount of useful pixels in ohdu 2 = " << amount_useful_pix_oh2 << endl;
    cout << "amount of useful pixels in ohdu 3 = " << amount_useful_pix << endl;
            
    //double readout_time_per_pxl_per_sample = 2*image_duration_day*86400/(1024*(6144/20+12.8)*atoi(&NSAMP[0])/2);//cálculo errado porque no contiene overscan en y
    //si factor = 2, multiplico por 2 para cancelar la división por 2 en imagedurationdays (que viene de que hay limpieza entre lecturas sucesivas)
    double readout_time_per_pxl_per_sample = image_duration_day/atoi(&NCOL[0])/atoi(&NROW[0])/atoi(&NSAMP[0]);// 51.6667 microsec
    //exposure_active_area = readout_time_per_pxl_per_sample*atoi(&NSAMP[0])*atoi(&NCOL[0])*atoi(&CCDNROW[0])/2/2;
    exposure_active_area = readout_time_per_pxl_per_sample*atoi(&NSAMP[0])*atoi(&NCOL[0])*atoi(&NROW[0])/2;
    double t_media_ov_px = atoi(&NSAMP[0])*SR*readout_time_per_pxl_per_sample;//;//51.6667 microsec /86400
    double t_expo_overscan_y = exposure_active_area/atoi(&NROW[0])+exposure_active_area;//

    //double t_expo_overscan_y = atoi(&NSAMP[0])*readout_time_per_pxl_per_sample*(filas_aa - 1)*atoi(&NCOL[0])/2;//0.017329 day
    //double t_expo_overscan_y = atoi(&NSAMP[0])*readout_time_per_pxl_per_sample*(SR+1)/2 + exposure_active_area;//0.017396 day
    //double t_expo_overscan_y = atoi(&NSAMP[0])*readout_time_per_pxl_per_sample*(filas_aa - 1)/2 + exposure_active_area;//0.017417 day //5.4153e-05 day + image_duration_day/2;
    //double t_expo_overscan_y = atoi(&NSAMP[0])*readout_time_per_pxl_per_sample*atoi(&NCOL[0])*(atoi(&NCOL[0])+1)/2 + exposure_active_area;//0.028249 day
    
    cout << "readout_time_per_pxl_per_sample = " << 1e6*readout_time_per_pxl_per_sample*24*60*60 << " microsec" << endl;
    cout << "exposure_active_area = " << exposure_active_area << " day " << endl;
    cout << "imageduration = " << setprecision(4)<< image_duration_day*24*60 << " minutes" << '\n';
    cout << "tiempo de exposición promedio de un pixel x-overscan  = " << std::setprecision(5) << t_media_ov_px << " day " << endl; //" sec"  
    cout << "tiempo de exposición promedio de un pixel y-overscan  = " << std::setprecision(5) << t_expo_overscan_y <<  " day " << endl; //" sec"

    cout << "-----------------------------------------------------------------------------------" << endl;
    cout << " ················· Active area  - Clusters ··········· " << endl;
    cout << "ohdu 1 see normalization = " << amount_useful_pix_oh1*image_duration_day <<  " day*pxl " << endl;
    cout << "ohdu 2 see normalization = " << amount_useful_pix_oh2*image_duration_day <<  " day*pxl " << endl;
    cout << "ohdu 3 & 4 see normalization = " << amount_useful_pix*image_duration_day <<  " day*pxl " << endl;

    occ_ohd1 = static_cast<float>(occ_ohd1)/ static_cast<float>(amount_useful_pix_oh1);// / image_duration_day; // / (image_duration_day*24*60*60); ///;//imagedurationindays(&DATESTART[0],&DATEEND[0],factor)/
    occ_ohd2 = static_cast<float>(occ_ohd2)/ static_cast<float>(amount_useful_pix_oh2);//  / image_duration_day; // / (image_duration_day*24*60*60);//occ_ohd2;///amount_useful_pix_oh2;
    occ_ohd3 = static_cast<float>(occ_ohd3)/ static_cast<float>(amount_useful_pix);//  / image_duration_day; // / (image_duration_day*24*60*60);//occ_ohd3;///amount_useful_pix;
    occ_ohd4 = static_cast<float>(occ_ohd4)/ static_cast<float>(amount_useful_pix);//  / image_duration_day; // / (image_duration_day*24*60*60); //occ_ohd4;///amount_useful_pix;
    see_oh1 = clusters_1e_oh1 / amount_useful_pix_oh1;//  /image_duration_day; // 
    see_oh2 = clusters_1e_oh2 / amount_useful_pix_oh2;//  /image_duration_day; // 
    see_oh3 = clusters_1e_oh3 / amount_useful_pix;//  /image_duration_day; // 
    see_oh4 = clusters_1e_oh4 / amount_useful_pix;//  /image_duration_day; // 

    cout << "see_oh1 = " << see_oh1 << " ---- " << "occ_oh1 = " << occ_ohd1 <<  " e/pix " << endl;///day
    cout << "see_oh2 = " << see_oh2 << " ---- " << "occ_oh2 = " << occ_ohd2 <<  " e/pix " << endl;
    cout << "see_oh3 = " << see_oh3 << " ---- " << "occ_oh3 = " << occ_ohd3 <<  " e/pix " << endl;
    cout << "see_oh4 = " << see_oh4 << " ---- " << "occ_oh4 = " << occ_ohd4 <<  " e/pix " << endl;

    // Resultados a archivo ···············
	// Verifico si el archivo ya existe
    bool ofile_already_exists = fileExists(ofilename);// sí = 0, no = 0
    ofstream ofile(ofilename, std::ios::app);//creo el archivo
	// Si el archivo no existía, escribir el encabezado
    if (ofile_already_exists!=0) {
        ofile << "RUNID " << '\t' << "DATESTART"<< '\t' << "see_oh1 [e/pxl]" << '\t' << "occ_oh1 [e/pxl]" << '\t' << "see_oh2 [e/pxl]" << '\t' << "occ_oh2 [e/pxl]" << '\t' << "see_oh3 [e/pxl]" << '\t' << "occ_oh3 [e/pxl]" << '\t' << "see_oh4 [e/pxl]" << '\t' << "occ_oh4 [e/pxl]" << '\t' << "lamb_1_ovx [e/pxl]" << '\t' << "err_lamb_1_ovx [e/pxl] " << '\t' << "lamb_2_ovx [e/pxl]" << '\t' << "err_lamb_2_ovx [e/pxl] " << '\t' << "RN_1_ovx [e]" << '\t' << "err_RN_1_ovx [e]" << '\t' << "RN_2_ovx [e] " << '\t' << "err_RN_2_ovx [e] " << '\t' << "lamb_1_ovy [e/pxl]" << '\t' << "err_lamb_1_ovy [e/pxl] " << '\t' << "lamb_2_ovy [e/pxl]" << '\t' << "err_lamb_2_ovy [e/pxl] " << "image_duration [day] " << "exposure_active_area [day] " << "exposure OVX [day] "<< "exposure OVY [day] " << '\n';
    }

	//cout << ofile_already_exists << endl;
	ofile << RUNID_head << '\t' << &DATESTART[0] << '\t' << see_oh1 << '\t' << occ_ohd1 << '\t' << see_oh2 << '\t' << occ_ohd2 << '\t' << see_oh3 << '\t' << occ_ohd3 << '\t' << see_oh4 << '\t' << occ_ohd4 << '\t'; 
	    
	// CalPixTree ---------------------------------------------------------------
	for(int j{0}; j < Entries_cptree; ++j){
		cptree->GetEntry(j);
		//Overscan x ·····
		if(x > OVX_INIT_COL && ePix < 2 && ePix > -1){//
			switch(ohdu){
			case 1:
				epix_ohdu1_ovx_hist->Fill(ePix);
				break;

			case 2:
				epix_ohdu2_ovx_hist->Fill(ePix);
				break;
				
			}
		}

		//SER ····
		isnthotcol = true;
		//Overscan y ····
		if(x <= OVX_INIT_COL  && x > 7 && y > atoi(&CCDNROW [0]) && ePix < 2 && ePix > -1){
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
					epix_ohdu1_ovy_hist->Fill(ePix);
					break;
				case 2:
					epix_ohdu2_ovy_hist->Fill(ePix);
					break;
				
	    		}

    	    }
		}

	}

	cout << '\n' << endl;
	cout << "·············calPixTree ············" << endl;
    
	//Normalizo histogramas ····
	epix_ohdu1_ovy_hist->Scale(1./epix_ohdu1_ovy_hist->Integral(),"WIDTH");//histogram describe a probability density function. 
	epix_ohdu1_ovx_hist->Scale(1./epix_ohdu1_ovx_hist->Integral(), "WIDTH");
	epix_ohdu2_ovy_hist->Scale(1./epix_ohdu2_ovy_hist->Integral(),"WIDTH");
	epix_ohdu2_ovx_hist->Scale(1./epix_ohdu2_ovx_hist->Integral(), "WIDTH");
	
	//··················
	//Ajustes 		····
	//··················
	// SER + Ruido de lectura: Poisson convolucionada con gaussiana ·······················
	TF1 *pgc = new TF1("poisson_gauss_conv", poisson_gauss_conv, -1, 2,4);//(,, x_lim_inf, x_lim_sup, #parameters)
	Nmax= 1;
	pgc->SetParameters(0.02,0.03, 1.0, 0.24); //lambda, mu, normalización, sigma. Inicializo.
	pgc->SetParNames("lambda", "mu", "norm", "RN");
	pgc->FixParameter(1,0);
	pgc->SetLineColor(kMagenta-4);
	
	//OVERSCAN ··············
	cout << '\n' << endl;
	cout << " --------- OVERSCAN X --------- " << endl;
	// --------- Convolution Poisson-Gauss---------;

	cout << " OHDU 1 " << endl;
	//“Q” Quiet mode (minimum printing); “S” The result of the fit is returned in the TFitResultPtr; "R" to restrict the fit to the range specified in the TF1 constructor
	auto resultovx = epix_ohdu1_ovx_hist->Fit("poisson_gauss_conv", "R+S");//+Q
	double pvalueovx = ROOT::Math::chisquared_cdf_c(resultovx->Chi2(), resultovx->Ndf()); //
	lamb_1_ovx = resultovx -> Parameter(0); err_lamb_1_ovx = resultovx -> ParError(0); RN_1_ovx = resultovx -> Parameter(3); err_RN_1_ovx = resultovx -> ParError(3);

	cout << "pvalue = " << pvalueovx << endl;
	//cout << "SER = " << std::setprecision(9) << lamb_1_ovx / t_media_ov_px << " +- " << err_lamb_1_ovx / t_media_ov_px << " e/pxl/day" << endl;
	
	cout << '\n' << endl;
	cout << " OHDU 2 " << endl;
	resultovx = epix_ohdu2_ovx_hist->Fit("poisson_gauss_conv", "R+S");//+Q
	pvalueovx = ROOT::Math::chisquared_cdf_c(resultovx->Chi2(), resultovx->Ndf());
	lamb_2_ovx = resultovx -> Parameter(0); err_lamb_2_ovx = resultovx -> ParError(0); RN_2_ovx = resultovx -> Parameter(3);	err_RN_2_ovx = resultovx -> ParError(3);
	
	cout << "pvalue = " << pvalueovx << endl;
	//cout << "SER = " << std::setprecision(9) << lamb_2_ovx/ t_media_ov_px  << " +- " << err_lamb_2_ovx/ t_media_ov_px  << " e/pxl/day" << endl;
	cout << '\n' << endl;
	// Resultados a archivo ····
	ofile << lamb_1_ovx << '\t' << err_lamb_1_ovx << '\t'  << lamb_2_ovx << '\t' << err_lamb_2_ovx << '\t' << RN_1_ovx << '\t' << err_RN_1_ovx << '\t' << RN_2_ovx << '\t' << err_RN_2_ovx << '\t';
	
	//Grafico OVERSCAN X ···
	//Poisson convolucionada con Gauss ·····························
	//gStyle->SetOptFit(0111);
    TCanvas * c1 = new TCanvas("c1", "OVX", 1200, 800);
	gStyle->SetOptStat(11);
	gStyle->SetOptFit(0111);
	c1 -> SetLogy();
	c1-> Divide(2,1);
	//ohdu 1 ····
	c1->cd(1);
	gPad->SetLogy();
    gPad->SetGrid();
	SetOhdu1Style(epix_ohdu1_ovx_hist);
	epix_ohdu1_ovx_hist->Draw();
	//ohdu 2 ····
	c1->cd(2);
	gPad->SetLogy();
    gPad->SetGrid();
	SetOhdu1Style(epix_ohdu2_ovx_hist);
	epix_ohdu2_ovx_hist->Draw();

	TObjArray *tokens = TString(filename).Tokenize("_, .");
	TString OVX_out = ((TObjString*)tokens->At(6))->GetString()+ "_" + ((TObjString*)tokens->At(18))->GetString() + "_OVX" + ".png";
	TString OVY_out = ((TObjString*)tokens->At(4))->GetString() + "_" + ((TObjString*)tokens->At(5))->GetString()+ "_" + ((TObjString*)tokens->At(12))->GetString() + "_OVY" + ".png";
	//std::cout << OVX_out << std::endl;
	delete tokens;
	
	//c1->SaveAs(OVX_out.Data());

	c1->SaveAs(Form("%s/%s", cfulloutdirname.Data(), OVX_out.Data()));
/*
	cout << " --------- OVERSCAN Y --------- " << endl;
	cout << " OHDU 1 " << endl;
	
	//Grafico y ajusto OVERSCAN Y ···
	//Poisson convolucionada con Gauss ·····························
	//gStyle->SetOptFit(0111);
    TCanvas * c2 = new TCanvas("c2", "OVY", 1200, 800);
	gStyle->SetOptStat(11);
	gStyle->SetOptFit(0111);
	c2-> Divide(2,1);
	c2->cd(1);
	//ohdu 1 ····
	auto resultovy = epix_ohdu1_ovy_hist->Fit("poisson_gauss_conv", "R+S");
	double pvalueovy = ROOT::Math::chisquared_cdf_c(resultovy->Chi2(), resultovy->Ndf());
	lamb_1_ovy = resultovy -> Parameter(0); err_lamb_1_ovy = resultovy -> ParError(0); RN_1_ovy = resultovy -> Parameter(3); err_RN_1_ovy = resultovy -> ParError(3);

	cout << "pvalue = " << pvalueovy << endl;
	//cout << "SER = " << std::setprecision(9) << lamb_1_ovy / t_expo_overscan_y << " +- " << err_lamb_1_ovy / t_expo_overscan_y << " e/pxl/day" << endl;
	gPad->SetLogy();
    gPad->SetGrid();
	SetOhdu1Style(epix_ohdu1_ovy_hist);
	epix_ohdu1_ovy_hist->Draw();
		
	//ohdu 2 ····
	c2->cd(2);
	gPad->SetLogy();
    gPad->SetGrid();
	
	cout << " OHDU 2 " << endl;
	resultovy = epix_ohdu2_ovy_hist->Fit("poisson_gauss_conv", "R+S");
	pvalueovy = ROOT::Math::chisquared_cdf_c(resultovy->Chi2(), resultovy->Ndf());
	lamb_2_ovy = resultovy -> Parameter(0); err_lamb_2_ovy = resultovy -> ParError(0); RN_2_ovy = resultovy -> Parameter(3); err_RN_2_ovy = resultovy -> ParError(3);
	
	cout << "pvalue = " << pvalueovy << endl;
	//cout << "SER = " << std::setprecision(9) << lamb_2_ovy / t_expo_overscan_y << " +- " << err_lamb_2_ovy / t_expo_overscan_y << " e/pxl/day" << endl;
	SetOhdu1Style(epix_ohdu2_ovy_hist);
	epix_ohdu2_ovy_hist->Draw();*/
	//c2->SaveAs(Form("%s/%s", cfulloutdirname.Data(), OVY_out.Data()));
	//c2->SaveAs(OVY_out.Data());
	// Resultados a archivo ····
	ofile  << lamb_1_ovy << '\t' << err_lamb_1_ovy << '\t' << lamb_2_ovy << '\t' << err_lamb_2_ovy  << '\t' << image_duration_day << '\t' << exposure_active_area << '\t' << t_media_ov_px << '\t' << t_expo_overscan_y << '\n';
	
	ofile.close();//cierro salida
	//file->close();//cierro *.root que contiene los trees usados para calcular
	
}

//··················
// Rutine		····
//··················
//to run over a full folder.

void RN_SER_Varios_calculos_app() {

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

delete files;//esto es importante para liberar memoria porque TObjArray 

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
    tree->SetBranchStatus("CCDNROW",1);
    tree->SetBranchStatus("NSAMP",1);

    tree->SetBranchAddress ("DATESTART",&DATESTART);//save branch content into variable
    tree->SetBranchAddress ("DATEEND",&DATEEND);
    tree->SetBranchAddress ("RUNID",&RUNID_head);
    tree->SetBranchAddress ("NROW",&NROW);
    tree->SetBranchAddress ("NCOL",&NCOL);
    tree->SetBranchAddress ("CCDNCOL",&CCDNCOL);
    tree->SetBranchAddress ("CCDNROW",&CCDNROW);
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

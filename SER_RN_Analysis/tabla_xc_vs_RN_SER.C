/*Autores: E.Depaoli, D.Rodrigues Fecha: 29/04/25
Este código sirve para generar una tabla de la función  xc = f(RN,SER)
[RN] = e- , [SER] = e-/pix , [xc] = e-
xc es el límite superior de integración de una distribución N(0,RN) y el límite inferior de integración de una distribución N(1,RN) tal que  N(0,RN) = SER * N(1,RN)
Esta tabla seré utilizada luego para obtener el valor xc que corresponde a un dado RN (img), SER (img)
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

// Function declaration ·················
double SER_2(double xc, double sigma, double bi, double bf);
double SER(double xc, double sigma); // Umbral epix, ruido de lectura
vector <double> make_vector(double beg, double end, double step);
bool fileExists(const std::string& filename);

//Configuracion global para los gráficos --------------------------
void SetGlobalStyle(){
	TStyle *st1 = new TStyle("st1","my style");
	st1->SetTitleSize(0.05,"XYZ");
	st1->SetLabelSize(0.05, "XYZ");
	st1->SetPadBorderMode(0);
	st1->SetLegendTextSize(0.03);
	st1->SetLegendBorderSize(-1);
	st1->SetGridColor(0);
	//st1->SetLegendFillColor(0);
	gStyle->SetLegendFillColor(0);
	gROOT->SetStyle("st1");//
	gROOT->ForceStyle();//esto no funciona: el estilo se aplica a todos los objetos creados después
	st1->cd ();//This is now the current style gStyle
}

//Notación ···························
//double ROOT::Math::normal_cdf (double	x, double sigma = 1, double x0 = 0) x0 es el centroide, x es el límite de integración
// _cdf Cumulative distribution function 

//Global variables ·····················
//Archivo de salida
TString cfulloutdirname = TString(gSystem->pwd()) + "/";
string ofilename = "tabla_RN_SER_xcePix_Atucha.txt";

//los usados en ruido_de_lectura_atucha.C
double Bi{-0.5};//{-10.5};//
double Bf{1.5};//{0.0};//

//Rutine ······························
void tabla_xc_vs_RN_SER(){

	SetGlobalStyle();

	auto xC_vec = make_vector(0.540, 0.802, 0.002);//(0.5, 0.8, 0.1);
	auto sigma_vec = make_vector(0.201, 0.242, 0.002);

	//cout << "Cantidad de valores xC_epix = " << xC_vec.size() << endl;
	//cout << "Cantidad de ruidos de lectura = " << sigma_vec.size() << endl;
	
	double lambda;
	
	// Verifico si el archivo ya existe
    bool ofile_already_exists = fileExists(ofilename);// sí = 0, no = 0
    ofstream ofile(ofilename, std::ios::app);//creo el archivo
	// Si el archivo no existía, escribir el encabezado
    if (ofile_already_exists!=0) {
        ofile << "Sigma Gauss" << '\t' << "Esperanza Poisson" << '\t' << "Límite de integración ePix " << '\n';
        ofile << " [e]" << '\t' << " [e/pxl]" << '\t' << " [e] " << '\n';
    }


	for(int j{0}; j < sigma_vec.size(); ++j){
	
		for(int k{0}; k < xC_vec.size(); ++k){
			//lambda = SER_2(xC_vec.at(k), sigma_vec.at(j), Bi, Bf);
			lambda = SER(xC_vec.at(k), sigma_vec.at(j));
			//cout << "RN =  " << sigma_vec.at(j) << endl;
			//cout << "Esperanza Poisson = " << lambda << endl;
			//cout << "xc = " << xC_vec.at(k) << endl;
			
			ofile << sigma_vec.at(j) << '\t' << lambda << '\t' << xC_vec.at(k) << '\n';
		}

	}
	
	ofile.close();//cierro salida

}

// Function definition ·····················
//El comportamiento de esta función fue evaluado y produce resultados razonables ·····························
double SER_2(double xc, double sigma, double bi, double bf){// Umbral epix, ruido de lectura, primer intervalo en el pico de 0 e, último intervalo en el pico de 1 e
	double num = ROOT::Math::normal_cdf(bf/sigma, 1,0)- ROOT::Math::normal_cdf(xc/sigma, 1,0);
	double denom = ROOT::Math::normal_cdf((xc-1)/sigma, 1,0)- ROOT::Math::normal_cdf((bi-1)/sigma, 1,0);
	double lambda = num/denom;
	return lambda;
}

double SER(double xc, double sigma){// Umbral epix, ruido de lectura
	//double num = ROOT::Math::normal_cdf_c(xc/sigma, 1,0);//1 - ROOT::Math::normal_cdf(xc/sigma, 1,0);//
	//double denom = ROOT::Math::normal_cdf((xc-1)/sigma, 1,0);
	double lambda = ROOT::Math::normal_cdf(-xc/sigma, 1,0)/ROOT::Math::normal_cdf((xc-1)/sigma, 1,0);
	return lambda;
}

vector <double> make_vector(double beg, double end, double step) 
{
    vector<double> vec;
    vec.reserve((end - beg)/step + 1);//
    
    while (beg < end + step) {
        vec.push_back(beg);
        beg += step;
    }

    return vec;
}

// Función para verificar si un archivo existe ·······························································
bool fileExists(const std::string& filename) {
    struct stat buffer;
    return (stat(filename.c_str(), &buffer));
}

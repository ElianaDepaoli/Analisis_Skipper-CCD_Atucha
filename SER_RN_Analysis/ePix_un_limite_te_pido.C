/* Autora: Eliana Depaoli Fecha: 17/07/25
Este código sirve para determinar el límite de integración de la variable ePix usando la corriente oscura y el 
ruido de lectura reales de un run, y el valor estimado de dicho límite para el cual se igualan las probabilidades de determinar
erróneamente que un pixel está prendido o apagado.
Para su ejecución, hay que dar la ruta completa a dos archivos que se usarán para el procesamiento:
-> RN_SER_Atucha_rXX.txt 
este contiene la corriente oscura (mu) y el ruido de lectura (RN) de cada imagen (runID) de una corrida (run), y también la exposición
-> tabla_RN_SER_xcePix_Atucha.txt
este contiene el valor calculado de mu para un conjunto de valores hipotéticos para RN y xC

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
using std::cout; using std::vector; using std::string; using std :: copy;

//Function declaration ·················································
void archivo_a_vectores(vector<double>& x_vec, vector<double>& y_vec, const char*filename, const int colx, const int coly);
void archivo_a_string_vectores(vector<string>& x_vec, vector<double>& y_vec, const char*filename, const int colx, const int coly);
double date_to_seconds (string DATESTART);
void archivo_a_vectores(vector<double>& v1_vec, vector<double>& v2_vec, vector<double>& v3_vec, const char*filename);
bool fileExists(const std::string& filename);

//·······································

double RN_on_table{0.0};
double mu_on_table{0.0};
double xc_epix_oh1_fi{0.0};


double xc_epix_val(vector <double> v1_vec, vector <double> v2_vec, vector <double> v3_vec, double RN_measured, double mu_measured ){
        
    int jj{0};
        for(int j{0}; j < v1_vec.size();++j){
            if(v1_vec.at(j) < RN_measured){
                ++jj;
            }else{
                    RN_on_table = v1_vec.at(jj);//to compare predicted value (in fixed table) vs measured value (in run's table)
                    if(v2_vec.at(jj) > mu_measured){
                        ++jj;
                    }else{
                        mu_on_table = v2_vec.at(jj);//to compare predicted value (in fixed table) vs measured value (in run's table)
                        break;
                    }
                }
        }
        cout << "jj" << jj << endl;

    return v3_vec.at(jj);

}

//Rutine ·······························································
const char* filename = "/home/eliana/Documentos/Atucha/Scripts/SEE_RN/RN_SER_Atucha_r41.txt";
const char* filename2 = "/home/eliana/Documentos/Atucha/Scripts/SEE_RN/tabla_RN_SER_xcePix_Atucha.txt";
string ofilename = "ePix_te_puse_un_limite_r41.txt";
// Mediciones ····
vector<vector<double>> RN_vec(2);
vector<vector<double>> err_RN_vec(2);
vector<vector<double>> lamb_ovx_vec(2);
vector<vector<double>> err_lamb_ovx_vec(2);
vector<vector<double>> lamb_ovy_vec(2);
vector<vector<double>> err_lamb_ovy_vec(2);
vector<vector<double>> see_vec(2);//vector de vectores
vector<string> date_start_vec;//fecha y horario de inicio de toma de la imagen
vector<double> RUNID_vec;
vector<double> timestamp_date_start_vec;//timestamp
vector<double> exposure_vec;//[day] exposure OVX
// Predicciones ····
vector<double> sigma_vec;
vector<double> mu_vec;
vector<double> xc_epix_vec;
// Valores resultantes ····
vector<double> xc_epix_oh1_vec;

//Salidas ····
TString run = "RUN_41";
TString cfulloutdirname = TString(gSystem->pwd());
TString label_Yaxis_units = "[e/pxl]";
TString label_occ = "Occupancy ";
TString label_see = "SEE ";
TString label_OVX = "OVX ";
TString label_OVY = "OVY ";
TString label_clusters = "Clusters ";
TString label_y_occ = label_occ + label_Yaxis_units;
TString label_y_ovy = label_see + label_OVY + label_Yaxis_units;
TString label_y_ovx = label_see + label_OVX + label_Yaxis_units;
TString label_y_clus = label_see + label_clusters + label_Yaxis_units;

void ePix_un_limite_te_pido(){

	//Mediciones ·····
    for(int i{0}; i < 2; ++i) archivo_a_vectores(RN_vec[i], see_vec[i], filename, 2*i+14, 2*i+2);
    for(int i{0}; i < 2; ++i) archivo_a_vectores(lamb_ovx_vec[i], lamb_ovy_vec[i], filename, 2*i+10, 2*i+18);
    archivo_a_vectores(RUNID_vec, exposure_vec, filename, 0, 24);
    // Ver vectores en pantalla ·····
    //for (int i{0}; i < lamb_ovy_vec[0].size(); ++i){cout << lamb_ovy_vec[0].at(i) << endl;}
    //for (int i{0}; i < lamb_ovx_vec[1].size(); ++i){cout << lamb_ovx_vec[1].at(i) << endl;}
    //for (int i{0}; i < RN_vec[0].size(); ++i){cout << RN_vec[0].at(i) << endl;}
    //for (int i{0}; i < see_vec[0].size(); ++i){cout << see_vec[0].at(i) << endl;}
    //for (int i{0}; i < exposure_vec.size(); ++i){cout << exposure_vec.at(i) << endl;}
    //Predicciones ·····
    archivo_a_vectores(sigma_vec, mu_vec, xc_epix_vec, filename2);
    // Ver vectores en pantalla ·····
	//for (int i{0}; i < sigma_vec.size(); ++i){cout << xc_epix_vec.at(i) << endl;}
    
    // Salida ·········
    // Verifico si el archivo ya existe
    bool ofile_already_exists = fileExists(ofilename);// sí = 0, no = 0
    ofstream ofile(ofilename, std::ios::app);//creo el archivo
    // Si el archivo no existía, escribir el encabezado
    if (ofile_already_exists!=0) {
        ofile << "RunID" << '\t' << "ohdu" << '\t' << "ePix_UL" << '\t' << '\t' << "Measured_RN" << '\t' << "Tabulated_RN" << '\t' << "Measured_mu" << '\t' << "Tabulated_mu" << '\n';
    }

    


    //Busco el valor xc_epix que le corresponde a un cuadrante
    /*cout << "measured RN = " << RN_vec[0].at(0) << "------------" << "measured mu = " << lamb_ovx_vec[0].at(0) << endl; 
    xc_epix_oh1_fi = xc_epix_val(sigma_vec, mu_vec, xc_epix_vec, RN_vec[0].at(0), lamb_ovx_vec[0].at(0));
    cout << "RN_on_table = " << RN_on_table << endl;
    cout << "mu_detected = " << setprecision(6) << mu_on_table << endl;
    cout << "xc_epix_oh1_fi = " << setprecision(6) << xc_epix_oh1_fi << endl;
    */
    cout << RN_vec[0].size() << endl;
    cout << RN_vec.size() << endl;

    for(int i{0}; i < RN_vec[0].size() ; ++i){//
        for(int k{0}; k < RN_vec.size(); ++k){
            ofile << RUNID_vec.at(i) << '\t' << '\t' << k+1 << '\t' << '\t' << xc_epix_val(sigma_vec, mu_vec, xc_epix_vec, RN_vec[k].at(i), lamb_ovx_vec[k].at(i)) << '\t' << '\t' << RN_vec[k].at(i) << '\t' << '\t' << RN_on_table << '\t' << '\t' << lamb_ovx_vec[k].at(i) << '\t' << '\t' << mu_on_table << endl;

        }
    }
}


//Function definition ·················································
void archivo_a_vectores(vector<double>& x_vec, vector<double>& y_vec, const char*filename, const int colx, const int coly) {
    
    // Abrir el archivo
    std::ifstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Error: No se pudo abrir el archivo " << filename << endl;
        return;
    }else{cout << "······" << " Abri " << filename << "······" << endl;}

    // Leer e ignorar el encabezado
    std::string header;
    std::getline(file, header);

    // Leer el archivo línea por línea
    std::string line;
    while (std::getline(file, line, '\n')) {
        
        //cout << line << endl;
        // Ignorar líneas vacías o comentarios
        if (line.empty() || line[0] == '#') continue;

        // Reemplazar tabulaciones por espacios para facilitar la lectura
        std::replace(line.begin(), line.end(), '\t', ' ');

        // Usar un stringstream para dividir la línea en columnas
        std::stringstream ss(line);//esta clase opera sobre strings
        string value;
        int column = 0;

        // Leer cada valor de la línea
        while (ss >> value) {
            //cout << column << endl;
            if (column == colx) x_vec.push_back(stod(value)); // Guardar columna colx
            if (column == coly) y_vec.push_back(stod(value)); // Guardar columna coly
            column++;
        }
    }

    // Cerrar el archivo
    file.close();

    // Verificar que se hayan leído datos
    /*if (x_vec.empty() || y_vec.empty()) {
        std::cerr << "Error: No se encontraron datos válidos" << std::endl;
        return;
    }else{ cout << "Vectores llenos" << endl;}*/

}

// ·········
void archivo_a_vectores(vector<double>& v1_vec, vector<double>& v2_vec, vector<double>& v3_vec, const char*filename) {
    
    // Abrir el archivo
    std::ifstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Error: No se pudo abrir el archivo " << filename << endl;
        return;
    }else{cout << "······" << " Abri " << filename << "······" << endl;}

    // Leer e ignorar el encabezado
    std::string header;
    std::getline(file, header);
    std::getline(file, header);

    // Leer el archivo línea por línea
    std::string line;
    while (std::getline(file, line, '\n')) {
        
        //cout << line << endl;
        // Ignorar líneas vacías o comentarios
        if (line.empty() || line[0] == '#') continue;

        // Reemplazar tabulaciones por espacios para facilitar la lectura
        std::replace(line.begin(), line.end(), '\t', ' ');

        // Usar un stringstream para dividir la línea en columnas
        std::stringstream ss(line);//esta clase opera sobre strings
        string value;
        int column = 0;

        // Leer cada valor de la línea
        while (ss >> value) {
            //cout << column << endl;
            if (column == 0) v1_vec.push_back(stod(value)); // Guardar columna colx
            if (column == 1) v2_vec.push_back(stod(value)); // Guardar columna coly
            if (column == 2) v3_vec.push_back(stod(value));

            column++;
        }
    }

    // Cerrar el archivo
    file.close();

    // Verificar que se hayan leído datos
    if (v1_vec.empty() || v2_vec.empty() || v3_vec.empty()) {
        std::cerr << "Error: No se encontraron datos válidos" << std::endl;
        return;
    }else{ cout << "Vectores llenos" << endl;}

}




// ········
void archivo_a_string_vectores(vector<string>& x_vec, vector<double>& y_vec, const char*filename, const int colx, const int coly) {
    
    // Abrir el archivo
    std::ifstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Error: No se pudo abrir el archivo " << filename << endl;
        return;
    }else{cout << "······" << " Abri " << filename << "······" << endl;}

    // Leer y ignorar la primera línea (encabezado)
    std::string header;
    std::getline(file, header);

    // Leer el archivo línea por línea
    std::string line;
    while (std::getline(file, line, '\n')) {
        
        //cout << line << endl;
        // Ignorar líneas vacías o comentarios (si los hay)
        if (line.empty() || line[0] == '#') continue;

        // Reemplazar tabulaciones por espacios para facilitar la lectura
        std::replace(line.begin(), line.end(), '\t', ' ');

        // Usar un stringstream para dividir la línea en columnas
        std::stringstream ss(line);//esta clase opera sobre strings
        string value;
        int column = 0;

        // Leer cada valor de la línea
        while (ss >> value) {
            //cout << column << endl;
            if (column == colx) x_vec.push_back(value); // Guardar columna colx
            if (column == coly) y_vec.push_back(stod(value)); // Guardar columna coly
            column++;// supuestamente incrementa column en uno luego de ejecutar lo que sigue pero no es cierto
        }
    }

    // Cerrar el archivo
    file.close();

    // Verificar que se hayan leído datos
    /*if (x_vec.empty() || y_vec.empty()) {
        std::cerr << "Error: No se encontraron datos válidos en las columnas 4 y 9." << std::endl;
        return;
    }else{ cout << "Vectores llenos" << endl;}*/
}

// Función para verificar si un archivo existe ·······························································
bool fileExists(const std::string& filename) {
    struct stat buffer;
    return (stat(filename.c_str(), &buffer));
}


/*
for(int j{0}; j < sigma_vec.size();++j){
            if(sigma_vec.at(j) < RN_ficcion){
                ++jj;
            }else{
                    RN_on_table = sigma_vec.at(jj);
                    if(mu_vec.at(jj) > mu_measured){
                        //cout << "jj = " << jj << endl;
                        ++jj;
                    }else{
                        mu_detected = mu_vec.at(jj);
                        xc_epix_oh1_fi = xc_epix_vec.at(jj);
                        ofile << "run bla" << '\t' << "ohdu bla" << '\t' << xc_epix_vec.at(jj) << endl;
                        break;
                    }
                }
        }
*/        
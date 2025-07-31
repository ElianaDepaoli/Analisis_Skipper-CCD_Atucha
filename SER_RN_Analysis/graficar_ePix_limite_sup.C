/* Autora: Eliana Depaoli Fecha: 23/07/25
Este código sirve para graficar el nuevo límite de integración de la variable ePix calculado con RN_SER_Varios_calculos_app.C, 
tabla_xc_vs_RN_SER.C y ePix_un_limite_te_pido.C 
Para su ejecución, hay que dar la ruta completa al archivo que se leerá:
-> ePix_te_puse_un_limite_rXX.txt 
este contiene runID, ohdu, ePix limite superior, RN medido y el usado el tabla_xxx.C, la corriente oscura (mu) medido y la
usada en tabla_xxx.C

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
void archivo_a_vectores(vector<double>& v1_vec, vector<double>& v2_vec, vector<double>& v3_vec, vector<double>& v4_vec,
 vector<double>& v5_vec, const char ohdu_value , const char*filename);

//Rutine ·······························································
const char* filename = "/home/eliana/Documentos/Atucha/Scripts/SEE_RN/ePix_te_puse_un_limite_r31.txt";
vector<vector<double>> RUNID_vec(2);
vector<vector<double>> ohdu_vec(2);
vector<vector<double>> ePix_vec(2);
vector<vector<double>> RN_vec(2);
vector<vector<double>> mu_vec(2);

//Plots ············
char color[4] = {kP10Yellow, kP10Violet, kP10Red,kP10Blue}; 
TString run = "RUN 31";
TString xaxis_label = "Run ID";
TString yaxis_ePix_label = "ePix [e-/pxl]";
TString yaxis_RN_label = "Readout Noise [e-]";
TString yaxis_mu_label = "Dark Currente [e-/pxl]";

void graficar_ePix_limite_sup(){
	archivo_a_vectores(RUNID_vec[0],ohdu_vec[0],ePix_vec[0],RN_vec[0],mu_vec[0], '2' ,filename);
	archivo_a_vectores(RUNID_vec[1],ohdu_vec[1],ePix_vec[1],RN_vec[1],mu_vec[1], '1' ,filename);

	// Ver vectores en pantalla ·····
    //for (int i{0}; i < RUNID_vec[0].size(); ++i){cout << RUNID_vec[0].at(i) << endl;}

    //Gráficos ···········
    //ePix ···
    TCanvas * c1 = new TCanvas("c1", "ePix Upper Limits", 1200, 800);
    gPad->SetGrid();
    //gStyle->SetOptStat(0);
    TGraph *go[2];
    for(int i=0; i < 2; ++i){
    	//cout << i << endl;
        //for (int j{0}; j < ohdu_vec[i].size(); ++j){cout << ohdu_vec[i].at(j) << endl;}
        go[i] = new TGraph(RUNID_vec[i].size(), &RUNID_vec[i][0], &ePix_vec[i][0]);//
        go[i]->SetMarkerColor(color[i]);//
        go[i]->SetMarkerStyle(20+i); //
    
    }
    
    go[0]->SetTitle(run);
    go[0]->GetYaxis()->SetTitle(yaxis_ePix_label);
    go[0]->GetXaxis()-> SetRangeUser(RUNID_vec[0].at(0)-RUNID_vec[0].at(0), RUNID_vec[0].at(RUNID_vec[0].size()-1)+13);
    go[0]->GetXaxis()->SetTitle(xaxis_label);
    go[0]-> Draw("AP");
    for(int i=1; i < 2; ++i) go[i]->Draw("sameP");
    auto l0 = new TLegend(0.75,0.15,0.90,0.30);
    
    for (int i=0; i<2; i++) {
      l0->AddEntry(go[i],TString:: Format("Ohdu %d", i+1), "p");
    }

    l0->Draw();

	//Readout Noise ···
    TCanvas * c2 = new TCanvas("c2", "Measured Readout Noise", 1200, 800);
    gPad->SetGrid();
    //gStyle->SetOptStat(0);
    TGraph *g1[2];
    for(int i=0; i < 2; ++i){
    	//cout << i << endl;
        //for (int j{0}; j < ohdu_vec[i].size(); ++j){cout << ohdu_vec[i].at(j) << endl;}
        g1[i] = new TGraph(RUNID_vec[i].size(), &RUNID_vec[i][0], &RN_vec[i][0]);//
        g1[i]->SetMarkerColor(color[i]);//
        g1[i]->SetMarkerStyle(20+i); //
    
    }
    
    g1[0]->SetTitle(run);
    g1[0]->GetYaxis()->SetTitle(yaxis_RN_label);
    g1[0]->GetYaxis()-> SetRangeUser(0.18, 0.26);
    g1[0]->GetXaxis()-> SetRangeUser(RUNID_vec[0].at(0)-RUNID_vec[0].at(0), RUNID_vec[0].at(RUNID_vec[0].size()-1)+13);
    g1[0]->GetXaxis()->SetTitle(xaxis_label);
    g1[0]-> Draw("AP");
    for(int i=1; i < 2; ++i) g1[i]->Draw("sameP");
    
    l0->Draw();
	
	//Dark Current ···
    TCanvas * c3 = new TCanvas("c3", "Measured Dark Current", 1200, 800);
    gPad->SetGrid();
    //gStyle->SetOptStat(0);
    TGraph *g2[2];
    for(int i=0; i < 2; ++i){
    	//cout << i << endl;
        //for (int j{0}; j < ohdu_vec[i].size(); ++j){cout << ohdu_vec[i].at(j) << endl;}
        g2[i] = new TGraph(RUNID_vec[i].size(), &RUNID_vec[i][0], &mu_vec[i][0]);//
        g2[i]->SetMarkerColor(color[i]);//
        g2[i]->SetMarkerStyle(20+i); //
    
    }
    
    g2[0]->SetTitle(run);
    g2[0]->GetYaxis()->SetTitle(yaxis_mu_label);
    g2[0]->GetYaxis()-> SetRangeUser(-0.01, 0.02);
    g2[0]->GetXaxis()-> SetRangeUser(RUNID_vec[0].at(0)-RUNID_vec[0].at(0), RUNID_vec[0].at(RUNID_vec[0].size()-1)+13);//0,1050);//
    g2[0]->GetXaxis()->SetTitle(xaxis_label);
    g2[0]->GetYaxis()->SetTitleOffset(1.2); // Desplazamiento de la etiqueta del eje Y
    g2[0]-> Draw("AP");
    for(int i=1; i < 2; ++i) g2[i]->Draw("sameP");
    
    l0->Draw();

	c1->SaveAs("ePix_UpperLimit_r31.png");
	c2->SaveAs("RN_measured_r31.png");
	c3->SaveAs("DarkCurrent_measured_r31.png");

}

//Function definition ·················································

// ·········
void archivo_a_vectores(vector<double>& v1_vec, vector<double>& v2_vec, vector<double>& v3_vec, vector<double>& v4_vec,
 vector<double>& v5_vec, const char ohdu_value , const char*filename) {
    //Saltea las filas cuya columna [6] es igual a ohdu_value

    // Abrir el archivo
    std::ifstream file(filename);
    //Si no pudo abrirlo
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
        
        //cout << " #4: " << line[4] << " #5: " << line[5] << " #6: " << line[6] << " #7: " << line[7] << endl;
        // Ignorar líneas vacías o comentarios
        if (line.empty() || line[0] == '#' || line[6]== ohdu_value || line[5]== ohdu_value || line[4]== ohdu_value) continue;

        // Reemplazar tabulaciones por espacios para facilitar la lectura
        std::replace(line.begin(), line.end(), '\t', ' ');

        // Usar un stringstream para dividir la línea en columnas
        std::stringstream ss(line);//esta clase opera sobre strings
        string value;
        int column = 0;

        // Leer cada valor de la línea
        while (ss >> value) {
            //cout << column << endl;
            if (column == 0) v1_vec.push_back(stod(value)); 
            if (column == 1) v2_vec.push_back(stod(value));
            if (column == 2) v3_vec.push_back(stod(value));
            if (column == 3) v4_vec.push_back(stod(value));
            if (column == 5) v5_vec.push_back(stod(value));
            column++;
        }
    }

    // Cerrar el archivo
    file.close();

    // Verificar que se hayan leído datos
    if (v1_vec.empty() || v2_vec.empty() || v3_vec.empty()) {
        std::cerr << "Error: No se encontraron datos válidos" << std::endl;
        return;}//else{ cout << "Vectores llenos" << endl;}

}

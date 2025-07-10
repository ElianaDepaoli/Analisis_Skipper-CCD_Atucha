/* Autora: Eliana Depaoli Fecha: 09/07/25
Este código sirve para hacer gráficos de ruido de lectura en función de la DC (SEE).
Para su ejecución, el archivo que contiene los datos debe estar en la misma carpeta que este script.
Los datos fuente fueron generados con RN_SER_Varios_calculos_app.C
Modificación:
09/07/25 Grafico SEE vs RunId y Ocupancy vs RunID de las 4 extensiones.
20/03/25 Implementé una solución para ordenar las fechas de las imágenes y los 
vectores enlazados a ella (como el que almacena el SEE y la ocupancia). 
Solución: creo un vector donde se guardan las fechas en formato numérico (unix
time stamp) y otro que guarda el valor de sus índices. Luego ordeno el vector de índices
según valores crecientes del tiempo en segundos. Finalmente, uso el vector de índices
ordenados para hacer una copia ordenada de los vectores que contienen las fechas y 
el SEE.    

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
void SEE_evolution_plot( vector<string> vecx, vector<Double_t> vecy,const int paso, const char *canv_name,const char *canv_title);
double date_to_seconds (string DATESTART);

//Rutine ·······························································
const char* filename = "/home/eliana/Documentos/Atucha/Scripts/SEE_RN/RN_SER_Atucha_r44.txt";
vector<vector<double>> RN_vec(2);
vector<vector<double>> err_RN_vec(2);
vector<vector<double>> lamb_ovx_vec(2);
vector<vector<double>> err_lamb_ovx_vec(2);
vector<vector<double>> lamb_ovy_vec(2);
vector<vector<double>> err_lamb_ovy_vec(2);
vector<vector<double>> see_vec(4);//vector de vectores
vector<vector<double>> occ_vec(4);
vector<string> date_start_vec;//fecha y horario de inicio de toma de la imagen
vector<double> RUNID_vec;
vector<double> timestamp_date_start_vec;//timestamp

TString run = "RUN_44";
TString cfulloutdirname = TString(gSystem->pwd());

double t_media_ov_px = 6.9942e-05 ;//[day]
double t_expo_overscan_y = 0.018687 ;//[day]
double exposure_oh1 = 2851;//[pxl*day]
double exposure_oh2 = 2942 ;//[pxl*day]
double exposure_oh3_4 = 3053;//[pxl*day]
double readout_time_per_pxl_per_sample = 65.61;// [microsec].q

void graf_SEE_Noise(){

    for(int i{0}; i < 4; ++i) archivo_a_vectores(see_vec[i], occ_vec[i], filename, 2*i+2, 2*i+3);
    for(int i{0}; i < 2; ++i) archivo_a_vectores(RN_vec[i], err_RN_vec[i], filename, 2*i+14, 2*i+15);
    for(int i{0}; i < 2; ++i) archivo_a_vectores(lamb_ovx_vec[i], err_lamb_ovx_vec[i], filename, 2*i+10, 2*i+11);
    for(int i{0}; i < 2; ++i) archivo_a_vectores(lamb_ovy_vec[i], err_lamb_ovy_vec[i], filename, 2*i+18, 2*i+19);
    archivo_a_string_vectores(date_start_vec,RUNID_vec,filename, 1, 0);
    
    // Ver vectores en pantalla ·····
    //for (int i{0}; i < lamb_ovy_vec[0].size(); ++i){cout << lamb_ovy_vec[0].at(i) << endl;}
    
    //Plots ············
    char color[4] = {kP10Yellow, kP10Violet, kP10Red,kP10Blue}; 
    
    //SER Overscan y ···
    TCanvas * c5 = new TCanvas("c5", "lambda OVY", 1200, 800);
    gPad->SetGrid();
    //gStyle->SetOptStat(0);
    TGraphErrors *glambOVY[2];
    for(int i=0; i < 2; ++i){
        glambOVY[i] = new TGraphErrors(RUNID_vec.size(), &RUNID_vec[0], &lamb_ovy_vec[i][0], 0, &err_lamb_ovy_vec[i][0]);
        glambOVY[i]->SetMarkerColor(color[i]);//
        glambOVY[i]->SetLineColor(color[i]);//
        glambOVY[i]->SetMarkerStyle(20+i); //
    
    }

    glambOVY[0]->SetTitle(run);
    glambOVY[0]->GetYaxis()->SetTitle("SEE Overscan Y [e-/pxl/day]");
    glambOVY[0]->GetYaxis()-> SetRangeUser(0.0, 1.0);
    glambOVY[0]->GetXaxis()->SetTitle("RUNID");
    glambOVY[0]-> Draw("AP");

    for(int i=1; i < 2; ++i) glambOVY[i]->Draw("sameP");
    
    auto l5 = new TLegend(0.75,0.75,0.90,0.90);;//(0.75,0.10,0.90,0.30);
    
    for (int i=0; i<2; i++) {
      l5->AddEntry(glambOVY[i],TString:: Format("Ohdu %d", i+1), "p");
    }

    l5->Draw();


    //SER Overscan x ···
    TCanvas * c4 = new TCanvas("c4", "lambda OVX", 1200, 800);
    gPad->SetGrid();
    //gStyle->SetOptStat(0);
    TGraphErrors *glambOVX[2];
    for(int i=0; i < 2; ++i){
        glambOVX[i] = new TGraphErrors(RUNID_vec.size(), &RUNID_vec[0], &lamb_ovx_vec[i][0], 0, &err_lamb_ovx_vec[i][0]);
        glambOVX[i]->SetMarkerColor(color[i]);//
        glambOVX[i]->SetLineColor(color[i]);//
        glambOVX[i]->SetMarkerStyle(20+i); //
    
    }

    glambOVX[0]->SetTitle(run);
    glambOVX[0]->GetYaxis()->SetTitle("SEE Overscan x [e-/pxl/day]");
    glambOVX[0]->GetYaxis()-> SetRangeUser(0.0, 250.0);
    glambOVX[0]->GetXaxis()->SetTitle("RUNID");
    glambOVX[0]-> Draw("AP");

    for(int i=1; i < 2; ++i) glambOVX[i]->Draw("sameP");
    
    auto l4 = new TLegend(0.75,0.75,0.90,0.90);
    
    for (int i=0; i<2; i++) {
      l4->AddEntry(glambOVX[i],TString:: Format("Ohdu %d", i+1), "p");
    }

    l4->Draw();

    //Readout Noise ···
    TCanvas * c3 = new TCanvas("c3", "Readout Noise", 1200, 500);
    gPad->SetGrid();
    //gStyle->SetOptStat(0);
    TGraphErrors *gRN[2];
    for(int i=0; i < 2; ++i){
        gRN[i] = new TGraphErrors(RUNID_vec.size(), &RUNID_vec[0], &RN_vec[i][0], 0, &err_RN_vec[i][0]);
        gRN[i]->SetMarkerColor(color[i]);//
        gRN[i]->SetLineColor(color[i]);//
        gRN[i]->SetMarkerStyle(20+i); //
    
    }

    gRN[0]->SetTitle(run);
    gRN[0]->GetYaxis()->SetTitle("Readout Noise [e-]");
    gRN[0]->GetYaxis()-> SetRangeUser(0.180, 0.30);
    gRN[0]->GetXaxis()->SetTitle("RUNID");
    gRN[0]-> Draw("AP");

    for(int i=1; i < 2; ++i) gRN[i]->Draw("sameP");
    
    auto l3 = new TLegend(0.75,0.75,0.90,0.90);//(0.75,0.10,0.90,0.35);
    
    for (int i=0; i<2; i++) {
      l3->AddEntry(gRN[i],TString:: Format("Ohdu %d", i+1), "p");
    }

   l3->Draw();

    //OCC ···
    TCanvas * c2 = new TCanvas("c2", "occupancy", 1200, 800);
    gPad->SetGrid();
    //gStyle->SetOptStat(0);
    TGraph *go[4];
    for(int i=0; i < 4; ++i){
        go[i] = new TGraph(RUNID_vec.size(), &RUNID_vec[0], &occ_vec[i][0]);
        go[i]->SetMarkerColor(color[i]);//
        go[i]->SetMarkerStyle(20+i); //
    
    }

    go[0]->SetTitle(run);
    go[0]->GetYaxis()->SetTitle("Occupancy [e-/pxl/day]");
    go[0]->GetYaxis()-> SetRangeUser(0.0, 1.0);
    go[0]->GetXaxis()->SetTitle("RUNID");
    go[0]-> Draw("AP");

    for(int i=1; i < 4; ++i) go[i]->Draw("sameP");
    
    auto lo = new TLegend(0.75,0.75,0.90,0.90);
    
    for (int i=0; i<4; i++) {
      lo->AddEntry(go[i],TString:: Format("Ohdu %d", i+1), "p");
    }

   lo->Draw();


   //SEE ···
    TCanvas * c1 = new TCanvas("c1", "SEE Clusters", 1200, 800);
    gPad->SetGrid();
    //gStyle->SetOptStat(0);
    TGraph *g[4];
    for(int i=0; i < 4; ++i){
        g[i] = new TGraph(RUNID_vec.size(), &RUNID_vec[0], &see_vec[i][0]);
        g[i]->SetMarkerColor(color[i]);//
        g[i]->SetMarkerStyle(20+i); //
    
    }

    g[0]->SetTitle(run);
    g[0]->GetYaxis()->SetTitle("SEE [e-/pxl/day]");
    g[0]->GetYaxis()-> SetRangeUser(0.0, 0.50);
    g[0]->GetXaxis()->SetTitle("RUNID");
    g[0]-> Draw("AP");

    for(int i=1; i < 4; ++i) g[i]->Draw("sameP");
    
    auto l1 = new TLegend(0.75,0.75,0.90,0.90);
    
    for (int i=0; i<4; i++) {
      l1->AddEntry(g[i],TString:: Format("Ohdu %d", i+1), "p");
    }

   l1->Draw();

    //SEE_evolution_plot(date_start_vec, see_vec[0], 10, "SEE evolution","SEE evolution");
    //SEE_evolution_plot(date_start_vec, RN_oh1_vec, 10, "RN","RN evolution");
    TString OVY_out = "SEE_OVY_" + run + ".png";
    TString OVX_out = "SEE_OVX_" + run + ".png";
    TString RN_out = "RN_" + run + ".png";
    TString occ_out = "occupancy_" + run + ".png";
    TString see_out = "active_area_see_" + run + ".png";
    c5->SaveAs(Form("%s/%s", cfulloutdirname.Data(), OVY_out.Data()));
    c4->SaveAs(Form("%s/%s", cfulloutdirname.Data(), OVX_out.Data()));
    c3->SaveAs(Form("%s/%s", cfulloutdirname.Data(), RN_out.Data()));
    c2->SaveAs(Form("%s/%s", cfulloutdirname.Data(), occ_out.Data()));
    c1->SaveAs(Form("%s/%s", cfulloutdirname.Data(), see_out.Data()));

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

void SEE_evolution_plot( vector<string> vecx, vector<Double_t> vecy,const int paso, const char *canv_name,const char *canv_title)
{
    TCanvas * c = new TCanvas(canv_name,canv_title, 1000, 750);
    c->cd();
    c->SetGrid(2);
    //fill the TGraph
    int n = TMath::Min(vecy.size(), vecx.size());
    TGraph* tg = new TGraph(n);
    for (Int_t i = 0; i <n; i++) tg->SetPoint(i, i + 1., vecy[i]);//Set x and y values for point number i.
    //Set x axis label
    auto h = new TH1F("h","h",n,0,n);
    tg->SetHistogram(h);
    int k{0};
    while(k<n)
    {
        h->GetXaxis()->SetBinLabel(k + 1, vecx[k].c_str());
        k +=paso;
    }

    h->GetXaxis()-> LabelsOption("v");
    h->SetStats(0);
    h->Draw("AXIS");
    //tg->GetXaxis()->SetTimeDisplay(1); not working 'cause must put labels in seconds. Later on.
    tg->SetMarkerColor(kOrange+1);//(kAzure-1);//kRed);
    tg->SetMarkerStyle(33);
    //tg->SetLineColor(kRed);
    tg->SetTitle("Atucha");//not working don't know why
    tg->GetYaxis()->SetTitle("#SEE [1e events/day/pix]");
    //tg->GetYaxis()->SetRangeUser(0., 0.2);
    tg->Draw("P same");
    gPad->RedrawAxis("g");
}

double date_to_seconds (string DATESTART)//Char_t DATESTART[0]
{
    struct tm tm {};
    strptime(DATESTART.c_str(),"%Y-%m-%dT%H:%M:%S",&tm);//converting string to date/time
    time_t date = mktime(&tm);//converting tm structure to time_t format 
    return date;
}



//······························································
    // Ordeno vectores por fecha creciente                      ····
    //······························································
    //fechas en formato unix time stamp ··
    //for(int i{0}; i < date_start_vec.size();++i) timestamp_date_start_vec.push_back(date_to_seconds(date_start_vec.at(i)));
    //Indexar ··
    //vector<int> indices_dates_vec(date_start_vec.size());//creo un vector donde guardar los índices reales del vector de fechas
    //std::iota(indices_dates_vec.begin(), indices_dates_vec.end(),0);//lleno vector de índices
    //reordeno los índices en orden creciente de tiempo usando las fechas como referencia ··
    //std::sort(indices_dates_vec.begin(),indices_dates_vec.end(), [&](int A, int B)-> bool {return timestamp_date_start_vec[A] < timestamp_date_start_vec[B];});
    //genero nuevo vector con fechas y see ordenados ··
    //for(int i{0}; i < indices_dates_vec.size(); ++i){
    //    ordered_SEE_vec.push_back(SEE_vec.at(indices_dates_vec.at(i)));
    //    ordered_date_start_vec.push_back(date_start_vec.at(indices_dates_vec.at(i)));
    //}

    //Pantalla ····································
    //for(int i{0}; i < 10 ;++i) cout << indices_dates_vec.at(i) << endl;
    //for(int i{0}; i < 10 ;++i) cout << timestamp_date_start_vec.at(i) << endl;
    //for(int i{0}; i < 10 ;++i) cout << ordered_SEE_vec.at(i) << endl;
    //for(int i{0}; i < 10 ;++i) cout << ordered_date_start_vec.at(i) << endl;
    //······························································
    


    // Graficos ·······················································
    // Ocupancia, RN ···
    /*TCanvas* c1 = new TCanvas("canvas", "Occ & Readout noise", 1200, 800);
    c1->SetGrid();
    c1->SetLeftMargin(0.12);   // Margen izquierdo (ajusta según sea necesario)
    c1->SetRightMargin(0.04);  // Margen derecho
    c1->SetBottomMargin(0.12); // Margen inferior
    c1->SetTopMargin(0.08);    // Margen superior
    TGraph* gr1 = new TGraph(RN_oh1_vec.size(), &occ_vec[0], &RN_oh1_vec[0]);
    TGraph* gr2 = new TGraph(RN_oh2_vec.size(), &occ_vec[0], &RN_oh2_vec[0]);
    // Ajustar los límites del eje
    double xmin = *std::min_element(occ_vec.begin(), occ_vec.end());
    double xmax = *std::max_element(occ_vec.begin(), occ_vec.end());
    double ymin = std::min(*std::min_element(RN_oh1_vec.begin(), RN_oh1_vec.end()),*std::min_element(RN_oh2_vec.begin(), RN_oh2_vec.end()));
    double ymax = std::max(*std::max_element(RN_oh1_vec.begin(), RN_oh1_vec.end()),*std::max_element(RN_oh2_vec.begin(), RN_oh2_vec.end()));

    gr1->GetXaxis()->SetLimits(0.5, 1.3);//xmin, xmax 
    gr1->GetYaxis()->SetRangeUser(ymin-0.001, ymax+0.001);
    // Configurar el gráfico
    gr1->SetTitle("RUN 31 - Reactor ON");//Occupancy & Readout noise - 
    gr1->SetMarkerStyle(20); //
    gr1->SetMarkerColor(kBlue);
    gr1->GetXaxis()->SetTitle("Occupancy (e-/pxls/day)");
    gr1->GetYaxis()->SetTitle("RN [e-]");
    gr1->GetXaxis()->SetTitleSize(0.05);
    gr1->GetYaxis()->SetTitleSize(0.05);
    gr1->GetXaxis()->SetTitleOffset(1.0); // Desplazamiento de la etiqueta del eje X
    gr1->GetYaxis()->SetTitleOffset(1.0); // Desplazamiento de la etiqueta del eje Y
    gr1->Draw("AP"); // "A" para dibujar ejes, "P" para dibujar puntos
    gr2->SetMarkerStyle(22); //
    gr2->SetMarkerColor(kRed); //
    gr2->Draw("P same");

    auto legend = new TLegend(0.15,0.45,0.25,0.65); 
    legend->AddEntry(gr1,"Ohdu 1","p");
    legend->AddEntry(gr2,"Ohdu 2","p");
    legend->SetTextSize(0.04); 
    legend->SetBorderSize(0);
    legend->Draw();
    //c1->SaveAs("Noise_occupancy_run31.png");
*/
    // SEE, RN ···
   /* TCanvas* c2 = new TCanvas("canvas", "DC & Readout noise", 1200, 800);
    c2->SetGrid();
    c2->SetLeftMargin(0.12);   // Margen izquierdo (ajusta según sea necesario)
    c2->SetRightMargin(0.04);  // Margen derecho
    c2->SetBottomMargin(0.12); // Margen inferior
    c2->SetTopMargin(0.08);    // Margen superior
    TGraph* gr3 = new TGraph(RN_oh1_vec.size(), &SEE_vec[0], &RN_oh1_vec[0]);
    TGraph* gr4 = new TGraph(RN_oh2_vec.size(), &SEE_vec[0], &RN_oh2_vec[0]);
    
    gr3->GetXaxis()->SetLimits(0.18, 0.42);//xmin, xmax 
    gr3->GetYaxis()->SetRangeUser(0.202, 0.214);//ymin-0.001, ymax+0.001
    // Configurar el gráfico
    gr3->SetTitle("DC & Readout noise - RUN 31 - Reactor ON");//Occupancy & Readout noise - 
    gr3->SetMarkerStyle(20); //
    gr3->SetMarkerColor(kBlue);
    gr3->GetXaxis()->SetTitle("SEE (e-/pxls/day)");
    gr3->GetYaxis()->SetTitle("RN [e-]");
    gr3->GetXaxis()->SetTitleSize(0.05);
    gr3->GetYaxis()->SetTitleSize(0.05);
    gr3->GetXaxis()->SetTitleOffset(1.0); // Desplazamiento de la etiqueta del eje X
    gr3->GetYaxis()->SetTitleOffset(1.0); // Desplazamiento de la etiqueta del eje Y
    gr3->Draw("AP"); // "A" para dibujar ejes, "P" para dibujar puntos
    gr4->SetMarkerStyle(22); //
    gr4->SetMarkerColor(kRed); //
    gr4->Draw("P same");

    auto l2 = new TLegend(0.75,0.20,0.95,0.30); 
    l2->AddEntry(gr3,"Ohdu 1","p");
    l2->AddEntry(gr4,"Ohdu 2","p");
    l2->SetTextSize(0.04); 
    l2->SetBorderSize(0);
    l2->Draw();
    //c2->SaveAs("Noise_SEE_run31.png");
*/
    
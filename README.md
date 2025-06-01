# Analisis_Skipper-CCD_Atucha
Códigos para recortar catálogos, calcular variables de interés (RN, SER, Hot Columns, etc) , obtener espectros.
Se listan en el orden en el que deben ser utilizados.

#### cutter.C ####
Aplica cortes de calidad a cada archivo root. 

Cómo usarlo:
	Copiar este script dentro de el directorio "run_**/", los archivos a procesar están en el subdirectorio hits/. 
	Dar un nombre al sudirectorio donde se almacenarán los archivos procesados. Típicamente "Cutted_cortes_aplicados/" Ejemplo: "/Cutted_calPixTree_HOT_COL_new_EDGES_xVarMin_0_xVarMax_4_yVarMin_0yVarMax_4"
  Este subdirectorio será generado por el código.

Outputs:	
	Un nuevo root file por cada archivo leído, con el mismo nombre pero con el prefijo Cutted_ agregado y ubicado dentro de "Cutted_cortes_aplicados/".

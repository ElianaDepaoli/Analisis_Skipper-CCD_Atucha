# Analisis_Skipper-CCD_Atucha
Códigos para recortar catálogos, calcular variables de interés (RN, SER, Hot Columns, etc) , obtener espectros.
Se listan en el orden en el que deben ser utilizados.

### cutter.C ###
Aplica cortes de calidad a cada archivo root. 

#### Cómo usarlo:
	Copiar este script dentro de el directorio "run_**/", los archivos a procesar están en el subdirectorio hits/. 
	Dar un nombre al sudirectorio donde se almacenarán los archivos procesados. Típicamente "Cutted_cortes_aplicados/" Ejemplo: "/Cutted_calPixTree_HOT_COL_new_EDGES_xVarMin_0_xVarMax_4_yVarMin_0yVarMax_4"

#### Outputs:	
	Un nuevo root file por cada archivo leído, con el mismo nombre pero con el prefijo Cutted_ agregado y ubicado dentro de "Cutted_cortes_aplicados/".
	El subdirectorio de destino es generado por el código.

### addtreetoimage.cpp ####

Incorpora a los catálogos ya cortados un nuevo tree 'calc' que contiene las siguientes nuevas branches: 
-> SEE, DEE, ocupancia sobre los 2 ohdu. Unidades: [e/pix/day]
-> Probabilidad de cada evento de ser neutrino en base a la pdf construida con las varianzas de una simulacion. 
-> SER [e/pix/day] , RN y sus errores para cada ohdu. Calculados por ajustes.
-> HOTCOL_amount, exposure time en días

No cambia el nombre de los root file.

#### Requisitos:
	Que exista atucha1.root en algún directorio, conteniendo las simulaciones MC necesarias para calcular la probabilidad de ser neutrino usando las varianzas.

#### Antes de usarlo:
	
	Solo si es necesario: editar en el código fuente *.cpp las siguientes variables: 
	-> hot_col_amount = cantidad de columnas brillantes recortadas con cutter.C; 
	-> factor = 0.5 si se lee entre dos imágenes sucesivas, factor = 1 si no
	-> filename_sim = {"/home/ruta_completa/atucha1.root"};
	Copiar todos los Cutted*.root del paso anterior dentro de una carpeta OFF y otra ON segun corresponda. No es obligatorio para hacer un análisis que no es ON vs OFF.
	Compilarlo con la siguiente linea: g++ -o addtreetoimage.exe addtreetoimage.cpp `root-config --cflags --glibs`
	Para que la compilación ocurra ROOT debe estar instalado.
	
#### Cómo usarlo: se corre el ejecutable generado con la compilación
	
	./addtreetoimage.exe /mnt/Data/Atucha_Data/OFF1_2022
	./addtreetoimage.exe /mnt/Data/Atucha_Data/ON1_2022
	./addtreetoimage.exe /mnt/Data/Atucha_Data/OFF2_2022
	./addtreetoimage.exe /mnt/Data/Atucha_Data/OFF_2023
	./addbranchestoroot.exe /mnt/Data/Atucha_Data/2022/ON1/Cutted_sin_extra_branches/ /mnt/Data/Atucha_Data/ON1_2022/
	

Una copia de atucha1.root en Disco Atucha #1: /0e0c1913-ca43-4252-b2bb-783703f57ea5/proc
### Generar una lista ordenada de los catálogos a graficar ###
ls -1A | sort -t'_' -k6,6n -k12.4,12n > OFF_R12_2023.txt
genera una lista donde los nombres de los catálogos aparecen ordenados por Número de RUN (campo 6) y luego por número de imágen (campo 12, a partir del 4º carácter)


### spectra.C ###




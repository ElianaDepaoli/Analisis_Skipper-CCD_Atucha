RN_SER_Varios_calculos_app.C
Sirve para calcular la corriente oscura y el ruido de lectura en cada uno de los 4 cuadrantes del Skipper-CCD. 
*Guarda los resultados en un archivo que se puede graficar con graf_SEE_Noise.C
*El archivo txt que produce sirve de entrada a ePix_un_limite_te_pido.C
tabla_xc_vs_RN_SER.C
*Produce tabla_RN_SER_xcePix_Atucha.txt
este contiene el valor calculado de mu (DC) para un conjunto de valores hipotéticos para RN y xC (límite que iguala las probabilidades de determinar erróneamente que un pixel está prendido o apagado)
ePix_un_limite_te_pido.C
Sirve para determinar el límite de integración de la variable ePix usando la corriente oscura y el ruido de lectura reales de un run
*Requiere como entradas las salidas de RN_SER_Varios_calculos_app.C y de tabla_xc_vs_RN_SER.C
graficar_ePix_limite_sup.C
Sirve para graficar el nuevo límite de integración de la variable ePix calculado con RN_SER_Varios_calculos_app.C, 
tabla_xc_vs_RN_SER.C y ePix_un_limite_te_pido.C 
*Para su ejecución, hay que dar la ruta completa al archivo que se leerá: ePix_te_puse_un_limite_rXX.txt 

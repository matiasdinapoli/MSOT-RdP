# MSOT-RdP
Modelo para la Simulación de Ondas de Tormenta - Río de la Plata y su Plataforma Continental adyacente

  El "Modelo para la Simulación de Ondas de Tormenta - Río de la Plata y Plataforma Continental adyacente" (MSOT-RdP) es un sistema de modelación de ondas de tormenta (ODT) establecido por el Dr. Dinápoli como parte de su Tesis Doctoral en Cs. de la Atmósfera y los Océanos en la Universidad de Buenos Aires. En particular, el Dr. Dinápoli realizó su investigación en el grupo "Procesos y Modelado de la Plataforma Continental Argentina y el Río de la Plata" (http://www.cima.fcen.uba.ar/pplata.php#) del Centro de Investigaciones del Mar y la Atmósfera (CIMA/CONICET-UBA). El sistema MSOT-RdP es una adecuación barotrópica-bidimensión (Dinápoli et. al, 2020a) del modelo oceánico CROCO (https://www.croco-ocean.org, Debreu et al., 2012) que trabaja en un esquema de anidamiento unidireccional y que es forzado superficialmente con viento y presión atmosférica en superfice y lateralmente con descarga continental y marea. Esta configuración es la óptima para una adecuada representación de la dinámica barotrópica de la región (Dinápoli et al., 2020b) e incluso su componente no lineal (Dinápoli et al., 2020c,d). Cabe mencionar que MSOT-RdP también funciona en modo "ensambles" (Dinápoli et al., 2021) pero de manera pre-operativa con el forzante HRES (https://www.ecmwf.int/en/forecasts/datasets/set-i).

  En este repositorio se presenta la configuración de MSOT-RdP para ser ejecutado en modo pronóstico determinístico forzado con las soluciones del modelo GFS (https://nomads.ncep.noaa.gov/), en particular con el subcojunto "GFS sflux". Tal adecuación consta de dos archivos: 
  (1) run_msot_forecast.py que tiene la configuracion del modelo, parámetros, los path hacia los archivos principales, locaciones, etc.
  (2) msot_forecast.py que tiene el pre y post procesamiento de todos los archivos que necesita y devuelve (respectivamente) el modelo a lo largo su ejecución.
 El sistema también necesita este conjunto de archivos para funcionar:
  (i) Dominios numéricos: dominio0.nc que cubre toda la Plataforma Continental Argentina para realizar el anidamiento dinámico de los modelo globales al modelo regional y dominio1.nc que es anidamiento hacia el Río de la Plata y su Plataforma Continental adyacente.
  (ii) Forzantes laterales: datos_marea.nc contiene un recorte del modelo de marea TXPO9 (https://www.tpxo.net/global/tpxo9-atlas) con los valores de amplitud y fase de las principales componentes de marea para la región; DescargaContinental.dat contiene las observaciones de descarga continental del Río de la Plata que son provistas por el Instituto Nacional del Agua.
  (iii) Modelo compilado: croco10nuc0 y croco10nuc1 son la adecuación del modelo CROCO para los dos dominios, el 10 indica la canidad de núcleos que utiliza para su ejecución.
  
 Para realizar un pronóstico sólo es necesario tener todos estos archivos juntos y ejecutar:
  $ ./run_msot_forecast.py
El programa según se lo indique descargará el forzante meteorológico o lo leerá de algún path, prepará el forzante ejecutará el modelo y devolverá el pronóstico en:

  path/prono{fecha}/{fecha}T{hora_de_inicio}_{dominio}_his.nc
  
donde "hora_de_inicio" indica alguno de los 4 ciclos diarios de inicializacion y domino si se refiere a la Plataforma Continental Argentina (0 ó dominio0) o al Río de la Plata y su Plataforma Continental adyacente (1 ó dominio1). En run_msot_forcast.py se puede explicitar si la fecha incial y final o si se quiere un pronostico operativo. 

Requerimientos:
  1) Para la simulación: FORTRAN
  2) Pre y post procesamientos: Python 3 y las librerías xarray y uptides.

Referencias

Debreu, L., P. Marchesiello, P. Penven, y G. Cambon, 2012: Two-way nesting in split-explicit ocean models: algorithms, implementation and validation. Ocean Modelling, 49-50, 1-21.
Dinápoli, M.G., Simionato, C.G. y Moreira, D., 2020a: Development and validation of a storm surge forecasting/hindcasting modelling system for the extensive Río de la Plata Estuary and its adjacent Continental Shelf. Natural Hazards. https://doi.org/10.1007/s11069-020-04079-5.
Dinápoli, M.G., Simionato, C.G. y Moreira, D., 2020b: Model Sensitivity during Extreme Positive and Negative Surges in the Río de la Plata Estuary: Highlighting the Need for an Appropriate Hindcast/Forecast System. Weather and Forecasting, Volume 35, pp. 1097-1112.. DOI: 10.1175/WAF-D-19-0171.1.
Dinápoli, M.G., Simionato, C.G. y Moreira, D., 2020c: Nonlinear tide-surge interactions in the Río de la Plata Estuary. Estuarine, Coastal and Shelf Science, Volume 241, 106834. ISSN 0272-7714. https://doi.org/10.1016/j.ecss.2020.106834.
Dinápoli, M.G., Simionato, C.G. y Moreira, D. 2020d: Nonlinear interaction between the tide and the storm surge with the current due to the flow of the tributary rivers in the tidal freshwater zone of the Río de la Plata. Estuaries and Coasts. https://doi.org/10.1007/s12237-020-00844-8.
Dinápoli, M.G., Simionato, C.G. y Moreira, D., 2021a: Development and evaluation of an ensemble forecast/hindcast system for storm surges at the Río de la Plata Estuary. QJR Meteorol Soc. 2021; 147: 557-572. https://doi.org/10.1002/qj.3933.

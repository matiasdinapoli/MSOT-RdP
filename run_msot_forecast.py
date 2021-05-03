#!/usr/bin/env python

""" 
En este repositorio se presenta la configuración de MSOT-RdP para ser
ejecutado en modo pronóstico determinístico forzado con las soluciones del
modelo GFS (https://nomads.ncep.noaa.gov/), en particular con el subcojunto
"GFS sflux". Tal adecuación consta de dos archivos: 
	(1) run_msot_forecast.py
que tiene la configuracion del modelo, parámetros, los path hacia los archivos
principales, locaciones, etc. 
	(2) msot_forecast.py que tiene el pre y post
procesamiento de todos los archivos que necesita y devuelve (respectivamente)
el modelo a lo largo su ejecución. 
El sistema también necesita este conjunto de archivos para funcionar: 
	(i) Dominios numéricos: dominio0.nc que cubre toda la Plataforma Continental Argentina 
para realizar el anidamiento dinámico de los modelo globales al modelo regional y 
dominio1.nc que es anidamiento hacia el Río de la Plata y su Plataforma Continental 
adyacente. 
	(ii) Forzantes laterales: datos_marea.nc contiene un recorte del modelo de marea TXPO9
(https://www.tpxo.net/global/tpxo9-atlas) con los valores de amplitud y fase
de las principales componentes de marea para la región;
DescargaContinental.dat contiene las observaciones de descarga continental del
Río de la Plata que son provistas por el Instituto Nacional del Agua. 
	(iii) Modelo compilado: croco10nuc0 y croco10nuc1 son la adecuación del modelo CROCO
para los dos dominios, el 10 indica la canidad de núcleos que utiliza para su
ejecución.
"""
##################################################################################
# Librerias
##################################################################################
from msot_forecast import MSOT
import pandas as pd

##################################################################################
# Archivos para el modelo
##################################################################################
outpath = "./"
dominio0 = outpath + "dominio0.nc"
dominio1 = outpath + "dominio1.nc"
marea = outpath + "datos_marea.nc"
rios = outpath + "DescargaContinental.dat"

##################################################################################
# Ejecucion
##################################################################################
# Si atmo dice "descargar" descargará los campos sino ira al link que le es indicado 

# En el caso que no se de una fecha final se corre indeterminadamente
fecha0 = pd.Timestamp.today().floor("D")
fecha1 = None

while fecha0 != fecha1:
	# Descripcion general del modelo
	sim = MSOT(name = fecha0.strftime("%Y%m%dT%H%M"),
	           date0 = fecha0,
	           domain = dominio0,
	           execut = "croco10nuc0",
	           threads = 10,
	           delta_t = 15,)
	# Modelo anidado
	sim.add_nested(domain = dominio1, 
				   execut = "croco10nuc1", 
				   threads = 10, 
				   delta_t = 5)
	sim.add_forcing(tide = [0], tide_file = marea,
	                atmo = "descargar", 
	                pressure = True)
	sim.add_runoff(1, rios, [[6,116,0,1],[11,121,1,-1],[12,121,1,-1]])
	# Parametros
	sim.add_parameters("rho0", 1025)
	sim.add_parameters("cl", 0.0e-4)
	sim.add_parameters("cd", 2.0e-3)
	sim.run()
	fecha0 += pd.Timedelta("6H")

#!/usr/bin/env python
################################################################################
# Liberias
################################################################################
import xarray as xr                     # Uso de campos
import pandas as pd                     # Uso de dataframes y fechas
import numpy as np                      # Operaciones matematicas
import uptide as tides                  # Analisis de marea, correciones nodales
import os                               # Instrucciones a linux
import multiprocessing                  # Paralelizar procesos
from itertools import product
from scipy.signal import hilbert        # Transformada de Hilbert
import scipy.cluster.hierarchy as hac   # Clusterizacion

################################################################################

################################################################################
# Definciones generales
################################################################################
"""
Posibilitar la conexion con otras redes
"""
os.environ["http_proxy"] = "http://proxy.fcen.uba.ar:8080/"

################################################################################
# Modelo para la Simulacion de Ondas de Tormenta : MSOT
################################################################################
class MSOT(object):
    """
    Descripcion general del modelo: nombre de la simulacion, fechas con un mes de spin up y
    5 dias de pronostico, dominio general, modelo compilado, resolucion temporal, 
    cantidad de nucleos y calculo de incerteza
    """
    def __init__(self, name, date0, domain, execut, delta_t, 
                 threads = 10):
        # Nombre del modelo
        self.name = name
        # |-----------||-------------|
        # d0 SpinUp   df    Prono    d1
        self.date0 = pd.to_datetime(date0) - pd.DateOffset(days = 3)
        self.datef = pd.to_datetime(date0)
        self.date1 = pd.to_datetime(date0) + pd.DateOffset(days = 10)
        self.dom0 = domain
        self.exect0 = execut
        self.threads_0 = threads
        self.dt = delta_t
        # Locacion de los archivos outputs
        self.locacion = "prono{}/".format(self.datef.strftime("%Y%m%d"))
         # Creo la carpeta contenedora en el caso que no exista
        if not os.path.isdir(self.locacion):
            os.system("mkdir {}".format(self.locacion))
        return

    """
    Agregar modelo anidado
    """
    def add_nested(self, domain, execut, delta_t, threads):
        level = str(len([x[-1] for x in self.__dict__.keys() if x[0:3]=="dom"]))
        setattr(self, "dom" + level, domain)
        setattr(self, "exect" + level, execut)
        lst = []; lst.append(self.dt); lst.append(delta_t)
        setattr(self, "dt", lst)
        setattr(self, "threads_" + level, threads)
        return

    """
    Forzantes
    """
    def add_forcing(self, **keys):
        if "atmo" in keys:
            setattr(self, "atmo_forcing", keys["atmo"])
        if "tide" in keys:
            setattr(self, "tide_forcing", keys["tide"])
            if "tide_file" in keys:
                setattr(self, "tide_file", keys["tide_file"])
        if "pressure" in keys:
            setattr(self, "Pair", keys["pressure"])
        else:
            setattr(self, "Pair", False)
        return

    """
    Descarga continental
    """
    def add_runoff(self, level, runoff_data, quantity):
        setattr(self, "runoff" + str(level), [runoff_data, quantity])
        return

    """
    Fijar parametros
    """
    def add_parameters(self, name, value):
        setattr(self, name, value)
        return

    """
    Correr el modelo
    """
    def run(self):
        # Fijo la cantidad de nucleos a utilizar
        os.environ["OMP_NUM_THREADS"] = str(self.threads_0)
        # Cento cuantos anidados ("dom") argegue para realizar la iteracion por dominio
        numero_de_modelos = len([x[-1] for x in self.__dict__.keys() if x[0:3]=="dom"])
        """
        Modelo dinamico: dominio grande y solucion de control
        """
        for i in range(numero_de_modelos):
            # Archivo input
            print("Creando el archivo input para el dominio {}".format(i))
            input_file = self.make_input(i)
            # Condicion inicial
            print("Creando la condicion inicial para el dominio {}".format(i))
            self.make_init(i)
            # Forzante
            print("Creando el forzante para el dominio {}".format(i))
            ## Creo el achivo forzante
            self.make_forcing(i)
            # Descarga continental
            if "runoff{}".format(i) in self.__dict__.keys():
                print("Creando la descarga continental para el dominio {}".format(i))
                self.make_runoff(i)
            # Condiciones de contorno para dominios anidados
            if i > 0:
                print("Creando condiciones de contorno para el dominio {}".format(i))
                self.make_boundary(i)
            # Run
            print("Corriendo el modelo {} nivel {}".format(self.name, i))
            sentencia = "./{} {} > log{}_{} && pid=$! && wait $pid".format(eval("self.exect{}".format(i)),
                                                                           input_file, i, self.name)
            ## En el caso que sea un dominio anidado recargo el numero de nucleos
            if i > 0:
                os.environ["OMP_NUM_THREADS"] = str(eval("self.threads_{}".format(i)))
            os.system(sentencia)


        # Ajuste de coordenadas del netcdf final
        print("Ajuste de coordenadas")
        self.coord_adjust("{}{}_{}_his.nc".format(self.locacion, self.name, i))

        """
        Elimino los archivos de input y log
        """
        os.system("rm *.in log*")       
        return

    ################################################################################
    #    Funciones del metodo run()
    ################################################################################

    """
    Archivo input
    """
    def make_input(self, n, **argkeys):
        # Creo el archivo input dependiendo si es un modelo deterministico o por ensambles
        txt = open("{}_{}.in".format(self.name, n),'w')
        txt.write("title:\n")
        txt.write(self.name.upper() + "\n")
        txt.write("start_date:\n")
        txt.write(self.date0.strftime("%Y-%b-%d" + " 00:00:00\n"))
        txt.write("time_stepping: NTIMES   dt[sec]  NDTFAST  NINFO\n")
        if len(self.dt) > 1:
            txt.write("{} {} 1 {} \n".format(int((self.date1-self.date0).total_seconds()/self.dt[n]), 
                                             int(self.dt[n]), int(21600/self.dt[n])))
        else:
            txt.write("{} {} 1 {} \n".format(int((self.date1-self.date0).total_seconds()/self.dt), 
                                             int(self.dt), int(21600/self.dt)))
        txt.write("grid:  filename\n")
        txt.write(eval("self.dom" + str(n)) + "\n")
        txt.write("forcing: filename\n")
        txt.write("{}{}_{}_frc.nc \n".format(self.locacion, self.name, n))
        if n != 0:
            txt.write("boundary: filename \n")
            txt.write("{}{}_{}_bry.nc \n".format(self.locacion, self.name, n))
        txt.write("initial: NRREC / filename\n")
        txt.write("0\n")
        txt.write("{}{}_{}_ini.nc \n".format(self.locacion, self.name, n))
        txt.write("restart:          NRST, NRPFRST / filename\n")
        if len(self.dt) > 1:
            txt.write("{} 0 \n".format(int(86400/self.dt[n])))
        else:
            txt.write("{} 0 \n".format(int(86400/self.dt)))
        txt.write(self.locacion + "rst_roms.nc\n")
        txt.write("history: LDEFHIS, NWRT, NRPFHIS / filename\n")
        if len(self.dt) > 1:
            txt.write("T {} 0 \n".format(int(3600/self.dt[n])))
        else:
            txt.write("T {} 0 \n".format(int(3600/self.dt)))
        txt.write("{}{}_{}_his.nc \n".format(self.locacion, self.name, n))
        txt.write("primary_history_fields: zeta UBAR VBAR  U  V   wrtT(1:NT)\n")
        txt.write("T    T   T   F  F    30*F\n")
        txt.write("rho0:\n")
        txt.write("{} \n".format(self.rho0))
        txt.write("lateral_visc:   VISC2,    VISC4    [m^2/sec for all] \n")
        txt.write("0.       0.\n")
        txt.write("bottom_drag:     RDRG [m/s],  RDRG2,  Zob [m],  Cdb_min, Cdb_max\n")
        txt.write("{} {} 1.d-3     1.d-4    1.d-1\n".format(self.cl, self.cd))
        txt.write("gamma2:\n")
        txt.write("1.\n")
        txt.write("nudg_cof:    TauT_in, TauT_out, TauM_in, TauM_out  [days for all]\n")
        txt.write("1.       360.      3.      360.\n")
        if "runoff{}".format(n) in self.__dict__.keys():
            txt.write("psource_ncfile:   Nsrc  Isrc  Jsrc  Dsrc qbardir  Lsrc  Tsrc   runoff file name\n")
            txt.write("{}{}_{}_runoff.nc \n".format(self.locacion, self.name, n))
            txt.write("{} \n".format(len(eval("self.runoff{}[1]".format(n)))))
            mouth = eval("self.runoff{}[1]".format(n))
            for s in range(len(eval("self.runoff{}[1]".format(n)))):
                txt.write("{} {} {} {} 30*T 5. 0.\n".format(mouth[s][0], mouth[s][1], mouth[s][2], mouth[s][3]))
        txt.close()
        return "{}_{}.in".format(self.name, n)

    """        
    Condicion inicial
    """
    def make_init(self, n):
        grd = xr.open_dataset(eval("self.dom{}".format(n)))
        # Leo la condicion inicial pasada
        nombre_del_modelo_previo = (pd.to_datetime(self.locacion[5:]) - pd.DateOffset(days = 3)).strftime("%Y%m%d%H")
        locacion_previa = "prono{}/".format(nombre_del_modelo_previo[:-2])
        solucion_previa = "{}_{}_his.nc".format(nombre_del_modelo_previo, n)
        # Cargo las variables generales
        variables = {"spherical": ("one", [1]), "Vtransform": ("one", [1]),
                         "Vstretching": ("one", [1]), "theta_s": ("one", [6.]),
                         "theta_b": ("one", [0.]), "Tcline": ("one", [0.5]),
                         "hc": ("one", [0.5]), "sc_r": ("one", [-0.5]),
                         "Cs_r": ("one", [-0.0497]), "tstart": ("one", [0]),
                         "tend": ("one", [0]), "ocean_time": ("one", [0]),
                         "scrum_time": ("one", [0]),}
        try:
            inicio = xr.open_dataset(locacion_previa + solucion_previa).sel(time = nombre_del_modelo_previo + "00")
            variables.update({"ubar": (["time", "eta_u", "xi_u"], inicio.ubar.values[np.newaxis, :, :]),
                              "vbar": (["time", "eta_v", "xi_v"], inicio.vbar.values[np.newaxis, :, :]),
                              "zeta": (["time", "eta_rho", "xi_rho"], inicio.zeta.values[np.newaxis, :, :]),})
        except:
            variables.update({"ubar": (["time", "eta_u", "xi_u"], 0 * grd.mask_u.values[np.newaxis, :, :]),
                              "vbar": (["time", "eta_v", "xi_v"], 0 * grd.mask_v.values[np.newaxis, :, :]),
                              "zeta": (["time", "eta_rho", "xi_rho"], 0 * grd.mask_rho.values[np.newaxis, :, :]),})
        ini = xr.Dataset(variables, coords = grd.coords)
        ini["tstart"].attrs = {"long_name": "start processing day","units": "day"}
        ini["tend"].attrs = {"long_name": "end processing day","units": "day"}
        ini["ocean_time"].attrs = {"long_name": "time since initialization","units": "second"}
        ini["scrum_time"].attrs = {"long_name": "time since initialization","units": "second"}
        ini.to_netcdf("{}{}_{}_ini.nc".format(self.locacion, self.name, n), 'w', format = 'NETCDF3_CLASSIC')
        ini.close()
        return

    """
    Descarga continental
    """
    def make_runoff(self, n):
        # Verifico si el valor de la descaga es un numero o un archivo
        if str(eval("self.runoff{}[0]".format(n))).isdigit():
            Qbar = eval("self.runoff[0]".format(n)) * np.ones(((self.date1 - self.date0).days, len(eval("self.runoff{}[1]".format(n)))))
            for i in range(Qbar.shape[0]):
                Qbar[i,0] *= (0.24 * np.tanh(i/12.))
                Qbar[i,1] *= (0.56 * np.tanh(i/12.))
                Qbar[i,2] *= (0.20 * np.tanh(i/12.))
        else:
            tabla = pd.read_csv(eval("self.runoff{}[0]".format(n)), index_col = "TIEMPO")
            tabla.index = pd.to_datetime(tabla.index)
            sub = tabla.loc[self.date0 : self.date1]
            if len(sub) == 0:
                descarga_mensual_climatologica = {1 : 22668.0, 2 : 24295.0, 3 : 27650.0, 
                                                  4 : 29190.0, 5 : 30327.0, 6 : 29718.0, 
                                                  7 : 28614.0, 8 : 25428.0, 9 : 22380.0, 
                                                  10: 25140.0, 11: 26903.0, 12: 24324.0,}
                Qbar = descarga_mensual_climatologica[self.date1.month] * np.ones(((self.date1 - self.date0).days, len(eval("self.runoff{}[1]".format(n)))))
                Qbar[:,0] *= 0.24
                Qbar[:,1] *= 0.56
                Qbar[:,2] *= 0.20
            else:
                # Serie temporal de descarga
                Qbar = sub.copy()
        # Time vector
        qbar_time = np.arange(len(Qbar))
    	# Create netCDF
        variables = {"qbar_time": ("qbar_time", qbar_time),
                     "runoff_position": (["n_qbar", "two"], np.array(eval("self.runoff{}[1]".format(n)))[:,:2]),
                     "runoff_direction": (["n_qbar", "two"], np.array(eval("self.runoff{}[1]".format(n)))[:,2:]),
                     "Qbar": (["n_qbar", "qbar_time"], Qbar.T) }
        rio = xr.Dataset(variables)
        rio.qbar_time.attrs = {"long_name": "runoff time", "units":"days", "cycle_length": float(10950)}
        rio.runoff_position.attrs = {"long_name": "position of the runoff (by line) in the ROMS grid"}
        rio.runoff_direction.attrs = {"long_name": "direction/sense of the runoff (by line) in the ROMS grid"}
        rio.Qbar.attrs = {"long_name": "runoff discharge", "units": "m3 s-1"}
        rio.to_netcdf("{}{}_{}_runoff.nc".format(self.locacion, self.name, n), 'w', format='NETCDF3_CLASSIC')
        rio.close()
        return

    """
    Compilacion de los campos de GFS 
    """
    def compilacion_forzantes(self, fecha = None):
        if fecha == None:
            fecha = self.datef
            locacion = self.locacion
        else:
            locacion = "prono{}/".format(fecha[:-2])
            fecha = pd.to_datetime(fecha + "00")
        # Dominio de referencia
        grd = xr.open_dataset(self.dom0)
        # Nombre de los archivos
        archivos = []
        for hora_pronostico in range(0, 243, 3):
            archivos.append("gfs.t{hora}z.sfluxgrbf{hora_pronostico:03d}.grib2".format(hora = fecha.strftime("%H"),
                                                                                       hora_pronostico = hora_pronostico))
        # Composicion de los campos
        for a, arch in enumerate(archivos):
            try:
                if "descargar" in self.atmo_forcing:
                    dataset = xr.open_dataset(locacion + arch, engine = "cfgrib")
                else:
                    viento = xr.open_dataset(self.atmo_forcing + arch, engine = "cfgrib",
                                             backend_kwargs = {'filter_by_keys': {'layer': "heightAboveGround"}})["u10", "v10"]
                    presion = xr.open_dataset(self.atmo_forcing + arch, engine = "cfgrib",
                                             backend_kwargs = {'filter_by_keys': {'layer': "surface"}})["sp"]
                    dataset = xr.merge([viento, presion])
                    # Acomodo las longitudes, recorto en el area, quito las coordenadas no utilizadas, 
                    # renombre la coordenada tiempo e interpolo
                aux = dataset.assign_coords(longitude = dataset.longitude - 360) \
                             .sel(longitude = slice(grd.lon_rho[0,0], grd.lon_rho[0,-1]),
                                  latitude = slice(grd.lat_rho[0,0], grd.lat_rho[-1,0])).drop(["time","step"])
                if a == 0:
                    frc_atmo = aux.copy(deep = True)
                else:
                    frc_atmo = xr.concat([frc_atmo, aux], dim = "valid_time")
            except:
                pass
        frc_atmo.rename({"valid_time": "time"}).to_netcdf(locacion + "{}_atmo.nc".format(fecha.strftime("%Y%m%d%H")))
        os.system("rm {locacion}*.grib2 {locacion}*.idx".format(locacion = locacion))
        return

    """
    Descarga y compilacion de los forzantes
    """
    def descarga_gfs(self, fecha = None):
        if fecha == None:
            fecha = self.datef
            locacion = self.locacion
        else:
            locacion = "prono{}/".format(fecha[:-2])
            fecha = pd.to_datetime(fecha + "00")
        # Listo de archivos a descargar
        descargas, nombres = [], []
        for hora_pronostico in range(0, 243, 3):
            nombre = "gfs.t{hora}z.sfluxgrbf{hora_pronostico:03d}.grib2".format(hora = fecha.strftime("%H"),
                                                                                hora_pronostico = hora_pronostico)
            nombres.append(nombre)
            descargas.append("curl -s 'https://nomads.ncep.noaa.gov/cgi-bin/filter_gfs_sflux.pl?file={name}".format(name = nombre) + \
                             "&lev_10_m_above_ground=on&lev_surface=on&var_PRES=on&var_UGRD=on&var_VGRD=on" + \
                             "&subregion=&leftlon=-70&rightlon=-40&toplat=-18&bottomlat=-62" + \
                             "&dir=%2Fgfs.{fecha}%2F{hora}%2Fatmos'".format(fecha = fecha.strftime("%Y%m%d"),
                                                                            hora = fecha.strftime("%H")))
        ## Descargo de a 10 archivos a la vez
        pool = multiprocessing.Pool(processes = 20)
        pool.map(os.system, [desc + " --output {loc}{name}".format(loc = locacion, name = nomb) for desc, nomb in zip(descargas, nombres)])
        return

    """
    Archivo forzante
    """
    def make_forcing(self, n, **kargs):
        # Descarga GFS
        if n == 0:
            if "descargar" in self.atmo_forcing:
                # Descargo los campos
                self.descarga_gfs()
            # Compilo los campos
            self.compilacion_forzantes()
        # Pregunto si ya existen, sino las descargo
        rango_fechas_previas = pd.date_range(self.date0, self.datef, freq = "6H")
        pronosticos_previos = ["prono{}/{}_atmo.nc".format(fechas_previas.strftime("%Y%m%d"),
                                                           fechas_previas.strftime("%Y%m%d%H")) 
                               for fechas_previas in rango_fechas_previas]
        pronosticos_faltantes = [previo for previo in pronosticos_previos if ~os.path.isdir(previo)]
        print(pronosticos_faltantes)
        for previos in pronosticos_faltantes:
            os.system("mkdir {}".format(previos.split("/")[0]))
            self.descarga_gfs(fecha = previos.split("/")[1][:10])
            self.compilacion_forzantes(fecha = previos.split("/")[1][:10])
        # Compilacion de estados previos (control = 0) hasta el dia del pronostico
        frc_atmo = xr.open_dataset("prono{}/{}_atmo.nc".format(rango_fechas_previas[0].strftime("%Y%m%d"),
                                                               rango_fechas_previas[0].strftime("%Y%m%d%H"),)) \
                     .sel(time = slice(rango_fechas_previas[0],
                                       rango_fechas_previas[0] + pd.DateOffset(hours = 3)))
        for fechas_previas in rango_fechas_previas[1:-1]:
            frc_atmo = xr.concat([frc_atmo,
                                  xr.open_dataset("prono{}/{}_atmo.nc".format(fechas_previas.strftime("%Y%m%d"),
                                                                              fechas_previas.strftime("%Y%m%d%H"))) \
                                  .sel(time = slice(fechas_previas,
                                                    fechas_previas + pd.DateOffset(hours = 3)))],
                                  dim = "time")
        # Le agrego el forzante del dia del pronostico
        frc_atmo = xr.concat([frc_atmo,
                              xr.open_dataset("prono{}/{}_atmo.nc".format(self.datef.strftime("%Y%m%d"),
                                                                          self.datef.strftime("%Y%m%d%H")))],
                                  dim = "time")
        # Interpolacion al dominio
        grd = xr.open_dataset(eval("self.dom{}".format(n)))
        frc_atmo_interp = frc_atmo.interp(longitude = grd.lon_rho.values[0,:]) \
                                  .interp(latitude = grd.lat_rho.values[:,0])
        # Variables
        w10 = xr.ufuncs.hypot(frc_atmo_interp.u10, frc_atmo_interp.v10)
        cd = np.empty_like(w10.values)
        cd[w10.values < 5.] = 1.1e-3
        cd[w10.values >= 5] = 1.1e-3 + 0.06e-3 * w10.values[w10.values >= 5]
        sustr = (cd * 1.20 * w10 * frc_atmo_interp.u10).values
        svstr = (cd * 1.20 * w10 * frc_atmo_interp.v10).values
        frc_atmo_dt = (frc_atmo.time.diff(dim = 'time')//(3600e9)).astype('float').values
        frc_atmo_dt = np.append(frc_atmo_dt, frc_atmo_dt[-1])
        vector_tiempo = frc_atmo_dt.cumsum() - frc_atmo_dt/2
        vector_tiempo[0] *= -1
        variables = {"sms_time": ("sms_time", vector_tiempo/24.),
                     "sustr": (["sms_time", "eta_u", "xi_u"], 0.5 * (sustr[:,:,0:-1] + sustr[:,:,1:])),
                     "svstr": (["sms_time", "eta_v", "xi_v"], 0.5 * (svstr[:,0:-1,:] + svstr[:,1:,:])),
                     "Pair": (["sms_time", "eta_rho", "xi_rho"], frc_atmo_interp.sp)}
        # Dataset
        frc = xr.Dataset(variables, coords = grd.coords)
        # Attributes
        frc.sms_time.attrs = {"long_name": 'surface momentum stress time',
                              "units": 'days', "cycle_length": float(10950),}
        frc.sustr.attrs = {"long_name": 'surface u-momentum stress',
                           "units": 'Newton meter-2',}
        frc.svstr.attrs = {"long_name": 'surface v-momentum stress',
                           "units": 'Newton meter-2',}
        frc.Pair.attrs = {"long_name": 'sea level pressure',
                           "units": 'Pascal',}
        # Tides
        if hasattr(self, "tide_forcing"):
            if n in self.tide_forcing:
                tide_file = xr.open_dataset(self.tide_file)
                sshR = tide_file.ssh_r.interp(lon = grd.lon_rho, lat = grd.lat_rho).drop(["lon", "lat"])
                sshI = tide_file.ssh_i.interp(lon = grd.lon_rho, lat = grd.lat_rho).drop(["lon", "lat"])
                nodal_factors = tides.Tides(["m2", "s2", "n2", "k2", "k1", "o1", "p1", "q1"])
                nodal_factors.set_initial_time(self.date0)
                Eamp = xr.ufuncs.hypot(sshR, sshI)
                Ephase = np.mod(xr.ufuncs.arctan2(-sshI, sshR) * 180/np.pi, 360)
                for p in range(Eamp.shape[0]):
                    Eamp[p,:,:] *= nodal_factors.f[p]
                    Ephase[p,:,:] = np.mod((Ephase[p,:,:] - 180/np.pi*(nodal_factors.phi + nodal_factors.u)[p]), 360)
                frc["tide_period"] = ('tide_period', tide_file.periods)
                frc.tide_period.attrs = {"long_name": 'Tide angular period',
                                         "units": 'hours'}
                frc["tide_Ephase"] = (['tide_period','eta_rho','xi_rho'], Ephase)
                frc.tide_Ephase.attrs = {"lon_name": 'Tidal elevation phase angle',
                                         "units": 'degrees'}
                frc["tide_Eamp"] = (['tide_period','eta_rho','xi_rho'], Eamp)
                frc.tide_Eamp.attrs = {"lon_name": 'Tidal elevation amplitude',
                                       "units": 'meters'}
        # Save netCDF
        frc.to_netcdf("{}{}_{}_frc.nc".format(self.locacion, self.name, n), 'w', format='NETCDF3_64BIT')
        return

    """
    Condiciones de contorno
    """
    def make_boundary(self, n):
        grd = xr.open_dataset(eval("self.dom{}".format(n)))
        his = xr.open_dataset("{}{}_{}_his.nc".format(self.locacion, self.name, n-1))
        # Reacomodo de dimensiones
        for vari in ["rho", "u", "v"]:
            his["lat_{}".format(vari)] = ("lat_{}".format(vari), eval("his.lat_{}[:,0]".format(vari)))
            his["lon_{}".format(vari)] = ("lon_{}".format(vari), eval("his.lon_{}[0,:]".format(vari)))
        variables = {"spherical": ("one", [1]), "Vtransform": ("one", [1]),
                     "Vstretching": ("one", [1]), "theta_s": ("one", [6.]),
                     "theta_b": ("one", [0.]), "Tcline": ("one", [0.5]),
                     "hc": ("one", [0.5]), "sc_r": ("one", [-0.5]),
                     "sc_w": ("two", [-1, 0]), "Cs_r": ("two", [-1, 0]),
                     "Cs_r": ("one", [-0.0497]), "tstart": ("one", [0]),
                     "tend": ("one", [his.time.shape[0]/4]),
                     "bry_time": ("bry_time", np.arange(his.time.shape[0])/(24.)),
                     "zeta_time": ("zeta_time", np.arange(his.time.shape[0])/(24.)),
                     "v2d_time": ("v2d_time", np.arange(his.time.shape[0])/(24.))}
        bry = xr.Dataset(variables)
    	# Interpolation of the solutions
        his["zeta"] = (["time", "lat_rho", "lon_rho"], his.zeta)
        zeta = his.zeta.interp(lon_rho = grd.lon_rho, lat_rho = grd.lat_rho)
        his["ubar"] = (["time", "lat_u", "lon_u"], his.ubar)
        ubar = his.ubar.interp(lon_u = grd.lon_u, lat_u = grd.lat_u)
        his["vbar"] = (["time", "lat_v", "lon_v"], his.vbar)
        vbar = his.vbar.interp(lon_v = grd.lon_v, lat_v = grd.lat_v)
        his.close()
    	# Iteracion por cada direccion
        direction = {"north": ["xi", "-1,:"], "south": ["xi", "0,:"],
                     "west": ["eta", ":,0"], "east": ["eta", ":,-1"]}
        for var in ["zeta", "ubar", "vbar"]:
            for dire in direction:
                if "bar" in var:
                    bry["{}_{}".format(var, dire)] = (["v2d_time", "{}_{}".format(direction[dire][0], var[0])], 
                                                                    eval("{}[:,{}]".format(var, direction[dire][1])))
                else:
                    bry["{}_{}".format(var, dire)] = (["zeta_time", "{}_rho".format(direction[dire][0])], 
                                                                    eval("{}[:,{}]".format(var, direction[dire][1])))
        bry.to_netcdf("{}{}_{}_bry.nc".format(self.locacion, self.name, n), 'w', format='NETCDF3_CLASSIC')
        bry.close()
        grd.close()
        return

    """
    Ajuste de coordenadas
    """
    def coord_adjust(self, file):
        his = xr.open_dataset(file)
        # Time
        his["time"] = (["time"], pd.date_range(self.date0, periods = his.time.size, freq = "H"))
        # Quito el nivel medio
        # his["zeta"] -= his.zeta.mean("time")
        # Remuevo el spin up
        his = his.sel(time = slice(self.datef, self.date1))
        # Coords
        his["lat_rho"] = ("lat_rho", his.lat_rho[:,0])
        his["lon_rho"] = ("lon_rho", his.lon_rho[0,:])
        his["lon_u"] = ("lon_u", his.lon_u[0,:])
        his["lat_u"] = ("lat_u", his.lat_u[:,0])
        his["lon_v"] = ("lon_v", his.lon_v[0,:])
        his["lat_v"] = ("lat_v", his.lat_v[:,0])
        # Variables y reinterpolacion
        his["zeta"] = (["time", "lat_rho", "lon_rho"], his.zeta)
        his["ubar"] = (["time", "lat_u", "lon_u"], his.ubar)
        his["vbar"] = (["time", "lat_v", "lon_v"], his.vbar)
        # Batimetria
        his["h"] = (["lat_rho", "lon_rho"], his.h)
        # Limpio variables
        coordenadas = ["xi_rho", "xi_u", "eta_rho", "eta_v"]
        variables = ["spherical", "xl", "el", "f", "pm", "pn", "angle", "time_step", "scrum_time"]
        his.drop(coordenadas + variables).to_netcdf(file)
        his.close()
        return

    

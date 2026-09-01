================================================================================
 N_EventosCEvNS_NSIAr  --  pipeline de ARGON (clon de N_EventosCEvNS_NSIXe)
================================================================================

Proposito
---------
Comparacion Xe vs Ar de la sensibilidad a NSI con parametros de RED-100
(arXiv:2411.18641). Clon del pipeline de Xe. Todo lo que NO cambia entre
blancos se copio sin tocar; los cambios estan acotados a:

  >>> Referencia metodologica completa (matematica, pipeline etapa por
      etapa, formalismo NSI, registro de cambios Ar vs Xe):
      ../../doc/metodologia.tex   (compilar con  pdflatex metodologia.tex).


  constants.f90        A,Z,N,M de Ar-40; masa activa 62 kg (ref.[46]);
                       EEE = 0.99 (ref.[46] + ReD + DarkSide-50); eff_ROI=1;
                       actividad 39Ar (atm 1 Bq/kg, UAr 7.3e-4), Q_beta=565 keV,
                       Z_hija=19, ancla SB_ref46=4 (para chi2_bkg_nest).
  Tnr_to_e.f90         lee  datos/nest_Ar_218V_dense.txt  (yield NR, LArNEST
                       ANCLADO a ReD 2025 en 2.4-7.6 keV; T_nr<2 keV = extrapol.)
  chi2_ideal_nest.f90  *** COMPARACION IDEAL SIMETRICA Xe vs Ar ***  (ver abajo)
                       copia identica al de la carpeta de Xe salvo TAG='Ar'.
                       T_nr_max 6 -> 3.5 keV (ref.[46] + cinematica E_nu=8 MeV).
  chi2red100_nest.f90  variante "ventana fisica de Ar" (N_e = 1..5, §VII):
                       sensibilidad Asimov de conteo puro. T_nr_max = 3.5 keV.
  chi2_bkg_nest.f90    *** NUEVO: RED-100 §V para Ar CON FONDO DE 39Ar SIMULADO ***
                       senal CEvNS simulada + fondo beta de 39Ar simulado
                       (isotopo conocido, actividad publicada) + Asimov.
                       3 escenarios: (0) sin fondo, (1) 39Ar UAr, (2) 39Ar atm.
                       Salida: datos/sensib_bkg_Ar.dat. Todo en N_e (SIN PE).
  red100_nest.f90      solo rutas de salida  -> *_Ar.dat
  mainred100_nest.f90  solo rutas de salida + etiquetas
  red100PE.f90         BORRADO. Era sobrante del clon de Xe: el espectro en PE
  mod_detector.f90     BORRADO. solo existe porque RED-100 publico su residuo
                       ON-OFF de Xe en unidades de PE, y chi2.f90 (Xe) lo
                       necesita para comparar. Para Ar no hay dato de RED-100
                       que comparar, y SEG/sig1 serian de LXe (Ar centellea en
                       128 nm, necesita TPB, RED-100 no publica el SEG). Todo
                       el analisis de Ar vive en N_e.

  flux.f90, xsections_nest.f90, mod_stats.f90  =  copia exacta de Xe

  Compilar Ar (SIN mod_detector):
    MODS="constants.f90 flux.f90 xsections_nest.f90 mod_stats.f90 Tnr_to_e.f90"

  chi2.f90  NO se clono. Es RED-100 §VI (limite OBSERVADO): compara la
            simulacion en PE contra el residuo ON-OFF REAL de RED-100
            (A_90 ~ 111 xSM). Solo Xe; requiere el dato publicado. En esa misma
            corrida chi2.f90 (Xe) imprime ahora tambien RED-100 §V
            (sensibilidad esperada): A_90_esp = 1 + sqrt(2.706/S2) ~ 126,
            /sqrt(3) ~ 73 (reproduce la Tabla I de RED-100: SM2018 58, KI 90,
            DB 56, INR 64) -> valida el metodo §V que chi2_bkg_nest.f90 aplica
            a Ar con fondo simulado.


METODOLOGIA: COMPARACION IDEAL SIMETRICA  (chi2_ideal_nest.f90)
================================================================================
Objetivo: comparar Xe y Ar en igualdad de condiciones, no "Xe de laboratorio
contra Ar ideal". El corte N_e >= 4 de RED-100-Xe NO es sensibilidad de los
aparatos (RED-100 ve un electron unico de sobra, ganancia 27 PE/e-); es el
fondo de electron unico (SE), ruido instrumental. Un analisis ideal (sin fondo)
no tiene por que imponerlo.

Reglas de la comparacion:
  - MISMO CODIGO para los dos blancos. chi2_ideal_nest.f90 es fisica agnostica
    al blanco (usa N_Ge/Z_Ge/M_Ge/A_Ge/EEE de `constants`, lee la tabla NEST via
    mod_tnr_to_e, calcula el espectro CEvNS internamente). La copia de Xe y la
    de Ar difieren en UNA linea: el parametro TAG.
  - ESTADISTICA: Asimov de conteo puro.  dN_k = R_k ,  sigma_k = sqrt(N_k) ,
    chi2(A) = (1 - A)^2 * N_tot_ROI ,  A_best = 1 ,
    A_90 (1 gdl, Dchi2 = 2.706) = 1 + sqrt(2.706 / N_tot_ROI).
    SIN fondo, SIN sistematicos, SIN eficiencia de seleccion (eff_ROI = 1).
    Mide la RESPUESTA INTRINSECA del blanco (sigma_CEvNS . Q_W^2 . yield .
    umbral), NO la sensibilidad experimental alcanzable.
  - VENTANA SIMETRICA:
      * umbral inferior: se BARRE  N_e >= {1, 2, 3, 4}.  N_e >= 1 es el limite
        fundamental (no se detecta cero cuantos; bajo ~0.15 keV el yield de NEST
        cae a 0 solo por quenching de Lindhard). N_e >= 4 es el umbral
        conservador comun (el que un detector real puede defender del SE).
      * SIN corte superior:  NE_HI = 30.  Los retrocesos de Ar son mas
        energeticos (T_max ~3.4 keV a E_nu = 8 MeV vs ~1 keV en Xe) y su
        espectro de ionizacion llega a N_e ~ 30-40; cortar en 7 "por simetria"
        tiraria senal real de Ar.
  - FLUCTUACION F(T):  se corre con la F de NEST propia de cada blanco (columna
    3 de la tabla) y con F = 1 (Poisson) como sistematico.
  - NSI 2D (ipar = 5: eps_ee^dV vs eps_emu^dV):  la NSI entra solo por
      q_eff2(eps) = (Q_W + q_ee)^2 + q_emu^2 + q_etau^2 ,
      A_amp(eps)  = q_eff2 / Q_W^2 ,  Dchi2 = (1 - A_amp)^2 * N_tot_ROI.
    Se hace para N_e >= 1 y para N_e >= 4 (F de NEST).
  - El limite OBSERVADO real de Xe (rama chi2.f90, datos ON-OFF 2024:
    A_90 ~ 111 xSM con nuisance, ~107 sin el) se reporta APARTE, en su propia
    caja, NO en la misma tabla que la comparacion ideal. La distancia
    ideal -> real (~x235 en A_90 - 1) es el coste de fondos + cortes +
    sistematicos + ajuste a 1 histograma, y le aplicaria igual a un Ar real.
  - Ar solo es fisico asumiendo ARGON DEPLETADO (UAr): el 39Ar atmosferico
    (~1 Bq/kg) domina el ROI por ~1e6-1e7.

Resultado (F de NEST, 192 kg*dia)   [R_tot en ev/(kg dia), A_90 en xSM]
  N_e>=   R_tot Xe   R_tot Ar   Ar/Xe   A_90 Xe   A_90 Ar
    1      11.03      38.00      3.4     1.036     1.019
    2       2.14      27.17     12.7     1.081     1.023
    3       0.372     20.27     54       1.195     1.026
    4       0.063     15.37    243       1.472     1.030

  - El "x470" de comparar Xe@4..7 vs Ar@1..5 era ARTEFACTO de ventana.
  - La coherencia sola (Q_W^2 x atomos/kg ~3.8x a favor de Xe) NO manda: el Qy
    de Xe es ~0 en el umbral y sus retrocesos son blandos, asi que N_e >= 1 ya
    recorta ~medio espectro de Xe; el Ar pierde poco. Es el argumento del §VII
    a favor de Ar, cuantificado.
  - Xe vuelca ~80% de su senal en N_e = 1 y cae en picado; el Ar decae suave con
    cola hasta N_e ~ 20. Subir el umbral hunde a Xe (N_e>=4: conserva 0.6%)
    mucho mas que a Ar (~40%).
  - F de NEST vs F = 1: A_90 cambia < 0.5 % en AMBOS blancos.
  - NSI 2D: SENSIBILIDAD PROYECTADA (no exclusion). Con la estadistica Asimov y
    sin fondo/sistematicos, la proyeccion alcanza el ~99.8 % (Xe) / ~99.9 % (Ar)
    del plano |eps| <= 1: solo el anillo "punto ciego" (q_eff^2 = Q_W^2) queda
    dentro del margen +/- del SM. Motivo: la NSI de quark down entra x(Z+2N)
    (~208 Xe, ~62 Ar) frente a Q_W ~ -37/-11, asi que |eps| ~ 0.05 ya cambia la
    tasa > 2x. Es un LIMITE OPTIMISTA: el analisis REAL de Xe (A_90 ~ 111 xSM,
    y A_amp max en |eps|<=1 es ~74 < 111) NO distingue del SM NINGUN punto de la
    caja. El alcance real esta mucho mas cerca de "nada".


================================================================================
RED-100 §V PARA Ar: SENAL SIMULADA + FONDO DE 39Ar SIMULADO  (chi2_bkg_nest.f90)
================================================================================
RED-100 hace dos analisis: §VI (limite OBSERVADO) compara la simulacion contra el
residuo ON-OFF REAL -> A_90 ~ 111 xSM (solo Xe, `chi2.f90`); y §V (sensibilidad
ESPERADA), ANTES de usar el reactor ON: fondo + senal CEvNS simulada -> Asimov ->
A_90 en Dchi2 = 2.706 (su columna "sensibilidad" de la Tabla I: SM2018 58, KI 90,
DB 56, INR 64 xSM).

  Para Xe: §V sale del S2 que `chi2.f90` ya calcula.
           A_90_esp = 1 + sqrt(2.706/S2) ~ 126 ; /sqrt(3) ~ 73  (RED-100 ajusta
           3 histogramas, aqui 1). Del orden de su Tabla I -> VALIDA el metodo.
  Para Ar: `chi2_bkg_nest.f90`, mismo metodo §V, pero el fondo es SIMULADO
           (no hay medida de RED-100-Ar). Fondo dominante y especifico de Ar =
           decaimiento beta del 39Ar (39Ar -> 39K + e- + nubar; Q_beta 565 keV,
           beta permitido, actividad publicada). Todo en N_e (SIN PE).
             senal S(N_e): flujo x seccion eficaz x binomial, ROI 1..5.
             fondo B(N_e): dN/dT ~ F(Z,T) p_e E_e (Q-T)^2, plegado por el yield
                           ER (datos/nest_Ar_ER_218V.txt) + binomial + EEE,
                           normalizado por actividad x masa x tiempo.
             Dchi2(A) = sum_k (1-A)^2 S_k^2/(S_k + B_k + se_floor*E)
             A_90 = 1 + sqrt(2.706 / sum_k S_k^2/(S_k+B_k+...))
           3 escenarios: (0) sin fondo [= respuesta intrinseca, reproduce
           chi2red100_nest], (1) 39Ar UAr (7.3e-4 Bq/kg), (2) 39Ar atmosferico
           (1 Bq/kg).  Salida: datos/sensib_bkg_Ar.dat.

RESULTADO (192 kg*dia, N_e = 1..5, F de NEST):
  escenario           A_90 (xSM)     S/sqrt(B) a 62 kg*dia
  sin fondo            1.022          -
  39Ar UAr            1.022          ~ 2300   (39Ar despreciable en N_e<=5)
  39Ar atmosferico    1.025          ~ 60

VALIDACION vs ref.[46]: ref.[46] declara S/sqrt(B) ~ 4 a 62 kg*dia. Nuestro
39Ar por SOLAPE ESPECTRAL PURO da ~2300 (UAr) -> MUCHO mas optimista que ref.[46].
Motivo: llegar a N_e <= 5 con un retroceso ELECTRONICO pide E_er <~ 0.1 keV, la
cola extrema del espectro beta (+ extrapolacion del yield ER sub-100 eV). El
"~4" de ref.[46] refleja probablemente el fondo de APILAMIENTO DE ELECTRON UNICO
por encima de 4 e-, que ref.[46] deja EXPLICITAMENTE SIN RESOLVER ("requires
special experimental study"). `chi2_bkg_nest.f90` imprime el se_floor
(~660 ev/(kg dia) por bin) que reconcilia el modelo con el "~4" -> es el tamano
implicito de ese fondo no resuelto. Para una proyeccion conservadora tipo
ref.[46], correr con se_floor a ese valor.

QUE NO SE HACE: §VI para Ar (ajuste al residuo ON-OFF real) -> necesita reactor
+ detector. Se dice explicito en la tesis.


PARAMETROS DE Ar (ref.[46] = Physics 5, 492 (2023); ReD 2025 = arXiv:2510.16404)
------------------------------------------------------------------------------
 masa activa    62 kg (no 126) -- ref.[46], proximo montaje.
 EEE            0.99 -- ref.[46]: umbral de emision ~0.2 kV/cm en Ar vs ~1.8 en
                Xe -> extraccion ~100 %. ReD: 3.8 kV/cm -> "100 %". DS-50 >99.9 %.
 T_nr_max       3.5 keV (era 6.0) -- cinematica E_nu=8 MeV (3.44 keV) + ref.[46].
                Trunca cola E_nu 8-10 MeV (<1 %). Efecto en A_90 de la ventana
                fisica: N_tot_ROI -1.9 %, A_90 sin cambio (1.022).
 ROI en N_e     1..5 primario -- §VII de arXiv:2411.18641 "below five ionization
                electrons". INCONSISTENCIA: ref.[46] dice "less than four".
                Se documenta; se usa 1..5.
 Qy (NR)        datos/nest_Ar_218V_dense.txt: LArNEST ANCLADO a los 5 puntos
                MEDIDOS de ReD 2025 Tabla 1 (2.4-7.6 keV; cociente ReD/LArNEST
                1.00-1.08, media 1.04). Por debajo de ~2 keV (el ROI CEvNS de
                RED-100: N_e<=5 <-> T_nr ~ 0.1-1 keV) NO HAY MEDIDA: es LArNEST
                reescalado = EXTRAPOLACION DE MODELO. ReD a 200 V/cm (RED-100
                218; <1 %). Bondar et al. 2017 (JINST 12, C05010; 80/233 keV)
                citado por consistencia con ref.[46], NO usado.
 F(T) (NR)      LArNEST (sub-Poissoniano, F~0.11-0.15 en el ROI) baseline;
                F=1 (binomial pura) como sistematico. ReD: la fluctuacion de
                ionizacion de NR en Ar "remains poorly characterized" -> SIN
                medida; ReD bracketea entre "sin fluctuaciones" y "binomial pura".
 39Ar           ~1 Bq/kg (atmosferico, ref.[46]); ~7.3e-4 Bq/kg (UAr, DarkSide-50,
                ~x1400 menos). Q_beta = 565 keV, Z_hija = 19.
 SEG_Ar         SIN valor (RED-100 planea TPB ~0.1 mg/cm2 pero no publica la
                ganancia). No se usa: Ar no tiene rama PE (ver arriba).


Compilar (no hay Makefile; orden de modulos obligatorio)
-------------------------------------------------------
  MODS="constants.f90 flux.f90 xsections_nest.f90 mod_stats.f90 Tnr_to_e.f90"
  gfortran -O2 -o chi2_ideal        $MODS chi2_ideal_nest.f90     # comparacion ideal
  gfortran -O2 -o chi2nsi_ar        $MODS chi2red100_nest.f90     # variante N_e=1..5
  gfortran -O2 -o chi2_bkg          $MODS chi2_bkg_nest.f90       # §V con fondo 39Ar
  gfortran -O2 -o mainred100_ar     $MODS mainred100_nest.f90
  gfortran -O2 -o red100_ar         $MODS red100_nest.f90
  (red100PE.f90 y mod_detector.f90 BORRADOS: Ar no tiene rama PE)

  (idem en ../N_EventosCEvNS_NSIXe/ para chi2_ideal con TAG='Xe')

Post-proceso:
  python/nest_Ar.py               genera datos/nest_Ar_218V_dense.txt (yield NR
                                  anclado a ReD) y datos/nest_Ar_ER_218V.txt
                                  (yield ER, para el fondo de 39Ar)
  python/compare_ideal_XeAr.py    tabla + chequeo analitico del alcance NSI +
                                  6 figuras


Salidas en datos/
-----------------
  --- comparacion ideal simetrica (chi2_ideal_nest.f90) ---
  espectro_Ne_ideal_Ar.dat      N_e  R_bin_Fnest  R_bin_F1   [ev/(kg dia)]
  espectro_Ne_ideal_Xe.dat      idem para Xe
  sensib_ideal_Ar.dat           NE_LO F_mode R_tot frac Ntot@192 A_90(x1..x335)
  sensib_ideal_Xe.dat           idem para Xe
  chi2_nsi_2DAr_ideal.dat       eps_x eps_y Dchi2   (N_e>=1, F de NEST)
  chi2_nsi_2DAr_ideal_ne4.dat   eps_x eps_y Dchi2   (N_e>=4, F de NEST)
  chi2_nsi_2DXe_ideal[_ne4].dat idem para Xe
  fig_ideal_1_espectro_Ne.png            espectro R(N_e) de Xe y Ar
  fig_ideal_2_retencion_umbral.png       % de senal retenida vs umbral
  fig_ideal_3_sensibilidad_exposicion.png A_90 proyectado vs exposicion + Xe real
  fig_ideal_4_figura_merito.png          cociente de eventos Ar/Xe vs umbral
  fig_ideal_5_plano_NSI.png              plano NSI 2D (sensibilidad proyectada)
  fig_ideal_XeAr.png                     resumen 2x2

  --- variante "ventana fisica de Ar" (chi2red100_nest.f90) ---
  chi2_nsi_2DAr.dat        eps_x  eps_y  Delta_chi2   (barrido NSI 2D, N_e=1..5)
  nsi_configAr.txt         etiquetas de los ejes
  sensib_expo_Ar.dat       multiplicador  exposicion_kgd  A_90   (Asimov)

  --- §V con fondo de 39Ar (chi2_bkg_nest.f90) ---
  sensib_bkg_Ar.dat        escenario  exposicion_kgd  A_90  S_tot  B_tot  S/sqrtB
                           (escenario: sin_fondo | Ar39_UAr | Ar39_atmosferico)
  nest_Ar_ER_218V.txt      E_er[keV]  Qy_er[e-/keV]  F_er   (yield ER, para 39Ar)

  --- diagnostico ---
  espectro_continuoAr.dat  T_nr  Kop  Mue  Comb
  ionization_electrones_Ar.dat   N_e  creados  extraidos
  eventos_sm_ar.dat        N_e  eventos SM Asimov absolutos
================================================================================

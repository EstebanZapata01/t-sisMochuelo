#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Exporta a datos_tutor/ los datos de Xe para que un codigo independiente
pueda reproducir el espectro y comparar.

  datos_tutor/espectros_Xe.csv      tabla principal (Ne = 0..15)
  datos_tutor/insumos/*.csv         flujo, espectro de retroceso, tabla NEST

Todo sale de archivos del pipeline (nada tecleado a mano):
  simulacion : datos/ionization_electrones.dat            (red100_nest.f90)
  RED-100 Fig.3 (creados/extraidos): datos/red100_fig3_digitalizado.csv
  RED-100 Fig.6 (senal antes/despues de cortes): datos/validacion_fig3_fig6_Xe.dat
               (red100PE.f90; f6_b/f6_a digitalizados de la figura)
  eff_ROI = sim_despues/sim_antes de ese mismo volcado (el de constants.f90)
El flujo se evalua compilando un driver minimo contra flux.f90 sin modificarlo.
"""
import csv
import os
import subprocess
import tempfile
import numpy as np

ROOT = "/home/oem/Desktop/Unipamplona/Trabajo de grado/Códigos"
D = f"{ROOT}/datos"
OUT = f"{ROOT}/datos_tutor"
os.makedirs(f"{OUT}/insumos", exist_ok=True)

fmt = lambda v: "" if v is None else f"{v:.9e}"

# ------------------------------------------------------------------ tabla principal
sim = np.loadtxt(f"{D}/ionization_electrones.dat", comments="#")       # Ne, creados, extraidos
f3 = np.loadtxt(f"{D}/red100_fig3_digitalizado.csv", delimiter=",", skiprows=1)
roi = np.loadtxt(f"{D}/validacion_fig3_fig6_Xe.dat", comments="#")     # Ne sim_a sim_d fig3 f6_a f6_d

sim_cre = {int(r[0]): r[1] for r in sim}
sim_ext = {int(r[0]): r[2] for r in sim}
r_cre = {int(r[0]): r[1] for r in f3}
r_ext = {int(r[0]): r[2] for r in f3}
eff = {int(r[0]): round(r[2] / r[1], 3) for r in roi}   # valores de constants.f90 (el volcado tiene 7 cifras)
f6_antes = {int(r[0]): r[4] for r in roi}
f6_desp = {int(r[0]): r[5] for r in roi}

cols = ["Ne", "sim_creados", "sim_extraidos", "red100_creados", "red100_extraidos",
        "sim_extraidos_x_effROI", "red100_senal_antes_cortes", "red100_senal_despues_cortes",
        "effROI"]
with open(f"{OUT}/espectros_Xe.csv", "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f)
    w.writerow(cols)
    for n in range(0, 16):
        w.writerow([n,
                    fmt(sim_cre.get(n)), fmt(sim_ext.get(n)),
                    fmt(r_cre.get(n)), fmt(r_ext.get(n)),
                    fmt(sim_ext[n] * eff[n]) if n in eff else "",
                    fmt(f6_antes.get(n)), fmt(f6_desp.get(n)),
                    fmt(eff.get(n))])

# ------------------------------------------------------------------ insumos
# espectro de retroceso: T_nr, Kopeikin, Mueller, combinado (hibrido)
rec = np.loadtxt(f"{D}/espectro_continuoXe.dat", comments="#")
with open(f"{OUT}/insumos/retroceso_Xe.csv", "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f)
    w.writerow(["T_nr_keV", "dRdT_Kopeikin", "dRdT_Mueller", "dRdT_hibrido"])
    for r in rec:
        w.writerow([fmt(x) for x in r[:4]])

# tabla NEST: T_nr, Qy, F
nest = np.loadtxt(f"{D}/nest_218V_dense.txt", comments="#")
with open(f"{OUT}/insumos/nest_Xe_Qy_F.csv", "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f)
    w.writerow(["T_nr_keV", "Qy_e_por_keV", "F_Fano"])
    for r in nest:
        w.writerow([fmt(x) for x in r[:3]])

# flujo hibrido dPhi/dE_nu, evaluado con flux.f90 sin modificarlo
DRIVER = """
program dump_flux
  use constants, only: dp
  use flux, only: flujo_diferencial
  implicit none
  integer :: u, k
  open(newunit=u, file='__OUT__', status='replace')
  do k = 0, 2000
     write(u,'(F8.4, ES20.10)') 0.005_dp*k, flujo_diferencial(0.005_dp*k)
  end do
  close(u)
end program dump_flux
"""
xe = f"{ROOT}/FORTRAN90/N_EventosCEvNS_NSIXe"
with tempfile.TemporaryDirectory() as tmp:
    out = os.path.join(tmp, "flux.dat")
    src = os.path.join(tmp, "dump_flux.f90")
    open(src, "w").write(DRIVER.replace("__OUT__", out))
    exe = os.path.join(tmp, "dump_flux")
    subprocess.run(["gfortran", "-O2", "-ffree-line-length-none", "-o", exe,
                    f"{xe}/constants.f90", f"{xe}/flux.f90", src],
                   check=True, cwd=tmp, capture_output=True)
    subprocess.run([exe], check=True, cwd=tmp, capture_output=True)
    flux = np.loadtxt(out)
with open(f"{OUT}/insumos/flujo_nu_hibrido.csv", "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f)
    w.writerow(["E_nu_MeV", "dPhi_dE_por_cm2_s_MeV"])
    for r in flux:
        w.writerow([f"{r[0]:.4f}", fmt(r[1])])

print(f"  -> {OUT}/espectros_Xe.csv  y  {OUT}/insumos/")

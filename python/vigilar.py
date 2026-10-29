"""vigilar.py -- regenera las figuras de la tesis al guardar.

Cada segundo mira las fechas de python/*.py y datos/*.dat:
  - si cambia un script de figuras_tesis.txt, corre ese script;
  - si cambia un modulo (estilo_tesis, nsi_core...), corre los de la lista que lo importan;
  - si cambia un .dat (por ejemplo, tras correr un binario Fortran), corre toda la lista.
Copia cada PNG a la carpeta de figuras de la tesis y guarda lo que imprime
el script en datos/logs/<script>.txt.      Uso:  python3 vigilar.py
"""
import glob, os, shutil, subprocess, time

AQUI = os.path.dirname(os.path.abspath(__file__))
DATOS = os.path.join(AQUI, "..", "datos")
TESIS = os.path.join(AQUI, "..", "..", "Proyecto_LaTeX", "figuras_metodologia")
LISTA = [l.split() for l in open(os.path.join(AQUI, "figuras_tesis.txt"))
         if l.strip() and not l.startswith("#")]

def fechas():
    archivos = glob.glob(os.path.join(AQUI, "*.py")) + glob.glob(os.path.join(DATOS, "*.dat"))
    return {f: os.path.getmtime(f) for f in archivos}

def correr(script, figura):
    os.makedirs(os.path.join(DATOS, "logs"), exist_ok=True)
    with open(os.path.join(DATOS, "logs", script + ".txt"), "w") as log:
        r = subprocess.run(["python3", script + ".py"], cwd=AQUI, stdout=log,
                           stderr=subprocess.STDOUT, env=dict(os.environ, MPLBACKEND="Agg"))
    if r.returncode == 0:
        shutil.copy2(os.path.join(DATOS, figura), TESIS)
    print(time.strftime("%H:%M:%S"), "ok   " if r.returncode == 0 else "ERROR", script, flush=True)

antes = fechas()
print("Vigilando", len(LISTA), "figuras de la tesis... (Ctrl+C para salir)", flush=True)
while True:
    time.sleep(1)
    ahora = fechas()
    cambiados = [os.path.basename(f) for f in ahora if ahora[f] != antes.get(f)]
    if not cambiados:
        continue
    time.sleep(2)                          # deja terminar a quien esté escribiendo
    for script, figura in LISTA:
        codigo = open(os.path.join(AQUI, script + ".py")).read()
        if (script + ".py" in cambiados
                or any(c.endswith(".dat") for c in cambiados)
                or any(c.endswith(".py") and "from " + c[:-3] + " import" in codigo for c in cambiados)):
            correr(script, figura)
    antes = fechas()                       # sus propias salidas no lo vuelven a disparar

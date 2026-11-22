# Ricostruisce index.html unendo la pagina, il motore e la busta sigillata.
# Uso: python3 strumenti/costruisci.py  (dalla cartella principale del repository)
import os
radice = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
def leggi(p): return open(os.path.join(radice, p), encoding='utf-8').read()
pagina, motore, busta = leggi('sorgenti/pagina.html'), leggi('sorgenti/engine.js'), leggi('sorgenti/busta.js')
assert '</script' not in motore and '</script' not in busta
open(os.path.join(radice, 'index.html'), 'w', encoding='utf-8').write(pagina.replace('/*ENGINE*/', motore).replace('/*BUSTA*/', busta))
print('index.html ricostruito')

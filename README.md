# Quantum Lab: laboratorio vergine

Un laboratorio di chimica che parte da zero. Il motore conosce solo il numero di protoni di un atomo; tutto il resto lo deve scoprire con le leggi della fisica. I risultati vengono confrontati con una "busta sigillata" di valori misurati, che il motore non vede mai.

L'app gira interamente sul telefono, anche senza rete, e si installa sulla schermata Home.

## Livello 0: l'atomo

Dato in ingresso: Z. Il motore calcola:

- la configurazione elettronica di energia più bassa (quali gusci e sottolivelli sono occupati, con quale spin);
- l'energia di prima ionizzazione, come differenza di energia tra lo ione e l'atomo neutro;
- la distribuzione degli elettroni attorno al nucleo, dove i gusci compaiono da soli.

### Leggi impostate

- Equazione di Schrödinger per gli elettroni, nella forma di Kohn-Sham (unità atomiche, simmetria sferica).
- Attrazione di Coulomb tra nucleo ed elettroni e repulsione tra elettroni (potenziale di Hartree).
- Principio di Pauli: al massimo 2l+1 elettroni per sottolivello e per spin.
- Scambio e correlazione nell'approssimazione LSDA: scambio esatto del gas di elettroni e correlazione PW92, i cui parametri vengono da simulazioni quantistiche del gas di elettroni, non da misure su atomi.

Il motore non contiene tabelle di configurazioni, la regola di riempimento di Madelung o valori misurati. Non conosce ancora la relatività.

### Come cerca la configurazione

1. Parte da uno schermo statistico (Thomas-Fermi) ricavato solo da Z.
2. Riempimento autoconsistente: a ogni passo gli elettroni occupano i livelli più bassi del potenziale del momento.
3. Prova varie disposizioni intere vicine e, per ciascuna, risolve il problema fino all'autoconsistenza.
4. Sposta un elettrone alla volta (anche invertendone lo spin) e tiene la disposizione di energia più bassa, finché nessuno spostamento la migliora.

### Risultati contro la busta (NIST ASD, elementi 1–103)

| Zona | Configurazioni uguali | Errore medio sull'energia di ionizzazione |
|---|---|---|
| Periodi 1–4 | 32 su 36 | 0,27 eV |
| Periodo 5 | 15 su 18 | 0,28 eV |
| Periodo 6 | 19 su 32 | 0,61 eV |
| Periodo 7 | 3 su 17 | 1,52 eV |
| Totale | 69 su 103 | 0,58 eV |

Il motore ritrova da solo i gusci, la periodicità e anomalie reali come Cr 3d⁵ 4s¹ e Cu 3d¹⁰ 4s¹. Gli errori sono sistematici: nei metalli di transizione l'approssimazione LSDA favorisce un elettrone d in più, e negli elementi pesanti l'errore cresce perché manca la relatività.

Validazione del metodo: con la stessa approssimazione, le energie totali di He, Ne e Ar coincidono con i valori di riferimento LDA del NIST entro 0,2 mHa, una volta tolta la differenza nota tra le parametrizzazioni di correlazione PW92 e VWN.

## File

- `index.html`, `manifest.webmanifest`, `sw.js`, icone: l'app pubblicata.
- `sorgenti/engine.js`: il motore (unico ingresso: Z).
- `sorgenti/busta.js`: i valori misurati, usati solo dal valutatore.
- `sorgenti/pagina.html`: l'interfaccia.
- `strumenti/costruisci.py`: ricostruisce `index.html` dai sorgenti.
- `strumenti/verifica_riferimenti.js`: confronta le energie totali con i riferimenti LDA.
- `strumenti/confronta_busta.js`: calcola gli elementi da riga di comando e li confronta con la busta (es. `node strumenti/confronta_busta.js 1 36`).

## Prossimi passi

- Livello 0.2: correzioni relativistiche, per vedere se gli errori dei periodi 6 e 7 si riducono.
- Livello 1: molecole piccole del mondo H-C-N-O.

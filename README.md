# Quantum Lab: laboratorio vergine

Un laboratorio di chimica che parte da zero. Il motore conosce solo il numero di protoni di un atomo e deve scoprire tutto il resto con le leggi della fisica. I risultati vengono confrontati con una "busta sigillata" di valori misurati, che il motore non vede mai.

L'app gira interamente sul telefono, anche senza rete, e si installa sulla schermata Home.

## Livello 0: l'atomo (motore 0.2)

Dato in ingresso: Z. Il motore calcola:

- la configurazione elettronica di energia più bassa (quali gusci e sottolivelli sono occupati, con quale spin);
- l'energia di prima ionizzazione, come differenza di energia tra lo ione e l'atomo neutro;
- la distribuzione degli elettroni attorno al nucleo, dove i gusci compaiono da soli.

### Leggi impostate

- Equazione di Schrödinger per gli elettroni, nella forma di Kohn-Sham (unità atomiche, simmetria sferica).
- Attrazione di Coulomb tra nucleo ed elettroni e repulsione tra elettroni (potenziale di Hartree).
- Principio di Pauli: al massimo 2l+1 elettroni per sottolivello e per spin.
- Scambio e correlazione PBE (gradiente generalizzato, Perdew-Burke-Ernzerhof 1996): costruito solo da vincoli fisici esatti. La sua base LSDA (scambio esatto del gas di elettroni e correlazione PW92) viene da simulazioni quantistiche del gas di elettroni, non da misure su atomi.
- Relatività scalare (Koelling-Harmon): massa dell'elettrone che cresce con la velocità e termine di Darwin, senza accoppiamento spin-orbita. Unica costante in più: la velocità della luce, c = 137,036 in unità atomiche. I termini relativistici usano il potenziale del nucleo più quello medio degli elettroni; lo scambio-correlazione resta nell'equazione ma non nella massa relativistica, così il rumore numerico del gradiente non li disturba.

Il motore non contiene tabelle di configurazioni, la regola di riempimento di Madelung o valori misurati.

### Come cerca la configurazione

1. Parte da uno schermo statistico (Thomas-Fermi) ricavato solo da Z.
2. Riempimento autoconsistente: a ogni passo gli elettroni occupano i livelli più bassi del potenziale del momento.
3. Prova varie disposizioni intere vicine e, per ciascuna, risolve il problema fino all'autoconsistenza.
4. Sposta un elettrone alla volta (anche invertendone lo spin) e tiene la disposizione di energia più bassa, finché nessuno spostamento la migliora.

### Risultati contro la busta (NIST ASD, elementi 1–103)

| Zona | Motore 0.1 (LSDA, senza relatività) | Motore 0.2 (PBE + relatività) |
|---|---|---|
| Periodi 1–4 | 32 su 36, errore medio 0,27 eV | 32 su 36, errore medio 0,19 eV |
| Periodo 5 | 15 su 18, 0,28 eV | 15 su 18, 0,24 eV |
| Periodo 6 | 19 su 32, 0,61 eV | 27 su 32, 0,23 eV |
| Periodo 7 | 3 su 17, 1,52 eV | 13 su 17, 0,17 eV |
| Totale | 69 su 103, 0,58 eV | 87 su 103, 0,21 eV |

"Su" = configurazioni elettroniche uguali a quelle misurate; l'errore è sull'energia di prima ionizzazione. Con LSDA più relatività le configurazioni giuste sarebbero 89 e l'errore 0,32 eV: PBE perde lo zirconio e il lantanio ma dimezza quasi l'errore sulle energie.

Errori che restano, con la causa probabile:

- Ti, V, Co, Ni, Zr, Tc, Rh, W, Ir: il motore sposta un elettrone dal livello s al livello d. È un difetto noto delle approssimazioni locali e semilocali (l'elettrone interagisce un po' con se stesso); conta anche che il motore confronta energie di configurazione mediate, mentre la misura indica il singolo livello più basso.
- La, Ce, Gd, U, Np, Cm: competizione tra livelli f e d molto vicini in energia.
- Lr, e gli errori più grandi sulle energie di Bi, Rn, Tl, At: serve l'accoppiamento spin-orbita.

### Validazione del metodo

- Energie totali LDA di He, Ne e Ar: coincidono con i riferimenti NIST entro 0,2 mHa, tolta la differenza nota tra le parametrizzazioni di correlazione PW92 e VWN.
- Energie totali PBE di H, He, Ne e Ar: coincidono con i valori di letteratura entro 0,04 mHa.
- Ione idrogenoide relativistico: l'energia 1s coincide con la soluzione esatta di Dirac (Z = 40 e Z = 80).

## File

- `index.html`, `manifest.webmanifest`, `sw.js`, icone: l'app pubblicata.
- `sorgenti/engine.js`: il motore (unico ingresso: Z).
- `sorgenti/busta.js`: i valori misurati, usati solo dal valutatore.
- `sorgenti/pagina.html`: l'interfaccia.
- `strumenti/costruisci.py`: ricostruisce `index.html` dai sorgenti.
- `strumenti/verifica_riferimenti.js`: confronta le energie totali con i riferimenti LDA.
- `strumenti/confronta_busta.js`: calcola gli elementi da riga di comando e li confronta con la busta. Esempi: `node strumenti/confronta_busta.js 1 36` (motore 0.2), `node strumenti/confronta_busta.js 1 103 --xc lsda --norel` (motore 0.1).

## Prossimi passi

- Accoppiamento spin-orbita (equazione di Dirac completa per gli orbitali).
- Livello 1: molecole piccole del mondo H-C-N-O.

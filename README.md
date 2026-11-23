# Quantum Lab: laboratorio vergine

Un laboratorio di chimica che parte da zero. Il motore conosce solo il numero di protoni di un atomo e deve scoprire tutto il resto con le leggi della fisica. I risultati vengono confrontati con una "busta sigillata" di valori misurati, che il motore non vede mai.

L'app gira interamente sul telefono, anche senza rete, e si installa sulla schermata Home.

## Livello 0: l'atomo (motore 0.3)

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

- Livelli singoli (dalla 0.3): il calcolo sferico dà l'energia media di una configurazione, la misura registra il livello più basso. Nei gusci aperti il motore risolve esattamente (interazione di configurazioni) la repulsione tra gli elettroni e l'accoppiamento spin-orbita, con integrali di Slater e costante spin-orbita calcolati dai propri orbitali, senza fattori di scala. La correzione è la differenza tra il livello più basso e la media delle disposizioni con gli stessi spin del calcolo sferico. Per ogni configurazione si usa come riferimento la disposizione degli spin che il calcolo sferico trova più bassa.

Il motore non contiene tabelle di configurazioni, la regola di riempimento di Madelung o valori misurati.

### Come cerca la configurazione

1. Parte da uno schermo statistico (Thomas-Fermi) ricavato solo da Z.
2. Riempimento autoconsistente: a ogni passo gli elettroni occupano i livelli più bassi del potenziale del momento.
3. Prova varie disposizioni intere vicine e, per ciascuna, risolve il problema fino all'autoconsistenza.
4. Sposta un elettrone alla volta (anche invertendone lo spin) e tiene la disposizione di energia più bassa, finché nessuno spostamento la migliora.

### Risultati contro la busta (NIST ASD, elementi 1–103)

Configurazioni uguali a quelle misurate ed errore medio sull'energia di prima ionizzazione:

| Zona | 0.1: LSDA | 0.2: PBE + relatività | 0.3: + multipletti e spin-orbita |
|---|---|---|---|
| Periodi 1–4 | 32/36, 0,27 eV | 32/36, 0,19 eV | 31/36, 0,19 eV |
| Periodo 5 | 15/18, 0,28 eV | 15/18, 0,24 eV | 16/18, 0,19 eV |
| Periodo 6 | 19/32, 0,61 eV | 27/32, 0,23 eV | 28/32, 0,13 eV |
| Periodo 7 | 3/17, 1,52 eV | 13/17, 0,17 eV | 14/17, 0,10 eV |
| Totale | 69/103, 0,58 eV | 87/103, 0,21 eV | 89/103, 0,16 eV |

Con la 0.3 gli errori più grandi della 0.2 spariscono: bismuto da +1,24 a −0,23 eV, radon da +0,83 a −0,24 eV, tallio da −0,79 a −0,24 eV; rodio, iridio, uranio e nettunio passano alla configurazione giusta. Il ferro e il torio invece passano a quella sbagliata.

Variante provata e scartata: solo spin-orbita, senza la struttura dei multipletti (`--soonly`). Indovina 91 configurazioni ma con errore medio 0,18 eV e sbaglia di 1,2 eV il polonio: tratta gli elementi p intermedi come se lo spin-orbita dominasse, il che è fisicamente sbagliato.

Errori che restano, con la causa probabile:

- Ti, V, Fe, Co, Ni, Zr, Tc, W e gli errori più grandi sulle energie (Cu, Cr, Ni, Pt, Pd, circa 0,4–0,65 eV): il motore sposta un elettrone dal livello s al livello d. È il difetto delle approssimazioni semilocali per cui ogni elettrone interagisce un po' con se stesso.
- La, Ce, Gd, Th, Cm: competizione tra livelli f e d molto vicini.
- Lr: 6d e 7p quasi alla pari. L'energia di ionizzazione è giusta (4,94 contro 4,96 eV), la configurazione no: servirebbe l'equazione di Dirac completa, con orbitali diversi per p₁/₂ e p₃/₂.

### Validazione del metodo

- Energie totali LDA di He, Ne e Ar: coincidono con i riferimenti NIST entro 0,2 mHa, tolta la differenza nota tra le parametrizzazioni di correlazione PW92 e VWN.
- Energie totali PBE di H, He, Ne e Ar: coincidono con i valori di letteratura entro 0,04 mHa.
- Ione idrogenoide relativistico: l'energia 1s coincide con la soluzione esatta di Dirac (Z = 40 e Z = 80).
- Multipletti: con integrali radiali dati, il calcolo riproduce le energie dei termini di p², p³, d² e d³ dei manuali (formule di Slater-Condon e Racah) e le separazioni spin-orbita di p¹, d¹ e p⁵.

## File

- `index.html`, `manifest.webmanifest`, `sw.js`, icone: l'app pubblicata.
- `sorgenti/engine.js`: il motore (unico ingresso: Z).
- `sorgenti/busta.js`: i valori misurati, usati solo dal valutatore.
- `sorgenti/pagina.html`: l'interfaccia.
- `strumenti/costruisci.py`: ricostruisce `index.html` dai sorgenti.
- `strumenti/verifica_riferimenti.js`: confronta le energie totali con i riferimenti LDA.
- `strumenti/confronta_busta.js`: calcola gli elementi da riga di comando e li confronta con la busta. Esempi: `node strumenti/confronta_busta.js 1 36` (motore 0.3), `--nolevels` (motore 0.2), `--xc lsda --norel --nolevels` (motore 0.1), `--soonly` (variante con solo spin-orbita).

## Prossimi passi

- Correzione dell'interazione di ogni elettrone con se stesso (Perdew-Zunger) o funzionale ibrido PBE0, per i metalli di transizione.
- Equazione di Dirac completa (orbitali p₁/₂ e p₃/₂ distinti), per il laurenzio.
- Livello 1: molecole piccole del mondo H-C-N-O.

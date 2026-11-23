# Quantum Lab: laboratorio vergine

Un laboratorio di chimica che parte da zero. Il motore conosce solo il numero di protoni di un atomo e deve scoprire tutto il resto con le leggi della fisica. I risultati vengono confrontati con una "busta sigillata" di valori misurati, che il motore non vede mai.

L'app gira interamente sul telefono, anche senza rete, e si installa sulla schermata Home.

## Livello 0: l'atomo (motore 0.4)

Dato in ingresso: Z. Il motore calcola:

- la configurazione elettronica di energia più bassa (quali gusci e sottolivelli sono occupati, con quale spin);
- l'energia di prima ionizzazione, come differenza di energia tra lo ione e l'atomo neutro;
- la distribuzione degli elettroni attorno al nucleo, dove i gusci compaiono da soli.

### Leggi impostate

- Equazione di Schrödinger per gli elettroni, nella forma di Kohn-Sham (unità atomiche, simmetria sferica).
- Attrazione di Coulomb tra nucleo ed elettroni e repulsione tra elettroni (potenziale di Hartree).
- Principio di Pauli: al massimo 2l+1 elettroni per sottolivello e per spin.
- Scambio e correlazione PBE (gradiente generalizzato, Perdew-Burke-Ernzerhof 1996): costruito solo da vincoli fisici esatti. La sua base LSDA (scambio esatto del gas di elettroni e correlazione PW92) viene da simulazioni quantistiche del gas di elettroni, non da misure su atomi.
- Scambio esatto (dalla 0.4): ibrido PBE0, cioè E = E_PBE + ¼ (Ex_esatto − Ex_PBE). Nello scambio esatto nessun elettrone interagisce con se stesso; il quarto viene da un argomento di teoria delle perturbazioni (Perdew, Ernzerhof, Burke 1996), non da un adattamento ai dati. Si valuta sugli orbitali PBE di ogni configurazione (non autoconsistente). Nei gusci aperti lo scambio esatto è la media vera sulle disposizioni degli elettroni, coerente con il riferimento dei multipletti: la formula ingenua con occupazioni frazionarie penalizzava i gusci aperti e dava al nichel la configurazione 3d¹⁰. Scartata invece la correzione di Perdew-Zunger: la letteratura recente (Maniar et al., PNAS 2025) mostra che sbilancia gli elettroni s e d di circa 2 eV.
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

| Zona | 0.1: LSDA | 0.2: PBE + relatività | 0.3: + multipletti e spin-orbita | 0.4: + scambio esatto |
|---|---|---|---|---|
| Periodi 1–4 | 32/36, 0,27 eV | 32/36, 0,19 eV | 31/36, 0,19 eV | 35/36, 0,14 eV |
| Periodo 5 | 15/18, 0,28 eV | 15/18, 0,24 eV | 16/18, 0,19 eV | 16/18, 0,12 eV |
| Periodo 6 | 19/32, 0,61 eV | 27/32, 0,23 eV | 28/32, 0,13 eV | 29/32, 0,17 eV |
| Periodo 7 | 3/17, 1,52 eV | 13/17, 0,17 eV | 14/17, 0,10 eV | 14/17, 0,18 eV |
| Totale | 69/103, 0,58 eV | 87/103, 0,21 eV | 89/103, 0,16 eV | 94/103, 0,15 eV |
| Errore massimo | | 1,24 eV (Bi) | 0,65 eV (Cu) | 0,44 eV (Ta) |

Cosa ha cambiato ogni versione:

- 0.2, relatività: gli elementi pesanti (periodo 7 da 3 a 13 configurazioni giuste).
- 0.3, multipletti e spin-orbita: gli elementi 6p (bismuto da +1,24 a −0,23 eV, radon da +0,83 a −0,24 eV, tallio da −0,79 a −0,24 eV). Variante scartata: solo spin-orbita senza multipletti (`--soonly`), che sbaglia il polonio di 1,2 eV.
- 0.4, scambio esatto: i metalli di transizione. Titanio, ferro, cobalto e nichel passano alla configurazione giusta; rame da +0,65 a +0,22 eV, cromo da +0,63 a +0,27 eV, palladio da +0,43 a −0,01 eV. Negli elementi 5d e negli attinidi l'energia peggiora un po' (tantalio −0,44, renio −0,41 eV): lo scambio esatto non è ancora autoconsistente.

Errori che restano:

- V, Zr, Tc, W: due configurazioni vicine in energia (dᶰ⁺¹s¹ contro dᶰs²).
- La, Ce, U, Cm: competizione tra livelli f e d.
- Lr: 6d contro 7p; servirebbe l'equazione di Dirac completa, con orbitali p₁/₂ e p₃/₂ distinti.

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
- `strumenti/confronta_busta.js`: calcola gli elementi da riga di comando e li confronta con la busta. Esempi: `node strumenti/confronta_busta.js 1 36` (motore 0.4), `--nohybrid` (motore 0.3), `--nohybrid --nolevels` (motore 0.2), `--xc lsda --norel --nolevels --nohybrid` (motore 0.1), `--soonly` (variante con solo spin-orbita).

## Prossimi passi

- Scambio esatto autoconsistente (orbitali ottimizzati con il funzionale ibrido), per gli elementi 5d e gli attinidi.
- Equazione di Dirac completa (orbitali p₁/₂ e p₃/₂ distinti), per il laurenzio.
- Livello 1: molecole piccole del mondo H-C-N-O.

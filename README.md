# Quantum Lab: laboratorio vergine

Un laboratorio di chimica che parte da zero. Il motore conosce solo il numero di protoni di un atomo e deve scoprire tutto il resto con le leggi della fisica. I risultati vengono confrontati con una "busta sigillata" di valori misurati, che il motore non vede mai.

L'app gira interamente sul telefono, anche senza rete, e si installa sulla schermata Home.

## Livello 0: l'atomo (motore 0.5)

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
- Correlazione della coppia ns² (dalla 0.5): una coppia ns² si mescola con np² perché i livelli s e p sono vicini (correlazione "statica"). PBE la include in modo implicito, la parte di scambio sostituita da quello esatto la perde. Il motore la calcola dai propri orbitali con un'interazione di configurazioni tra ns² npᵏ e npᵏ⁺² (energie delle due configurazioni da calcoli autoconsistenti, repulsione esatta tra gli elettroni) e ne aggiunge la stessa frazione dello scambio esatto, un quarto. Verifica: per il berillio il mescolamento vale −0,038 Ha, contro circa −0,044 Ha del calcolo di riferimento a due configurazioni. Gli altri gusci aperti non partecipano a questo calcolo.
- Margine (dalla 0.5): per ogni atomo il motore riporta la distanza in energia tra la configurazione scelta e la seconda. Sotto 0,2 eV (circa il doppio dell'errore medio, soglia fissata prima di vedere i risultati) il verdetto è segnalato come incerto.
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

| Zona | 0.1 | 0.2 | 0.3 | 0.4 | 0.5 |
|---|---|---|---|---|---|
| Periodi 1–4 | 32/36, 0,27 | 32/36, 0,19 | 31/36, 0,19 | 35/36, 0,14 | 35/36, 0,11 |
| Periodo 5 | 15/18, 0,28 | 15/18, 0,24 | 16/18, 0,19 | 16/18, 0,12 | 18/18, 0,10 |
| Periodo 6 | 19/32, 0,61 | 27/32, 0,23 | 28/32, 0,13 | 29/32, 0,17 | 28/32, 0,12 |
| Periodo 7 | 3/17, 1,52 | 13/17, 0,17 | 14/17, 0,10 | 14/17, 0,18 | 15/17, 0,11 |
| Totale | 69/103, 0,58 | 87/103, 0,21 | 89/103, 0,16 | 94/103, 0,15 | 96/103, 0,11 |
| Errore massimo (eV) | | 1,24 (Bi) | 0,65 (Cu) | 0,44 (Ta) | 0,42 (Ni) |

Versioni: 0.1 LSDA senza relatività; 0.2 PBE e relatività scalare; 0.3 multipletti e spin-orbita; 0.4 scambio esatto (PBE0); 0.5 correlazione della coppia s².

Cosa ha cambiato ogni versione:

- 0.2, relatività: gli elementi pesanti (periodo 7 da 3 a 13 configurazioni giuste).
- 0.3, multipletti e spin-orbita: gli elementi 6p (bismuto da +1,24 a −0,23 eV). Variante scartata: solo spin-orbita senza multipletti (`--soonly`), che sbaglia il polonio di 1,2 eV.
- 0.4, scambio esatto: titanio, ferro, cobalto e nichel passano alla configurazione giusta; rame da +0,65 a +0,22 eV. Scartata la correzione di Perdew-Zunger, che secondo la letteratura sbilancia s e d di circa 2 eV.
- 0.5, coppia s²: zirconio, tecnezio e uranio passano alla configurazione giusta; berillio da −0,35 a −0,09 eV, manganese da −0,35 a −0,18, tantalio da −0,44 a −0,32. Il platino invece passa a quella sbagliata. Una prima versione applicava la correzione solo quando il guscio p era vuoto, e peggiorava alluminio, gallio, indio e tallio (correzione sullo ione ma non sull'atomo); la versione finale tratta allo stesso modo anche ns² npᵏ.

Errori che restano (7) e margine con cui il motore li ha scelti:

- Segnalati come incerti (margine sotto 0,2 eV): V (0,16), La (0,17), Ce (0,00), W (0,07), Pt (0,04). Sono casi in cui due configurazioni sono quasi alla pari.
- Sbagliati con margine ampio: Cm (0,28 eV, competizione 5f contro 6d) e Lr (0,40 eV, 6d contro 7p: servirebbe l'equazione di Dirac completa).

Tra le 96 configurazioni giuste, 4 sono segnalate come incerte (Zr, Tc, Th, U). In tutto quindi 9 verdetti incerti, di cui 5 sbagliati; dei 94 verdetti "sicuri", 92 sono giusti.

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
- `strumenti/confronta_busta.js`: calcola gli elementi da riga di comando e li confronta con la busta. Esempi: `node strumenti/confronta_busta.js 1 36` (motore 0.5), `--nopair` (motore 0.4), `--nopair --nohybrid` (motore 0.3), `--nopair --nohybrid --nolevels` (motore 0.2), `--xc lsda --norel --nolevels --nohybrid --nopair` (motore 0.1), `--soonly` (variante con solo spin-orbita).

## Prossimi passi

- Scambio esatto autoconsistente (orbitali ottimizzati con il funzionale ibrido).
- Correlazione della coppia s² insieme agli altri gusci aperti, per i casi quasi alla pari (V, W, Pt, La, Ce).
- Equazione di Dirac completa (orbitali p₁/₂ e p₃/₂ distinti), per il laurenzio.
- Livello 1: molecole piccole del mondo H-C-N-O.

# FrankWolfeSolver — Design document

Porting in SMS++ del framework Julia [FrankWolfe.jl](../FrankWolfe.jl) per gli
algoritmi tipo Frank-Wolfe (conditional gradient). L'obiettivo è implementare
quanto più possibile *verbatim* gli algoritmi di FrankWolfe.jl, adattati alle
strutture di SMS++ (Block tree, `C05Function`, `Solver`, `Solution`,
`Modification`).

Testo in italiano, identificatori tecnici in inglese, come da convenzioni del
primer (`SMS++-CONTEXT.md`, §13).

Stato: **design approvato sui punti principali; v1 = vanilla Frank-Wolfe.**

---

## 1. Cosa è e dove si attacca

`FrankWolfeSolver` è un `CDASolver` che si registra a un Block "padre" con
questa struttura:

- una `FRealObjective` la cui `Function` è una `C05Function` (la funzione di
  *linking*, di cui useremo la *diagonal linearization* = gradiente);
- un numero arbitrario di **sub-Block**, che contengono *tutte* le variabili e
  i vincoli (il padre non ha né variabili né vincoli propri);
- ogni sub-Block ha una `FRealObjective` con dentro una `LinearFunction`,
  `DQuadFunction` o `QuadFunction`.

A ciascun sub-Block `FrankWolfeSolver` registra un `:Solver` (l'**LMO** di quel
Block). I nomi dei solver vengono dai parametri di `Configuration`, esattamente
come in `LagrangianDualSolver` (LDS): se un sub-Block ha già un `:Solver`
registrato lo si usa così com'è.

La regione ammissibile complessiva è il **prodotto** delle regioni dei figli:

```
C = conv(C_1) x conv(C_2) x ... x conv(C_k)
```

dove `C_j` è la regione ammissibile del sub-Block `j` (anche con vincoli di
integralità: `FrankWolfeSolver` non li guarda; combinando soluzioni intere come
combinazione convessa minimizza di fatto sull'inviluppo convesso).

### Analogia con LagrangianDualSolver

| Aspetto | LagrangianDualSolver | FrankWolfeSolver |
|---|---|---|
| Linking nel padre | **vincoli** lineari | **obiettivo** `C05Function` |
| Cosa fa col linking | dualizza (moltiplicatori) | linearizza (gradiente) |
| Solver per i figli | uno per sub-Block | uno per sub-Block (= LMO) |
| Direzione "verbatim" da copiare | registrazione solver, cache `Configuration`, parametri BSC per-figlio, gestione `Modification` | idem |

Gran parte del *plumbing* (registrazione/deregistrazione dei solver dei figli,
cache di `Configuration`, parametri `str`/`vint` per i `BlockSolverConfig`
per-figlio, snapshot/ripristino della f.o. dei figli, gestione `Modification`)
viene copiato il più possibile verbatim da LDS.

---

## 2. L'obiettivo composito (unificazione delle tre semantiche)

Le tre semantiche possibili per come la f.o. dei figli entra nel problema sono
casi particolari di un unico schema di **Frank-Wolfe generalizzato / conditional
gradient composito**, in cui l'oracolo (= il solver del figlio) tratta
*esattamente* una parte separabile dell'obiettivo.

Obiettivo totale minimizzato:

```
F(x) = f_father(x) + Σ_j h_j(x_j),     h_j(y) = α·⟨c_j, y⟩ + β·q_j(y)
```

dove `c_j`, `q_j` sono la parte lineare e quadratica **originali** della f.o.
del figlio `j`, e `g = ∇f_father(x)`, ristretto al figlio `j` come `g_j`. La
parte smooth `f_father` viene linearizzata; le `h_j` restano esatte
nell'oracolo. Le tre opzioni sono i valori di `(α, β)`, selezionati dal
parametro intero `intLMOObj`:

| `intLMOObj` | nome | `(α,β)` | l'oracolo del figlio `j` risolve | tipo LMO |
|---|---|---|---|---|
| `0` | `LMOLinear` (default) | `(0,0)` | `min ⟨g_j, v_j⟩` | LP/MILP |
| `1` | `LMOQuad` | `(0,1)` | `min ⟨g_j, v_j⟩ + q_j(v_j)` | QP/MIQP |
| `2` | `LMOFull` | `(1,1)` | `min ⟨c_j + g_j, v_j⟩ + q_j(v_j)` | QP/MIQP |

(Il quarto combo `(1,0)` non è richiesto; lo si aggiunge solo se servisse.)

**Nota (semantica Block-tree)**: l'obiettivo totale di un Block tree in SMS++ è
la *somma* di tutte le f.o. dell'albero; è esattamente ciò che un `:MILPSolver`
sul padre minimizza (f.o. padre + Σ f.o. figli, ricorsivamente). Quindi
**`LMOFull` è il modo che coincide con l'obiettivo "vero" del Block tree** e con
ciò che calcola il `:MILPSolver` monolitico. `LMOLinear` (figli ignorati)
coincide col monolitico solo se i figli hanno f.o. nulla; `LMOQuad` solo se
hanno f.o. puramente quadratica (parte lineare nulla).

### Meccanica (snapshot + scatter + restore)

Sullo stile di `LagBFunction::cleanup_inner_objective` di LDS:

- in `set_Block`: **snapshot** della f.o. originale di ogni figlio (`c_j`, e la
  parte quadratica `q_j` se presente);
- ogni iterazione: scrivo i coefficienti lineari della f.o. del figlio a
  `α·c_j + g_j` via `LinearFunction::modify_coefficients` /
  `DQuadFunction::modify_linear_coefficients`. La parte quadratica non viene mai
  toccata: resta se `β=1`, va azzerata una volta sola in `set_Block` se `β=0`;
- in `set_Block(nullptr)`/destructor: **ripristino** la f.o. originale, così il
  Block resta pulito.

### Il gap si calcola allo stesso modo in tutti i modi

Il **gap di Frank-Wolfe generalizzato** collassa in una forma unica:

```
gap(x) = Σ_j [ M_j(x_j) − M_j(v_j) ]   ≥   F(x) − F*
```

dove `M_j` è la f.o. del figlio *così come è stata modificata*
(`⟨α·c_j+g_j, ·⟩ + β·q_j`): `M_j(v_j)` è semplicemente **il valore ottimo che il
solver del figlio restituisce**, e `M_j(x_j)` è la stessa f.o. valutata
nell'iterate corrente `x`. In modo `LMOLinear` si riduce esattamente a
`⟨g,x⟩ − ⟨g,v⟩`. Un'unica routine per il criterio d'arresto, indipendente dal
modo.

### Sottigliezza per Away-step / Blended Pairwise (post-v1)

Nei modi `LMOQuad`/`LMOFull` l'oracolo è un QP: la direzione *away* e
l'`active_set_argminmax` vanno definiti rispetto alla sola parte smooth
`∇f_father`, con `h_j` che entra solo nell'oracolo e nella line search. Il caso
classico pulito è `LMOLinear`. La line search esatta resta in forma chiusa se
*sia* il padre *sia* le `q_j` sono quadratiche (l'Hessiano totale è la somma).

### 2.bis I due problemi: valore all'iterata vs combinazione convessa (`intCvxComb`)

Quando i figli hanno una f.o. **non lineare** (`h_j` con `β=1`, modo `LMOFull`),
lo stesso schema FW può calcolare due quantità diverse — *due problemi diversi* —
a seconda di come la parte `h_j` entra nel **valore** e nel **gap**. È una
**feature qualificante** del solver, controllata da `intCvxComb`:

- **(P1) — valore all'iterata** (`intCvxComb = eObjAtX = 0`, default):
  ```
  F_P1(x) = f_father(x) + Σ_j h_j(x_j)
  ```
  si minimizza la f.o. composita *vera* sull'inviluppo convesso `conv(C)`: `h_j`
  è valutata **nell'iterate corrente** `x_j` (che è un punto interno di
  `conv(C_j)`, non un vertice). È il bound che si otterrebbe ottimizzando
  esattamente `F` sulla chiusura convessa della regione prodotto.

- **(P2) — combinazione convessa / Dantzig-Wolfe** (`intCvxComb = eObjCvxComb = 1`):
  ```
  F_P2(x) = f_father(x) + Σ_j ( Σ_k λ_{jk} h_j(v_{jk}) )
  ```
  `h_j` è la **combinazione convessa dei costi nei vertici** `v_{jk}` dell'active
  set (i `λ_{jk}` sono i pesi FW), cioè l'**inviluppo convesso** di `h_j` ristretto
  ai vertici generati. È esattamente il bound di **decomposizione di
  Dantzig-Wolfe / simplicial decomposition**, e — quando `conv(C_j)` è l'inviluppo
  convesso intero del figlio — coincide con il **bound di perspective/Perspective-
  Cut** (P/C) che la formulazione DP rilassata-continua cut-separata calcola.

**Relazione.** Per la disuguaglianza di Jensen (`h_j` convessa, `x_j = Σ_k λ_{jk}
v_{jk}`): `F_P2(x) ≥ F_P1(x)`, con **uguaglianza per figli lineari** (`q_j=0`):
in quel caso `h_j` è lineare, la combinazione convessa dei costi-vertice *è*
`h_j(x_j)`, e i due problemi coincidono. La differenza è puramente nel termine
non-lineare dei figli.

**Cosa cambia nell'algoritmo** (tutto il resto è identico — **stessa LMO**, stessi
vertici, stessa direzione):
- **bookkeeping del costo `f_ci`**: ogni atomo memorizza `c_i = Σ_j h_j(atom_j)`,
  il costo-figlio *nel vertice*; il termine-figlio del valore corrente è allora
  `cbar = Σ_i λ_i c_i` (P2) invece di `Σ_j h_j(x_j)` (P1). `cost_x = ⟨g..⟩-pezzo +
  (cvx ? cbar : Σ h_j(x_j))`.
- **line search**: in P2 il termine-figlio è **lineare in γ** lungo `x+γ(v−x)`
  (interpola `c_i` tra i due atomi), in P1 è **quadratico** (è `h_j` valutata in
  un punto che si muove). Quindi la line search esatta e il valore lungo il
  segmento differiscono, **e di conseguenza differiscono le iterate** generate.
- **gap / criterio d'arresto**: usano `cost_x` coerente col modo; in P2 il gate
  della line search esatta è rilassato (il modello è lineare-nei-pesi).

In sintesi: **stesso oracolo, stesso schema, due funzioni-valore diverse** → P2 ≥
P1, e P2 è il ponte formale tra Frank-Wolfe e Dantzig-Wolfe / Perspective-Cut.
Questo è ciò che ha permesso di validare `FrankWolfeSolver` (modo P2/`LMOFull`)
contro `MILPSolver`+DPForm+P/C su `ThermalUnitBlock`: stesso bound a meno di
tolleranza, anche nel caso di commitment frazionario in cui `P1 ≠ P2`.

---

## 3. Estrazione del gradiente e scatter sui figli

Il padre `C05Function` ha come *active variables* le `ColVariable` fisicamente
residenti nei sub-Block. In `set_Block` costruisco una mappa

```
gradient position p  →  (sub-Block j, coefficient index i in B_j's objective)
```

iterando le active variables del padre e localizzandone il Block proprietario
(`ColVariable::get_Block()` / risalita al sub-Block diretto del padre) e la
posizione nella f.o. lineare del figlio. Questa mappa è statica finché la
struttura non cambia (cfr. §8 Modification).

Ogni iterazione (nel **main thread**, serialmente, prima del fan-out):

1. scrivo l'iterate `x` nelle variabili del Block (`Solution::write`);
2. sul padre: `compute()` → `compute_new_linearization(true)` →
   `get_linearization_coefficients(buf)` (assumiamo che il padre non produca mai
   linearizzazioni *verticali*);
3. *scatter* di `buf` nei coefficienti lineari dei figli via la mappa e
   `modify_coefficients`, applicando il modo `(α,β)`;
4. (parallelo) eseguo gli LMO, raccolgo vertici `v_j` e valori `M_j(v_j)`.

Gestione `min`/`max`: il *sense* del padre deve coincidere con quello dei figli
(check in `set_Block`, come LDS). `FrankWolfeSolver` minimizza internamente; per
`eMax` si nega il gradiente / si gestisce il segno coerentemente con il sense
dei solver dei figli.

---

## 4. Lo schema di Frank-Wolfe (v1: vanilla)

Iterate `x` mantenuto come `Solution` del Block padre, ottenuta via
`Block::get_Solution()` (il tipo concreto dipende dal Block: non ci serve
saperlo). Atomi e active set arrivano con Away-step/BPCG (post-v1); per il
vanilla basta l'iterate corrente.

```
x ← punto iniziale (da un primo LMO con gradiente in un punto di partenza)
for t = 0, 1, 2, ... :
    scrivi x nel Block; compute() padre; g ← ∇f_father(x)          # §3
    scatter g sui figli (modo α,β)                                  # §2,§3
    parallel for j in 1..k:  v_j, M_j(v_j) ← LMO_j.compute()       # §5
    gap ← Σ_j [ M_j(x_j) − M_j(v_j) ]                               # §2
    if gap ≤ ε  →  STOP (kOK)
    d ← x − v            # v = (v_1,...,v_k) come Solution del Block padre
    γ ← line_search(...) ∈ [0,1]                                    # §6
    x ← (1−γ)·x + γ·v      # Solution::scale + Solution::sum
```

La combinazione `(1−γ)x + γv` sfrutta il supporto **nativo** di `Solution` alla
combinazione (lineare/convessa) — `scale`, `sum` — che ricorre automaticamente
nei sub-Block. Tutte le `Solution` hanno la stessa interfaccia e qui si
combinano sempre e solo `Solution` provenienti dallo **stesso** Block (il
padre), quindi compatibili tra loro qualunque sia il tipo concreto generato dal
Block: `FrankWolfeSolver` non assume mai un tipo concreto. È anche la base
"gratuita" dell'Away-step (post-v1): l'active set sarà un insieme di `Solution`
del padre con pesi `λ_i`.

---

### 4.bis Varianti con active set: Away-step e Blended Pairwise (fatto)

Parametro `intAlgorithm` ∈ {`AlgVanilla`=0, `AlgAwayStep`=1, `AlgBPCG`=2}.
Le varianti con active set mantengono l'insieme degli **atomi** (vertici =
`Solution` del padre) con pesi `λ_i`, `x = Σ λ_i a_i`. Per scegliere l'**away
vertex** serve `⟨∇F, a_i⟩` per ogni atomo: si scrive l'atomo nel Block e si
valuta `Σ_j M_j` (le f.o. modificate dei figli) — `eval_modified_objective()`.

- **Away vertex**: `a = argmax_i ⟨∇F,a_i⟩` (argmin se `eMax`), peso `λ_a`.
- **FW vertex**: `v` dall'LMO. `fw_gap = ⟨∇F, x−v⟩`, `away_gap = ⟨∇F, a−x⟩`
  (entrambi ≥0; `fw_gap` certifica l'ottimalità → criterio d'arresto).
- **Away-step**: se `fw_gap ≥ away_gap` → FW step (`d=x−v`, `γ_max=1`); altrimenti
  away step (`d=a−x`, `γ_max=λ_a/(1−λ_a)`); aggiornamento `x ← x − γd`. Drop step
  se `γ=γ_max` (l'atomo away esce). Pesi: FW `λ_i←(1−γ)λ_i`, `λ_v+=γ`; away
  `λ_i←(1+γ)λ_i`, `λ_a−=γ`.
- **BPCG (pairwise)**: `d=a−v`, `γ_max=λ_a`; `x ← x+γ(v−a)`; `λ_a−=γ`, `λ_v+=γ`.
- **Line search**: stessa esatta, `γ*=clamp(gd/(2Q), 0, γ_max)` con `gd=⟨∇F,d⟩`
  (= il gap della direzione) e `Q=½⟨d,Ad⟩` da `quad_form` (riusa la struttura
  quadratica cachata). Agnostic `min(2/(t+2), γ_max)` altrimenti.
- **Dedup**: ogni nuovo vertice `v` viene cercato fra gli atomi (`find_atom`, per
  valori delle variabili del padre, tol relativa); se già presente, si somma il
  peso invece di creare un duplicato.

**Active set limitato (`intMaxAtoms`)**: per problemi grossi, dove l'active set può
crescere e l'argmax `O(|aset|)` diventa il collo di bottiglia, si tiene un tetto.
`intMaxAtoms` (0 = illimitato). Quando si sfora, gli atomi **meno attivi**
(contatore `f_count` di iterazioni consecutive con peso > 0 minore; a parità, peso
minore) vengono **aggregati** in un unico atomo: `ā = Σ_{i∈E}(w_i/W) a_i`, peso
`W=Σ w_i`. Questo **preserva esattamente `x`** (analogo all'aggregazione dei metodi
bundle). L'aggregato è **indistinguibile** da un vertice (solo, non è estremo): può
essere a sua volta away-vertex, droppato, o ri-aggregato (un aggregato recente può
contenerne di precedenti), e conta nel budget. **Trade-off**: gli aggregati non
sono vertici → si indebolisce la convergenza lineare dell'away-step (regime
"bundle aggregato"), in cambio di memoria/costo-per-iterazione limitati.

**Argmax cachato (figli lineari)**: `⟨∇F(x), a_i⟩ = ⟨g(x), a_i⟩ + c_i`, con
`c_i = ⟨c, a_i⟩` **indipendente da x** (memorizzato per atomo, `Atom::f_ci`,
calcolato all'inserimento come `M_j(v)−⟨g,v⟩`, aggregato linearmente). Così
l'argmax è un **prodotto scalare** `⟨f_grad, a_i.f_val⟩ + c_i` — niente write nel
Block né `Function::compute` per atomo. Con figli quadratici (`!all_lin`) si torna a
scrivere+valutare. `a_val` (valori dell'away atom) si legge da `f_aset[a_idx].f_val`
(già memorizzato), niente capture.

Stato: implementato e validato (away-step + BPCG, dedup, active set limitato con
aggregazione, argmax cachato). 72 combo base + 32 con tetto (incl. `intMaxAtoms=5`)
+ 32 col caching: `|aset|` al cap, `x` esatto, converge. TODO perf (rimandato): gram
matrix completa (`active_set_quadratic`) → argmax `O(|aset|²)` invece di
`O(|aset|·G)`, utile a G grande se l'argmax domina.

---

## 5. LMO via i solver dei figli + parallelismo

### 5.1 LMO

**Acquisizione del `:Solver` LMO (scelta v1)**: per ogni sub-Block si usa il
`:Solver` **già registrato** ad esso (registrato dalla `BlockSolverConfig`
ricorsiva del tester, `-S`), all'indice `intLMOSlvr` (default 0) in
`get_registered_solvers()`. L'acquisizione è **lazy**, alla prima `compute()`
(non in `set_Block`), perché l'ordine di registrazione dei `:Solver` ai vari
Block non è garantito quando `FrankWolfeSolver` riceve `set_Block`. È
responsabilità della config che un `:Solver` adeguato sia registrato a ogni
figlio (altrimenti `compute()` tira eccezione). La self-registration dai
parametri propri di `FrankWolfeSolver` à la LDS (`str`/`vstr`/`vint` per
`BlockSolverConfig` per-figlio) è un'estensione v1.1.

Ogni sub-Block `j` ha un `:Solver` registrato (l'LMO). Eseguire l'LMO `j`:

1. (già fatto serialmente: scatter del gradiente nella f.o. del figlio);
2. `solver_j->compute()`;
3. `M_j(v_j) ← solver_j->get_var_value()` (valore ottimo della f.o. modificata);
4. il vertice `v_j` resta scritto nelle variabili del sub-Block; lo si legge in
   una `Solution` (via `Block::get_Solution`) quando serve (per `d`, per la
   combinazione). Le `Solution` dei singoli figli compongono la `Solution` del
   padre, che ricorre nei sub-Block.

### 5.2 Parallelismo: thread pool persistente bulk-synchronous

Il pattern è bulk-synchronous (barriera dopo i `k` LMO), più regolare di quello
on-demand di `ParallelBundleSolver`. Schema proposto (più efficiente di
`std::async` per-iterazione + active-wait di PBS):

- **Pool persistente** di `min(intMaxThread, k) − 1` worker thread, creato lazy
  alla prima `compute()`, **riusato** per tutte le iterazioni e tutte le
  `compute()`, distrutto in `set_Block(nullptr)`/destructor.
- **Dispatch ad atomic counter**: `std::atomic<Index> next; while( (i =
  next.fetch_add(1)) < k ) run_lmo(i);` → bilanciamento dinamico del carico
  (LMO con tempi diversi), senza partizionamento statico.
- **Il main thread partecipa** al lavoro. Con `intMaxThread = 1` (default) zero
  thread extra → path puramente sequenziale (debug, e nessun costo quando non
  serve).
- **Sincronizzazione con `condition_variable` + generation counter** (no
  `sleep`/busy-wait): barriera di dispatch (sveglia i worker) e di join (sveglia
  il main a completamento).
- **Risultati in array per-indice** (`v_value[j]`, `v_status[j]`,
  `v_exception[j]`): nessun lock; ogni worker scrive solo il suo slot; ogni LMO
  lavora su un Block distinto → nessuna contesa sui lock dei Block.
- **Eccezioni**: catturate per-slot come `std::exception_ptr` e ri-lanciate nel
  main dopo la barriera.

Precondizioni (documentate): ogni `:Solver` dei figli su un Block distinto;
lo scatter del gradiente (che emette `Modification`) avviene nel main *prima*
del fan-out; nessuno stato mutabile condiviso tra worker.

Parametri: `intMaxThread` (default 1).

**Implementato (versione semplice)**: visto che per i nostri casi `k` (= numero di
sub-Block/LMO) è piccolo (2–3), il costo di creare i thread per-chiamata è
trascurabile rispetto alla risoluzione degli LMO. Quindi `run_LMOs` parallelo è
load-balanced con **atomic counter** (`next.fetch_add`) su `min(intMaxThread, k)`
thread (il main partecipa), niente pool persistente. Ogni LMO gira sul suo Block
figlio **distinto** con l'identità `f_id` prestata al `:Solver` (`set_id(f_id)` in
`acquire_LMOs`), quindi `lock(f_id)` sul figlio già posseduto da `f_id` non
contende; nessuno stato mutabile condiviso (ogni LMO scrive solo nel suo figlio e
nel suo slot `SubBlockData`); eccezioni catturate per-slot (`SubBlockData::excp`) e
ri-lanciate nel main dopo il join. Validato: 48 run paralleli (`intMaxThread=4`,
tutti gli algoritmi, con ripetizioni) identici al seriale, nessuna race.

**Pool persistente (rimandato)**: lo schema con thread pool + condition_variable
sopra descritto serve solo se `k` è grande e gli LMO sono minuscoli (creazione
thread non più trascurabile). Da fare se mai servisse; idem eventuale
armonizzazione di `ParallelBundleSolver`.

---

## 6. Line search

Da `FrankWolfe.jl/src/linesearch.jl`, per v1:

- **Agnostic** (`2/(t+2)`, ovvero `l/(t+l)`): nessuna valutazione extra, default
  sicuro; convergenza `O(1/t)`.
- **Esatta per obiettivo quadratico**: se `f_father` è `QuadFunction`/
  `DQuadFunction` (più, nei modi `LMOQuad`/`LMOFull`, le `q_j` quadratiche),
  `F` è quadratica lungo il segmento e

  ```
  γ* = clamp( −⟨∇F(x), d⟩ / ⟨d, A d⟩ , 0, γ_max )
  ```

  con `A` Hessiano totale (`Quad/DQuadFunction` del padre + Σ `q_j`). Per
  `DQuadFunction` `A` è diagonale; per `QuadFunction` si usa la matrice
  (Eigen sparse, `get_matrix`).

Post-v1: `Adaptive`/`Backtracking`, `Shortstep`, `Goldenratio`.

`FrankWolfeSolver` controlla in `set_Block` se la f.o. del padre è
`Quad`/`DQuadFunction` per abilitare il path esatto (parametro
`intLineSearch` per la scelta; default = auto: esatta se quadratico,
altrimenti Agnostic).

---

## 7. Interfaccia Solver / CDASolver

Base class: **`CDASolver`** (come LDS).

- `compute(bool)`: il loop di §4. Lock del Block padre (pattern LDS:
  `is_owned_by` / `lock(f_id)`), `process_outstanding_Modification()`, loop,
  unlock. Return code: `kOK` (gap ≤ ε), `kStopIter`, `kStopTime`,
  `kLowPrecision`, `kInfeasible`/`kUnbounded` (propagati/tradotti dagli LMO).
- `get_var_solution(Configuration*)`: scrive l'iterate `x` (la combinazione
  convessa mantenuta) nelle variabili del Block padre.
- `has_dual_solution()`: `true` **iff** ogni LMO dei figli è un `CDASolver` e il
  suo `has_dual_solution()` è `true`. Altrimenti `false`.
- `get_dual_solution(Configuration*)`: **puro forwarding** — per ogni figlio,
  `dynamic_cast<CDASolver*>(solver_j)->get_dual_solution(...)`. A terminazione
  (idealmente esatta) i duali ottimi degli LMO sono moltiplicatori di Lagrange
  validi per il problema complessivo.
- `get_lb()`/`get_ub()`: per minimizzazione, `ub = F(x)` (primale ammissibile),
  `lb = F(x) − gap` (best visto); a terminazione `gap → 0` e coincidono. Segno
  invertito per massimizzazione.
- Eventi e parametri base di `ThinComputeInterface` gestiti come di norma.

---

### 7.bis Variabili dei figli: assunzione "leaf" e TODO nipoti

La mappa gradiente→figli e lo scatter assumono che ogni variabile *attiva* della
f.o. del padre appartenga a un **figlio diretto** del padre. In generale le
variabili di un sub-Block "appartengono" ricorsivamente ai suoi antenati, quindi
una f.o. del padre potrebbe dipendere da variabili definite in un **nipote** (un
figlio di un figlio). v1 **assume il caso semplice**: o i figli sono *leaf* (niente
sub-Block), oppure la f.o. del padre dipende solo dalle variabili dei figli
diretti.

**TODO** (non neutro): se la variabile è definita in un nipote, modificarne il
coefficiente nella f.o. di quel discendente (non del figlio diretto) potrebbe non
essere "gradito" al figlio intermedio. Serve, in analogia a
`LagBFunction::PushCostToOwner`, la possibilità di modificare il coefficiente nel
Block che *definisce* la variabile quando non è un figlio diretto. Rimandato.

---

### 7.ter Convergenza: garanzia solo nel caso smooth (f.o. padre)

`FrankWolfeSolver` linearizza la f.o. del padre. La **convergenza globale al
minimo è garantita solo se la f.o. del padre è differenziabile** (smooth con
gradiente Lipschitz): è il caso `[D]QuadFunction`. Se la f.o. del padre è
**nondifferenziabile** (es. una funzione poliedrale/piecewise-linear), Frank-Wolfe
con un (sub)gradiente **non ha garanzia di convergere all'ottimo** — è un risultato
noto: F-W non è il metodo del subgradiente (si "impegna" sul vertice LMO di *quel*
subgradiente, direzione potenzialmente di salita), e ridurre il passo a 0 (DSS)
**non** lo risolve. I metodi nonsmooth convergenti richiedono **smoothing**
(Moreau-Yosida / Nesterov), rate $O(1/\sqrt k)$ (tight); vedi i paper in
`papers/` (White 1993; survey Braun-Pokutta 2022 §3.2.3; Thekumparampil et al.
NeurIPS 2020 — *"FW fails to converge if subgradients are used instead of
gradients"*).

**Decisione**: NON implementiamo un metodo nonsmooth dedicato (nessuna
applicazione attuale, fix non minore). `FrankWolfeSolver` **non fa alcun
riferimento diretto** a `PolyhedralFunction` o ad altri tipi nonsmooth: tratta la
f.o. del padre genericamente come `C05Function` e ne usa la *diagonal
linearization*. Si **documenta esplicitamente** che la convergenza globale è
garantita solo per il caso smooth; con f.o. padre nonsmooth l'algoritmo gira
comunque (può oscillare/stallare, regime Kelley/cutting-plane).

Due cose si tengono comunque, gratis e corrette:
- il **gap di F-W resta un bound valido** anche nonsmooth:
  $\text{gap}_t=\langle g_t,x_t-v_t\rangle\ge f(x_t)-f^*$ per ogni convessa e ogni
  subgradiente $g_t$; quindi `get_lb()`/`get_ub()` sono sempre validi — solo che il
  bracket **può non chiudersi** se la f.o. padre è nonsmooth;
- **hook progettuale** (rimandato): un "father-gradient provider" intercambiabile
  + parametro di smoothing $\mu$ permetterebbe di aggiungere lo smoothing (per una
  poliedrale: softmax delle righe attive) se mai servisse, senza che il resto dello
  schema cambi.

---

## 8. Modification e warm-start (implementato)

`FrankWolfeSolver` gestisce **lazy** le `Modification` provenienti dagli inner
Block, mantenendo coerenti le strutture cache e **warm-startando** active set /
iterata tra una `compute()` e l'altra. La logica ricalca `LagBFunction`, ma è più
semplice perché `FrankWolfeSolver` è in cima alla catena: **consuma** le
Modification, non le ri-traduce in `FunctionMod` per un solver esterno.

**Meccanica.**
- Niente `inhibit` globale fra due `compute()`: le Modification esterne si
  **accodano** in `v_mod`.
- `process_modifications()`, a inizio `compute()`, drena la coda e la categorizza.
- `inhibit_Modification(true)` riusato come `play_dumb` **solo attorno
  all'algoritmo e al restore**: lo scatter (e il restore finale) cambiano le f.o.
  dei figli e generano Modification "self-inflicted" da ignorare; l'inhibit scarta
  solo quelle nuove in arrivo, **non** la coda esterna già accumulata.
- **Restore-at-end**: a fine `compute()` le f.o. dei figli tornano al `c0`
  originale (`restore_objectives`), così un cambio esterno fra due solve agisce su
  uno stato pulito e il re-snapshot di `c0` è corretto. A teardown il restore usa
  `eNoMod` (gli altri `:Solver` potrebbero non esserci più — vedi nota sul crash
  di teardown).

**Categorizzazione** (`guts_of_process_modifications` → bit-mask), che ricalca la
*relax test* di LagBFunction (si ricontrolla la fattibilità solo quando il cambio
può *restringere* la regione):
- **f.o. padre** (bit 1) → re-cache della struttura quadratica (`analyze_father`);
  atomi intatti (valori e costi non cambiano).
- **f.o. figlio** (bit 2) → re-snapshot `c0` (`snapshot_c0`) + ricalcolo dei costi
  `f_ci` degli atomi (`recompute_atom_costs`), che vengono tenuti.
- **struttura** (bit 4) — `C05FunctionModVars` sull'obiettivo, `NBModification`,
  *fissaggio* di variabile → re-analisi completa + drop dell'active set.
- **regione ammissibile** (bit 8) — `RowConstraintMod` (LHS/RHS) → check;
  `VariableMod` (segno/integralità), `ConstraintMod` (enforce), `BlockModAD`
  (del-var / add-constraint) → `0` se *rilassa*, `8`/`4` se *restringe*; `BlockMod`
  generico → check.
- **Fissaggio di variabili**: instradato al **reset** (bit 4), non a 8.
  `ColVariable::is_feasible()` verifica segno/integralità ma **non** il fissaggio
  (uno *stato*, non un vincolo), quindi un feasibility-check non scarterebbe gli
  atomi che lo violano → serve il reset.
- *Nota* (come nel TODO di LagBFunction): la direzione di un `RowConstraintMod`
  (es. RHS che cresce su `≤` → rilassa) **non** è deducibile senza il valore
  precedente, quindi quel caso si ricontrolla sempre.

**Parametro `intHandleMod`** (`eModReset` default / `eModFine`): `eModReset` fa
sempre il caso pessimo (re-analisi completa + drop dell'active set su qualunque
cambio, nessun warm-start); `eModFine` aggiorna solo la cache interessata e
**conserva il warm-start**.

**Warm-start** (in `compute_active_set` e `compute_vanilla`): se una `compute()`
precedente ha lasciato active set / iterata ancora validi (cosa che
`process_modifications` garantisce sui cambi di sola f.o.), si riparte da lì.
Vanilla non ha atomi (non può ricostruire `cbar` esatto), ma riusa comunque
l'iterata `f_x` reinizializzando `cbar = C_lit(f_x)` (esatto per figli lineari,
"lavato via" dai pesi open-loop altrimenti).

**`feasibility_check`** (eModFine, cambio di regione): scrive ogni atomo, chiama
`f_Block->is_feasible()`, scarta gli inammissibili e ricostruisce `f_x` dai
sopravvissuti (combinazione convessa di punti ammissibili → ammissibile); per
vanilla verifica direttamente `f_x` e lo tiene se ancora ammissibile, altrimenti
cold restart. **Guard sugli aggregati**: il `f_ci` di un atomo aggregato è la
combinazione convessa dei costi dei vertici fusi (persi), non `h()` nel punto; con
active set limitato (`intMaxAtoms>0`) e figli quadratici un cambio-f.o.-figlio non
lo può ricalcolare esattamente → in quel caso si scarta il warm-start.

**Propagazione dell'infeasibility**: se una LMO ritorna `kInfeasible` la regione
prodotto è vuota → `compute()` ritorna `kInfeasible` (`run_LMOs` setta
`f_lmo_infeas`, controllato nel loop e all'init).

---

## 9. Struttura del modulo e build

Layout che ricalca `LagrangianDualSolver/`:

```
FrankWolfeSolver/
├── include/FrankWolfeSolver.h
├── src/FrankWolfeSolver.cpp
├── obj/
├── makefile        # FWSlvOBJ/INC/H/LIB
├── makefile-s      # + deps (SMS++ escluso)
├── makefile-c      # + SMS++/lib/makefile-c
├── CMakeLists.txt
└── README.md  CHANGELOG.md  LICENSE
```

I tester **non** stanno in `FrankWolfeSolver/test`: stanno in `tests/` (vedi
§10), perché dipendono dai Block "base" e dai loro `:Solver`, che non sono
dipendenze strette di `FrankWolfeSolver`.

Factory: `SMSpp_insert_in_factory_h` (header), `SMSpp_insert_in_factory_cpp_0(
FrankWolfeSolver )` (cpp). Prefisso macro makefile: `FWSlv` (sullo stile
`StcBlk`, `SDDPBk`, `MILP`).

### Parametri (sistema a 8 punti, §7 del primer)

I nomi non hanno l'infix `FWSlv`: essendo nel ComputeConfig di un
`FrankWolfeSolver` è ovvio che siano suoi (a differenza di LDS, dove i nomi
disambiguano il redirect verso il suo unico inner Solver — meccanismo che qui
non si applica).

- `int` (v1, no forwarding): `intLMOObj` (0/1/2), `intLineSearch`
  (auto/agnostic/exact), `intLMOSlvr` (indice del `:Solver` già registrato a
  ogni figlio da usare come LMO, default 0). Più gli ereditati `intMaxThread`
  (default 1), `intMaxIter`.
- `dbl`: ereditati `dblRelAcc`/`dblAbsAcc` (tolleranza sul gap), `dblMaxTime`.
- v1.1: `intLMOSlvr` → **`vint_LMOSlvr`** (un indice per LMO; vettore più corto
  del numero di LMO → i restanti al default 0; vettore vuoto [default] → tutti a
  0). Più, à la LDS, `str`/`vstr`/`vint` per i `BlockSolverConfig` per-figlio e
  la cache di `Configuration`, per la self-registration degli LMO.

---

## 10. Test

**Posizione**: un'unica directory `tests/FrankWolfeSolver/` (scaffolding condiviso
in `fw_test_common.h`, vedi "Tester implementati"), ricalcando
`tests/LagrangianDualSolver_{MMCF,UC,Box}` — da cui si copia a man bassa
(`test.cpp`, `makefile`, config `*.txt`, batch). Sta in `tests/` perché dipende
dai Block base e dai loro `:Solver`, non dipendenze strette di
`FrankWolfeSolver`.

Ogni tester costruisce un Block padre con f.o. semplice
(`QuadFunction`/`DQuadFunction` o `PolyhedralFunction`) sopra uno o più Block
"base" dello stesso tipo:

- `MCFBlock`, `ThermalUnitBlock` (`UCBlock`), `BinaryKnapsackBlock`.

Così lo stesso Block è risolvibile da più `:Solver` e li si confronta.

**Validità del cross-check (importante)**: F-W minimizza `F` sull'**inviluppo
convesso** delle regioni dei figli; un `:MILPSolver` monolitico minimizza sulla
regione **vera**. I due ottimi coincidono solo quando i due insiemi danno lo
stesso valore:

| f.o. padre | figli | FW vs MILP monolitico | `SolverReading` |
|---|---|---|---|
| lineare / `PolyhedralFunction` | anche interi | uguali (vertici dell'inviluppo interi) | `Exact` (ma vedi caveat sotto) |
| quadratica (`DQuad`/`Quad`) | **continui** (es. `MCFBlock`) | uguali (conv = regione stessa) | `Exact` |
| quadratica | interi generici (Knapsack) | FW = rilassamento ≤ ottimo intero | `LowerBound` (bound, non uguaglianza) |
| quadratica | `ThermalUnitBlock` form. DP (P/C) | **uguali**: la DP descrive l'inviluppo convesso delle soluzioni intere | `Exact` |

Ordine dei test: **(1) `MCFBlock` continui + padre `DQuadFunction`** → modo
`LMOLinear`, line search esatta, cross-check `Exact` (primo test, fatto); (2)
`PolyhedralFunction` padre + figli interi — **caveat sciolto** (vedi §7.ter): `F`
nondifferenziabile ⇒ **nessuna garanzia di convergenza globale**, F-W gira ma può
non raggiungere l'ottimo; quindi il test (2) non è un cross-check `Exact` ma una
verifica empirica (e `get_lb`/`get_ub` danno un bracket valido che può non
chiudersi). Niente smoothing implementato. (3) `ThermalUnitBlock` con formulazione
DP (P/C) + padre quadratico → `Exact`, perché la DP caratterizza l'inviluppo
convesso delle soluzioni intere. Un passo per volta.

**Solver-agnostico, tutto via configurazione**: il codice C++ del tester **non
fa nessuna assunzione** su quali `:Solver` sono attaccati — né ai sub-Block, né
al Block padre. Si usa l'infrastruttura di `tests/common_utils.h`:

- `process_args` (CLI standard: instance, `-B` BlockConfig, `-S`
  BlockSolverConfig, `-c`/`-p` prefissi); `require_block_config` /
  `require_solver_config`;
- `b_config_Block` applica il `BlockConfig` (-B); `s_config_Block` applica il
  `BlockSolverConfig` (-S) che **registra i `:Solver`** (al padre e ai figli) e
  li libera a fine corsa;
- `SolveAll( block, classify, ref, tol )` gira **tutti** i `:Solver` registrati
  al padre, stampa la riga uniforme (`print_instance_line`) e fa il
  **cross-check** (`cross_check`). Per v1: classifier `exact_getter(
  ObjGetter::VarValue )` (sia `FrankWolfeSolver` sia il solver di confronto
  letti come ottimo Exact via `get_var_value()`).

I `:Solver` per i sub-Block (gli LMO) e il `:Solver` di confronto sul padre sono
scelti nei file di config. Il solver di confronto sul padre **deve** essere un
`:MILPSolver` (è l'unico che risolve il problema monolitico appiattito), ma è
responsabilità di chi scrive il file di config che il Solver sia corretto:
altrimenti si tira eccezione. `PolyhedralFunction` come f.o. del padre è un buon
caso non-quadratico (LMO puramente lineare, line search Agnostic).

### Tester implementati

Un'unica directory `tests/FrankWolfeSolver/`, con scaffolding condiviso in
`fw_test_common.h` (namespace `fwtest`: `build_father`, `make_father_objective`,
`generate_poly`, `collect_vars`, `rnd`/`pos`), e **due** tester che si dividono
sull'asse "indipendenza dal tipo di Block":

- **`test.cpp` (generico, Block-agnostico)** — costruisce un padre con f.o.
  `DQuad`/`Quad`/`Polyhedral` sopra `-k` figli (di qualunque tipo, da CLI/config)
  e cross-checka `FrankWolfeSolver` contro il solver di confronto. Con **`-M`**
  esegue round randomici di Modification **della sola f.o.** (padre, ed
  eventualmente del primo figlio), che si esprimono attraverso la `FRealObjective`
  astratta e quindi **non** richiedono di conoscere il tipo del Block. Esercita la
  branch *cambio-f.o.* di `process_modifications` (re-snapshot `c0`, ricalcolo
  `f_ci`) e il warm-start, sia vanilla che active-set, sia `eModReset` che
  `eModFine`.

- **`test-mcf.cpp` (MCF-specifico)** — esercita i cambi della **regione
  ammissibile**, che intrinsecamente richiedono di conoscere il tipo: costruisce
  `-k` `MCFBlock` sotto un padre `DQuad`, e in round randomici (stile
  `tests/MCF_MILP`) cambia **costi** (`chg_costs` → cambio f.o.), **capacità**
  (`chg_ucaps` → cambio regione) o **fissa/sfissa** un arco (chiusura/riapertura =
  `VariableMod` che restringe/rilassa), sempre con `eModBlck, eModBlck` e capacità
  intere pulite ≥ 1 (MCFSimplex è sensibile a capacità a molte cifre). Cross-check
  contro `:MILPSolver` a ogni round, inclusa la **propagazione dell'infeasibility**
  quando la chiusura di un arco rende l'MCF inammissibile.

Entrambi girano in `regression` (statico + round `-M`) e in CMake/ctest
(`FWS_test`, `FWS_mcf_test`); la matrice di validazione (varianti × `eModReset`/
`eModFine` × `DQuad`/`Quad`, e cost/cap/fix × seed × istanze) è quella riportata
nel README del tester.

### Note di implementazione / limiti v1 (compute)

Il loop vanilla e le varianti active-set sono implementati in `compute()` e
testati a runtime (vedi "Tester implementati"). Scelte/limiti v1:

- **Inizializzazione**: un primo LMO al punto corrente (variabili a 0 di
  default) dà il vertice iniziale `x0 = v0` (sempre ammissibile).
- **`f_value = F(x)`** calcolato come `f_father(x) + (mx_sum − ⟨g,x⟩)` (vale per
  tutti i modi; in `LMOLinear` si riduce a `f_father(x)`). `f_bound = F(x) ∓ gap`.
- **Modification**: gestione lazy completa con warm-start — vedi §8.
  L'`inhibit_Modification(true)` non è più globale in `set_Block` ma riusato come
  `play_dumb` solo attorno ad algoritmo e restore.
- **`LMOLinear` richiede figli `LinearFunction`** (niente parte quadratica da
  tenere nell'oracolo): se un figlio è DQuad/Quad in `LMOLinear`, `compute()`
  tira eccezione (usare `LMOQuad`/`LMOFull`). `LMOQuad`/`LMOFull` lasciano la
  parte quadratica del figlio intatta nell'oracolo.
- **Line search esatta** quando il padre è quadratico (`DQuadFunction` *o*
  `QuadFunction`) e i figli sono lineari (Hessiano totale = Hessiano del padre):
  `γ* = clamp(−gd/(2Q), 0, 1)` con `gd = mv_sum − mx_sum = ⟨∇F,d⟩` e
  `Q = ½⟨d,Ad⟩ = Σ_p a_p d_p² + Σ_(r,c) q_rc d_r d_c` (diagonale + off-diagonale).
  La struttura quadratica statica del padre (diagonale `a_p` + off-diagonale da
  `mat_nd`) è **cachata in `set_Block`**. Padre non-quadratico → `Agnostic`.
- **LMO sequenziali** (`run_LMOs`); il thread pool persistente (§5.2) è il passo
  successivo.
- Il vertice `v_j` viene materializzato nelle variabili del figlio via
  `lmo->get_var_solution()` dopo `lmo->compute()`.

---

## 11. Milestone

1. **v1 — vanilla Frank-Wolfe end-to-end** (questo doc):
   - scheletro modulo + factory + parametri base;
   - `set_Block`: validazione struttura, snapshot f.o. figli, mappa gradiente,
     registrazione LMO dei figli (verbatim LDS);
   - gradiente + scatter (modo `LMOLinear` prima, poi `LMOQuad`/`LMOFull`);
   - loop FW + gap + line search (Agnostic + esatta quadratica);
   - `get_var_solution` / `get_dual_solution` / `get_lb`-`get_ub`;
   - parallelismo (pool persistente; ma prima versione seriale per correttezza);
   - test in `tests/FrankWolfeSolver_<base>` (MCF/Thermal/Knapsack), config-driven,
     cross-check via `SolveAll` contro un `:MILPSolver` monolitico.
2. **v2** — active set + Away-step + Blended Pairwise (`afw.jl`,
   `blended_pairwise.jl`, `active_set.jl`, `active_set_quadratic.jl`).
3. **v3** — line search aggiuntive, LMO lazy/cached, trattamento incrementale
   fine delle `Modification`.

---

*Aggiornare man mano che le decisioni si consolidano.*

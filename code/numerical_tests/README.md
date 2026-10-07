# Numerical test on scaled analogues

This folder contains a **test, not part of the proof**. It checks that the links of the argument compose and that the
definitions are consistent, on small exact examples. It cannot check the constants of [Q] and [R], which are replaced
here by values measured on each packing. The test was run on 7 October 2026 on an earlier draft of the argument. Its
numbering differs from the paper:

| first draft | paper |
|---|---|
| Lemma 1 (tilt budget along one path) | Lemmas 3.1 and 3.2 |
| Lemma 2 (lateral drift) | Lemma 3.3 |
| Lemma 3 (deficit and quantization) | Lemma 3.5, with θ as in Definition 3.4 |
| Lemma 4 (unreachable squares) | Lemma 4.1 |
| Lemma 5 (many short chords) | Lemma 4.2 |
| Lemma 6 (shadow budget) | Lemma 4.3 |
| Lemma 7 (total of the terminations) | Lemma 4.4 (not tested here) |
| Theorem A | Theorem 5.1 (the draft stated (k − 8), the paper has the exact count (k − 6)) |
| Theorem B | Section 6, Theorems 1.1 and 1.2 (only the constants were recomputed here) |

The "items" of the programs are: 1 = Lemmas 1–2, 2 = Lemma 3, 3 = Lemma 4, 4 = Lemma 6, 5 = Theorem A with the chain of
Lemma 5, 6 = the exact closed-packing verifier, 7 = negative controls.

## The scaled analogue

**Flows.** The flow is the master flow F(θ, δ, k/2 − 1) of [Q], and the ceiling flow is the floor flow of the reflected
packing.
- k ∈ {14, 16}, δ ∈ {1/20, 1/10}, inclinations 0.01–0.34.
- θ is chosen so that the chain of Lemmas 1 and 3 guarantees deficit ≤ τ₁ (modes `chain` and `chainfp`). Mode `a+` puts
  θ just above the tilt of a family; it is a stress mode and a negative control.

**Arithmetic.** Every quantity of the flows is rational and computed exactly (`fractions.Fraction`): pieces, contact
points, gaps, heights, deficits, R1 comparisons, chords and line measures. Only inclinations (asin of a rational) use
mpmath, with 50 digits.

**Measured constants.**
- A_meas(β̄)·W := ∫_{ω<½} 2·min(maxincl(y), β̄) dy, the largest possible left side of [R, Lemma 4.13] on the packing.
- E is computed exactly on the lines of K.
- The losses are the measured losses of both flows.
- The bound 1.0001β of Lemma 4 is replaced by the exact w(β) = (1 − cos β) + sin β, since the angles here are large.

## Volume

- 2,560 closed packings in 11 families (`random`, `rows`, `column`, `deepcol`, `phase`, `lshape`, `rotgrid`, `rotband`,
  `stair`, `converge`, `wall`), each verified exactly in rational arithmetic.
- 6,480 (packing, θ-mode) records and 12,960 flows.
- 893,085 path pieces, 1,458,151 passes, 1,489,203 entry candidates and 311,040 exact line samples.
- Every process stayed below 40 MB; see `results/resource_log.txt`.

## Results

**No counterexample was found to any tested link.**

**Lemma 1.**
- The per-path inequality holds and reaches equality (ratio 1.0).
- The measured global bound is used up to at most 48 %.
- The pass intervals of one path were disjoint in all 1,458,151 cases.

**Lemma 3.**
- The identity q_y − m + D = Σ gᵢ cos φᵢ and the range q_y − m ∈ [−D, δ) hold exactly at all 1,489,203 entry candidates.
- D ≤ (θ/2)Σa is tight (ratio up to 0.99992).
- With the chain-consistent θ the largest deficit was 0.0505 = 0.34·τ₁.

**Lemma 4.**
- No entry on a type-Z square in any flow with deficit ≤ τ₁.
- When θ was deliberately too large (`stair` family, mode `a+`), entries on type-Z squares occurred, always with deficit > τ₁ (exact reproducer: `programs/tbx_repro.py`, output `results/repro_stair.json`).
- The geometric core of the proof holds in 45,935 instances, with margin down to 5.4·10⁻⁶.
- In the top-ramp case the term 1 − cos a is needed for large a. In the paper a ≤ 2·10⁻⁴, where (1 − cos a) + sin a ≤ 1.0001a.

**Lemma 6.**
- The global, integrated and pointwise forms all hold; the pointwise ratio reaches 0.98.
- The shadow identity of [Q, Theorem 6.9] holds exactly at all line samples.

**Theorem A (measured constants).**
- The inequality holds with large margins at this scale, because W/ν₀ > |𝓗|.
- Its last step (K ≠ ∅ with 1 − ν₀ − E > 0) could be exercised in only 12 packings (`rotband`).

**[R, Lemma 3.5]** holds exactly on all 83,068 elementary line intervals of 132 packings (`results/lemma35.json`).

**Negative controls.** Removing one hypothesis or one term made the corresponding check fail in every case where it can
act at this scale. For example:
- the pass intervals overlap if the vertical extent is used;
- Lemma 3 fails if D or the gap term is dropped;
- Lemma 6 fails if passed squares are included or Λ is dropped.

Two controls cannot act at this scale and were replaced; the details are in the summaries.

**Cross-checks of the tracer.** No disagreement with:
- an independent exact point tracer (`programs/tbx_point.py`: 6,946 points identical, 119 merge points);
- the earlier floating-point tracer of the k^{1/4} repository (termination measures agree to 5.2·10⁻¹³; `results/xcheck_summary.txt`).

**Constants of Theorem B.** c_∞ = 0.0628538… (with [R, Lemma 4.10]) and 0.0522975…, in interval arithmetic
(`results/constB.txt`).

## Limitations

- The analogues are small (k ≤ 16) and the angles are large (up to 0.34 instead of ≤ 1.5·10⁻⁶).
- The constants of [Q] and [R] are replaced by measured values, so the test checks the structure of the argument, not
  these constants.
- Several inequalities have structural slack at this scale.
- About 75 % of the floor measure dies at this scale (δ is large).
- Campaign A partly ran with an earlier version of the per-path statistics (0/0 recorded as "inf"; two counters added
  later). Campaigns B and C ran with the final code.

## Files

| Path | Content |
|---|---|
| `programs/tbx_core.py` | exact rational squares, the exact closed-packing verifier, line structure, the exact piecewise-affine flow tracer |
| `programs/tbx_point.py` | independent exact point tracer of free paths (for the cross-check) |
| `programs/tbx_gen.py` | generators of the 11 families of exact rational packings |
| `programs/tbx_checks.py` | items 1–5 and the negative controls |
| `programs/tbx_run.py` | driver of one campaign job (JSONL output) |
| `programs/run_campaign.ps1` | campaign runner (Windows PowerShell; at most 3 processes, memory guard) |
| `programs/tbx_summary.py` | aggregation of the campaign records |
| `programs/tbx_xcheck.py`, `programs/tbx_xsummary.py` | cross-checks against the point tracer and the earlier floating-point tracer |
| `programs/tbx_lemma35.py`, `programs/tbx_geo4.py` | exact check of [R, Lemma 3.5]; the geometric core of Lemma 4 |
| `programs/tbx_repro.py` | exact reproducer of the Lemma 4 control |
| `programs/tbx_constB.py` | the constants of Theorem B (interval arithmetic) |
| `programs/tbx_smoke.py`, `programs/tbx_probe_*.py` | smoke test and probes |
| `results/summary_*.txt` | summaries of campaigns A, B, C and of all three |
| `results/xcheck_summary.txt`, `results/crosscheck/` | cross-check summary and records |
| `results/lemma35.json`, `results/geo4.json`, `results/constB.txt`, `results/repro_stair.json` | side checks |
| `results/resource_log.txt` | memory and time log of the runs |
| `results/campaign_records.zip` | the raw campaign records (74 JSONL files) |

## Reproduction

Use Python 3.12 with mpmath; numpy is needed only for the comparison with the earlier floating-point tracer. Run from
`programs/` (outputs go to `programs/out/`).

```
python tbx_smoke.py stair 16 1 a+
python tbx_run.py out/one.jsonl deepcol 16 77 2 1/20 chainfp,chain,a+ 16 200
powershell -File run_campaign.ps1 -Tag B -MaxPar 3 -NPack 30 -Fams "random,rows,column,deepcol,phase,lshape,rotgrid,stair,converge,wall" -Thetas "chainfp,chain,a+" -SeedBase 2000
powershell -File run_campaign.ps1 -Tag C -MaxPar 2 -NPack 40 -Fams "rotband" -Ks "16" -Thetas "chainfp,a+" -SeedBase 3000
python tbx_summary.py "out/campA_*.jsonl" "out/campB_*.jsonl" "out/campC_*.jsonl"
python tbx_lemma35.py out/lemma35.json 6 5151
python tbx_geo4.py out/geo4.json 6 6161
python tbx_constB.py
python tbx_repro.py out/repro_stair.json stair 16 2031 1 1/20 a+ 3/10 15/100 0.08
```

- Campaign A used `run_campaign.ps1 -Tag A -NPack 40` with the default families and θ modes and seeds 1001–1032.
- The stair generator was fixed after campaign A, so campaign A's stair packings do not contain the target square.
- For `tbx_xcheck.py`, set the environment variable `ASMX_DIR` to the folder `code/numerical_tests/assembly/` of the
  k^{1/4} repository (it contains `asmx_core.py`).

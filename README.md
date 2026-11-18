# Kineshoe: simulated four-experiment protocol (ngspice)

This report runs the four experiments of the *Kineshoe* paper (Sec. 3.3) as circuit simulations from scratch. Each experiment's data comes from **ngspice 45.2** transient analyses of the full circuits in Fig. 1 and Fig. 2. Python only generates the gait-force stimuli, writes the netlists, runs ngspice in batch mode and post-processes the waveforms. The paper's protocol is kept throughout: the same variables (Table 4), the same procedures, three trials per point (mean ± SD), least-squares R², and the Eq. 6 deviation. The filled versions of the paper's Tables 6–9 are in Section 3.

> **These are simulated values, not measurements.** They can show what the physics and the component models predict. They should be used to plan the bench experiments, not to replace them. Every parameter that was assumed is listed in Section 2, and the main uncertainties are tested with sensitivity variants.

---

## 1. Key findings (short version)

| # | Paper assumption | Simulation result | Verdict |
|---|---|---|---|
| 1 | The optimal load lies in the 1–100 kΩ decade-box range (Exp. 1). | Energy per strike rises monotonically across 1–100 kΩ. The true optimum is **R_opt ≈ 3.45 MΩ**, giving **1.33 µJ/strike** (single disc, full body weight). This lies between \|Z_CP\| = 2.65 MΩ and π/(2ωC_P) = 4.17 MΩ. At 100 kΩ the disc delivers only 9 % of the maximum. | Sweep range must be extended to ~100 MΩ. The 10 MΩ scope probe also loads the disc. |
| 2 | Schottky bridges reduce the 2V_F dead zone, so they are the better choice (Eq. 3). | At piezo current levels (~µA), the 1N5819's **reverse leakage** (tens of µA, a ~1 kΩ shunt at low bias) dominates. With the vendor SPICE model, **η ≈ 0 %**, and it is still ≈ 0 % with leakage cut 10×. Only at 1/100 of the model leakage does η reach 44 %. A 1N4148 bridge reaches **η = 87 %**; a near-ideal rectifier reaches 96 %. | **Refuted** for a µA piezo source. Eq. 3 ignores leakage. |
| 3 | Per-disc bridges decouple the channels, so E_step ∝ N (Eq. 4). | Rectified E_step vs N is clearly **sub-linear**: E(8)/(8·E(1)) = 0.44–0.47 with R² ≈ 0.93. The cause is placement: forefoot and midfoot discs see less force. A small cross-coupling also appears (array/sum = 0.94–0.96 at N = 8), because the shared bus back-biases the weaker discs. | Partly supported, as the paper anticipated. |
| 4 | V_C rises linearly and is clamped at 5.1 V within 20 min (Exp. 3). | V_C is linear (R² = 1.000) but tiny. As designed (1N5819, 1 F): **V_C(20 min) ≈ 3 µV** at running. With a near-ideal rectifier: **1.9 mV**, which is **1.9 µJ** banked. The Zener never conducts. Reaching the converter's V_min = 1.8 V would take **~490 h of continuous running**. Run/walk slope ratio = 3.8 > 2, so force dominates over strike count, consistent with Eq. 2. | Linearity is supported. The charging target is **off by roughly 10⁴–10⁵×**. |
| 5 | Eq. 6 predicts the massage runtime (Exp. 4). | With a 1 Ω-ESR supercap, t_m = **11.3 / 27.0 / 45.5 s** from 3 / 4 / 5 V at P_M ≈ 208 mW. Eq. 6 is within **−5 % to −12 %**, and the single fitted **η_c = 87.6 %**. With a 30 Ω coin-type 1 F/5.5 V cell, the converter enters UVLO hiccup (thousands of restarts). The motor never starts from 3 V, and **about 3 V (≈ 4.5 J) stays stranded**. | Eq. 6 is supported if ESR is low. Coin cells break it. |

**Bottom line.** The downstream half of the power chain (supercap → buck-boost → ERM) works and agrees with Eq. 6. The upstream half (discs → bridges → 1 F) delivers about **0.2–2 µJ per 20-minute run**, while a 45 s massage needs **12.5 J**. The 1N5819 bridges make this worse by shunting the discs. Section 5 lists the changes that matter.

---

## 2. Simulation set-up

### 2.1 Toolchain
* **ngspice 45.2** (official 64-bit Windows build, `tools/ngspice`). It provides transient analysis with the Gear integrator, XSPICE `filesource` stimuli, behavioural B-sources and hysteretic switches.
* Python 3.15 with numpy, scipy and matplotlib, used for stimulus generation, netlist assembly and analysis.
* Run everything with `python sim/run_all.py`, or run `sim/exp1.py` … `sim/exp4.py` and then `sim/figures.py`. Total runtime is about 25 min on 14 cores.

### 2.2 Circuit models

| Block | ngspice model | Key parameters |
|---|---|---|
| 27 mm piezo disc (7BB-27-4 class) | Thevenin charge source: V = d₃₃F(t)/C_P behind C_P, with R_P across the terminals (this is Eq. 1 in circuit form) | d₃₃,eff = 350 pC/N, C_P = 20 nF, R_P = 100 MΩ |
| Bridge diodes | SPICE diode models | **1N5819**: Is = 31.7 µA, N = 1.373, Rs = 51 mΩ, Cjo = 190 pF (vendor/LTspice values). Variants with Is/10 and Is/100 cover the datasheet leakage spread (50 µA–1 mA max at rated V_R). **1N4148**: Is = 2.52 nA, N = 1.752. **Near-ideal** bound: Is = 1 pA, N = 0.25. |
| Zener D_Z | soft-knee breakdown model | BV = 5.1 V at 49 mA, NBV = 5 (≈10 µA at 4 V, matching the datasheet IR) |
| Supercap C_S | C + ESR + R_leak | 1 F, ESR = 1 Ω (radial EDLC) or 30 Ω (coin type); R_leak = 2.2 MΩ (≈2.3 µA at 5 V) |
| DMM | resistor | 10 MΩ, left connected during logging |
| Scope probe ×10 | resistor | 10 MΩ across the measured node |
| Buck-boost (Exp. 4A) | **switch-level** 4-switch synchronous buck-boost | 1 MHz, L = 10 µH (DCR 0.3 Ω), C_out = 22 µF, R_on = 0.2 Ω per switch, buck / buck-boost / boost modes. Duty is solved iteratively for V_out = 3.000 V. Gate-charge (0.5 nC/switch), overlap (t_r = t_f = 8 ns) and I_q = 50 µA losses are added analytically. |
| Buck-boost (Exp. 4B) | **averaged** model | P_in = P_out/η_c(V_in), using the η_c table from 4A. UVLO turns off at 1.8 V and on at 2.0 V. 10 µF input cap, 1.5 A input limit, 1 ms soft start, output disconnected when off. |
| ERM motor 3 V | electromechanical: R-L-back-EMF plus J, viscous b and Coulomb friction (an electrical analogue on node ω) | R = 24 Ω, L = 60 µH, K_e = K_t = 1.05 mV·s/rad, J = 1.5×10⁻⁹ kg·m², b = 4.65×10⁻⁸ N·m·s, T_c = 15 µN·m. Gives 70 mA / 12.5 krpm / 210 mW at 3 V. |
| Flyback diode D_F | 1N5819 | across the motor |

### 2.3 Mechanical stimulus (gait)

![gait forces](../figures/fig0_gait_forces.png)

* Subject: **65 kg** (the paper leaves this as [TODO]). Body weight = 638 N.
* **Plate rig** (Exp. 1, 2A): one foot steps on the disc at 1 Hz, with a 0.6 s stance and a walking double-hump vGRF (1.16 / 0.75 BW). The whole load passes through the disc, giving a peak of 780 ± 20 N.
* **In-shoe** (Exp. 2B, 3): the vGRF is walking (1.16 BW) or heel-strike running (impact peak + active peak: jog 2.3 BW, run 2.5 BW; contact 0.30 / 0.25 s). It is split over heel, midfoot and forefoot with time-varying weights (heel loads first, then forefoot). Each disc captures a fixed share of its region through its rigid spreader. Discs: 4 heel (H1–H4), 3 forefoot (F1–F3), 1 midfoot (M1). Peak per-disc forces while running are 65–323 N.
* **Cadence convention**: cadence f counts the steps of both feet, so **each shoe is struck at f/2** (0.75, 1.25, 1.5 Hz). This is a correction to the reading n = f·t in Eq. 4.
* **Variability**, per trial: body-load scale σ = 3 %. Per strike: amplitude σ = 5 %, contact time σ = 3 %, period σ = 1 %. Each trial uses a different random seed. Results are reported as the mean ± SD of the three trial means.

---

## 3. Results

### 3.1 Experiment 1: single-disc output and optimal load

![exp1](../figures/fig1_exp1_load_sweep.png)

#### Table 6 – Experiment 1: single-disc output vs. external load

(n = 3 trials × 6 strikes, mean ± SD; subject mass 65 kg; plate rig, 1 Hz; peak force 780 ± 20 N.  V_oc read through the ×10 probe = **8.54 ± 0.18 V** (unloaded: 9.45 V; Eq. 1 ideal: 13.04 V))

| R_L | V_oc (V) | V_R (V) | E_step (µJ) |
|---|---|---|---|
| 1 kΩ | 8.54 ± 0.18 | 0.0052 ± 0.0013 | 0.0013 ± 0.0001 |
| 10 kΩ | 8.54 ± 0.18 | 0.0326 ± 0.0026 | 0.0130 ± 0.0007 |
| 47 kΩ | 8.54 ± 0.18 | 0.1407 ± 0.0061 | 0.0605 ± 0.0031 |
| 100 kΩ | 8.54 ± 0.18 | 0.2937 ± 0.0112 | 0.1262 ± 0.0064 |
| 1.00 MΩ *(extension)* | 8.54 ± 0.18 | 2.5345 ± 0.0830 | 0.8979 ± 0.0431 |
| 3.16 MΩ *(extension)* | 8.54 ± 0.18 | 5.4317 ± 0.1527 | 1.3357 ± 0.0571 |
| 10.00 MΩ *(extension)* | 8.54 ± 0.18 | 7.5963 ± 0.1782 | 1.0081 ± 0.0397 |
| 100.00 MΩ *(extension)* | 8.54 ± 0.18 | 8.4756 ± 0.1776 | 0.1683 ± 0.0065 |

Simulated optimum: **R_opt = 3.45 MΩ, E_max = 1.33 µJ/strike**. Dominant frequency of the disc current 3.0 Hz → |Z_CP| = 1/(ωC_P) = 2.65 MΩ, rectified-source optimum π/(2ωC_P) = 4.17 MΩ.



**Discussion.** The curve has the expected peak for a capacitive source, but the peak sits two decades above the planned sweep. The strike is slow (about 3 Hz content), so |Z_CP| = 1/(2π·3 Hz·20 nF) ≈ 2.7 MΩ, and the measured optimum lies between the impedance-match and rectified-source predictions. Over 1–100 kΩ the energy scales almost exactly ∝ R_L (slope 1.0 on log–log), so that sweep would report no maximum at all.

Two measurement artefacts also appear:
1. The **10 MΩ probe** makes R_P∥R_probe·C_P = 0.18 s, so the "open-circuit" reading is 8.5 V instead of 9.4 V.
2. **Disc leakage** (R_P·C_P = 2 s) leaves a negative residual after each release. The sustained V_oc is therefore well below the Eq. 1 value of 13 V.

**Recommendation:** extend the decade box to 100 MΩ, or measure with a 1 GΩ electrometer or a high-impedance buffer.

### 3.2 Experiment 2: rectifier efficiency and array scaling

![exp2](../figures/fig2_exp2_rectifier_array.png)

#### Table 7a – Experiment 2A: bridge efficiency, single disc on plate at R_opt

E_raw = 1.470 ± 0.128 µJ per strike (no bridge, R_opt = 3.45 MΩ)

| Rectifier | E_rect (µJ) | η (%) | peak V_bus (V) |
|---|---|---|---|
| 1N5819 Schottky (paper BOM) | 2.171e-10 ± 4.40e-11 | 0.00 ± 0.00 | 0.00158 |
| 1N5819, leakage /10 | 2.661e-6 ± 5.98e-7 | 0.00 ± 0.00 | 0.0165 |
| 1N5819, leakage /100 | 0.6486 ± 0.0861 | 44.00 ± 2.01 | 4.9 |
| 1N4148 silicon | 1.286 ± 0.118 | 87.48 ± 0.42 | 6.15 |
| near-ideal rectifier (bound) | 1.406 ± 0.125 | 95.64 ± 0.16 | 6.29 |


#### Table 7 – Experiment 2B: rectification and number of active discs (in-shoe, walking 1.5 Hz)

**1N5819 Schottky (paper BOM)** — linear fit E_rect = -0.00·N +0.00 nJ, R² = 0.025; linearity index E(8)/(8·E(1)) = 0.14

| N | E_step,raw (nJ) | E_step,rect (nJ) | η (%) | E_rect(array)/ΣE_rect(single) |
|---|---|---|---|---|
| 1 | 45.83 ± 0.48 | 1.382e-10 ± 5.91e-12 | 0.00 ± 0.00 | n/a |
| 2 | 77.66 ± 0.82 | 1.983e-10 ± 8.49e-12 | 0.00 ± 0.00 | n/a |
| 4 | 107.02 ± 1.13 | 1.883e-10 ± 8.06e-12 | 0.00 ± 0.00 | n/a |
| 8 | 185.20 ± 2.03 | 1.499e-10 ± 6.65e-12 | 0.00 ± 0.00 | n/a |

**1N4148 silicon** — linear fit E_rect = 9.74·N +29.60 nJ, R² = 0.931; linearity index E(8)/(8·E(1)) = 0.44

| N | E_step,raw (nJ) | E_step,rect (nJ) | η (%) | E_rect(array)/ΣE_rect(single) |
|---|---|---|---|---|
| 1 | 45.83 ± 0.48 | 29.44 ± 0.358 | 64.25 ± 0.22 | 1.000 |
| 2 | 77.66 ± 0.82 | 57.11 ± 0.691 | 73.55 ± 0.31 | 1.181 |
| 4 | 107.02 ± 1.13 | 73.83 ± 0.926 | 68.99 ± 0.35 | 1.179 |
| 8 | 185.20 ± 2.03 | 104.1 ± 1.49 | 56.21 ± 0.22 | 0.941 |

**near-ideal rectifier (bound)** — linear fit E_rect = 14.19·N +38.41 nJ, R² = 0.935; linearity index E(8)/(8·E(1)) = 0.47

| N | E_step,raw (nJ) | E_step,rect (nJ) | η (%) | E_rect(array)/ΣE_rect(single) |
|---|---|---|---|---|
| 1 | 45.83 ± 0.48 | 38.80 ± 0.430 | 84.66 ± 0.09 | 1.000 |
| 2 | 77.66 ± 0.82 | 77.00 ± 0.826 | 99.16 ± 0.19 | 1.187 |
| 4 | 107.02 ± 1.13 | 104.0 ± 1.13 | 97.18 ± 0.25 | 1.198 |
| 8 | 185.20 ± 2.03 | 146.7 ± 1.87 | 79.21 ± 0.15 | 0.961 |



**Discussion.**
* **Rectifier.** Eq. 3 models the bridge as a 2V_F voltage drop. At the current a disc can supply (i = d₃₃·dF/dt ≈ 1–5 µA), the 1N5819 hardly drops any voltage at all. Its small-signal resistance at zero bias, nV_T/I_s ≈ 1.1 kΩ, together with its tens-of-µA saturation leakage, **shorts the disc**: the bus peaks at 1.6 mV. This "Schottky vs. silicon" trade-off reverses the paper's reasoning. For a µA source, leakage matters more than forward drop, and a 1N4148-class bridge (nA leakage) delivers 87 % even with a 0.5–0.6 V drop. Real 1N5819 leakage varies by maker and temperature (it roughly doubles every 10 °C), so this should be checked on the bench. The leakage/100 variant shows a 1N5819 would need roughly < 0.3 µA leakage to be usable.
* **Array scaling.** E_rect(N) is sub-linear (linearity index 0.44–0.47), as Sec. 4.2 of the paper predicted. Discs H3, H4, F3 and M1 sit in lower-pressure areas (bar chart), so this is a **placement problem**. At N = 8 the array delivers 4–6 % less than its channels would separately. While the heel discs keep the bus charged, the forefoot bridges stay reverse-biased. So the channels are not fully independent, as Eq. 4 assumes. The array/sum ratio > 1 at N = 2–4 is a load effect: R_opt/N happens to match the in-shoe waveform better.

### 3.3 Experiment 3: charging profile (20 min, S₁ open)

![exp3](../figures/fig3_exp3_charging.png)

#### Table 8 – Experiment 3: supercapacitor charging rate vs. cadence (20 min, n = 3)

**As designed (1N5819, 1 F)**

| Activity (cadence) | strikes | dV_C/dt (mV/s) | R² | V_C at 20 min (V) | E_C at 20 min (J) | flattening at | Zener peak current |
|---|---|---|---|---|---|---|---|
| walking 1.5 Hz | 900 | 5.792e-8 ± 4.80e-9 | 0.5468 | 1.242e-7 ± 8.99e-9 | 7.714e-15 | — | 2.4e-14 A |
| jogging 2.5 Hz | 1500 | 8.975e-7 ± 5.86e-8 | 0.5503 | 1.932e-6 ± 1.37e-7 | 1.866e-12 | — | 1.7e-12 A |
| running 3.0 Hz | 1800 | 1.455e-6 ± 8.90e-8 | 0.5248 | 3.198e-6 ± 8.88e-8 | 5.113e-12 | — | 1.6e-12 A |

run/walk slope ratio = 25.13 (cadence ratio 2.0)

**1N4148 bridges, 1 F**

| Activity (cadence) | strikes | dV_C/dt (mV/s) | R² | V_C at 20 min (V) | E_C at 20 min (J) | flattening at | Zener peak current |
|---|---|---|---|---|---|---|---|
| walking 1.5 Hz | 900 | 3.811e-4 ± 1.89e-5 | 1.0000 | 4.574e-4 ± 2.23e-5 | 1.046e-7 | — | 1.8e-11 A |
| jogging 2.5 Hz | 1500 | 1.136e-3 ± 4.63e-5 | 1.0000 | 1.363e-3 ± 5.61e-5 | 9.288e-7 | — | 5.3e-11 A |
| running 3.0 Hz | 1800 | 1.469e-3 ± 2.76e-5 | 1.0000 | 1.764e-3 ± 2.82e-5 | 1.555e-6 | — | 6.8e-11 A |

run/walk slope ratio = 3.85 (cadence ratio 2.0)

**Near-ideal rectifier, 1 F**

| Activity (cadence) | strikes | dV_C/dt (mV/s) | R² | V_C at 20 min (V) | E_C at 20 min (J) | flattening at | Zener peak current |
|---|---|---|---|---|---|---|---|
| walking 1.5 Hz | 900 | 4.241e-4 ± 1.96e-5 | 1.0000 | 5.089e-4 ± 2.30e-5 | 1.295e-7 | — | 2.0e-11 A |
| jogging 2.5 Hz | 1500 | 1.242e-3 ± 4.70e-5 | 1.0000 | 1.490e-3 ± 5.70e-5 | 1.110e-6 | — | 5.9e-11 A |
| running 3.0 Hz | 1800 | 1.604e-3 ± 2.79e-5 | 1.0000 | 1.926e-3 ± 2.86e-5 | 1.854e-6 | — | 7.4e-11 A |

run/walk slope ratio = 3.78 (cadence ratio 2.0)

**Near-ideal rectifier, 100 µF (accelerated)**

| Activity (cadence) | strikes | dV_C/dt (mV/s) | R² | V_C at 20 min (V) | E_C at 20 min (J) | flattening at | Zener peak current |
|---|---|---|---|---|---|---|---|
| walking 1.5 Hz | 900 | 0.5358 ± 0.0218 | 0.7418 | 0.8940 ± 0.0367 | 3.996e-5 | 15.2 min | 1.0e-09 A |
| jogging 2.5 Hz | 1500 | 0.9306 ± 0.0310 | 0.6469 | 1.809 ± 0.0679 | 1.637e-4 | 11.5 min | 1.0e-09 A |
| running 3.0 Hz | 1800 | 0.9617 ± 0.0155 | 0.6051 | 2.010 ± 0.0293 | 2.020e-4 | 10.8 min | 1.0e-09 A |

run/walk slope ratio = 1.80 (cadence ratio 2.0)

**1N4148, 100 µF (accelerated)**

| Activity (cadence) | strikes | dV_C/dt (mV/s) | R² | V_C at 20 min (V) | E_C at 20 min (J) | flattening at | Zener peak current |
|---|---|---|---|---|---|---|---|
| walking 1.5 Hz | 900 | 0.4315 ± 0.0202 | 0.7104 | 0.7456 ± 0.0346 | 2.779e-5 | 14.2 min | 1.0e-09 A |
| jogging 2.5 Hz | 1500 | 0.8153 ± 0.0303 | 0.6361 | 1.609 ± 0.0661 | 1.295e-4 | 11.5 min | 1.0e-09 A |
| running 3.0 Hz | 1800 | 0.8559 ± 0.0150 | 0.5984 | 1.807 ± 0.0285 | 1.633e-4 | 10.7 min | 1.0e-09 A |

run/walk slope ratio = 1.98 (cadence ratio 2.0)



**Discussion.**
* **As designed** (left panel), V_C settles within minutes at **0.1–3 µV**. Here the bridge leakage (16 diode strings, about 140 Ω at µV bias) balances the harvest current.
* **With low-leakage rectifiers** (centre), Eq. 4's linear law holds exactly (R² = 1.0000), but the slope is only 0.4–1.6 µV/s. One strike delivers about 1 µC (ΔV ≈ 1 µV on 1 F). The run/walk slope ratio of **3.8** is well above the cadence ratio of 2: the larger running forces add charge per strike (Eq. 1, Q ∝ F) on top of the strike rate. This is the "force-dominated" outcome described in the paper's Sec. 4.3.
* **Accelerated 100 µF test** (right). This shows the predicted flattening: V_C levels off after 11–15 min at 0.9–2.0 V. As V_C rises, only the part of each strike's V_oc above V_C + 2V_F transfers charge, while the DMM (10 MΩ) and leakage keep drawing current. **The 5.1 V clamp is never reached**, because the in-shoe disc V_oc peaks are only about 2–6 V. Even the small store barely reaches the converter's 1.8 V V_min.
* **Energy view.** Charging a 1 F store at millivolts converts almost none of the disc energy (E = Q·V_C). Over the same 20 min the 100 µF store banks **202 µJ**, against **1.9 µJ** for 1 F. Storage should be sized for voltage, not capacity.

### 3.4 Experiment 4: massage delivery

![exp4](../figures/fig4_exp4_massage.png)

![eq6](../figures/fig5_exp4_eq6.png)

#### Table 9a – Experiment 4A: switch-level buck-boost characterisation (V_out = 3 V, 70 mA)

| V_in (V) | mode | duty | P_out (mW) | conduction loss (mW) | gate (mW) | switching (mW) | η_c (%) |
|---|---|---|---|---|---|---|---|
| 1.8 | boost | 0.429 | 209.9 | 10.85 | 3.00 | 2.94 | 92.6 |
| 2.0 | boost | 0.359 | 209.9 | 8.73 | 3.00 | 2.62 | 93.6 |
| 2.2 | boost | 0.290 | 210.0 | 7.11 | 3.00 | 2.37 | 94.3 |
| 2.5 | boost | 0.186 | 209.5 | 5.34 | 3.00 | 2.06 | 95.2 |
| 2.8 | bb | 0.536 | 209.9 | 17.38 | 6.00 | 7.01 | 87.3 |
| 3.0 | bb | 0.517 | 209.9 | 16.28 | 6.00 | 6.98 | 87.7 |
| 3.2 | bb | 0.500 | 209.9 | 15.38 | 6.40 | 6.97 | 87.9 |
| 3.5 | buck | 0.871 | 210.0 | 3.54 | 3.50 | 1.96 | 95.8 |
| 4.0 | buck | 0.762 | 210.0 | 3.76 | 4.00 | 2.24 | 95.4 |
| 4.5 | buck | 0.678 | 209.9 | 3.94 | 4.50 | 2.52 | 94.9 |
| 5.0 | buck | 0.610 | 210.0 | 4.29 | 5.00 | 2.80 | 94.5 |
| 5.1 | buck | 0.598 | 210.0 | 4.34 | 5.10 | 2.86 | 94.4 |


#### Table 9 – Experiment 4B: measured (simulated) and predicted massage runtime (n = 3)

Eq. 6 evaluated with V_min = 1.8 V and η_c = 92.8 % (mean converter efficiency from 4A). Single-parameter fit over the three ESR = 1 Ω points: **η_c = 87.6 %**.

| Supercap | V_C at start (V) | E_C (J) | P_M (mW) | t_m measured (s) | t_m Eq. 6 (s) | Deviation (%) | η_c fitted (%) | V_C left at stop (V) | UVLO restarts |
|---|---|---|---|---|---|---|---|---|---|
| ESR 1 Ω | 3.0 | 4.50 | 208.5 ± 1.4 | 11.27 ± 0.09 | 12.82 | -12.0 | 81.6 | 1.93 | 0 |
| ESR 1 Ω | 4.0 | 8.00 | 205.4 ± 4.8 | 26.96 ± 0.67 | 28.83 | -6.5 | 86.8 | 1.92 | 0 |
| ESR 1 Ω | 5.0 | 12.50 | 210.9 ± 6.4 | 45.53 ± 1.42 | 47.89 | -4.9 | 88.2 | 1.93 | 0 |
| ESR 30 Ω | 3.0 | 4.50 | 118.4 ± 6.7 | 0.00 ± 0.00 | 22.62 | -100.0 | 0.0 | 3.00 | 4806 |
| ESR 30 Ω | 4.0 | 8.00 | 145.4 ± 1.7 | 20.53 ± 1.09 | 40.73 | -49.6 | 46.8 | 2.95 | 6067 |
| ESR 30 Ω | 5.0 | 12.50 | 153.7 ± 0.8 | 29.74 ± 1.27 | 65.67 | -54.7 | 42.0 | 3.04 | 8160 |



**Discussion.**
* With a low-ESR store, the stopwatch time follows Eq. 6 closely. The single fitted **η_c = 87.6 %** is about 5 points below the converter's mean 92.8 %. That gap includes the ESR loss, I²·ESR (1–2 % at 70–140 mA), and the poorer buck-boost-mode efficiency (87–88 % at V_in ≈ 3 V, where the discharge spends much of its time). The deviation is largest at 3 V (−12 %), because most of that run is spent near V_in ≈ V_out.
* The motor stops at V_C ≈ 1.93 V. At that point the terminal voltage V_C − I·ESR hits the 1.8 V UVLO, so Eq. 6 should use the **internal** cap voltage, V_min,eff = V_min + I_in·ESR.
* With a **30 Ω coin-type** 1 F/5.5 V cell, the store can deliver at most V_C²/(4·ESR) = 75–208 mW, which is barely the motor's demand. The converter hiccups between UVLO and restart thousands of times, and the motor runs at reduced, pulsating speed. From 3 V it never reaches perceptible vibration. Eq. 6 overestimates t_m by about 50–100 %, and about 3 V of charge stays unused. **Specify a low-ESR (≤ 1–2 Ω) EDLC.**
* At ≈ 210 mW, a 45 s massage needs the full 12.5 J of the 5 V store.

### 3.5 Energy autonomy (Sec. 4.5 figure of merit)

![autonomy](../figures/fig6_autonomy.png)

#### Synthesis – energy autonomy (running, 3.0 Hz)

* average harvest current into 1 F (near-ideal bridge): 1.60 µA; leakage + DMM = 1.80 MΩ → asymptote 2.89 V
* time to reach V_min = 1.8 V: 488 h of continuous running (Eq. 4 linear, no leakage: 312 h); to 5 V (linear): 866 h
* energy banked by a 20-min run: as designed 5.11e-12 J; near-ideal 1 F 1.85 µJ; near-ideal 100 µF 202 µJ

The figure of merit t_m / t_run is therefore **0** for any realistic run with the BOM as specified. Even with ideal rectifiers and no leakage, it would take **~310 h** of running to reach 1.8 V and **~870 h** to reach 5 V. With realistic leakage (2.2 MΩ supercap ∥ 10 MΩ DMM), the store levels off at ~2.9 V, and 1.8 V takes about **490 h**. For scale, Shenck & Paradiso's 8.4 mW PZT dimorph would bank about 10 J in 20 min. These 27 mm buzzer discs give about 0.5–1 µW average per shoe, roughly 10⁴ × less.

---

## 4. Answers to the paper's diagnostic questions (Table 5)

1. **Energy per strike and optimal load:** 1.33 µJ/strike at 3.45 MΩ on the plate rig with full body weight. In the shoe, each disc gives 4–45 nJ per walking strike at that load, which is placement-dependent.
2. **Bridge loss and array additivity:** the 1N5819 bridge loses ~100 % (leakage-dominated), 1N4148 loses 13 % and a near-ideal rectifier 4 %. The array is sub-linear (index ≈ 0.45), mainly because of placement plus slight bus back-bias coupling.
3. **Charging speed:** linear (R² = 1), at 0.42 / 1.24 / 1.60 µV/s (walk / jog / run) with good rectifiers and ~0 as designed. The clamp is never reached.
4. **Massage duration:** 11.3 / 27.0 / 45.5 s from 3 / 4 / 5 V. Eq. 6 holds with η_c = 87.6 % and a low-ESR store.

## 5. Implications for the "Future Work" section

1. **Replace the 1N5819 bridges** with low-leakage rectifiers: 1N4148/BAS116-class diodes, or better, an integrated piezo harvester IC (e.g. LTC3588-1 / LTC3331 / BQ25570-class, with nA-leakage active bridges, UVLO and a buck output). This is the largest single gain (from 0 % to 85–95 %).
2. **Shrink and decouple the storage.** First charge a small capacitor (10–100 µF) to its optimum voltage (≈ V_oc/2), then transfer to the store with a switched converter or MPPT stage, or use SSHI/SECE interfaces (4–10× gain reported in the literature). A 1 F store cannot be charged by µC strikes.
3. **Re-plan Experiment 1** with a 1 kΩ–100 MΩ sweep and a ≥ 1 GΩ measurement front end.
4. **Increase transducer energy.** Use larger or stacked discs, PZT or PVDF benders in the heel cup, or flexure amplification that turns the 1–5 mm sole deflection into bending. The energy budget needs a 10³–10⁴× improvement before a 30–45 s, 200 mW massage becomes possible per run. Alternatively, cut actuator power, e.g. a low-duty-cycle pulsed ERM or an LRA at ~20–50 mW.
5. **Use a low-ESR EDLC** (≤ 1–2 Ω), and set V_min in Eq. 6 to the internal cap voltage at UVLO.

## 6. Reproducibility and files

| Path | Content |
|---|---|
| `sim/kshoe.py` | Shared parameters, SPICE models, gait model, ngspice runner |
| `sim/exp1.py` … `sim/exp4.py` | One script per experiment (writes `results/expN.json`) |
| `sim/figures.py` | All figures and the filled tables (`results/tables.md`) |
| `sim/run_all.py` | Runs everything and rebuilds this report |
| `figures/*.png` | Figures shown above |
| `results/*.json` | Raw summary data, including all trial-level statistics |
| `work/` | Temporary force stimulus files (safe to delete; regenerated on each run) |
| `tools/ngspice/` | ngspice 45.2 binaries |

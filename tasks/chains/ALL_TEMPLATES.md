# All 10 Chain Templates — Fill-In Reference

Edit the `[FILL: ...]` sections below. Once done, copy each step's content into its corresponding JSON file in `tasks/chains/tasks/`.

---

## Chain 1: Paper to Experiment (4 steps)
**paper_to_experiment** — Paper finding → data interpretation → statistical test → hypothesis generation

### Step 1 — Extract key finding from literature
**File:** `chain01_step1.json`
- **Question:** [FILL: Describe a specific paper finding about a protein-protein interaction. Example: 'A 2023 study in Cell found that protein X binds receptor Y with Kd = 2.3 nM, but only when zinc ions are present. What is the role of zinc in this interaction?']
- **Ideal answer:** [FILL: The correct mechanistic answer. Example: 'Zinc acts as an allosteric cofactor that stabilizes the binding interface between X and Y by coordinating residues His142 and Cys205.']
- **Verification:** llm-judge
- **Cascade risk:** If the model misunderstands the mechanism, all downstream steps use wrong reasoning

### Step 2 — Quantitative data interpretation
**File:** `chain01_step2.json`
- **Question:** [FILL: Present data and ask the model to interpret it in light of step 1. Example: 'Given that zinc facilitates the X-Y interaction, the following dose-response data shows binding at zinc concentrations of 0, 0.1, 1, 10, 100 µM with percent binding of 2%, 8%, 31%, 78%, 95%. At what zinc concentration does binding reach half-maximal? Is this consistent with physiological zinc levels (10-100 µM)?']
- **Ideal answer:** [FILL: Quantitative answer plus interpretation. Example: 'Half-maximal binding occurs at approximately 3 µM zinc. This is below the physiological range, suggesting the interaction is saturated under normal conditions.']
- **Verification:** llm-judge
- **Cascade risk:** Wrong concentration estimate → wrong experimental design downstream

### Step 3 — Statistical test selection
**File:** `chain01_step3.json`
- **Question:** [FILL: Ask which statistical test for the follow-up experiment. Example: 'You want to test whether the zinc-dependent binding (half-max ~3 µM) differs significantly between wild-type protein X and a H142A mutant. You have 8 biological replicates per condition, continuous measurements. Which statistical test is appropriate and why?']
- **Ideal answer:** [FILL: Correct test with justification. Example: 'Unpaired two-sample t-test (or Welch's t-test). Two groups, continuous outcome, independent replicates. Mann-Whitney U as non-parametric alternative with n=8.']
- **Verification:** llm-judge
- **Cascade risk:** Wrong test → flawed experimental conclusions

### Step 4 — Hypothesis generation
**File:** `chain01_step4.json`
- **Question:** [FILL: Synthesize everything into hypotheses. Example: 'Based on your analysis — zinc mediates X-Y binding at physiological concentrations, half-max ~3 µM, testable via t-test with WT vs mutant — propose 3 testable hypotheses for how zinc dysregulation in Wilson disease might disrupt this pathway.']
- **Ideal answer:** [FILL: 3 specific, falsifiable hypotheses. Example: '1) Copper accumulation displaces zinc from His142/Cys205, reducing affinity >10-fold. 2) Hepatic zinc depletion drops below 3 µM threshold, abolishing signaling. 3) H142 undergoes copper-mediated oxidation, permanently disrupting the site.']
- **Verification:** llm-judge
- **Cascade risk:** Terminal step — quality depends on everything above

---

## Chain 2: Structure to Drug (4 steps)
**structure_to_drug** — Protein structure → binding mechanism → SAR prediction → experimental validation

### Step 1 — Identify structural features
**File:** `chain02_step1.json`
- **Question:** [FILL: Ask about a protein's structure from PDB data. Example: 'The crystal structure of human carbonic anhydrase II (PDB: 1CA2) contains a zinc ion in the active site. How many protein chains are in this structure, and what amino acids coordinate the catalytic zinc?']
- **Ideal answer:** [FILL: Structural facts. Example: 'One chain. Zinc coordinated by His94, His96, His119, plus a water/hydroxide molecule.']
- **Verification:** llm-judge
- **Cascade risk:** Wrong active site residues → wrong binding analysis → wrong drug design

### Step 2 — Explain binding mechanism
**File:** `chain02_step2.json`
- **Question:** [FILL: Ask about inhibitor binding given the structure. Example: 'Given the active site architecture you described, how do sulfonamide inhibitors like acetazolamide bind to carbonic anhydrase II? What is the key interaction with the zinc ion?']
- **Ideal answer:** [FILL: Binding mechanism. Example: 'Sulfonamide NH displaces zinc-bound water, forming direct N-Zn coordination. Sulfonamide oxygens H-bond with Thr199, aromatic ring makes hydrophobic contacts.']
- **Verification:** llm-judge
- **Cascade risk:** Wrong binding model → wrong SAR predictions

### Step 3 — SAR prediction
**File:** `chain02_step3.json`
- **Question:** [FILL: Multiple-choice SAR prediction. Example: 'Based on the binding mode, which modification would MOST improve affinity: A) Fluorine on ring, B) Replace sulfonamide with carboxylate, C) Extend ring into hydrophobic pocket, D) Add hydroxyl to ring']
- **Ideal answer:** [FILL: Letter answer. Example: 'C']
- **Verification:** programmatic (multiple_choice)
- **Cascade risk:** Wrong SAR → wrong lead optimization

### Step 4 — Experimental validation plan
**File:** `chain02_step4.json`
- **Question:** [FILL: Design experiment to validate prediction. Example: 'Design a 3-step experimental plan to validate your SAR prediction. Include: synthesis approach, binding assay, and selectivity assessment.']
- **Ideal answer:** [FILL: Concrete plan. Example: 'Suzuki coupling for synthesis, CO2 hydration stopped-flow assay for Ki, profile against CA isoform panel for selectivity.']
- **Verification:** llm-judge
- **Cascade risk:** Terminal step

---

## Chain 3: Stats Pipeline (3 steps)
**stats_pipeline** — Test selection → multiple testing correction → pathway interpretation

### Step 1 — Choose statistical framework
**File:** `chain03_step1.json`
- **Question:** [FILL: Present gene expression data, ask for test. Example: 'RNA-seq from 12 tumor + 12 matched normal. 20,000 genes, normalized counts. What statistical test for DE genes and why?']
- **Ideal answer:** [FILL: Example: 'DESeq2/edgeR with negative binomial GLM. Paired design for matched samples. T-test inappropriate for count data.']
- **Verification:** llm-judge
- **Cascade risk:** Wrong test → wrong p-values → wrong gene list → wrong pathways

### Step 2 — Multiple testing correction
**File:** `chain03_step2.json`
- **Question:** [FILL: Give p-value results, ask about MTC. Example: '1,247 genes at p<0.05 from 20,000. How to correct? How many survive BH at q<0.05? How many false positives expected?']
- **Ideal answer:** [FILL: Example: 'BH procedure. ~800 survive, expect ~40 FPs (5% of significant). Bonferroni too conservative for 20K tests.']
- **Verification:** llm-judge
- **Cascade risk:** Wrong correction → inflated/deflated gene lists

### Step 3 — Pathway interpretation
**File:** `chain03_step3.json`
- **Question:** [FILL: Present GO enrichment, ask for interpretation. Example: 'Enriched: inflammatory response (p=1e-12, 45 genes), cell cycle (p=3e-8, 38 genes), apoptosis (p=0.02, 12 genes). Colleague says all equally important. What is wrong?']
- **Ideal answer:** [FILL: Example: 'Confusing significance with importance. P-values differ by 10 orders of magnitude. Gene counts differ. Apoptosis marginal after GO-level MTC. Prioritize inflammatory response.']
- **Verification:** llm-judge
- **Cascade risk:** Terminal step

---

## Chain 4: Critical Appraisal (3 steps)
**critical_appraisal** — Evaluate weak evidence → conflicting data → definitive experiment

### Step 1 — Evaluate strength of evidence
**File:** `chain04_step1.json`
- **Question:** [FILL: Present a claim with weak evidence. Example: 'Preprint: drug X reduces tumor growth 40% in mice (n=5, p=0.048), single cell line (4T1), subcutaneous model. Sufficient to conclude efficacy?']
- **Ideal answer:** [FILL: Example: 'No. Marginal p with n=5, single cell line, subQ model not representative, no dose-response. Preliminary evidence only.']
- **Verification:** llm-judge
- **Cascade risk:** Accepting weak evidence → wrong conclusions downstream

### Step 2 — Integrate contradictory evidence
**File:** `chain04_step2.json`
- **Question:** [FILL: Present conflicting follow-up data. Example: 'New study, 3 cell lines, n=10: 4T1 15% reduction (p=0.21), B16 5% increase (p=0.73), LLC 32% reduction (p=0.01). How does this change your assessment?']
- **Ideal answer:** [FILL: Example: '4T1 did not replicate (40%→15%, now NS). Inconsistent across lines. Only LLC significant — may be cell-line-specific, not broad anti-cancer.']
- **Verification:** llm-judge
- **Cascade risk:** Failure to update → wrong experimental recommendation

### Step 3 — Design definitive experiment
**File:** `chain04_step3.json`
- **Question:** [FILL: Ask for experiment to resolve discrepancy. Example: 'Design a definitive experiment: cell lines, sample sizes, controls, endpoint, analysis plan.']
- **Ideal answer:** [FILL: Example: '6+ cell lines, 3 cancer types, n=15/group, vehicle + positive control, tumor volume day 21, dose-response, two-way ANOVA, pre-registered.']
- **Verification:** llm-judge
- **Cascade risk:** Terminal step

---

## Chain 5: Genetics to Therapy (3 steps)
**genetics_to_therapy** — Genetic finding → structural analysis → therapeutic strategy

### Step 1 — Gene function and disease mechanism
**File:** `chain05_step1.json`
- **Question:** [FILL: Present a genetic finding. Example: 'Family with AR early-onset Parkinson. Homozygous PINK1 p.Gly309Asp. What is PINK1 function, how might this cause disease?']
- **Ideal answer:** [FILL: Example: 'PINK1 is mitochondrial kinase for mitophagy via Parkin. G309D in kinase domain disrupts ATP binding, impairs mitophagy, ROS accumulation, dopaminergic neuron loss.']
- **Verification:** llm-judge
- **Cascade risk:** Wrong pathway → wrong therapeutic targets

### Step 2 — Structural consequence of mutation
**File:** `chain05_step2.json`
- **Question:** [FILL: Ask about structural impact. Example: 'G309 is in the activation loop. What does Gly→Asp do structurally? How to confirm experimentally?']
- **Ideal answer:** [FILL: Example: 'Asp introduces steric clash, disrupts DFG motif, prevents activation loop phosphorylation. Confirm: CD, thermal shift, in vitro kinase assay WT vs G309D.']
- **Verification:** llm-judge
- **Cascade risk:** Wrong structural model → wrong drug design

### Step 3 — Therapeutic strategy
**File:** `chain05_step3.json`
- **Question:** [FILL: Ask for two therapeutic approaches. Example: 'Propose one strategy targeting PINK1-Parkin pathway directly, one targeting downstream consequences. Name specific compounds.']
- **Ideal answer:** [FILL: Example: '1) Kinetin triphosphate as neo-substrate for partial PINK1 activity. 2) MitoQ/SS-31 as mitochondria-targeted antioxidant for ROS from damaged mitochondria.']
- **Verification:** llm-judge
- **Cascade risk:** Terminal step

---

## Chain 6: Protocol Troubleshoot (3 steps)
**protocol_troubleshoot** — Diagnose error → interpret corrected result → quantitative follow-up

### Step 1 — Diagnose protocol error
**File:** `chain06_step1.json`
- **Question:** [FILL: Describe a failed experiment. Example: 'Western blot for phospho-ERK1/2: RIPA lysis, 12% gel, PVDF, block 5% milk/TBST, anti-pERK 1:1000 O/N 4C. No signal but positive control works. What is wrong?']
- **Ideal answer:** [FILL: Example: 'Milk contains casein (a phosphoprotein) that competes with phospho-epitopes. Use 5% BSA for phospho-antibodies.']
- **Verification:** llm-judge
- **Cascade risk:** Wrong diagnosis → wrong fix → wasted experiment

### Step 2 — Interpret corrected result
**File:** `chain06_step2.json`
- **Question:** [FILL: Present results after fix. Example: 'After switching to BSA: strong band at ~42 kDa, weaker at ~44 kDa. Are these expected for pERK1/2? What do relative intensities mean?']
- **Ideal answer:** [FILL: Example: 'ERK2=42kDa, ERK1=44kDa. Stronger 42kDa = ERK2 dominant activation, typical for HeLa. Both bands confirm antibody specificity.']
- **Verification:** llm-judge
- **Cascade risk:** Wrong band assignment → wrong interpretation

### Step 3 — Quantitative follow-up
**File:** `chain06_step3.json`
- **Question:** [FILL: Ask for quantification plan. Example: 'Quantify pERK1/2 across EGF time course (0, 5, 15, 30, 60 min), 3 replicates. Describe quantification and statistical test for peak activation.']
- **Ideal answer:** [FILL: Example: 'Densitometry (ImageJ), normalize pERK to total ERK, one-way repeated measures ANOVA + Dunnett post-hoc vs t=0. Report mean ± SEM.']
- **Verification:** llm-judge
- **Cascade risk:** Terminal step

---

## Chain 7: Paradox Resolution (3 steps)
**paradox_resolution** — Explain paradox → discriminating experiment → synthesize conclusion

### Step 1 — Explain paradoxical result
**File:** `chain07_step1.json`
- **Question:** [FILL: Present a surprising finding. Example: 'CRISPR KO of ZEB1 (EMT promoter) in A549 cells INCREASED migration 2.1x (p<0.001). How do you explain this paradox?']
- **Ideal answer:** [FILL: Example: 'Possible: 1) ZEB2 paralog compensation/hyperactivation, 2) Collective migration phenotype (faster sheet migration via restored E-cadherin) vs. mesenchymal single-cell, 3) Wound healing assay measures 2D not 3D invasion.']
- **Verification:** llm-judge
- **Cascade risk:** Can't generate hypotheses → can't design discriminating experiments

### Step 2 — Design discriminating experiment
**File:** `chain07_step2.json`
- **Question:** [FILL: Ask for ONE experiment distinguishing top two hypotheses. Example: 'Design an experiment distinguishing ZEB2-compensation from collective-migration hypothesis. Specify assay, controls, expected results for each scenario.']
- **Ideal answer:** [FILL: Example: 'Transwell invasion through Matrigel. If ZEB2 compensation: increased invasion too. If collective migration: decreased invasion (lost mesenchymal phenotype). Controls: WT, ZEB1-KO, ZEB1/ZEB2 double KO. Western for ZEB2, E-cadherin, vimentin.']
- **Verification:** llm-judge
- **Cascade risk:** Bad experiment → can't resolve paradox

### Step 3 — Synthesize conclusion from results
**File:** `chain07_step3.json`
- **Question:** [FILL: Present results, ask for conclusion. Example: 'ZEB1-KO: 60% decreased Matrigel invasion (p<0.01). Double KO: 80% decreased. E-cadherin restored, vimentin reduced. Which hypothesis wins? What does this mean for wound healing assays in EMT?']
- **Ideal answer:** [FILL: Example: 'Supports collective migration. More migratory in 2D but less invasive in 3D = epithelial sheet migration, not mesenchymal invasion. ZEB2 compensation weakened (double KO further reduces). Wound healing assays misleading for EMT — use 3D invasion.']
- **Verification:** llm-judge
- **Cascade risk:** Terminal step

---

## Chain 8: Sequence to Function (3 steps)
**sequence_to_function** — Identify protein → predict adaptations → design validation

### Step 1 — Identify protein from sequence
**File:** `chain08_step1.json`
- **Question:** [FILL: Present a protein sequence, ask for identification. Example: 'Sequence from deep-sea organism: [paste ~100 AA]. What protein family? What conserved domains expected?']
- **Ideal answer:** [FILL: Example: 'Based on [motif], this is a [family] protein. Expected domains: [domain 1] residues X-Y, [domain 2] residues A-B.']
- **Verification:** llm-judge
- **Cascade risk:** Wrong family → wrong functional predictions → wrong experiment

### Step 2 — Predict functional consequences
**File:** `chain08_step2.json`
- **Question:** [FILL: Ask about specific mutations vs. human homolog. Example: 'Compared to human homolog, this has [2-3 mutations] in the active site. Organism lives at 4C, 200 atm. How might these be cold/pressure adaptations?']
- **Ideal answer:** [FILL: Example: 'Mutation 1 increases flexibility at low temp. Mutation 2 strengthens hydrophobic packing against pressure. Consistent with psychrophilic enzyme adaptations.']
- **Verification:** llm-judge
- **Cascade risk:** Wrong adaptation reasoning → wrong biochemical predictions

### Step 3 — Design validation experiment
**File:** `chain08_step3.json`
- **Question:** [FILL: Ask for experimental plan. Example: 'Design experiment to test cold/pressure tolerance. Compare deep-sea vs. human variant. Expression, purification, activity measurements.']
- **Ideal answer:** [FILL: Example: 'E. coli expression, Ni-NTA + SEC purification. Kinetics at 4/15/25/37C at 1 atm. High-pressure cell at 1/100/200/500 atm at 4C. DSF for Tm. Predict: deep-sea >50% activity at 4C/200atm where human <10%.']
- **Verification:** llm-judge
- **Cascade risk:** Terminal step

---

## Chain 9: Data to Mechanism (3 steps)
**data_to_mechanism** — Interpret ambiguous data → update with evidence → correct prior analysis

### Step 1 — Interpret ambiguous data
**File:** `chain09_step1.json`
- **Question:** [FILL: Present confusing experimental result. Example: 'Western for GAPDH (37 kDa loading control): Lane 1 (WT) strong 37kDa. Lane 2 (Treated) strong 37kDa + unexpected 25kDa band. Lane 3 (KO) faint 37kDa. What is the 25kDa band?']
- **Ideal answer:** [FILL: Example: 'Likely GAPDH degradation product. Caspase-3 cleavage during apoptosis (GAPDH is a known substrate). Or non-specific cross-reactivity (less likely, lane-specific). Or treatment-activated protease.']
- **Verification:** llm-judge
- **Cascade risk:** Wrong interpretation → wrong mechanistic model

### Step 2 — Update with new evidence
**File:** `chain09_step2.json`
- **Question:** [FILL: Present confirmatory data. Example: '45% apoptotic (Annexin V), 8x caspase-3 activity, Z-VAD-FMK eliminates 25kDa band. Which hypothesis confirmed? What does this mean for GAPDH as loading control?']
- **Ideal answer:** [FILL: Example: 'Confirms caspase-3 cleavage. Z-VAD rescue definitive. GAPDH is NOT valid loading control here — treatment degrades it. Switch to Ponceau S or Na/K-ATPase.']
- **Verification:** llm-judge
- **Cascade risk:** Wrong loading control → all quantification wrong

### Step 3 — Correct prior analysis
**File:** `chain09_step3.json`
- **Question:** [FILL: Ask model to re-analyze with corrected understanding. Example: 'Original paper: treatment reduces protein X by 70% (normalized to GAPDH). But GAPDH is degraded. Is the 70% real? How to re-analyze?']
- **Ideal answer:** [FILL: Example: '70% overestimated. Degraded GAPDH denominator artificially low → inflated ratio. If GAPDH ~50% degraded, true reduction ~40%. Re-probe with Ponceau S normalization.']
- **Verification:** llm-judge
- **Cascade risk:** Terminal step — tests backward error correction

---

## Chain 10: Evidence Synthesis (3 steps)
**evidence_synthesis** — Compare conflicting studies → meta-analytic reasoning → clinical recommendation

### Step 1 — Compare conflicting papers
**File:** `chain10_step1.json`
- **Question:** [FILL: Present two conflicting studies. Example: 'Paper A (Nature 2022, n=500): biomarker Z predicts 5yr survival, AUC=0.89. Paper B (Lancet Oncol 2024, n=1200): AUC=0.61 for same biomarker. Explain the discrepancy.']
- **Ideal answer:** [FILL: Example: 'Overfitting in smaller Paper A, population differences, assay differences, post-hoc optimization vs. pre-registered, publication bias. Larger/newer Paper B more credible.']
- **Verification:** llm-judge
- **Cascade risk:** Uncritical acceptance → wrong clinical recommendation

### Step 2 — Meta-analytic reasoning
**File:** `chain10_step2.json`
- **Question:** [FILL: Present 3 more studies, ask for synthesis. Example: 'Study C (n=300, AUC=0.78), D (n=150, AUC=0.85), E (n=800, AUC=0.65). Best estimate of true AUC? How to combine? What pattern emerges?']
- **Ideal answer:** [FILL: Example: 'Random-effects meta-analysis, pooled AUC ~0.70-0.72. Clear small-study effect: small n → high AUC. High I2 heterogeneity. Funnel plot asymmetry. True AUC likely 0.65-0.72.']
- **Verification:** llm-judge
- **Cascade risk:** Wrong pooled estimate → wrong clinical utility assessment

### Step 3 — Clinical recommendation
**File:** `chain10_step3.json`
- **Question:** [FILL: Ask for clinical adoption decision. Example: 'Based on likely AUC ~0.70, significant heterogeneity, small-study bias — should biomarker Z enter clinical practice? What next steps?']
- **Ideal answer:** [FILL: Example: 'Not yet. AUC 0.70 insufficient alone. Need: prospective multi-center validation (n>1000, pre-registered, locked assay), added value over staging, clinical utility demonstration, cost-effectiveness. Currently suitable as research stratification biomarker only.']
- **Verification:** llm-judge
- **Cascade risk:** Terminal step

# Transfer and revision record: Mathematics of Computation submission

## Revision of 1 October 2026 (desk-screening pass)

Basis: the transfer package of 29 September 2026 (`previous_versions/main_2026-09-29.tex`, SHA-256 recorded in `TRANSFER_VERIFICATION.json`). The target journal, title, author identity and metadata are unchanged.

### What changed in the manuscript

- **Abstract (274 words).** Now opens with the physical observation law, states the computational question in one sentence, and lists the four results in plain terms. All claims are restatements of Proposition 4.6, Theorems 4.2, 4.3, 5.2 and 6.2, Corollary 5.3 and Proposition 6.3.
- **Introduction.** New Section 1.1 states the application (impedance spectroscopy, distribution of relaxation times, calibration parameters), the normalized model, the information width and the numerical excess before any result is claimed. Section 1.2 keeps the previous summary of results, with "worst-class" replaced by "worst-case over the class" and one sentence pointing to the executed examples. Section 1.3 adds the classical Markov--Krein moment theory and the regularized DRT literature as context. Section 1.4 adds a notation paragraph that defines $\mathcal P_n$, the coefficient norm, the norms used, "raw feasible", and the constant convention, which were previously used without a definition.
- **Lead-ins.** Sections 2--6 open with a short paragraph saying what the section does and why; the baseline prediction distance used in Theorem 4.2 is now defined explicitly before the theorem. Subsection 4.3's title uses "worst-case".
- **Section 7.** Tables 1 (equal-effort efficiency, twelve certified enclosures), 2 (certified cost lower bounds) and 4 (guard-grid counts) and Figure 1 (exact response difference with certified maximizers and uniform bounds) are added. Every number is taken from the archived certificates or drawn from the closed form that those certificates verify; the figure script cross-checks its maximizers against the archived enclosures before drawing. The stopping-path table is unchanged and is now Table 3.
- **Conclusions.** A paragraph on limitations (known noise level and interval; class-tuned constants; synthetic data) and the scope of Proposition 4.4 is added.
- **Declarations.** The AI-use disclosure keeps the same scope and is reworded to state that the author verified all results.
- **References.** Five entries added (Karlin--Studden 1966; Krein--Nudel'man 1977; Orazem--Tribollet 2017; Ciucci--Chen 2015; Wan et al. 2015), all cited in Section 1; the Wasilkowski--Wo\'zniakowski entry now follows the format of the other journal entries. No existing entry was removed.

### What did not change

All 4 theorems, 3 lemmas, 5 propositions, 2 corollaries, 2 remarks and 14 proofs, all 60 inherited display blocks, Algorithm 1 and the stopping-path table are byte-identical to the previous version (`source_integrity.json`). No numerical value, seed, raw observation or error radius was changed. No new experiment or optimization run was performed.

### Supplement

The 210 code and data files under `code/` are byte-identical to the previous package. Only `SECTION_MAP.txt` (table numbers and the figure note) and `TRANSFER_NOTE.txt` (revision note) were updated, and `MANIFEST_SHA256.json` regenerated. The frozen evidence was replayed from the unzipped tree (`evidence_replay/`) and from a fresh extraction of the new ZIP (`zip_evidence_replay/`): 80 certificate tasks and 154 checks including 53 required rejections passed in both.

### Build and preflight

`tools/assemble_package.py` rebuilds the three PDFs with the standard `amsart` class, rebuilds them again in an isolated copy and confirms identical page text, checks embedded fonts, black text and page boxes, writes the integrity records, regenerates both manifests and runs `verify_delivery.py`. The compile logs contain one pre-existing harmless microtype font-expansion notice (typewriter font in URLs) and one underfull vertical box on a full bibliography page; neither affects the output.

### Remaining author decisions

Confirm no concurrent submission, the three-submissions-per-12-month rule, the accuracy of the AI-use disclosure, and the supplement upload route of the AMS portal. Journal policy statements were checked through search summaries because ams.org was not directly reachable from this environment; re-confirm on ams.org before filing.

## Transfer of 29 September 2026: SINUM final submission to Mathematics of Computation

### Source discipline

The source is the provided final SINUM package, not a newly inferred Science China version. The source manuscript title was "Frequency Saturation and Certified Recovery of Relaxation Functionals". The new manuscript title is "Frequency Saturation and Certified Computation of Relaxation Functionals". The signed-off author identity supplied in the prior materials is retained.

No theorem hypotheses, quantifiers, proof bodies, numerical values, seeds, raw observations, or error radii were changed. The source diff (`SINUM_to_MCOM.diff`) identifies the changes exactly. Algorithm 1 is a structured restatement of the existing finite-refinement argument. Its theoretical fallback and acceptance checks have the same meanings as Theorem 5.2; it is not described as a newly benchmarked solver or a runtime-optimal method.

### Editorial strategy

The first page starts from a concrete computational output: the full feasible range of a functional and its certified enclosure. The identifying threshold and finite/all-frequency modulus explain the information supplied to this computation. The Gaussian acquisition result is presented as an attainable measurement-effort law, with numerical excess controlled separately.

The title stays within relaxation functionals rather than asserting a sharp theorem for all positive-measure inverse operators. The existing sufficient operator proposition remains in the manuscript. Mathematical scope conditions stay in the model and theorem statements rather than being removed to make broader promotional claims.

### Literature update (29 September)

One addition: G. W. Wasilkowski and H. Wo\'zniakowski, "The power of standard information for multivariate approximation in the randomized setting," Math. Comp. 76 (2007), no. 258, 965--988, DOI 10.1090/S0025-5718-06-01944-2. The present RC problem is not claimed to be a corollary. All 24 pre-existing bibliography entries are preserved.

### Submission materials

The manuscript uses the installed, unmodified AMS `amsart` class. The letter is addressed to the Editors of Mathematics of Computation and suggests Markus Bachmayr on the basis of his current associate-editorship and research in approximation, adaptivity, inverse problems and computational complexity. Assignment remains with the journal.

No submission, cloud write or external message was made in either round. The two earlier editorial returns establish closure only of those submissions; they do not validate or invalidate the mathematical proofs.

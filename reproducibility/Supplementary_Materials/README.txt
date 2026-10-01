SUPPLEMENTARY MATERIALS
Frequency Saturation and Certified Computation of Relaxation Functionals
Haotian Zhong, Shanghai University
Contact: zhonghaotian1219@gmail.com

1. PURPOSE
This is the research-facing reproducibility snapshot accompanying the final
manuscript. It contains exact inputs, certificates, verification software,
optional candidate-generation software, recorded seeds, and table scripts.
It excludes historical manuscripts, cover letters, draft comparisons, and
internal editorial reports. The mathematical proofs are in the main article.

2. QUICK START (NO THIRD-PARTY PYTHON PACKAGES)
Extract the ZIP before running commands. From this directory, run:

    python3 run_all.py

On Windows, use "python" instead of "python3" where appropriate.
An alternative output directory can be selected with:

    python3 run_all.py --output /path/to/replay_results

Python 3.10 or later is required. Frozen verification uses the standard library.
The script checks the file manifest, verifies the required case inventory,
replays each suite in a disposable workspace, and regenerates the numerical
TeX tables. The frozen code and inputs remain unchanged. New reports appear in
replay_output/ by default. The final status is replay_output/SUMMARY.json.
A failed check or timeout returns a nonzero exit code rather than a success.

3. EXPECTED ACCOUNTING
80 independent certificate replay tasks:
  - 22 full-frequency positive pairs and 4 operator maps;
  - 18 stopped-cost certificates and 2 stopping-trajectory tasks;
  - 6 strict-feasibility anchors and 12 equal-effort design certificates;
  - 16 noiseless-identification certificates.
The 2 trajectory tasks contain 6 original stage data vectors and 12 continuous
endpoint certificates; these 12 endpoints are not counted as 12 extra datasets.

154 named algebra/mechanism/rejection checks, including 53 intentional invalid
input or certificate rejections. A separate optional SymPy program reconstructs
the 16 identification examples without importing their production checker.
The counts identify finite computational checks, not independent mathematical
proofs, new experiments, or Monte Carlo estimates of confidence coverage.

4. MAP FROM THE PAPER TO THE FILES
See SECTION_MAP.txt. INDEX.txt lists every shipped file. PROVENANCE.json maps
retained source/data files to the author-supplied v1.4 research archive and
records their SHA-256 hashes. Directory names inside code/ are retained to
preserve the original import paths and permit byte-identical verification code.

5. RESULTS AND REGENERATION
verification/ contains the completed packaging replay reports. tables/ contains
TeX/CSV/JSON regenerated from the exact results. A fresh run produces new reports
in replay_output/ and can differ in execution-time fields while agreeing in
mathematical quantities. Source-truth records in the sequential provenance/
directory are used only for evaluating generated synthetic cases; the isolated
endpoint verifier does not copy or load that directory.

For OPTIONAL regeneration of synthetic inputs and optimization candidates, see
REGENERATION.txt and the recorded environment below code/inherited/fs2/EIS_PEM_Math_FS2_Sequential_Cost/checks/.
Do this in a separate writable copy. Generation may use NumPy/SciPy and is more
expensive than frozen verification. Those packages are not needed for run_all.py.

6. VERIFICATION SCOPE
The full-frequency response bound is checked through an exact formula and a
rational enclosure of its unique maximizer. Discrete frequency samples are
additional checks, not the continuum proof. A continuous dual certificate is
checked in the original observation coordinates. Resource-capped candidate
search can return UNRESOLVED; only verified certificates are accepted.

7. RIGHTS AND DEPENDENCIES
See RIGHTS_AND_DEPENDENCIES.txt. This assembly does not assign a new open-source
license, alter existing ownership, publish a repository, or create a DOI.

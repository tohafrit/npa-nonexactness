# Copy-paste kit for the arXiv submissions

## PAPER 1 (submit first): the proof paper

**File to upload:** paper/npa_nonexactness_proof.tex

**Title:**
No finite level of the NPA hierarchy is exact for the doubly-tilted CHSH functional near the critical tilt

**Authors:** Anton Pakhunov

**Abstract (plain text, fits the 1920-char limit):**
Gigena, Panwar, Scala, Araujo, Farkas and Chaturvedi [npj Quantum Inf. 11, 82 (2025)] determined the quantum maximum of the doubly-tilted CHSH functionals, observed that the Navascues-Pironio-Acin level needed for exactness grows without evident bound toward the critical tilt, and asked whether any finite level suffices. We answer this in the negative. For the symmetric critical family $B_s=(1-s/2)(\langle A_0\rangle+\langle B_0\rangle)+\mathrm{CHSH}$ we prove: for every NPA level $k\ge 2$ there are an explicit rational $g_k>0$ and an $s^*_k>0$ with $c_k(s)\ge 4-s+g_k s^2$ on $(0,s^*_k]$; since the quantum value leaves the local bound only cubically, every finite level strictly overshoots on an interval: no finite level is exact on any neighbourhood of the critical point. Unconditionally $a_2>1/39$, $a_3>1/188$, $a_4>1/641$. The proof is a primal construction: an exactly feasible moment curve at each level, built from level-uniform structural laws and one level-independent signed witness -- a closed-form class function $y^*$ with $N_k^T\Gamma(y^*)N_k=u_k u_k^T$ at every level. The mechanism forces the sign: for $k\ge 3$ no quantum state and no smooth curve of quantum models can realize the gain direction, so the overshoot lives strictly in the non-quantum part of the NPA tangent cone. The proof is computer-assisted in the strict sense: finite exact-integer verifications with proven degree bounds are constituent parts of the argument; the chain has been re-verified against independent implementations, including a symbolic per-regime proof of the witness identity and a clean-room implementation written from the paper text alone. The one external input is the published quantum value of Gigena et al., cross-checked to twelve digits.

**Comments field:**
21 pages. Verification code (12 exact-arithmetic verifiers, one-command reproduction) at https://github.com/tohafrit/npa-nonexactness

**Primary category:** quant-ph
**Cross-list:** math.OC
**License:** arXiv.org perpetual, non-exclusive license (the default)

---

## PAPER 2 (submit right after): the companion

**File to upload:** companion/npa_overshoot_note.tex

**Title:**
A phase transition in the exactness of the NPA hierarchy at the critical doubly-tilted CHSH functional

**Authors:** Anton Pakhunov

**Abstract:** use the .tex abstract, flattened; if over 1920 chars, cut the certificate-method sentence.

**Comments field:**
Companion to "No finite level of the NPA hierarchy is exact for the doubly-tilted CHSH functional near the critical tilt". Certificates and verification code at https://github.com/tohafrit/npa-nonexactness

**Primary category:** quant-ph

---

## After both are announced
1. Note the two arXiv IDs (arXiv:26MM.NNNNN).
2. Put each ID into the other paper's \bibitem (proof paper's \bibitem{overshoot}; companion's \bibitem{proof}), recompile, submit as replacement (v2). Also add the IDs to the repo README.

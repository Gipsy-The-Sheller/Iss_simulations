"""issbench -- simulation and closed-form tooling for the Iss audit project.

Layout
------
entropy   h, H_FSS under both conventions (A = Xia Eq. 3, B = Part I/II), exact
          multinomial entropy, the convention ratio c
models    reversible Q matrices normalised to unit rate; spectral utilities
gamma     Yang (1994) discrete Gamma (mean + median) and the continuous
          Laplace transform of Part I Lemma 6
trees     Node/Newick, quartet/star builders, pairwise distances, the
          Felsenstein->Farris transfer path
exact     exact per-site law of H_i and exact E[I_ss] (closed form on star
          trees, pattern enumeration for small topologies)
design    uniform-Iss sampling density (Item 0)
alisim    AliSim driver
stats     Iss, p-distance, C-value, lambda
infer     IQ-TREE inference and the P_true criterion
"""
__version__ = "0.1.0"

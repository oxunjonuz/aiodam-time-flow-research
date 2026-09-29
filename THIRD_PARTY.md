# Third-party materials

The authors' CC BY 4.0 and MIT licenses do not relicense third-party materials.

## Gravitational-wave data

The strain files in `data/` and the published GW150914 sky map in `work/sky/`
originate from the LIGO/Virgo public data releases, not from AIODAM.
GWOSC states that its data are released under CC BY 4.0:
<https://gwosc.org/about/>. Please follow its acknowledgement guidelines:
<https://gwosc.org/acknowledgement/>.

- Event and strain release: <https://gwosc.org/eventapi/html/GWTC-1-confident/GW150914/v3/>
- Sky-map release: <https://dcc.ligo.org/LIGO-P1500227/public>
- Strain provenance and SHA-256 values: [DATA_PROVENANCE.md](data/DATA_PROVENANCE.md).

The redistributed strain and sky-map input files are unmodified research inputs.
Derived figures and measurement artifacts are identified separately in the paper.

## Papers and video

Third-party papers, their extracted full text and the full video transcript are
not redistributed in this package. References and retrieval links are provided
in [sources/README.md](sources/README.md). They remain in the local research
workspace; this packaging decision does not delete or alter those originals.

Dependencies such as Mathlib, NumPy and phenomxpy are not bundled. Their own
licenses apply when installed separately.

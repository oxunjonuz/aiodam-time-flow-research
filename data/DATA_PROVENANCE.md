# GW150914 strain files: download provenance

Downloaded on 2026-09-27 for the time-flow research task. The official
GWOSC GW150914-v3 catalog lists both 32-second, 16384-Hz HDF5 strain files:
https://gwosc.org/eventapi/json/GWTC-1-confident/GW150914/v3/

Direct connections from this host to `gwosc.org` timed out. The two files were
therefore downloaded from an institutional mirror hosted at the Albert Einstein
Institute (AEI), Hannover:
https://www.atlas.aei.uni-hannover.de/work/yifan.wang/ringdown/GW150914/maxisi-ringdown/

| File | Bytes | SHA-256 |
| --- | ---: | --- |
| `H-H1_GWOSC_16KHZ_R1-1126259447-32.hdf5` | 4068797 | `81040e1ecfaf40ffe15a5efc59dbc3a888653162f1613425f1e68d0828dd1b97` |
| `L-L1_GWOSC_16KHZ_R1-1126259447-32.hdf5` | 3919118 | `23a207023ef8b49cf0b70b78d7ba80635881799324e50bfb2165c60039e41b1f` |

Both transfers completed successfully; byte lengths match the AEI server's
Content-Length headers. Both files have the HDF5 signature, and both are visible
inside the AIODAM container under `/work/time_flow_research_20260927/data/`.
Byte-for-byte equality with files served directly by GWOSC has **not** been
checked because this host could not connect to GWOSC. Do not treat the mirror
copy as independently authenticated by GWOSC until that comparison is possible.
These are publicly released calibrated strain time series, not raw detector
electronics data.

---

## AMENDMENT, 2026-09-27 (turn 220): the GWOSC caveat is LIFTED

The owner downloaded the same two files directly from GWOSC on their own host and
compared hashes. They match one for one. Independently of that, this agent
re-measured the local bytes with `file_fingerprint`:

| File | Bytes | SHA-256 (local) | SHA-256 (GWOSC, owner) |
| --- | ---: | --- | --- |
| `H-H1_GWOSC_16KHZ_R1-1126259447-32.hdf5` | 4068797 | `81040e1ecfaf40ffe15a5efc59dbc3a888653162f1613425f1e68d0828dd1b97` | `81040e1e…dd1b97` |
| `L-L1_GWOSC_16KHZ_R1-1126259447-32.hdf5` | 3919118 | `23a207023ef8b49cf0b70b78d7ba80635881799324e50bfb2165c60039e41b1f` | `23a20702…e41b1f` |

All three values agree: the AEI mirror copy is byte-identical to the official
GWOSC files, and the sizes match (4068797 and 3919118 bytes). The earlier
caveat — "the mirror is not authenticated by GWOSC" — is therefore **withdrawn**.
It was a provenance/legal question, not a scientific one; the owner's point that
the data were already functionally authenticated by three independent physical
checks (measured ASD at 100 Hz matching published aLIGO O1, peak time matching
the published merger time, and H1-L1 agreement to 0.37 ms) stands on its own.

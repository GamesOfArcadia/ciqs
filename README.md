# CIQS: Causal Iteration Quantum Solution
# The 'million-qubit compiler' and The '1:5 QEC'

by [arcadialab.fr](https://arcadialab.fr)

---

**CIQS** is composed of: **CIQA** (error correction), **CIQM** (mapper & router), and **CIQO** (optimizer)

CIQA is optional in the pipeline, not all architectures need error correction ([See matter qubits Vs photon qubits](https://zenodo.org/records/19143454)).

---

```
1. Description
2. Files
3. Requirements
4. Usage
- 4.1 With CIQA
- 4.2 Without CIQA
- 4.3 Topologies
- 4.4 Get Results
5. CIQA External pipeline integration
6. Performances
- 6.1 CIQS - IBM Benchpress
- 6.2 CIQA - hardware validation
- 6.3 Hayden-Preskill - combined pipeline validation
- 6.4 Concatenation
7. Published work
8. Contact
9. License
```
 
---

## Description

**CIQA** is an analytic, hardware-agnostic quantum error correction code. Each logical qubit is encoded into 5 physical data qubits based on the [[5,1,3]] perfect quantum code (the distance parameter doesn't apply to CIQA), with 2 ancilla quDits per logical block for syndrome extraction. CIQA is not a stabilizer code. Its protection is not built from redundancy defeating random errors by majority, it comes from two analytically derived geometric constants which define closed-form error floors for two independent parity channels. An error either crosses a channel's geometric floor or it doesn't. There is no notion of weight, no minimum-distance decoding, and no threshold that scales with block size in the stabilizer sense.

**CIQM + CIQO** form a fully analytic transpilation pipeline. Every placement, routing, and optimization decision is determined by an exact criterion derived from the circuit structure. No heuristics, no stochastic components, no tunable parameters. The same input produces the same output on every run. Hardware-agnostic, supports qudits of any dimension natively.

CIQA encodes first. CIQO receives the encoded circuit, isolates the encoding prefix automatically, routes across the target topology, and CIQO applies optimization exclusively to the operational gates outside the protected prefix. A single pipeline from logical circuit to hardware-ready output. *See section 6 for a performance summary and section 7 for all data and published papers.*

Both tools are compiled binaries, free for non-commercial use. Contact for licensing options.

---

## Files

In order of execution:

| File | Role |
|---|---|
| `CIQA.pyd` | QEC encoder: 1 logical qubit to 5 physical + 2 ancilla |
| `CIQS_pipeline.cp312-win_amd64.pyd` | Pipeline entry point: runs CIQM and CIQO in sequence |
| `CIQM.cp312-win_amd64.pyd` | Mapper and router: qubit placement and SWAP insertion |
| `CIQO_v2.cp312-win_amd64.pyd` | Optimizer: gate removal on operational gates only |
|---|---|
| `CIQS_adapters.cp312-win_amd64.pyd` | Format adapters: Qiskit, OpenQASM 3, generic gate list |
| `CIQA_integration_guide.pdf` | Standalone adapter functions for external pipeline integration |

---

## Requirements

- Python 3.9+
- Qiskit >= 1.0 (for IBM hardware runs and the Qiskit adapter)
- pip install numpy scipy

Place all `.pyd` modules in your working directory or add them to your Python path.

---

## Usage

### With CIQA

```python
from CIQA import CIQA
from CIQS_pipeline import CIQS_pipeline, CIFPipelineResult

ciqa = CIQA(n=5, d=2)
encoded = ciqa.encode(logical_circuit)

result = CIQS_pipeline(encoded, CouplingMap.heavy_hex(127))
```

### Without CIQA

Everything goes through CIQS_pipeline. No need to call CIQM or CIQO directly.

CIQS_pipeline accepts circuits directly without QEC encoding.
*coupling_map is required for all input types except Perceval*.

```python
from CIQS_pipeline import CIQS_pipeline, CIFPipelineResult

# Qiskit
result = CIQS_pipeline(qiskit_circuit, CouplingMap.line(1000))

# OpenQASM 2.0
with open("circuit.qasm") as f:
    qasm_str = f.read()
result = CIQS_pipeline(qasm_str, CouplingMap.heavy_hex(127))

# Perceval (coupling map auto-derived from processor topology)
result = CIQS_pipeline(perceval_processor)

# Blackbird / Strawberry Fields
result = CIQS_pipeline(blackbird_program, CouplingMap.line(100))
```

---

### Topologies

```python
from CIQM import CouplingMap

CouplingMap.line(n)
CouplingMap.heavy_hex(n)
CouplingMap.grid(rows, cols)
CouplingMap.from_edges(edge_list)
```

### Get Results

```python
print(result.swap_count)        # SWAPs inserted by router
print(result.gates_removed)     # gates removed by optimizer
print(result.placement.valid)   # True = all qubit pairs placed natively

# Export
qc       = result.to_qiskit()
qasm     = result.to_openqasm()
processor = result.to_perceval()
```

For a complete end-to-end example including QPU submission, see `page_curve_ibm_experiment_v5.py`.

---

### CIQA External pipeline integration

`ciqa.encode()` returns a dict with four keys:

| Key | Type | Description |
|---|---|---|
| `gates` | `list[dict]` | Ordered list of physical gate operations |
| `n_logical` | `int` | Number of logical qubits in the input circuit |
| `n_physical` | `int` | Total physical qubits after encoding (`n_logical x 5` + ancillae) |
| `report` | `EncodingReport` | Encoding parameters for auditing |

Three gate types are produced: `Ry` (single-qubit rotation), `CX` (controlled-X), and `M` (measurement). Measurement gates carry `_ciqa_physical: True` to identify CIQA-inserted ancilla operations. Any pipeline supporting single-qubit rotations and two-qubit controlled gates can accept CIQA output directly.

Adapters for Qiskit, OpenQASM 3, and a generic gate list format are in `CIQA_integration_guide.pdf`.

**Coupling map sizing at n=5:**

| Logical qubits | Data qubits (x5) | Ancillae (x2) | Minimum total |
|---|---|---|---|
| 1 | 5 | 2 | 7 |
| 10 | 50 | 20 | 70 |
| 100 | 500 | 200 | 700 |
| 1,000 | 5,000 | 2,000 | 7,000 |

The downstream router may require additional qubits for SWAP insertion depending on hardware topology.

**Hard constraint:** CIQA encoding angles are exact fixed values. Gate fusion or small-angle removal applied to the encoding prefix breaks the error correction guarantee. Pass the encoded circuit to the router without a pre-routing simplification pass. Apply gate optimization only after routing and only to gates outside the encoding prefix. *When using CIQS, this is handled automatically*.

---

## Performance

### CIQS - IBM Benchpress

The full IBM Benchpress suite (892 tests, circuits up to 930 qubits) was run with CIQS on a Ryzen 5 / 32GB machine, against Qiskit's published results on a Ryzen 9 / 128GB machine.

| SDK | Passed / Attempted | Failed | Wall time |
|---|---|---|---|
| **CIQS** | **892 / 892** | **0** | **1h 15m** |
| Qiskit | 892 / 892 | 0 | 17h 15m |
| QTS | 877 / 892 | 15 | >= 35h 18m |
| Tket | 814 / 892 | 78 | >= 33h 11m |
| BQSKit | 716 / 892 | 176 | >= 64h 08m |

Per-test geometric mean: CIQS at 0.051s, Qiskit at 0.310s. In the 200-433 qubit range, CIQS is 14.7x faster than Qiskit.

**Scaling** (Heisenberg model, same hardware):

| Qubits | Topology | 2Q Gates | SWAPs | Gates removed | Runtime |
|---|---|---|---|---|---|
| 100,000 | linear | 499,995 | 0 | 299,997 | 9.6s |
| 100,000 | heavy-hex | 739,658 | 239,663 | 299,997 | 32.5s |
| 1,000,000 | linear | 4,999,995 | 0 | 2,999,997 | 93.4s |
| 1,000,000 | heavy-hex | 7,513,506 | 2,513,511 | 2,999,997 | 848.8s |

On the linear topology, runtime scales proportionally with qubit count across three orders of magnitude.

### CIQA - hardware validation

IBM ibm_fez (Heron r2, 156 qubits), April 3, 2026.
Three independent runs, 8,192 shots each, pinned layout, live calibration data.

| Circuit | Mean fidelity | SD |
|---|---|---|
| Bare qubit, no error | 0.994 | 0.002 |
| Bare qubit, X error, no correction | 0.015 | 0.002 |
| CIQA encoded, no error | 0.818 | 0.011 |
| CIQA, X error, corrected (q0 and q2) | 0.847 | 0.008 |
| Mean improvement over unprotected | **+0.832** | 0.008 |


### Hayden-Preskill - combined pipeline validation

The Hayden-Preskill circuit models black hole evaporation as a quantum system across multiple entangled registers over several steps. No experiment ever successfully ran this circuit, until now.
Paired with CIQA, CIQS compiled and routed the circuit onto IBM Heron r2, a heavy-hex topology that never was a design target for it.
The routing required 987 SWAP insertions and 120 gate removals across 300 circuits per run. The experiment ran across three IBM backends at two system sizes, with up to 307,200 shots per job.
The mutual information I₂(D:R) remained within ±0.004 of zero across all runs and all time steps, indicating no measurable noise introduced by the compilation.
CIQA prevented the decoherence to happen early on, allowing the very first observation of the black hole evaporation past the Page time.
The results have been reproduced a dozen times.

---

## Concatenation

Each physical qudit of a CIQA block can itself be encoded in another CIQA block, reducing the logical error rate at the cost of additional overhead.

| Level | Overhead | p_L at p = 0.1% |
|---|---|---|
| 1 | 1:5 | ~10^-5 |
| 2 | 1:25 | ~10^-9 |
| 3 | 1:125 | ~10^-17 |

Level 3 is sufficient to run Shor's algorithm. The surface code requires approximately 1:1,000 overhead to reach equivalent logical error rates at the same physical error rate.

---

## Published work

- **CIQS, architecture and benchmarks:** [10.5281/zenodo.19056796](https://zenodo.org/records/19056796)
- **CIQS vs Qiskit, full Benchpress comparison:** [10.5281/zenodo.19896662](https://zenodo.org/records/19896662)
- **CIQA, architecture and benchmarks:** [10.5281/zenodo.19405503](https://zenodo.org/records/19405503)
- **Hayden-Preskill black hole circuit, full dataset:** [10.5281/zenodo.21326186](https://zenodo.org/records/21326186)

**All benchmarks data, runners, QPU side data, are attached alongside their respective paper**

---

## Contact
If something is missing, not clear, or not working as described, please let us know.

For every enquiries including licensing options:

jessy.pensedent@gmail.com

---

## License

Free for non-commercial use. See `LICENSE` for terms.

(c) Jessy Pensédent - [arcadialab.fr](https://arcadialab.fr) 

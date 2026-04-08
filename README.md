# CIQS — Causal Iteration Quantum Solution

CIQS is an analytic quantum circuit pipeline including: mapper, router, and optimizer.

The pipeline natively handles quDits of any dimension ( qubit, qutrit, ququart, and above ) with no additional parameters required.
The optimizer scans the mapped and routed circuit once and applies three removal passes in sequence:

- Pass1 removes individual redundant gates,
- Pass2 removes gate sequences impossible to remove individually,
- Pass3 removes gates the target device cannot execute faithfully, at zero fidelity cost, based on vendor-supplied hardware calibration data.

CIQS handles the full pipeline from raw circuit, to hardware connection.

No heuristics, no pattern matching, no approximation, no simulation, no statevector.
Zero fidelity loss by construction.

CIQA (1:5 quantum error correction) is not included in this distribution. Contact for licensing.
---
---

## IBM Benchpress Results

All numbers below are from IBM Benchpress, made on an AMD Ryzen 5, 32Go RAM, no GPU.

Circuit family: Heisenberg model. Topologies: linear chain and IBM heavy-hex.

| Qubits    | Topology  | Time (s) | Gates Removed | SWAPs     |
|-----------|-----------|----------|---------------|-----------|
| 100,000   | linear    | ~9       | 299,997       | 0         |
| 100,000   | heavy-hex | ~31      | 299,997       | 239,663   |
| 500,000   | linear    | ~47      | 1,499,997     | 0         |
| 500,000   | heavy-hex | ~303     | 1,499,997     | 1,302,621 | 
| 1,000,000 | linear    | ~96      | 2,999,997     | 0         |
| 1,000,000 | heavy-hex | ~809     | 2,999,997     | 2,513,511 |

**Gate removal: ~30% consistently across all scales.**  
Linear scaling confirmed from 100,000 to 1,000,000 qubits.

---

## Published Results

- CIQS pipeline paper: [10.5281/zenodo.18750722](https://zenodo.org/records/19056796)
- Benchpress raw data included in this repository

---

## Requirements

```
pip install numpy scipy
```

Python 3.10, 3.11, or 3.12. Linux x86_64 or Windows x64.

---

## Usage

Everything goes through `CIQS_pipeline`. No need to call CIQM or CIQO directly.

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

# Blackbird Strawberry Fields (input only — to_blackbird not exposed. My bad)
result = CIQS_pipeline(blackbird_program, CouplingMap.line(100))


# Results
print(result.swap_count)        # SWAPs inserted by router
print(result.gates_removed)     # gates removed by optimizer
print(result.placement.valid)   # True = all qubit pairs placed natively


# Export
qc = result.to_qiskit()
qasm = result.to_openqasm()
processor = result.to_perceval()
```

`coupling_map` is required for all input types except Perceval.

---

## Supported Input Formats

- Qiskit `QuantumCircuit`
- OpenQASM 2.0 string
- Perceval `Processor` (coupling map auto-derived)
- Blackbird Strawberry Fields

---

## Topologies

```python
CouplingMap.line(n)
CouplingMap.heavy_hex(n)
CouplingMap.grid(rows, cols)
CouplingMap.from_edges(edge_list)
```

`coupling_map` is required for all input types except Perceval.

---

## License

*Free for academic and non-commercial research use with attribution.    
Commercial use requires a separate agreement.  
See LICENSE.txt for terms*

---

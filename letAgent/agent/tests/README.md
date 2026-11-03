# Testing the agent without MATLAB

`mock_engine.py` stands in for MATLAB/Octave. It answers `letAgent_prepareCase` and
`letAgent_runBatch` calls with a toy dose / LET / robustness trade-off model
(no physics, no patient data), so the search, rules, Pareto front and report can be
tested in seconds.

From `letAgent/agent`:

    python -m letagent search --case tests/CNS_mock_case.json --engine octave --engine-cmd "$(pwd)/tests/mock_engine.py" --budget 40 --workers 2

`mock_report_example.html` is the report this produced (synthetic numbers only).

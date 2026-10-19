# Regression baselines

Each `<case>.mat` file holds the outputs (`expected`) of one case in
`tests/helpers/regression_cases.m`, plus provenance (`meta`: commit,
MATLAB version, k-Wave path, date).

**Reference code version:** upstream `PROTEUS-SIM/PROTEUS` branch
`transducer-calibration` at `6a863b4`, i.e. the code base *before* the
revamp. Baselines must only be (re)recorded from that reference, never from
a version under test:

```matlab
% in a checkout/worktree of the reference commit with tests/ copied in
cd tests
generate_baseline()                       % all cases whose requirements are met
generate_baseline({'rf_hybrid_standard'}) % selected cases
```

Cases that need k-Wave and the Statistics and Machine Learning Toolbox
(`rf_*`, `pressure_field`) are skipped when those are unavailable; their
baseline files are added once recorded.

An intentional change of behaviour (a bug fix, a new default) must update
the affected baseline in the same commit, and the commit message must say
why.

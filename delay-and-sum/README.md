# Delay-and-sum reconstruction

The delay-and-sum (DAS) code is independent of the simulation pipeline: it
reads the RF data and the settings file written by PROTEUS and returns an
image. Any beamformer with the same interface can be used in its place.

## Interface

```matlab
[IMG, z, x] = beamformer(folder, filename_parameters)
```

- `folder`: folder with the RF data (`load_RF_data`).
- `filename_parameters`: settings file (`Acquisition`, `SimulationParameters`,
  `Transmit`, `Transducer`, `Geometry`, `Medium`).
- `IMG`: log-compressed image, `length(x)`-by-`length(z)`-by-`Nframes` [dB].
- `x`, `z`: lateral and axial pixel coordinates [m].

## Beamformers

| Function | Arrays | Transmits | Status |
|---|---|---|---|
| `DAS_reconstruction` + `compute_das_matrix` | linear | plane wave (no steering), diverging | reference |
| `DAS_reconstruction_general` + `initialize_das_parameters` + `compute_das_matrix_general` | linear, curved | plane wave and diverging (steered), focused, curved diverging | under test against the reference (`tests/DelayAndSumTest.m`) |

`compute_das_matrix_general` uses `das_transmit_distance` for the transmit
path, which matches the PROTEUS transmit delays (`compute_delays`) for every
transmit type, and `get_element_geometry` (in `GUIfunctions`, add it to the
path) for the element positions.

With binary apodization and the same pixel grid and time axis,
`compute_das_matrix_general` reproduces `compute_das_matrix` exactly for a
plane wave, and within one sample for diverging waves (the reference
assumes an element at the array centre).

To add a beamformer, write a function with the interface above and pass it
to `example_scripts/example_reconstruction.m` (`beamformer = @my_function`).

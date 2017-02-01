# wfirst-phot-sim

*Photometric-only SN light curve simulation scipts.*

**Requirements:** numpy, scipy, matplotlib, sncosmo, tqdm

In the future, this might be merged into the `wfirst-sim` repository.
Currently there is just a single script: `phot-sim`.

Run simulation, output files to `lcs` directory (A running iteration
total is shown):

```
$ phot-sim lcs
7131it [02:22, 50.09it/s]
```

Note that the total number of light curves saved is less than the
number of iterations (SNe generated). This is because there is a
minimum signal-to-noise ratio required for a light curve to be saved:
we don't bother saving low S/N light curves that would never be
detected.

Plot bandpasses to a file, just to check:

```
phot-sim --plotbands bands.png
```


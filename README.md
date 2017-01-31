# wfirst-phot-sim

*Photometric-only SN light curve simulation scipts.*

**Requirements:** numpy, scipy, matplotlib, sncosmo, tqdm

In the future, this might be merged into the `wfirst-sim` repository.
Currently there is just a single script:

```
phot-sim lcs  # run simulation, output files to "lcs" directory
phot-sim --plotbands bands.png  # plot bandpasses to a file.
```

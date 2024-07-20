Shyam Kumar; INFN Bari, Italy; shyam055119@gmail.com
Method to produce the tracking performances with ePIC tracker
The scripts can be used to create the debug plots for the momentum resolutions.

To run a full simulation-reconstruction-analysis chain do:
```
snakemake -c2 results/tracking_performances/local
```
or, referencing it by the rule name:
```
snakemake --cores 1 tracking_performance_local
```

To process an individual campaign do:
```
snakemake -c2 results/tracking_performances/24.04.0
```

It will produce the results with truth/realistic seeding.

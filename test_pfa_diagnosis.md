# PFA Plugin Diagnosis

## Issue Summary
PFA plugins (pandora and arbor) load successfully but do not produce output collections.

## Test Configuration
- **EICrecon version:** v1.40.0-183-gb47dade12
- **Input file:** sim_dis_5x41_minQ2=1_epic_craterlake.edm4hep.root
- **Command:** `eicrecon input.root -Pplugins=pandora,arbor`

## Observations

### ✅ Plugins Load Successfully
```
19:50:06.536  [info] Loading plugin 'pandora' from '/opt/local/lib/EICrecon/plugins/pandora.so'
19:50:06.537  [info] Loading plugin 'arbor' from '/opt/local/lib/EICrecon/plugins/arbor.so'
```

### ✅ Plugins Added to Plugin List
```
value: log,dd4hep,...,pandora,arbor
```

### ✅ Required Input Collections Exist
- EcalBarrelImagingRecHits
- EcalBarrelScFiRecHits
- EcalEndcapNRecHits
- EcalEndcapPRecHits
- HcalBarrelRecHits
- HcalEndcapNRecHits
- CalorimeterTrackProjections

### ❌ Expected Output Collections NOT Created
- PandoraPFAParticles (MISSING)
- ArborPFAParticles (MISSING)

### ❌ No Factory Registration Messages
- No log entries mentioning "PandoraPFAParticles" or "ArborPFAParticles"
- No factory initialization messages
- No errors about missing factories

## Suspected Issues

1. **Factories not registered:** The plugins may load but fail to register their factories
2. **Missing XML configuration:** PFA algorithms may require PandoraPFASettings.xml / ArborPFASettings.xml
3. **Collection name mismatch:** The factory names may differ from expected collection names
4. **Silent initialization failure:** Factory constructors may be failing silently

## Next Steps for EICrecon Session

1. Verify factory registration in plugin source code
2. Check if XML configuration files are required and present
3. Add debug logging to confirm factory Init() is called
4. Verify output collection names match factory declarations
5. Check for any conditional compilation or runtime guards

## Files for Reference
- Plugin libraries: /opt/local/lib/EICrecon/plugins/pandora.so, arbor.so
- Test log: eicrecon_test3.log
- Test output: test_pfa_rec3.edm4eic.root (5.7 MB, no PFA collections)

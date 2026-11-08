MCsquare scanner consistent with matRad's pencil beam (matRad/hluts/matRad_default.hlut):
all voxels water, density = matRad relative stopping power. Used by UHN_CNS_07 (Nov 8 2026)
because MCsquare's matRad_default scanner (all water, older linear table: HU 1200 -> 1.83,
HU 1800 -> 2.32) gives ~10-30 % more stopping in bone than the pencil-beam HLUT
(HU 1200 -> 1.62, HU 1800 -> 1.80), so MCsquare and the pencil beam had different ranges
behind the skull / petrous bone. Interim: replace both with the UHN RayStation CT calibration.

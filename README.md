# acs_flowthrough
Matlab codes for processing data measured by underway AC-S flow-through system. Use these functions in the order of:

1) acs_rd_wetview.m
2) acs_rm_spikes.m
3) acs_filemerge_1minbin.m
4) acs_TempSalCorr.m
5) acs_calc_apcp.m
6) acs_ResidTempScatCorr.m or acs_ResidTempScatCorr_RR.m

References:
1. Liu et al. (2018)
Underway spectrophotometry in the Fram Strait (European Arctic Ocean): a highly resolved chlorophyll a data source for complementing satellite ocean color. Optics Express, 26(14), A678-A696.
2. Liu et al. (2019)
Retrieval of Phytoplankton pigments from underway spectrophotometry in the Fram Strait. Remote Sensing, 11(3), 318.
3. Literature cited by the above two papers.

Alternative codes can be found at https://github.com/OceanOptics/InLineAnalysis, https://github.com/OceanOptics/ACCode, https://github.com/OceanOptics/pyACS.

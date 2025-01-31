# EDDI (Evaporative Demand Drought Index)

This repository contains Python code for calculating the Evaporative Demand Drought Index (EDDI), a physically based drought monitoring and early warning guidance tool that measures drought dynamics based on atmospheric evaporative demand (E₀).

![eddi](/images/idn_cli_agera5_eddi_1month_20000630.png)

## Overview

EDDI examines the signal of actual or potential drought conditions by measuring the anomalous atmospheric evaporative demand using reference evapotranspiration (E₀). This approach complements traditional land-surface perspective drought indicators by providing an independent assessment of drought conditions from a purely atmospheric perspective.

## Input Data Requirements

The code expects the following input variables:

AgERA5 - https://cds.climate.copernicus.eu/datasets/sis-agrometeorological-indicators?tab=download
- 2m air temperature (minimum, maximum, mean)
- 2m dew point temperature
- 10m wind speed
- Surface solar radiation

SRTM
- Digital elevation model

All meteorological variables should be in NetCDF format with dimensions (time, lat, lon).

## Methodology

### 1. Reference Evapotranspiration (ET₀) Calculation

The FAO-56 Penman-Monteith equation is used to calculate daily reference evapotranspiration (ET₀). This equation combines both energy balance and mass transfer principles to estimate ET₀:

$ET_0 = \frac{0.408\Delta(R_n - G) + \gamma\frac{900}{T + 273}u_2(e_s - e_a)}{\Delta + \gamma(1 + 0.34u_2)}$

where:
- $ET_0$ = reference evapotranspiration [mm day⁻¹]
- $R_n$ = net radiation at the crop surface [MJ m⁻² day⁻¹]
- $G$ = soil heat flux density [MJ m⁻² day⁻¹]
- $T$ = mean daily air temperature at 2 m height [°C]
- $u_2$ = wind speed at 2 m height [m s⁻¹]
- $e_s$ = saturation vapor pressure [kPa]
- $e_a$ = actual vapor pressure [kPa]
- $e_s - e_a$ = saturation vapor pressure deficit [kPa]
- $\Delta$ = slope vapor pressure curve [kPa °C⁻¹]
- $\gamma$ = psychrometric constant [kPa °C⁻¹]

#### Net Radiation Calculation

Net radiation ($R_n$) is calculated as the difference between net shortwave and net longwave radiation:

$R_n = R_{ns} - R_{nl}$

where:
- $R_{ns}$ = net shortwave radiation = $(1 - \alpha)R_s$
- $\alpha$ = albedo (0.23 for reference crop)
- $R_s$ = incoming solar radiation [MJ m⁻² day⁻¹]
- $R_{nl}$ = net longwave radiation [MJ m⁻² day⁻¹]

Net longwave radiation is calculated using:

$R_{nl} = \sigma \left(\frac{T_{max,K}^4 + T_{min,K}^4}{2}\right)(0.34 - 0.14\sqrt{e_a})\left(1.35\frac{R_s}{R_{so}} - 0.35\right)$

where:
- $\sigma$ = Stefan-Boltzmann constant [4.903 × 10⁻⁹ MJ K⁻⁴ m⁻² day⁻¹]
- $T_{max,K}, T_{min,K}$ = maximum and minimum absolute temperature [K]
- $R_{so}$ = clear-sky radiation [MJ m⁻² day⁻¹]

### 2. EDDI Calculation Process

The EDDI calculation involves the following steps:

1. **Daily ET₀ Calculation**:
   - Calculate daily ET₀ using the Penman-Monteith equation
   - Process input variables: temperature, humidity, wind speed, solar radiation

2. **Temporal Aggregation**:
   - Sum ET₀ over specified time scales (e.g., 1, 3, 6 months)
   - Use dekad-based windows for consistent temporal aggregation

3. **Climatological Analysis**:
   - Compare current ET₀ sum against 1991-2020 climatology
   - Calculate percentiles using empirical distribution
   - Rank current conditions within historical context

4. **EDDI Classification**:
   The final EDDI value is expressed as a percentile, classified into categories:

| Percentile Range | Category | Color (HEX) | RGB |
|-----------------|----------|-------------|-----|
| 98-100 | ED4 (Exceptional Drought) | #67001F | rgb(103, 0, 31) |
| 95-98  | ED3 (Extreme Drought)     | #B2182B | rgb(178, 24, 43) |
| 90-95  | ED2 (Severe Drought)      | #D6604D | rgb(214, 96, 77) |
| 80-90  | ED1 (Moderate Drought)    | #F4A582 | rgb(244, 165, 130) |
| 70-80  | ED0 (Abnormally Dry)      | #FDDBC7 | rgb(253, 219, 199) |
| 30-70  | Normal                    | #F7F7F7 | rgb(247, 247, 247) |
| 20-30  | EW0 (Abnormally Wet)      | #D1E5F0 | rgb(209, 229, 240) |
| 10-20  | EW1 (Moderately Wet)      | #92C5DE | rgb(146, 197, 222) |
| 5-10   | EW2 (Severely Wet)        | #4393C3 | rgb(67, 147, 195) |
| 2-5    | EW3 (Extremely Wet)       | #2166AC | rgb(33, 102, 172) |
| 0-2    | EW4 (Exceptionally Wet)   | #053061 | rgb(5, 48, 97) |

### Code Output

The code produces:
1. NetCDF files containing EDDI values for each time scale
2. Visualization plots showing:
   - Input variables (temperature, wind speed, solar radiation)
   - Calculated ET₀
   - Final EDDI percentiles with drought/wet categories

## References

1. Hobbins, M. T., A. Wood, D. J. McEvoy, J. L. Huntington, C. Morton, M. Anderson, and C. Hain (2016): The Evaporative Demand Drought Index. Part I: Linking drought evolution to variations in evaporative demand. J. Hydrometeor., 17(6), 1745-1761, [doi:10.1175/JHM-D-15-0121.1](https://doi.org/10.1175/JHM-D-15-0121.1)

2. McEvoy, D. J., J. L. Huntington, M. T. Hobbins, A. Wood, C. Morton, M. Anderson, and C. Hain (2016): The Evaporative Demand Drought Index. Part II: CONUS-wide assessment against common drought indicators. J. Hydrometeor., 17(6), 1763-1779, [doi:10.1175/JHM-D-15-0122.1](https://doi.org/10.1175/JHM-D-15-0122.1)

3. NOAA Physical Sciences Laboratory EDDI Documentation:
   - [EDDI Homepage](https://psl.noaa.gov/eddi/)
   - [EDDI User Guide](https://psl.noaa.gov/eddi/documentation.html)

4. Allen, R. G., Pereira, L. S., Raes, D., & Smith, M. (1998). Crop evapotranspiration-Guidelines for computing crop water requirements. FAO Irrigation and drainage paper 56, FAO, Rome, 300(9), D05109.

## Acknowledgments

- NOAA Physical Sciences Laboratory for developing and documenting the EDDI methodology
- The World Bank GOST/DEC Data Group for supporting this implementation
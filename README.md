# EDDI (Evaporative Demand Drought Index)

This repository contains Python code for calculating the Evaporative Demand Drought Index (EDDI), a physically based drought monitoring and early warning guidance tool that measures drought dynamics based on atmospheric evaporative demand (E₀).

![eddi](/images/idn_cli_agera5_eddi_3month_20151231.png)
Figure 1. EDDI 3-month as of 31 December 2015 (Very Strong El-Niño)

## Overview

EDDI examines the signal of actual or potential drought conditions by measuring the anomalous atmospheric evaporative demand using reference evapotranspiration (E₀). This approach complements traditional land-surface perspective drought indicators by providing an independent assessment of drought conditions from a purely atmospheric perspective.

## Case study

[Sumba Island](https://en.wikipedia.org/wiki/Sumba), located in East Nusa Tenggara, Indonesia, is a captivating destination known for its unique blend of natural beauty and traditional culture. 

Despite its stunning landscapes of rolling hills and pristine beaches, the island faces significant climate challenges, particularly during the dry season from May to October. The region experiences a distinct monsoon climate with recurring drought events that affect local agriculture and water resources. The island's average annual rainfall varies significantly between the northern and southern regions, with some areas receiving less than 1000mm of rainfall annually. 

Despite these environmental challenges, Sumba has gained international recognition for sustainable luxury tourism, most notably through Nihi Sumba (formerly Nihiwatu), which was voted the best hotel in the world by Travel + Leisure magazine in 2016 and 2017. Located on the island's remote southwestern coast, the resort exemplifies how tourism can thrive while respecting both environmental constraints and local traditions. The island's combination of dramatic coastlines, traditional villages, and unique Marapu culture continues to draw visitors, even as it grapples with climate-related challenges.

## Input Data Requirements

> **Important**: This documentation focuses on the EDDI calculation methodology and implementation. The preliminary steps of downloading AgERA5 data from the Copernicus Climate Data Store (CDS), preprocessing individual variables, and merging them into the required format are not discussed here. These data preparation steps should be completed separately following CDS guidelines and data handling best practices. For information about data acquisition and preprocessing, please refer to the CDS documentation and API guidelines.

### AgERA5 Data

The code uses AgERA5 (Agrometeorological indicators from 1979 to present derived from reanalysis) available from the [Copernicus Climate Data Store](https://cds.climate.copernicus.eu/cdsapp#!/dataset/sis-agrometeorological-indicators):

Required variables:
- 2m air temperature (minimum, maximum, mean) [K]
- 2m dew point temperature [K]
- 10m wind speed [m s⁻¹]
- Solar radiation flux [J m⁻² day⁻¹]

### Additional Data
- Digital elevation model (e.g., SRTM) [m]

All meteorological variables should be in NetCDF format with dimensions (`time, lat, lon`).

## Methodology

### 1. Data Pre-processing and Adjustments for FAO-56

Several adjustments are required to use AgERA5 data in the FAO-56 Penman-Monteith equation:

1. Temperature Conversions:
   - AgERA5 provides all temperatures in Kelvin
   - Conversion to Celsius required for all temperature variables:
     ```python
     T[°C] = T[K] - 273.15
     ```
   - Applied to:
     * Mean air temperature
     * Maximum air temperature
     * Minimum air temperature
     * Dew point temperature

2. Wind Speed Height Correction:
   - AgERA5 provides wind speed at 10m height
   - FAO-56 requires wind speed at 2m height
   - Conversion using logarithmic wind profile:
     ```python
     z = 10.0  # measurement height in meters
     u2 = u10 * (4.87 / np.log(67.8 * z - 5.42))
     ```

3. Solar Radiation Processing:
   - AgERA5 provides Solar_Radiation_Flux in J m⁻² day⁻¹
   - Conversion to MJ m⁻² day⁻¹ required:
     ```python
     Rs = solar_rad_j * 1.0e-6  # [MJ m⁻² day⁻¹]
     ```

### 2. Reference Evapotranspiration (ET₀) Calculation with AgERA5

The FAO-56 Penman-Monteith equation using adjusted AgERA5 data:

$ET_0 = \frac{0.408\Delta(R_n - G) + \gamma\frac{900}{T + 273}u_2(e_s - e_a)}{\Delta + \gamma(1 + 0.34u_2)}$

Key components calculated from AgERA5:

1. Net Radiation ($R_n$) Components:
   - Not directly available in AgERA5
   - Calculated in two parts:

   1. Net shortwave radiation:
      ```python
      ALBEDO = 0.23  # reference crop
      R_ns = (1.0 - ALBEDO) * Rs  # [MJ m⁻² day⁻¹]
      ```

   2. Net longwave radiation:
      ```python
      R_nl = σ * ((T_kmax⁴ + T_kmin⁴)/2) * (0.34 - 0.14√e_a) * (1.35 * Rs/R_so - 0.35)
      ```
      where:
      - T_kmax, T_kmin are from AgERA5 max/min temperatures
      - e_a is calculated from AgERA5 dew point
      - R_so is clear-sky radiation estimated from elevation and extraterrestrial radiation

2. Vapor Pressure Components:
   - Actual vapor pressure (e_a) from AgERA5 dew point:
     ```python
     ea = 0.6108 * np.exp((17.27 * tdew_c) / (tdew_c + 237.3))  # [kPa]
     ```
   - Saturation vapor pressure (e_s) from AgERA5 max/min temperatures:
     ```python
     es_max = 0.6108 * np.exp((17.27 * tmax_c) / (tmax_c + 237.3))
     es_min = 0.6108 * np.exp((17.27 * tmin_c) / (tmin_c + 237.3))
     es = (es_max + es_min) / 2  # [kPa]
     ```

3. Other Components:
   - Soil heat flux (G) assumed zero for daily calculations
   - Psychrometric constant (γ) calculated using elevation-adjusted pressure
   - Temperature coefficient uses AgERA5 mean temperature converted to Kelvin


### 3. EDDI Calculation Process

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

| Percentile Range | Category | Color | HEX | RGB |
|-----------------|----------|-----|-------------|-----|
| 98-100 | ED4 (Exceptional Drought) | <span style="background-color: #67001F; padding: 0 10px;">&nbsp;&nbsp;&nbsp;</span> | #67001F | rgb(103, 0, 31) |
| 95-98  | ED3 (Extreme Drought)     | <span style="background-color: #B2182B; padding: 0 10px;">&nbsp;&nbsp;&nbsp;</span> | #B2182B | rgb(178, 24, 43) |
| 90-95  | ED2 (Severe Drought)      | <span style="background-color: #D6604D; padding: 0 10px;">&nbsp;&nbsp;&nbsp;</span> | #D6604D | rgb(214, 96, 77) |
| 80-90  | ED1 (Moderate Drought)    | <span style="background-color: #F4A582; padding: 0 10px;">&nbsp;&nbsp;&nbsp;</span> | #F4A582 | rgb(244, 165, 130) |
| 70-80  | ED0 (Abnormally Dry)      | <span style="background-color: #FDDBC7; padding: 0 10px;">&nbsp;&nbsp;&nbsp;</span> | #FDDBC7 | rgb(253, 219, 199) |
| 30-70  | Normal                    | <span style="background-color: #F7F7F7; padding: 0 10px;border: 1px solid #000;">&nbsp;&nbsp;&nbsp;</span> | #F7F7F7 | rgb(247, 247, 247) |
| 20-30  | EW0 (Abnormally Wet)      | <span style="background-color: #D1E5F0; padding: 0 10px;">&nbsp;&nbsp;&nbsp;</span> | #D1E5F0 | rgb(209, 229, 240) |
| 10-20  | EW1 (Moderately Wet)      | <span style="background-color: #92C5DE; padding: 0 10px;">&nbsp;&nbsp;&nbsp;</span> | #92C5DE | rgb(146, 197, 222) |
| 5-10   | EW2 (Severely Wet)        | <span style="background-color: #4393C3; padding: 0 10px;">&nbsp;&nbsp;&nbsp;</span> | #4393C3 | rgb(67, 147, 195) |
| 2-5    | EW3 (Extremely Wet)       | <span style="background-color: #2166AC; padding: 0 10px;">&nbsp;&nbsp;&nbsp;</span> | #2166AC | rgb(33, 102, 172) |
| 0-2    | EW4 (Exceptionally Wet)   | <span style="background-color: #053061; padding: 0 10px;">&nbsp;&nbsp;&nbsp;</span> | #053061 | rgb(5, 48, 97) |

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
   - [EDDI User Guide](https://psl.noaa.gov/eddi/#resources)

4. Allen, R. G., Pereira, L. S., Raes, D., & Smith, M. (1998). Crop evapotranspiration-Guidelines for computing crop water requirements. FAO Irrigation and drainage paper 56, FAO, Rome, 300(9), D05109. https://www.fao.org/4/x0490e/x0490e00.htm

## Acknowledgments

- NOAA Physical Sciences Laboratory for developing and documenting the EDDI methodology
- The World Bank GOST/DEC Data Group for supporting this implementation
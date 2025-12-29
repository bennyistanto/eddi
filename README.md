# EDDI (Evaporative Demand Drought Index)

This repository contains Python code for calculating the Evaporative Demand Drought Index (EDDI), a physically based drought monitoring and early warning guidance tool that measures drought dynamics based on atmospheric evaporative demand (E₀).

![eddi](/images/idn_cli_agera5_eddi_3month_20151231.png)
*Figure 1. EDDI 3-month as of 31 December 2015 (Very Strong El-Niño)*

## Overview

EDDI examines the signal of actual or potential drought conditions by measuring the anomalous atmospheric evaporative demand using reference evapotranspiration (ET₀). This approach complements traditional land-surface perspective drought indicators by providing an independent assessment of drought conditions from a purely atmospheric perspective.

### Key Features

- **Multiple ET₀ Methods**: Supports FAO-56 Penman-Monteith and ASCE Standardized equations
- **Flexible Time Scales**: Calculate EDDI for periods from 1 to 12 months
- **Robust Calculations**: Includes numerical safeguards for polar regions and edge cases
- **30-Year Climatology**: Reference period based on WMO-recommended 1991-2020 baseline

## Case Study

[Sumba Island](https://en.wikipedia.org/wiki/Sumba), located in East Nusa Tenggara, Indonesia, is a captivating destination known for its unique blend of natural beauty and traditional culture. 

Despite its stunning landscapes of rolling hills and pristine beaches, the island faces significant climate challenges, particularly during the dry season from May to October. The region experiences a distinct monsoon climate with recurring drought events that affect local agriculture and water resources. The island's average annual rainfall varies significantly between the northern and southern regions, with some areas receiving less than 1000mm of rainfall annually. 

Despite these environmental challenges, Sumba has gained international recognition for sustainable luxury tourism, most notably through Nihi Sumba (formerly Nihiwatu), which was voted the best hotel in the world by Travel + Leisure magazine in 2016 and 2017. Located on the island's remote southwestern coast, the resort exemplifies how tourism can thrive while respecting both environmental constraints and local traditions. The island's combination of dramatic coastlines, traditional villages, and unique Marapu culture continues to draw visitors, even as it grapples with climate-related challenges.

---

## Input Data Requirements

> [!IMPORTANT]
> This documentation focuses on the EDDI calculation methodology and implementation. The preliminary steps of downloading AgERA5 data from the Copernicus Climate Data Store (CDS), preprocessing individual variables, and merging them into the required format are not discussed here. These data preparation steps should be completed separately following CDS guidelines and data handling best practices. For information about data acquisition and preprocessing, please refer to the CDS documentation and API guidelines.

### AgERA5 Data

The code uses AgERA5 (Agrometeorological indicators from 1979 to present derived from reanalysis) available from the [Copernicus Climate Data Store](https://cds.climate.copernicus.eu/cdsapp#!/dataset/sis-agrometeorological-indicators).

#### Required Variables

| Variable | Band Name | Unit | Description |
|----------|-----------|------|-------------|
| 2m Air Temperature (Min) | `Temperature_Air_2m_Min_24h` | K | Daily minimum temperature at 2m height |
| 2m Air Temperature (Max) | `Temperature_Air_2m_Max_24h` | K | Daily maximum temperature at 2m height |
| 2m Air Temperature (Mean) | `Temperature_Air_2m_Mean_24h` | K | Daily mean temperature at 2m height |
| 2m Dew Point Temperature | `Dew_Point_Temperature_2m_Mean` | K | Daily mean dew point temperature |
| 10m Wind Speed | `Wind_Speed_10m_Mean` | m s⁻¹ | Daily mean wind speed at 10m height |
| Solar Radiation Flux | `Solar_Radiation_Flux` | J m⁻² day⁻¹ | Daily downward solar radiation |

### Elevation Data

A Digital Elevation Model (DEM) is required for calculating clear-sky radiation and the psychrometric constant:

- **Recommended**: SRTM (Shuttle Radar Topography Mission) at 90m resolution
- **Alternative**: Any DEM resampled to match AgERA5 resolution (~0.1°)

The elevation data should be provided as a NetCDF file with the same spatial extent and coordinate system as the AgERA5 data.

### File Structure

The code expects input files organized as follows:

```
main_folder/
├── input/
│   ├── idn_cli_agera5_sumba_tmin2m_1981_2024.nc
│   ├── idn_cli_agera5_sumba_tmax2m_1981_2024.nc
│   ├── idn_cli_agera5_sumba_tavg2m_1981_2024.nc
│   ├── idn_cli_agera5_sumba_tdew2m_1981_2024.nc
│   ├── idn_cli_agera5_sumba_ws10m_1981_2024.nc
│   ├── idn_cli_agera5_sumba_sr_1981_2024.nc
│   └── idn_cli_srtm_sumba_elev.nc
├── temp/
├── climatology/
├── output/
├── images/
└── bnd/
    └── ne_10m_admin_0_countries.shp
```

---

## Methodology

### 1. Reference Evapotranspiration (ET₀) Methods

The implementation supports three ET₀ calculation methods based on the Penman-Monteith equation:

| Option | Method | Reference Crop | Product ID | Description |
|:------:|--------|----------------|------------|-------------|
| 1 | ASCE Standardized | Tall (Alfalfa) | ETrs | ASCE-EWRI standardized for alfalfa reference |
| 2 | ASCE Standardized | Short (Grass) | ETos | ASCE-EWRI standardized for grass reference |
| 3 | FAO-56 | Short (Grass) | ETos | FAO Irrigation and Drainage Paper 56 |

#### Penman-Monteith Equation

The general form of the Penman-Monteith equation used:

$$ET_0 = \frac{0.408 \cdot \Delta \cdot (R_n - G) + \gamma \cdot \frac{C_n}{T + 273} \cdot u_2 \cdot (e_s - e_a)}{\Delta + \gamma \cdot (1 + C_d \cdot u_2)}$$

Where:
- $ET_0$ = reference evapotranspiration [mm day⁻¹]
- $\Delta$ = slope of saturation vapor pressure curve [kPa °C⁻¹]
- $R_n$ = net radiation [MJ m⁻² day⁻¹]
- $G$ = soil heat flux (assumed 0 for daily calculations) [MJ m⁻² day⁻¹]
- $\gamma$ = psychrometric constant [kPa °C⁻¹]
- $T$ = mean daily air temperature [°C]
- $u_2$ = wind speed at 2m height [m s⁻¹]
- $e_s$ = saturation vapor pressure [kPa]
- $e_a$ = actual vapor pressure [kPa]
- $C_n$, $C_d$ = method-specific constants

#### Method Coefficients

| Parameter | FAO-56 / ASCE Short | ASCE Tall |
|-----------|:-------------------:|:---------:|
| Cₙ (numerator constant) | 900 | 1600 |
| Cₐ (denominator constant) | 0.34 | 0.38 |

### 2. Data Pre-processing

Several adjustments are applied to the AgERA5 input data:

#### 2.1 Temperature Conversion

All temperature variables are converted from Kelvin to Celsius:

```python
T_celsius = T_kelvin - 273.15
```

Applied to: mean, maximum, minimum air temperatures, and dew point temperature.

#### 2.2 Wind Speed Height Adjustment

Wind speed is adjusted from 10m measurement height to the standard 2m reference height using the logarithmic wind profile:

$$u_2 = u_{10} \times \frac{4.87}{\ln(67.8 \times z - 5.42)}$$

Where:
- $u_{10}$ = wind speed at 10m height [m s⁻¹]
- $z$ = measurement height (10 m)
- $u_2$ = wind speed at 2m height [m s⁻¹]

#### 2.3 Solar Radiation Conversion

Solar radiation is converted from J m⁻² day⁻¹ to MJ m⁻² day⁻¹:

$$R_s = R_{s,J} \times 10^{-6}$$

### 3. Radiation Components

#### 3.1 Net Shortwave Radiation

$$R_{ns} = (1 - \alpha) \times R_s$$

Where $\alpha$ = 0.23 (FAO-56 reference crop albedo).

#### 3.2 Extraterrestrial Radiation

$$R_a = \frac{24 \times 60}{\pi} \times G_{sc} \times d_r \times [\omega_s \sin(\varphi) \sin(\delta) + \cos(\varphi) \cos(\delta) \sin(\omega_s)]$$

Where:
- $G_{sc}$ = 0.0820 MJ m⁻² min⁻¹ (solar constant)
- $d_r$ = $1 + 0.033 \cos(2\pi \times DOY / 365)$ (inverse relative distance Earth-Sun)
- $\delta$ = $0.409 \sin(2\pi \times DOY / 365 - 1.39)$ (solar declination)
- $\omega_s$ = $\arccos[-\tan(\varphi) \tan(\delta)]$ (sunset hour angle)
- $\varphi$ = latitude [radians]
- $DOY$ = day of year

> **Note**: The arccos argument is clamped to [-1, 1] to ensure numerical stability at polar regions.

#### 3.3 Clear-Sky Radiation

$$R_{so} = (0.75 + 2 \times 10^{-5} \times z) \times R_a$$

Where $z$ = elevation [m].

#### 3.4 Net Longwave Radiation

$$R_{nl} = \sigma \times \frac{T_{max,K}^4 + T_{min,K}^4}{2} \times (0.34 - 0.14\sqrt{e_a}) \times (1.35 \times \frac{R_s}{R_{so}} - 0.35)$$

Where:
- $\sigma$ = 4.903 × 10⁻⁹ MJ K⁻⁴ m⁻² day⁻¹ (Stefan-Boltzmann constant)
- $T_{max,K}$, $T_{min,K}$ = maximum and minimum absolute temperatures [K]
- $e_a$ = actual vapor pressure [kPa]

#### 3.5 Net Radiation

$$R_n = R_{ns} - R_{nl}$$

### 4. Vapor Pressure Components

#### 4.1 Saturation Vapor Pressure

$$e°(T) = 0.6108 \times \exp\left[\frac{17.27 \times T}{T + 237.3}\right]$$

$$e_s = \frac{e°(T_{max}) + e°(T_{min})}{2}$$

#### 4.2 Actual Vapor Pressure

Derived from dew point temperature:

$$e_a = 0.6108 \times \exp\left[\frac{17.27 \times T_{dew}}{T_{dew} + 237.3}\right]$$

#### 4.3 Vapor Pressure Deficit

$$VPD = e_s - e_a$$

### 5. Other Parameters

#### 5.1 Slope of Saturation Vapor Pressure Curve

$$\Delta = \frac{4098 \times e_s}{(T_{mean} + 237.3)^2}$$

#### 5.2 Psychrometric Constant

$$P = 101.3 \times \left[\frac{293 - 0.0065 \times z}{293}\right]^{5.26}$$

$$\gamma = 0.000665 \times P$$

Where $P$ = atmospheric pressure [kPa] and $z$ = elevation [m].

### 6. Numerical Safeguards

The implementation includes safeguards to prevent computational errors:

| Issue | Solution | Function |
|-------|----------|----------|
| Division by zero | Minimum value clamping (10⁻⁶) | `clamp_min_value()` |
| Invalid arccos domain | Input clamping to [-1, 1] | `safe_arccos()` |
| Clear-sky radiation = 0 | Minimum value clamping (10⁻⁴) | Applied in calculation |
| Negative ET₀ | Result clamped to ≥ 0 | Applied in calculation |

---

## EDDI Calculation Process

The EDDI calculation transforms reference evapotranspiration into a drought index through temporal aggregation and percentile ranking against climatology.

### Step 1: Daily ET₀ Calculation

For each day in the dataset:
1. Load AgERA5 meteorological variables
2. Apply unit conversions (K→°C, J→MJ)
3. Adjust wind speed from 10m to 2m
4. Calculate radiation components using DEM-derived elevation
5. Compute ET₀ using the selected Penman-Monteith method
6. Apply mask for missing input data

### Step 2: Temporal Aggregation

Sum ET₀ over the specified time scale using dekad-based windows:

| Dekad | Days | Description |
|:-----:|:----:|-------------|
| 1 | 1-10 | First 10 days of month |
| 2 | 11-20 | Middle 10 days of month |
| 3 | 21-end | Remaining days of month |

For multi-month time scales, the window extends backward from the end dekad:

$$ET_{0,sum} = \sum_{d=start}^{end} ET_0(d)$$

### Step 3: Climatology Distribution

Build the reference distribution by computing equivalent ET₀ sums for the same calendar window across all years in the reference period (1991-2020):

For each reference year $Y$ (1991 to 2020):
1. Define window with same start month/day and end month/day as analysis period
2. Handle year-crossing windows (e.g., Nov-Jan spans Y-1 to Y)
3. Sum daily ET₀ for the window
4. Add to climatology distribution

$$\text{Climatology} = \{ET_{0,sum}(Y) : Y \in [1991, 2020]\}$$

### Step 4: Percentile Calculation

EDDI percentile represents the rank of current ET₀ sum within the climatology:

$$EDDI_{percentile} = \frac{\text{count of climatology values} \leq \text{current } ET_{0,sum}}{n} \times 100$$

Where $n$ = number of valid climatology years (typically 30).

**Interpretation**: Higher percentile = Higher evaporative demand = Drier conditions

### Step 5: EDDI Classification

The EDDI percentile is classified into 11 categories following NOAA conventions:

| Percentile Range | Category Code | Category Name | HEX Color |
|-----------------:|---------------|---------------|-----------|
| 98 - 100 | ED4 | Exceptional Drought | `#730000` |
| 95 - 98 | ED3 | Extreme Drought | `#E60000` |
| 90 - 95 | ED2 | Severe Drought | `#FFAA00` |
| 80 - 90 | ED1 | Moderate Drought | `#FCD37F` |
| 70 - 80 | ED0 | Abnormally Dry | `#FFFF00` |
| 30 - 70 | — | Near Normal | `#FFFFFF` |
| 20 - 30 | EW0 | Abnormally Wet | `#8CCDEF` |
| 10 - 20 | EW1 | Moderately Wet | `#00BFFF` |
| 5 - 10 | EW2 | Severely Wet | `#1D90FF` |
| 2 - 5 | EW3 | Extremely Wet | `#4169E1` |
| 0 - 2 | EW4 | Exceptionally Wet | `#0000FF` |

### Workflow Diagram

```
┌─────────────────────────────────────────────────────────────────────┐
│                     EDDI CALCULATION WORKFLOW                        │
└─────────────────────────────────────────────────────────────────────┘

┌──────────────────┐    ┌──────────────────┐    ┌──────────────────┐
│     AgERA5       │    │       DEM        │    │   User Inputs    │
│   Daily Data     │    │   (Elevation)    │    │  (Date, Scale,   │
│   (6 variables)  │    │                  │    │   ET₀ Method)    │
└────────┬─────────┘    └────────┬─────────┘    └────────┬─────────┘
         │                       │                       │
         └───────────────────────┼───────────────────────┘
                                 │
                                 ▼
              ┌──────────────────────────────────────┐
              │     Daily ET₀ Calculation            │
              │  (FAO-56 / ASCE ETos / ASCE ETrs)    │
              └────────────────┬─────────────────────┘
                               │
              ┌────────────────┴────────────────┐
              │                                 │
              ▼                                 ▼
┌──────────────────────┐          ┌──────────────────────┐
│  Current Period      │          │  Reference Period    │
│  ET₀ Sum             │          │  (1991-2020)         │
│  (Analysis Window)   │          │  30 Annual ET₀ Sums  │
└──────────┬───────────┘          └──────────┬───────────┘
           │                                 │
           └─────────────┬───────────────────┘
                         │
                         ▼
              ┌──────────────────────────────────┐
              │     Percentile Calculation       │
              │  Rank current vs. climatology    │
              └────────────────┬─────────────────┘
                               │
                               ▼
              ┌──────────────────────────────────┐
              │     EDDI Classification          │
              │  (ED4 to EW4, 11 categories)     │
              └────────────────┬─────────────────┘
                               │
                               ▼
              ┌──────────────────────────────────┐
              │     Output: NetCDF & PNG         │
              └──────────────────────────────────┘
```

---

## Usage

### Basic Example

```python
from compute_eddi import process_climate_data, get_dekad_dates

# Configuration
main_folder = '/path/to/data'
year = 2015
month = 12
dekad = 3  # End of December
time_scales = [1, 3, 6, 12]  # 1, 3, 6, and 12 month EDDI

# Select ET₀ method:
# 1 = ASCE Standardized tall reference crop (alfalfa, ETrs)
# 2 = ASCE Standardized short reference crop (grass, ETos)
# 3 = FAO-56 Penman-Monteith (grass, ETos)
et_option = 1

# Process EDDI
ds_result = process_climate_data(
    main_folder, 
    year, 
    month, 
    dekad, 
    time_scales, 
    et_option
)
```

### Configuration Options

Modify the configuration at the top of `compute_eddi.py`:

```python
# Select ET₀ method (1, 2, or 3)
ET_OPTION = 1

# Or use get_et_config() for programmatic access
from compute_eddi import get_et_config
config = get_et_config(1)  # Returns ETConfig object
print(config.description)  # "ASCE Standardized tall reference crop (alfalfa)"
```

---

## Output

### NetCDF Files

The code produces NetCDF files containing EDDI values for each time scale:

- **File naming**: `idn_cli_agera5_eddi_{scale}month_{YYYYMMDD}.nc`
- **Variables**: `EDDI_{scale}month` (percentile values 0-100)
- **Attributes**: CF-1.8 compliant metadata including method information

### Visualization

The `plot_eddi_and_inputs()` function creates a 3×3 grid showing:

| Row 1 | Row 2 | Row 3 |
|-------|-------|-------|
| Tavg | Tdew | ET₀ Sum |
| Tmin | Wind Speed | Clim. ET₀ |
| Tmax | Solar Radiation | EDDI |

Output saved as PNG at 300 DPI.

---

## Dependencies

```
numpy
pandas
xarray
geopandas
matplotlib
netcdf4
```

Install with:
```bash
pip install numpy pandas xarray geopandas matplotlib netcdf4
```

---

## References

1. Hobbins, M. T., A. Wood, D. J. McEvoy, J. L. Huntington, C. Morton, M. Anderson, and C. Hain (2016): The Evaporative Demand Drought Index. Part I: Linking drought evolution to variations in evaporative demand. *J. Hydrometeor.*, 17(6), 1745-1761, [doi:10.1175/JHM-D-15-0121.1](https://doi.org/10.1175/JHM-D-15-0121.1)

2. McEvoy, D. J., J. L. Huntington, M. T. Hobbins, A. Wood, C. Morton, M. Anderson, and C. Hain (2016): The Evaporative Demand Drought Index. Part II: CONUS-wide assessment against common drought indicators. *J. Hydrometeor.*, 17(6), 1763-1779, [doi:10.1175/JHM-D-15-0122.1](https://doi.org/10.1175/JHM-D-15-0122.1)

3. Allen, R. G., Pereira, L. S., Raes, D., & Smith, M. (1998). Crop evapotranspiration - Guidelines for computing crop water requirements. *FAO Irrigation and Drainage Paper 56*. Rome: FAO. [https://www.fao.org/4/x0490e/x0490e00.htm](https://www.fao.org/4/x0490e/x0490e00.htm)

4. ASCE-EWRI (2005). The ASCE Standardized Reference Evapotranspiration Equation. *ASCE-EWRI Task Committee Report*. Reston, VA: ASCE.

5. NOAA Physical Sciences Laboratory EDDI Documentation:
   - [EDDI Homepage](https://psl.noaa.gov/eddi/)
   - [EDDI User Guide](https://psl.noaa.gov/eddi/#resources)

---

## Acknowledgments

- NOAA Physical Sciences Laboratory for developing and documenting the EDDI methodology
- The World Bank GOST/DEC Data Group for supporting this implementation
- ECMWF and Copernicus Climate Data Store for providing AgERA5 reanalysis data

---

## Contact

Benny Istanto  
Climate Geographer - GOST/DEC Data Group, The World Bank  
Email: bistanto@worldbank.org

---

## License

Public domain.

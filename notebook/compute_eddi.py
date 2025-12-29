"""
EDDI (Evaporative Demand Drought Index)

This script calculates the Evaporative Demand Drought Index (EDDI) using AgERA5 
meteorological data. EDDI is a drought monitoring tool that measures the anomalous 
atmospheric evaporative demand (E₀) using reference evapotranspiration.

Supported ET₀ Methods:
---------------------
1. ASCE Standardized tall reference crop (alfalfa, ETrs)
2. ASCE Standardized short reference crop (grass, ETos)
3. FAO-56 Penman-Monteith (grass, ETos)

Required Input Data:
-------------------
AgERA5 variables (NetCDF format):
    - Temperature_Air_2m_Mean_24h [K]
    - Temperature_Air_2m_Max_24h [K]
    - Temperature_Air_2m_Min_24h [K]
    - Dew_Point_Temperature_2m_Mean [K]
    - Wind_Speed_10m_Mean [m s⁻¹]
    - Solar_Radiation_Flux [J m⁻² day⁻¹]
    
Additional Data:
    - Digital elevation model (e.g., SRTM) [m]

Main Functions:
--------------
- calculate_reference_et(): Calculates daily reference ET using selected method
- compute_eddi_for_timescale(): Computes EDDI for a specific time window
- process_climate_data(): Main processing pipeline
- plot_eddi_and_inputs(): Creates visualizations

Output:
------
- NetCDF files containing EDDI values for each time scale
- Visualization plots showing input variables and EDDI results

References:
----------
1. Hobbins et al. (2016): doi:10.1175/JHM-D-15-0121.1
2. McEvoy et al. (2016): doi:10.1175/JHM-D-15-0122.1
3. Allen et al. (1998): https://www.fao.org/4/x0490e/x0490e00.htm
4. ASCE-EWRI (2005): The ASCE Standardized Reference Evapotranspiration Equation

Notes:
-----
- Uses 1991-2020 as the climatological reference period (WMO standard)
- Implements dekad-based temporal aggregation
- Supports multiple time scales (1, 2, 3, 6, 9, 12 months)
- All computations handle NaN values appropriately
- Includes numerical safeguards for polar regions and edge cases

Author: Benny Istanto, GOST/DEC Data Group/The World Bank
Version: 2.0.0
Last Updated: 2025
"""

# =======================
# Import Libraries
# =======================
import os
import sys
import calendar
import numpy as np
import pandas as pd
import geopandas as gpd
import xarray as xr
import datetime
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from dataclasses import dataclass
from typing import Tuple, List, Optional, Union


# =======================
# Configuration
# =======================
@dataclass
class ETConfig:
    """Configuration for reference ET calculation method."""
    method: str          # 'FAO56' or 'ASCE'
    ref_type: str        # 'short' or 'tall'
    product_id: str      # 'ETos' or 'ETrs'
    description: str     # Human-readable description
    cn: float            # Numerator constant
    cd: float            # Denominator constant


def get_et_config(option: int) -> ETConfig:
    """
    Get ET configuration based on selected option.
    
    Parameters
    ----------
    option : int
        1: ASCE Standardized tall reference crop (alfalfa, ETrs)
        2: ASCE Standardized short reference crop (grass, ETos)
        3: FAO-56 Penman-Monteith (grass, ETos)
    
    Returns
    -------
    ETConfig
        Configuration object with method parameters.
    """
    configs = {
        1: ETConfig(
            method='ASCE',
            ref_type='tall',
            product_id='ETrs',
            description='ASCE Standardized tall reference crop (alfalfa)',
            cn=1600.0,
            cd=0.38
        ),
        2: ETConfig(
            method='ASCE',
            ref_type='short',
            product_id='ETos',
            description='ASCE Standardized short reference crop (grass)',
            cn=900.0,
            cd=0.34
        ),
        3: ETConfig(
            method='FAO56',
            ref_type='short',
            product_id='ETos',
            description='FAO-56 Penman-Monteith (grass)',
            cn=900.0,
            cd=0.34
        )
    }
    
    if option not in configs:
        print(f"[WARNING] Invalid option {option}. Defaulting to ASCE tall reference crop.")
        return configs[1]
    
    return configs[option]


# Default configuration
ET_OPTION = 1  # Change this to select different methods (1, 2, or 3)
ET_CONFIG = get_et_config(ET_OPTION)


# =======================
# Utility Functions
# =======================
def clamp_min_value(data: Union[xr.DataArray, np.ndarray], 
                    min_value: float = 1e-6) -> Union[xr.DataArray, np.ndarray]:
    """
    Clamp values to ensure they are greater than or equal to a minimum value.
    This helps prevent division by zero and other numerical issues.
    
    Parameters
    ----------
    data : xarray.DataArray or numpy.ndarray
        Input data to clamp.
    min_value : float
        Minimum value to enforce (default: 1e-6).
    
    Returns
    -------
    Same type as input with values clamped.
    """
    if isinstance(data, xr.DataArray):
        return data.where(data >= min_value, min_value)
    return np.maximum(data, min_value)


def safe_arccos(x: Union[xr.DataArray, np.ndarray]) -> Union[xr.DataArray, np.ndarray]:
    """
    Safely calculate arccos by first clamping input values to [-1, 1].
    This prevents domain errors at polar regions where numerical issues 
    might push values slightly outside the valid range.
    
    Parameters
    ----------
    x : xarray.DataArray or numpy.ndarray
        Input values for arccos calculation.
    
    Returns
    -------
    arccos of the clamped input.
    """
    if isinstance(x, xr.DataArray):
        x_clipped = x.clip(-1, 1)
    else:
        x_clipped = np.clip(x, -1, 1)
    return np.arccos(x_clipped)


def kelvin_to_celsius(kelvin: Union[xr.DataArray, np.ndarray]) -> Union[xr.DataArray, np.ndarray]:
    """
    Convert temperature from Kelvin to Celsius.
    
    Parameters
    ----------
    kelvin : xarray.DataArray or numpy.ndarray
        Temperature in Kelvin.
    
    Returns
    -------
    Temperature in Celsius.
    """
    return kelvin - 273.15


def adjust_wind_speed_10m_to_2m(ws10m: Union[xr.DataArray, np.ndarray]) -> Union[xr.DataArray, np.ndarray]:
    """
    Adjust 10m wind speed to 2m using the FAO-56 recommended logarithmic profile.
    
    u2 = u10 * 4.87 / ln(67.8 * z - 5.42)
    
    where z = 10 m (measurement height).
    
    Parameters
    ----------
    ws10m : xarray.DataArray or numpy.ndarray
        Wind speed at 10 m height [m s⁻¹].
    
    Returns
    -------
    Wind speed at 2 m height [m s⁻¹].
    """
    z = 10.0  # measurement height in meters
    inside = 67.8 * z - 5.42
    inside_safe = clamp_min_value(inside, 1e-6)
    return ws10m * 4.87 / np.log(inside_safe)


def calc_saturation_vapor_pressure(temp_c: Union[xr.DataArray, np.ndarray]) -> Union[xr.DataArray, np.ndarray]:
    """
    Calculate saturation vapor pressure from temperature in Celsius.
    
    e°(T) = 0.6108 * exp[17.27 * T / (T + 237.3)]
    
    Parameters
    ----------
    temp_c : xarray.DataArray or numpy.ndarray
        Temperature in Celsius.
    
    Returns
    -------
    Saturation vapor pressure in kPa.
    """
    return 0.6108 * np.exp((17.27 * temp_c) / (temp_c + 237.3))


def calc_actual_vapor_pressure(tdew_c: Union[xr.DataArray, np.ndarray]) -> Union[xr.DataArray, np.ndarray]:
    """
    Calculate actual vapor pressure from dew point temperature in Celsius.
    
    Parameters
    ----------
    tdew_c : xarray.DataArray or numpy.ndarray
        Dew point temperature in Celsius.
    
    Returns
    -------
    Actual vapor pressure in kPa.
    """
    return calc_saturation_vapor_pressure(tdew_c)


def calc_psychrometric_constant(elevation: Union[xr.DataArray, np.ndarray]) -> Union[xr.DataArray, np.ndarray]:
    """
    Calculate psychrometric constant (gamma) using elevation.
    
    P = 101.3 * [(293 - 0.0065 * z) / 293]^5.26
    γ = 0.000665 * P
    
    Parameters
    ----------
    elevation : xarray.DataArray or numpy.ndarray
        Elevation in meters.
    
    Returns
    -------
    Psychrometric constant in kPa °C⁻¹.
    """
    # Atmospheric pressure based on elevation [kPa]
    p = 101.3 * ((293.0 - 0.0065 * elevation) / 293.0) ** 5.26
    # Psychrometric constant
    return 0.000665 * p


def calc_extraterrestrial_radiation(lat_rad: Union[xr.DataArray, np.ndarray],
                                    doy: Union[xr.DataArray, int]) -> Union[xr.DataArray, np.ndarray]:
    """
    Calculate extraterrestrial radiation (Ra) in MJ m⁻² day⁻¹.
    
    Ra = (24 * 60 / π) * Gsc * dr * [ωs * sin(φ) * sin(δ) + cos(φ) * cos(δ) * sin(ωs)]
    
    Parameters
    ----------
    lat_rad : xarray.DataArray or numpy.ndarray
        Latitude in radians.
    doy : xarray.DataArray or int
        Day of year (1-365/366).
    
    Returns
    -------
    Extraterrestrial radiation in MJ m⁻² day⁻¹.
    """
    Gsc = 0.0820  # Solar constant [MJ m⁻² min⁻¹]
    
    # Inverse relative distance Earth-Sun
    dr = 1.0 + 0.033 * np.cos(2.0 * np.pi * doy / 365.0)
    
    # Solar declination [radians]
    solar_dec = 0.409 * np.sin(2.0 * np.pi * doy / 365.0 - 1.39)
    
    # Sunset hour angle [radians]
    # Clamp the argument to [-1, 1] for numerical stability at polar regions
    x = -np.tan(lat_rad) * np.tan(solar_dec)
    omega_s = safe_arccos(x)
    
    # Extraterrestrial radiation
    Ra = (24.0 * 60.0 / np.pi) * Gsc * dr * (
        omega_s * np.sin(lat_rad) * np.sin(solar_dec) +
        np.cos(lat_rad) * np.cos(solar_dec) * np.sin(omega_s)
    )
    
    return Ra


# =======================
# Dekad Helper Functions
# =======================
def get_dekad_dates(year: int, month: int, dekad: int) -> Tuple[datetime.date, datetime.date]:
    """
    Get the start and end dates for a dekad.
    
    Dekad definitions:
        - Dekad 1: Day 1-10
        - Dekad 2: Day 11-20
        - Dekad 3: Day 21 to end of month
    
    Parameters
    ----------
    year : int
        Year.
    month : int
        Month (1-12).
    dekad : int
        Dekad number (1, 2, or 3).
    
    Returns
    -------
    Tuple[datetime.date, datetime.date]
        Start and end dates for the dekad.
    """
    if dekad == 1:
        start = datetime.date(year, month, 1)
        end = datetime.date(year, month, 10)
    elif dekad == 2:
        start = datetime.date(year, month, 11)
        end = datetime.date(year, month, 20)
    elif dekad == 3:
        # Last day of month
        last_day = calendar.monthrange(year, month)[1]
        start = datetime.date(year, month, 21)
        end = datetime.date(year, month, last_day)
    else:
        raise ValueError("dekad must be 1, 2, or 3.")
    return start, end


def shift_date_by_months(date_obj: datetime.date, n_months: int) -> datetime.date:
    """
    Shift a date backward by n_months.
    
    Parameters
    ----------
    date_obj : datetime.date
        Original date.
    n_months : int
        Number of months to shift backward.
    
    Returns
    -------
    datetime.date
        Shifted date.
    """
    month = date_obj.month - n_months
    year = date_obj.year
    day = date_obj.day
    
    while month <= 0:
        month += 12
        year -= 1
    
    # Clamp day if out of range for the resulting month
    last_day_new_month = calendar.monthrange(year, month)[1]
    if day > last_day_new_month:
        day = last_day_new_month
    
    return datetime.date(year, month, day)


def get_dekad_based_month_window(year: int, month: int, dekad: int, 
                                  scale_in_months: int) -> Tuple[datetime.date, datetime.date]:
    """
    Return a (start_date, end_date) that covers 'scale_in_months' of
    discrete, dekad-based months, ending on the actual final date for
    (year, month, dekad).
    
    Parameters
    ----------
    year : int
        Year of the ending dekad.
    month : int
        Month of the ending dekad.
    dekad : int
        Dekad number (1, 2, or 3).
    scale_in_months : int
        Number of months to aggregate.
    
    Returns
    -------
    Tuple[datetime.date, datetime.date]
        Start and end dates for the aggregation window.
    """
    # Final date from (year, month, dekad)
    _, final_end = get_dekad_dates(year, month, dekad)
    end_date = final_end
    
    def step_back_one_month(yr, mo):
        mo_new = mo - 1
        yr_new = yr
        if mo_new < 1:
            mo_new += 12
            yr_new -= 1
        return yr_new, mo_new
    
    # Start from the month/year of the final date
    start_yr, start_mo = year, month
    
    for _ in range(scale_in_months - 1):
        start_yr, start_mo = step_back_one_month(start_yr, start_mo)
    
    # Determine start day based on dekad and scale
    if scale_in_months == 1:
        if dekad == 1:
            start_yr, start_mo = step_back_one_month(start_yr, start_mo)
            start_day = 11
        elif dekad == 2:
            start_yr, start_mo = step_back_one_month(start_yr, start_mo)
            start_day = 21
        else:  # dekad == 3
            start_day = 1
    else:
        # For multi-month scales, start from day 1
        start_day = 1
    
    # Clamp start_day if needed
    last_day_of_startmo = calendar.monthrange(start_yr, start_mo)[1]
    if start_day > last_day_of_startmo:
        start_day = last_day_of_startmo
    
    start_date = datetime.date(start_yr, start_mo, start_day)
    
    return start_date, end_date


# =======================
# Reference ET Calculation
# =======================
def calculate_reference_et(
        tmin: xr.DataArray,
        tmax: xr.DataArray,
        tmean: xr.DataArray,
        tdew: xr.DataArray,
        ws10m: xr.DataArray,
        srad: xr.DataArray,
        elevation: xr.DataArray,
        lat: xr.DataArray,
        doy: xr.DataArray,
        config: Optional[ETConfig] = None
    ) -> xr.DataArray:
    """
    Calculate reference evapotranspiration (ET₀) using the Penman-Monteith equation.
    
    Supports three methods:
        - FAO-56 Penman-Monteith (grass reference)
        - ASCE Standardized short reference crop (grass, ETos)
        - ASCE Standardized tall reference crop (alfalfa, ETrs)
    
    Parameters
    ----------
    tmin : xarray.DataArray
        Daily minimum temperature [K].
    tmax : xarray.DataArray
        Daily maximum temperature [K].
    tmean : xarray.DataArray
        Daily mean temperature [K].
    tdew : xarray.DataArray
        Dew point temperature [K].
    ws10m : xarray.DataArray
        Wind speed at 10m [m s⁻¹].
    srad : xarray.DataArray
        Downward solar radiation [J m⁻² day⁻¹].
    elevation : xarray.DataArray
        Elevation [m].
    lat : xarray.DataArray
        Latitude [degrees].
    doy : xarray.DataArray
        Day of year (1-365/366).
    config : ETConfig, optional
        Configuration object specifying the ET method. 
        Defaults to global ET_CONFIG.
    
    Returns
    -------
    xarray.DataArray
        Reference evapotranspiration [mm day⁻¹].
    
    Notes
    -----
    The Penman-Monteith equation:
    
    ET₀ = [0.408 * Δ * (Rn - G) + γ * (Cn / (T + 273)) * u2 * (es - ea)] /
          [Δ + γ * (1 + Cd * u2)]
    
    Where coefficients (Cn, Cd) vary by method:
        - FAO-56 / ASCE Short: Cn=900, Cd=0.34
        - ASCE Tall: Cn=1600, Cd=0.38
    """
    if config is None:
        config = ET_CONFIG
    
    print(f"[INFO] Calculating ET₀ using {config.description}")
    
    # Convert temperatures from Kelvin to Celsius
    tmin_c = kelvin_to_celsius(tmin)
    tmax_c = kelvin_to_celsius(tmax)
    tmean_c = kelvin_to_celsius(tmean)
    tdew_c = kelvin_to_celsius(tdew)
    
    # Adjust wind speed from 10m to 2m
    ws2m = adjust_wind_speed_10m_to_2m(ws10m)
    
    # Convert solar radiation from J/m²/day to MJ/m²/day
    Rs = srad * 1e-6
    
    # Use constant albedo of 0.23 as per FAO-56 reference crop
    ALBEDO = 0.23
    
    # Net shortwave radiation: Rns = (1 - albedo) * Rs
    Rns = (1.0 - ALBEDO) * Rs
    
    # Convert latitude from degrees to radians
    lat_rad = np.deg2rad(lat)
    
    # Calculate extraterrestrial radiation (Ra)
    Ra = calc_extraterrestrial_radiation(lat_rad, doy)
    
    # Clear-sky radiation: Rso = (0.75 + 2e-5 * elevation) * Ra
    Rso_raw = (0.75 + 2e-5 * elevation) * Ra
    Rso = clamp_min_value(Rso_raw, 1e-4)
    
    # Net longwave radiation (Rnl)
    # Average of tmin^4 and tmax^4 (using original Kelvin values)
    avg_T4 = 0.5 * (tmin**4 + tmax**4)
    sigma = 4.903e-9  # Stefan-Boltzmann constant [MJ m⁻² day⁻¹ K⁻⁴]
    
    # Actual vapor pressure from dew point
    ea = calc_actual_vapor_pressure(tdew_c)
    
    # Handle division by zero for Rs/Rso ratio
    ratio_Rs_Rso = xr.where(Rso > 0, Rs / Rso, np.nan)
    
    Rnl = sigma * avg_T4 * (0.34 - 0.14 * np.sqrt(ea)) * (1.35 * ratio_Rs_Rso - 0.35)
    
    # Net radiation: Rn = Rns - Rnl
    Rn = Rns - Rnl
    
    # Saturation vapor pressure (average of es at tmin and tmax)
    es_tmin = calc_saturation_vapor_pressure(tmin_c)
    es_tmax = calc_saturation_vapor_pressure(tmax_c)
    es = 0.5 * (es_tmax + es_tmin)
    
    # Vapor pressure deficit
    vpd = es - ea
    
    # Slope of saturation vapor pressure curve (delta)
    delta = (4098.0 * es) / ((tmean_c + 237.3) ** 2)
    
    # Psychrometric constant (gamma)
    gamma = calc_psychrometric_constant(elevation)
    
    # Soil heat flux (G) - assumed zero for daily calculations
    G = 0.0
    
    # Get coefficients from configuration
    Cn = config.cn
    Cd = config.cd
    
    # Penman-Monteith equation
    numerator = (0.408 * delta * (Rn - G)) + (gamma * (Cn / (tmean_c + 273)) * ws2m * vpd)
    denominator_raw = delta + gamma * (1.0 + Cd * ws2m)
    denominator = clamp_min_value(denominator_raw, 1e-6)
    
    et0 = numerator / denominator
    
    # Ensure non-negative values
    et0 = et0.where(et0 > 0, 0.0)
    
    return et0


# =======================
# EDDI Calculation
# =======================
def compute_eddi_for_timescale(
        ds_full_et0: xr.DataArray,
        climatology_et0: xr.DataArray,
        scale_in_months: int,
        end_date: datetime.date
    ) -> xr.DataArray:
    """
    Compute EDDI percentile for a specific time scale.
    
    EDDI is calculated by:
    1. Summing ET₀ over the analysis window
    2. Comparing to the distribution of sums from 1991-2020 climatology
    3. Computing the percentile rank
    
    Higher percentile = Higher evaporative demand = Drier conditions
    
    Parameters
    ----------
    ds_full_et0 : xarray.DataArray
        Daily ET₀ for the full dataset (e.g., 1981-present).
    climatology_et0 : xarray.DataArray
        Daily ET₀ for the climatology period (1991-2020).
    scale_in_months : int
        Number of months for aggregation.
    end_date : datetime.date
        End date of the analysis window.
    
    Returns
    -------
    xarray.DataArray
        EDDI percentile (0-100).
    """
    print(f"[INFO] Computing EDDI for {scale_in_months}-month window ending on {end_date}")
    
    # Determine dekad from end date
    day = end_date.day
    if day <= 10:
        dekad = 1
    elif day <= 20:
        dekad = 2
    else:
        dekad = 3
    
    year = end_date.year
    month = end_date.month
    
    # Get the analysis window dates
    start_date_window, actual_end_date = get_dekad_based_month_window(
        year, month, dekad, scale_in_months
    )
    
    print(f"[INFO] --> Window: {start_date_window} to {actual_end_date}")
    
    # Sum ET₀ for the current window
    current_sel = ds_full_et0.sel(time=slice(str(start_date_window), str(actual_end_date)))
    current_sum = current_sel.sum(dim='time', skipna=False)
    
    # Extract month/day for climatology window construction
    start_mo = start_date_window.month
    start_day = start_date_window.day
    end_mo = actual_end_date.month
    end_day = actual_end_date.day
    
    # Build climatology distribution from 1991-2020
    all_clim_sums = []
    years_in_clim = range(1991, 2021)
    valid_years = []
    
    print("[INFO] --> Building climatology distribution for 1991-2020...")
    
    for y_clim in years_in_clim:
        # Handle year-crossing windows (e.g., Nov-Jan)
        if (start_mo > end_mo) or (start_mo == end_mo and start_day > end_day):
            baseline_start_yr = y_clim - 1
        else:
            baseline_start_yr = y_clim
        
        baseline_end_yr = y_clim
        
        try:
            # Handle leap year edge cases
            try:
                baseline_start = datetime.date(baseline_start_yr, start_mo, start_day)
            except ValueError:
                # Feb 29 in non-leap year -> use Feb 28
                baseline_start = datetime.date(baseline_start_yr, start_mo, 28)
            
            try:
                baseline_end = datetime.date(baseline_end_yr, end_mo, end_day)
            except ValueError:
                # Feb 29 in non-leap year -> use Feb 28
                baseline_end = datetime.date(baseline_end_yr, end_mo, 28)
            
            sel_da = climatology_et0.sel(time=slice(str(baseline_start), str(baseline_end)))
            
            # Check if we have data for this window
            if sel_da.time.size > 0:
                sum_da = sel_da.sum(dim='time', skipna=False)
                all_clim_sums.append(sum_da)
                valid_years.append(y_clim)
        except Exception as e:
            print(f"[WARNING] Skipping year {y_clim}: {e}")
            continue
    
    if len(all_clim_sums) == 0:
        raise ValueError("[ERROR] No valid climatology sums found. Check data coverage.")
    
    print(f"[INFO] --> Using {len(valid_years)} years for climatology")
    
    # Stack climatology sums
    clim_sums_da = xr.concat(all_clim_sums, dim=pd.Index(valid_years, name="year"))
    
    # Calculate percentile rank
    # EDDI percentile = (count of climatology values <= current sum) / n * 100
    n_years = clim_sums_da["year"].size
    count_le = (clim_sums_da <= current_sum).sum(dim="year", skipna=False)
    eddi_percentile = (count_le / n_years) * 100.0
    
    # Mask where current_sum is NaN
    eddi_percentile = eddi_percentile.where(current_sum.notnull())
    eddi_percentile.name = f"EDDI_{scale_in_months}month"
    
    return eddi_percentile


# =======================
# Main Processing Pipeline
# =======================
def process_climate_data(
        main_folder: str,
        year: int,
        month: int,
        dekad: int,
        time_scales: List[int],
        et_option: int = 1
    ) -> xr.Dataset:
    """
    Process climate data from NetCDF files and produce EDDI output.
    
    Parameters
    ----------
    main_folder : str
        Path to the main project folder.
    year : int
        Year of the analysis end date.
    month : int
        Month of the analysis end date.
    dekad : int
        Dekad of the analysis end date (1, 2, or 3).
    time_scales : list
        List of time scales in months (e.g., [1, 3, 6, 12]).
    et_option : int
        ET₀ method option (1=ASCE tall, 2=ASCE short, 3=FAO-56).
    
    Returns
    -------
    xarray.Dataset
        Dataset containing EDDI variables for each time scale.
    """
    print("[INFO] Starting EDDI processing...")
    print(f"[INFO] ET₀ Method: {get_et_config(et_option).description}")
    
    # Get ET configuration
    config = get_et_config(et_option)
    
    user_start_date, user_end_date = get_dekad_dates(year, month, dekad)
    
    # Output file naming
    dekad_str = {1: "01", 2: "11", 3: "21"}[dekad]
    end_date_str = f"{year}{month:02d}{dekad_str}"
    
    # File paths
    merged_path = f"{main_folder}/temp/idn_cli_agera5_sumba_full_1981_2024.nc"
    clim_path = f"{main_folder}/climatology/idn_cli_agera5_sumba_clim_1991_2020.nc"
    
    if os.path.exists(merged_path) and os.path.exists(clim_path):
        print("[INFO] Loading existing merged data and climatology...")
        ds = xr.open_dataset(merged_path)
        ds_clim = xr.open_dataset(clim_path)
    else:
        # Load input files
        files = {
            'tmin2m': f'{main_folder}/input/idn_cli_agera5_sumba_tmin2m_1981_2024.nc',
            'tmax2m': f'{main_folder}/input/idn_cli_agera5_sumba_tmax2m_1981_2024.nc',
            'tavg2m': f'{main_folder}/input/idn_cli_agera5_sumba_tavg2m_1981_2024.nc',
            'td2m': f'{main_folder}/input/idn_cli_agera5_sumba_tdew2m_1981_2024.nc',
            'wind_speed': f'{main_folder}/input/idn_cli_agera5_sumba_ws10m_1981_2024.nc',
            'solar_radiation': f'{main_folder}/input/idn_cli_agera5_sumba_sr_1981_2024.nc',
            'elevation': f'{main_folder}/input/idn_cli_srtm_sumba_elev.nc'
        }
        
        loaded_datasets = []
        for key, filepath in files.items():
            print(f"[INFO] Loading {key} from {filepath}")
            try:
                ds_i = xr.open_dataset(filepath, decode_cf=True)
                loaded_datasets.append(ds_i)
            except FileNotFoundError:
                print(f"[ERROR] File not found: {filepath}")
                sys.exit(1)
            except Exception as e:
                print(f"[ERROR] Could not open {filepath}: {e}")
                sys.exit(1)
        
        print("[INFO] Merging input datasets...")
        ds = xr.merge(loaded_datasets)
        ds = ds.sel(time=slice("1981-01-01", "2024-12-31"))
        
        # Create climatology subset
        ds_clim = ds.sel(time=slice("1991-01-01", "2020-12-31"))
        
        # Compute daily ET₀
        print("[INFO] Computing daily ET₀...")
        
        ds["et0"] = calculate_reference_et(
            tmin=ds["Temperature_Air_2m_Min_24h"],
            tmax=ds["Temperature_Air_2m_Max_24h"],
            tmean=ds["Temperature_Air_2m_Mean_24h"],
            tdew=ds["Dew_Point_Temperature_2m_Mean"],
            ws10m=ds["Wind_Speed_10m_Mean"],
            srad=ds["Solar_Radiation_Flux"],
            elevation=ds["elev"],
            lat=ds.coords["lat"],
            doy=ds["time"].dt.dayofyear,
            config=config
        )
        
        # Apply mask for missing input data
        mask = (
            ds["Temperature_Air_2m_Min_24h"].notnull() &
            ds["Temperature_Air_2m_Max_24h"].notnull() &
            ds["Temperature_Air_2m_Mean_24h"].notnull() &
            ds["Dew_Point_Temperature_2m_Mean"].notnull() &
            ds["Wind_Speed_10m_Mean"].notnull() &
            ds["Solar_Radiation_Flux"].notnull() &
            ds["elev"].notnull()
        )
        ds["et0"] = ds["et0"].where(mask)
        
        # Add ET₀ to climatology
        ds_clim["et0"] = ds["et0"].sel(time=ds_clim.time)
        
        # Save datasets
        print("[INFO] Saving merged datasets...")
        ds.to_netcdf(merged_path)
        ds_clim.to_netcdf(clim_path)
    
    # Calculate EDDI for each time scale
    ds_out = xr.Dataset()
    for scale in time_scales:
        eddi_map = compute_eddi_for_timescale(
            ds_full_et0=ds["et0"],
            climatology_et0=ds_clim["et0"],
            scale_in_months=scale,
            end_date=user_end_date,
        )
        var_name = f"EDDI_{scale}month"
        ds_out[var_name] = eddi_map
    
    # Add metadata
    ds_out.attrs.update({
        "Conventions": "CF-1.8",
        "title": "Evaporative Demand Drought Index (EDDI)",
        "institution": "GOST/DEC Data Group/The World Bank",
        "source": f"{config.description} Reference ET₀ from AgERA5 data",
        "history": f"Created on {datetime.datetime.now(datetime.UTC).isoformat()}",
        "references": "https://psl.noaa.gov/eddi/",
        "et0_method": config.method,
        "et0_reference_type": config.ref_type,
        "et0_product": config.product_id
    })
    
    for scale in time_scales:
        var_name = f"EDDI_{scale}month"
        ds_out[var_name].attrs.update({
            "long_name": f"EDDI {scale}-month percentile",
            "standard_name": "evaporative_demand_drought_index",
            "units": "%",
            "description": f"Rank-based percentile of {scale}-month aggregated ET₀ vs. 1991-2020 climatology"
        })
    
    # Save output files
    for scale in time_scales:
        var_name = f"EDDI_{scale}month"
        output_path = f"{main_folder}/output/idn_cli_agera5_eddi_{scale}month_{end_date_str}.nc"
        print(f"[INFO] Saving {var_name} to {output_path}")
        ds_out[[var_name]].to_netcdf(output_path, format='NETCDF4', engine='netcdf4')
    
    print("[INFO] EDDI processing complete.")
    return ds_out


# =======================
# Visualization
# =======================
def plot_eddi_and_inputs(
        ds_full: xr.Dataset,
        ds_clim: xr.Dataset,
        ds_eddi: xr.Dataset,
        scale: int,
        user_start_date: datetime.date,
        user_end_date: datetime.date,
        main_folder: str
    ) -> None:
    """
    Create a 3x3 grid of maps visualizing EDDI and input variables.
    
    Subplot layout:
        Row 1: Tavg, Tmin, Tmax
        Row 2: Tdew, Wind Speed, Solar Radiation
        Row 3: ET₀ Sum, Climatology ET₀, EDDI
    
    Parameters
    ----------
    ds_full : xarray.Dataset
        Full dataset with meteorological variables and ET₀.
    ds_clim : xarray.Dataset
        Climatology dataset (1991-2020).
    ds_eddi : xarray.Dataset
        Dataset containing EDDI variables.
    scale : int
        Time scale in months.
    user_start_date : datetime.date
        Start date of the analysis window.
    user_end_date : datetime.date
        End date of the analysis window.
    main_folder : str
        Path to the main project folder.
    """
    # EDDI classification bins and labels
    my_bins = [0, 2, 5, 10, 20, 30, 70, 80, 90, 95, 98, 101]
    my_labels = [
        "EW4 (Exceptionally Wet)",
        "EW3 (Extremely Wet)",
        "EW2 (Severely Wet)",
        "EW1 (Moderately Wet)",
        "EW0 (Abnormally Wet)",
        "Normal",
        "ED0 (Abnormally Dry)",
        "ED1 (Moderate Drought)",
        "ED2 (Severe Drought)",
        "ED3 (Extreme Drought)",
        "ED4 (Exceptional Drought)"
    ]
    
    ncat = len(my_bins) - 1
    cmap = plt.get_cmap("RdBu_r", ncat)
    norm = mcolors.BoundaryNorm(my_bins, ncat)
    
    ymd_str = f"{user_end_date.year}{user_end_date.month:02d}{user_end_date.day:02d}"
    in_slice = slice(str(user_start_date), str(user_end_date))
    
    # Calculate means over the window
    tavg_mean = (ds_full["Temperature_Air_2m_Mean_24h"].sel(time=in_slice) - 273.15).mean("time")
    tmin_mean = (ds_full["Temperature_Air_2m_Min_24h"].sel(time=in_slice) - 273.15).mean("time")
    tmax_mean = (ds_full["Temperature_Air_2m_Max_24h"].sel(time=in_slice) - 273.15).mean("time")
    tdew_mean = (ds_full["Dew_Point_Temperature_2m_Mean"].sel(time=in_slice) - 273.15).mean("time")
    wind_mean = ds_full["Wind_Speed_10m_Mean"].sel(time=in_slice).mean("time")
    srad_mean = ds_full["Solar_Radiation_Flux"].sel(time=in_slice).mean("time")
    
    # ET₀ sum
    et0_sum = ds_full["et0"].sel(time=in_slice).sum("time")
    
    # Climatology ET₀
    def sum_et0_one_year(ds_sub, yr, sd, ed):
        st = sd.replace(year=yr)
        en = ed.replace(year=yr)
        sel_da = ds_sub["et0"].sel(time=slice(str(st), str(en)))
        return sel_da.sum("time")
    
    years_in_clim = range(1991, 2021)
    clim_sums = []
    for y in years_in_clim:
        try:
            clim_sums.append(sum_et0_one_year(ds_clim, y, user_start_date, user_end_date))
        except ValueError:
            pass
    
    if len(clim_sums) > 0:
        climatology_et0 = xr.concat(clim_sums, dim="year").mean("year")
    else:
        climatology_et0 = xr.zeros_like(et0_sum) * np.nan
    
    # EDDI map
    eddi_var = f"EDDI_{scale}month"
    eddi_map = ds_eddi[eddi_var]
    
    # Create figure
    fig, axs = plt.subplots(nrows=3, ncols=3, figsize=(20, 14))
    axlist = axs.ravel()
    
    def plot_map(data, ax, title, cmap_name="viridis", units=None, vmin=None, vmax=None,
                 is_discrete=False, disc_norm=None, disc_labels=None):
        if is_discrete and disc_norm is not None:
            img = data.plot.imshow(
                ax=ax, cmap=cmap_name, norm=disc_norm,
                add_colorbar=False, interpolation="bilinear"
            )
            cb = fig.colorbar(img, ax=ax, spacing="proportional", shrink=0.8)
            bin_edges = disc_norm.boundaries
            tick_locs = [(bin_edges[i] + bin_edges[i+1]) / 2 for i in range(len(bin_edges)-1)]
            cb.set_ticks(tick_locs)
            if disc_labels:
                cb.set_ticklabels(disc_labels)
        else:
            img = data.plot.imshow(
                ax=ax, cmap=cmap_name, add_colorbar=False,
                interpolation="bilinear", vmin=vmin, vmax=vmax
            )
            cb = fig.colorbar(img, ax=ax, orientation="vertical", shrink=0.8)
            if units:
                cb.set_label(units)
        ax.set_title(title, fontsize=11)
    
    # Temperature range
    t_vmin, t_vmax = 15.0, 35.0
    
    # Plot all maps
    plot_map(tavg_mean, axlist[0], "Tavg (°C)", cmap_name="RdYlBu_r", units="°C", vmin=t_vmin, vmax=t_vmax)
    plot_map(tmin_mean, axlist[1], "Tmin (°C)", cmap_name="RdYlBu_r", units="°C", vmin=t_vmin, vmax=t_vmax)
    plot_map(tmax_mean, axlist[2], "Tmax (°C)", cmap_name="RdYlBu_r", units="°C", vmin=t_vmin, vmax=t_vmax)
    plot_map(tdew_mean, axlist[3], "Tdew (°C)", cmap_name="RdYlBu_r", units="°C", vmin=t_vmin, vmax=t_vmax)
    plot_map(wind_mean, axlist[4], "Wind Speed (m/s)", cmap_name="YlGnBu", units="m/s")
    plot_map(srad_mean, axlist[5], "Solar Rad (J m⁻² day⁻¹)", cmap_name="inferno", units="J/m²/day")
    plot_map(et0_sum, axlist[6], "ET₀ Sum (mm)", cmap_name="YlOrBr", units="mm")
    plot_map(climatology_et0, axlist[7], "Climatology ET₀ (mm)", cmap_name="YlOrBr", units="mm")
    plot_map(eddi_map, axlist[8], f"EDDI {scale}-month", cmap_name=cmap, units=None,
             is_discrete=True, disc_norm=norm, disc_labels=my_labels)
    
    # Overlay boundaries
    boundary_path = f"{main_folder}/bnd/ne_10m_admin_0_countries.shp"
    try:
        world = gpd.read_file(boundary_path)
        for ax_i, data_arr in zip(axlist, [tavg_mean, tmin_mean, tmax_mean, tdew_mean,
                                            wind_mean, srad_mean, et0_sum, climatology_et0, eddi_map]):
            lon_min = float(data_arr.coords['lon'].min())
            lon_max = float(data_arr.coords['lon'].max())
            lat_min = float(data_arr.coords['lat'].min())
            lat_max = float(data_arr.coords['lat'].max())
            ax_i.set_xlim(lon_min, lon_max)
            ax_i.set_ylim(lat_min, lat_max)
            world.plot(ax=ax_i, color='none', edgecolor='black', linewidth=0.8)
    except Exception as e:
        print(f"[WARNING] Could not load boundary shapefile: {e}")
    
    fig.suptitle(f"EDDI and Input Variables - {scale}-month as of {ymd_str}", fontsize=16, y=0.93)
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    
    out_png = f"{main_folder}/images/idn_cli_agera5_eddi_{scale}month_{ymd_str}.png"
    print(f"[INFO] Saving figure to {out_png}")
    plt.savefig(out_png, dpi=300)
    plt.show()


# =======================
# Main Entry Point
# =======================
if __name__ == "__main__":
    # Configuration
    main_folder = '/mnt/d/temp/eddi'
    year = 2000
    month = 6
    dekad = 3
    time_scales = [1, 2, 3, 6, 9, 12]
    
    # Select ET₀ method:
    # 1 = ASCE Standardized tall reference crop (alfalfa, ETrs)
    # 2 = ASCE Standardized short reference crop (grass, ETos)
    # 3 = FAO-56 Penman-Monteith (grass, ETos)
    et_option = 1
    
    # Process EDDI
    ds_result = process_climate_data(main_folder, year, month, dekad, time_scales, et_option)
    print("[INFO] Final EDDI dataset:")
    print(ds_result)
    
    # Get dates for visualization
    user_start_date, user_end_date = get_dekad_dates(year, month, dekad)
    
    # Load datasets for visualization
    ds_full_path = f"{main_folder}/temp/idn_cli_agera5_sumba_full_1981_2024.nc"
    ds_clim_path = f"{main_folder}/climatology/idn_cli_agera5_sumba_clim_1991_2020.nc"
    ds_full = xr.open_dataset(ds_full_path)
    ds_clim = xr.open_dataset(ds_clim_path)
    
    # Create plots for each time scale
    for scale in time_scales:
        this_start_date = shift_date_by_months(user_end_date, scale)
        plot_eddi_and_inputs(
            ds_full=ds_full,
            ds_clim=ds_clim,
            ds_eddi=ds_result,
            scale=scale,
            user_start_date=this_start_date,
            user_end_date=user_end_date,
            main_folder=main_folder
        )

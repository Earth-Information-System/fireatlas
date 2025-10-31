import pytest 
from fireatlas.FireTime import aprox_local_datetime
import datetime as dt
import pandas as pd 


"""
Test cases: 

Single values 
Pandas series 

On UTC line (no change)
Positive change 
negative change 

After that: run Ofunato again to check? 

"""

def test_aprox_local_scalar_values_positive_longitude():
        """Test with single datetime and positive longitude."""
        utc_time = dt.datetime(2025, 1, 1, 12, 0, 0)  # Noon UTC
        lon = 15.0  # Should add 1 hour
        
        result = aprox_local_datetime(utc_time, lon)
        expected = dt.datetime(2025, 1, 1, 13, 0, 0)
        
        assert result == expected

def test_aprox_local_scalar_values_negative_longitude():
    """Test with single datetime and negative longitude (western hemisphere)."""
    utc_time = dt.datetime(2025, 1, 1, 12, 0, 0)  # Noon UTC
    lon = -75.0  # Should subtract 5 hours
    
    result = aprox_local_datetime(utc_time, lon)
    expected = dt.datetime(2025, 1, 1, 7, 0, 0)
    
    assert result == expected

def test_aprox_local_scalar_values_zero_longitude():
        """Test with longitude = 0 (Prime Meridian)."""
        utc_time = dt.datetime(2025, 1, 1, 12, 0, 0)
        lon = 0.0
        
        result = aprox_local_datetime(utc_time, lon)
        
        assert result == utc_time

def test_aprox_local_crossing_dateline_forward():
       """Test with different expected day"""
       utc = dt.datetime(2025, 1, 1, 22, 0, 0)
       lon = 180.0
       res = aprox_local_datetime(utc, lon)
       expected = dt.datetime(2025, 1, 2, 10, 0, 0)
       assert res == expected

def test_aprox_local_crossing_dateline_backward():
       """Test with different expected day"""
       utc = dt.datetime(2025, 1, 1, 10, 0, 0)
       lon = -180.0
       res = aprox_local_datetime(utc, lon)
       expected = dt.datetime(2024, 12, 31, 22, 0, 0)
       assert res == expected

def test_aprox_local_series():
       """Test with series inputs"""
       utc = pd.Series([
              dt.datetime(2025, 1, 1, 12, 0, 0), 
              dt.datetime(2025, 1, 1, 12, 0, 0),
              dt.datetime(2025, 1, 1, 12, 0, 0)
       ])

       lons = pd.Series([15.0, -30.0, 45.0])

       res = aprox_local_datetime(utc, lons)
       expected = pd.Series([
              dt.datetime(2025, 1, 1, 13, 0 ,0), # +1 hour 
              dt.datetime(2025, 1, 1, 10, 0, 0), # -2 hours
              dt.datetime(2025, 1, 1, 15, 0, 0) # + 3
       ])

       pd.testing.assert_series_equal(res, expected)


import pytest 
from fireatlas.FireTime import dt2t
import datetime as dt

def test_dt2t_simple():
       """Test with old version (00:00:00 and 12:00:00 only)"""
       am = dt.datetime(2025, 1, 1, 0, 0, 0)
       pm = dt.datetime(2025, 1, 1, 12, 0, 0)

       assert dt2t(am) == [2025, 1, 1, "AM"]
       assert dt2t(pm) == [2025, 1, 1, "PM"]

def test_dt2t():
       
       assert dt2t(dt.datetime(2025, 1, 1, 7, 0, 0)) == [2025, 1, 1, "PM"]
       assert dt2t(dt.datetime(2025, 1, 1, 17, 59, 59)) == [2025, 1, 1, "PM"]
       assert dt2t(dt.datetime(2025, 1, 1, 18, 0, 0 )) == [2025, 1, 2, "AM"]
       assert dt2t(dt.datetime(2025, 1, 2, 6, 59, 59)) == [2025, 1, 2, "AM"]
       assert dt2t(dt.datetime(2025, 1, 2, 7, 0, 0)) == [2025, 1, 2, "PM"]
       assert dt2t(dt.datetime(2025, 1, 2, 17, 59, 59)) == [2025, 1, 2, "PM"]
       assert dt2t(dt.datetime(2025, 1, 2, 18, 0, 0)) == [2025, 1, 3, "AM"]
       
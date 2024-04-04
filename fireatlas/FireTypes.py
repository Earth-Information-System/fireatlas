from typing import Dict, Tuple, Literal

Bbox = Tuple[float, float, float, float]
Region = Tuple[str, Bbox]
AMPM = Literal["AM", "PM"]
TimeStep = Tuple[int, int, int, AMPM]
Location = Literal["s3", "local"]

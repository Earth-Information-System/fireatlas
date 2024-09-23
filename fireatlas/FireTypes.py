from typing import Dict, Tuple, Literal, Optional

Bbox = Tuple[float, float, float, float]
Region = Tuple[str, Bbox]
AMPM = Literal["AM", "PM"]
TimeStep = Tuple[int, int, int, AMPM]
Location = Optional[Literal["s3", "local"]]
Timeout_param = Optional[Literal[int, None]]

from typing import Optional

import pydantic


class RuntimeSettings(pydantic.BaseSettings):

    number_of_multi_proc_workers: int = 4

    epsg: int = 9311

    debug: bool = False

    ftyp_opt: int = 1

    cont_opt: int = 1

    optional_example_variable: Optional[str] = None

    class Config:
        env_prefix = "FIRE_"
        env_file = ".env"
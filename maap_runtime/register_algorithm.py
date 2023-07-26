#!/usr/bin/env python
from maap.maap import MAAP
import argparse
from yaml import load as yaml_load, dump as yaml_dump

try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper
import json

KEY_MAP = {
    "algo_name": "algorithm_name",
    "version": "code_version",
    "environment": "environment_name",
    "description": "algorithm_description",
    "docker_url": "docker_container_url",
    "inputs": "algorithm_params",
    "run_command": "script_command",
    "repository_url": "repo_url",
}


def register(config, maap_api_host):
    maap = MAAP(maap_host=maap_api_host)
    resp = maap.registerAlgorithm(json.dumps(config))
    print(resp.text)
    pass


def main(algo_config, maap_api_host):
    with open(algo_config, "r") as fr:
        config = yaml_load(fr, Loader=Loader)
        output_config = {}
        for key, value in config.items():
            if key in KEY_MAP:
                if key == "inputs":
                    inputs = []
                    for argument in value:
                        inputs.append(
                            {
                                "field": argument.get("name"),
                                "type": argument.get("type"),
                            }
                        )
                    output_config.update({"algorithm_params": inputs})
                else:
                    output_config.update({KEY_MAP.get(key): value})
            else:
                output_config.update({key: value})
        print(json.dumps(output_config, indent=2))
        register(output_config, maap_api_host)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Register Algorithm to MAS")
    parser.add_argument(
        "config_file",
        metavar="config",
        nargs=1,
        help="Path to ./maap_runtime/<algorithm-dir-name>/algorithm_config.yaml",
    )
    parser.add_argument(
        "--maap_api_host",
        help="Optional maap api host url",
        default="api.maap-project.org",
    )
    args = parser.parse_args()
    main(args.config_file[0], args.maap_api_host)

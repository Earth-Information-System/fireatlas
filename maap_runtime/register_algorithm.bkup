from maap.maap import MAAP
import yaml
import json

maap = MAAP(maap_host="api.ops.maap-project.org")

with open("./maap_registry_yamls/algorithm_config.yaml", "r") as f:
    yml = yaml.safe_load(f)

info = {
    "script_command" : yml["run_command"],
    "algorithm_name" : yml["algo_name"],
    "code_version": yml["version"],
    "algorithm_description" : yml["description"],
    "environment_name": yml["environment"],
    "docker_container_url": yml["docker_url"],
    "repo_url" : yml["repository_url"],
    "disk_space": yml["disk_space"],
    "queue": yml["queue"],
    "algorithm_params" : [
        yml["inputs"]
    ]
}

version = yml["version"]
alg_id = f"{yml['algo_name']}_{yml['environment']}:{yml['version']}"
adel = maap.deleteAlgorithm(alg_id)
print(adel.text)

response = maap.registerAlgorithm(json.dumps(info))
print(response.text)

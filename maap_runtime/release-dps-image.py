from maap.maap import MAAP
import os

maap = MAAP(maap_host='api.maap-project.org')
current_dir = os.getcwd()
coordinator_config = f'{current_dir}/maap_runtime/coordinator/algorithm_config.yaml'
response = maap.register_algorithm_from_yaml_file(coordinator_config)
archive_config = f'{current_dir}/maap_runtime/archive/algorithm_config.yaml'
response = maap.register_algorithm_from_yaml_file(archive_config)

print(response.text)
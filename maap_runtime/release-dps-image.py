from maap.maap import MAAP
import os

maap = MAAP(maap_host='api.maap-project.org')
current_dir = os.getcwd()
target_config = f'{current_dir}/maap_runtime/coordinator/algorithm_config.yaml'
response = maap.register_algorithm_from_yaml_file(target_config)
print(response.text)
import yaml

def get_yaml_data(yaml_file: str) -> object:
    with open(yaml_file, 'r') as f:
        parsed_yaml = yaml.safe_load(f)
    return parsed_yaml


def write_yaml_data(data: object, out_yaml: str, sort_keys:bool = True) -> None:
    with open(out_yaml, 'w') as yml:
        yaml.dump(data, yml, default_flow_style=False, sort_keys=sort_keys)
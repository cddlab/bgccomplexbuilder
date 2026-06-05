import argparse
import json
from pathlib import Path

METRIC_KEYS = {
    "ipSAE": ("max", "ipSAE"),
    "ipSAE_min": ("min", "ipSAE"),
    "ipTM": ("max", "ipTM_af"),
    "pDockQ": ("max", "pDockQ"),
    "pDockQ2": ("max", "pDockQ2"),
    "LIS": ("max", "LIS"),
}


def extract_metrics(json_file: Path) -> tuple[str, dict[str, int | float]]:
    with open(json_file) as f:
        data = json.load(f)

    model_name = next(iter(data))
    ab = data[model_name]["A-B"]

    return model_name, {
        key: ab[section][field] for key, (section, field) in METRIC_KEYS.items()
    }


def main():
    parser = argparse.ArgumentParser(
        description="Extract complex metrics from AlphaFold2 JSON files."
    )
    parser.add_argument(
        "-i", "--input", required=True, help="Input directory (e.g. BGC0000003)"
    )
    args = parser.parse_args()

    input_dir = Path(args.input)
    if not input_dir.is_dir():
        raise SystemExit(f"Error: '{input_dir}' is not a directory.")

    results: dict[str, list[dict]] = {key: [] for key in METRIC_KEYS}

    for json_file in sorted(input_dir.rglob("*_model_10_10.json")):
        model_name, values = extract_metrics(json_file)
        for key, value in values.items():
            results[key].append({model_name: value})

    for key in results:
        results[key].sort(key=lambda x: list(x.values())[0], reverse=True)

    output_path = input_dir / "complexmetrics.json"
    with open(output_path, "w") as f:
        json.dump(results, f, indent=2)

    print(
        f"Written: {output_path}  ({sum(len(v) for v in results.values()) // len(results)} entries per metric)"
    )


if __name__ == "__main__":
    main()

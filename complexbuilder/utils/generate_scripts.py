import complexbuilder.utils.machinetype as mt


def generate_colabfold_search_runner(type: str, inputname: str, outputdir: str) -> str:
    """Generate bash scripts to submit jobs to the cluster."""
    if type == "yayoi":
        return mt.yayoi.colabfold_search(inputname, outputdir)
    elif type == "flow":
        return mt.flow.colabfold_search(inputname, outputdir)
    else:
        raise ValueError(f"Unknown machine type: {type}")

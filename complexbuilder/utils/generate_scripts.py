import complexbuilder.utils.machinetype as mt


def generate_scripts(type: str, inputname: str, outputdir: str) -> str:
    """Generate bash scripts to submit jobs to the cluster."""
    if type == "yayoi":
        return mt.yayoi.format(inputname=inputname, outputdir=outputdir)
    elif type == "foodin":
        return mt.foodin.format(inputname=inputname, outputdir=outputdir)
    elif type == "flow":
        return mt.flow.format(inputname=inputname, outputdir=outputdir)
    else:
        raise ValueError(f"Unknown machine type: {type}")

from pathlib import Path

import fancylog

import oscar_colony

output_dir = "tmp/logs/"


def init_logging(output_dir="tmp/logs/", verbose=True):
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    fancylog.start_logging(
        output_dir,
        oscar_colony,
        verbose=verbose,
        timestamp=True,
    )

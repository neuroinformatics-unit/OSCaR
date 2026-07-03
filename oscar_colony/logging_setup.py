from pathlib import Path

import fancylog

import oscar_colony

output_dir = "tmp/logs/"


def init_logging(
    output_dir="tmp/logs/", filename=None, verbose=True, log_to_console=False
):
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    fancylog.start_logging(
        output_dir,
        oscar_colony,
        verbose=verbose,
        timestamp=True,
        filename=filename,
        log_to_console=log_to_console,
    )

from pathlib import Path

import fancylog

import oscar_colony  # your package, used for git-hash/version logging


def init_logging(output_dir="tmp/logs/", verbose=True):
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    fancylog.start_logging(
        output_dir,
        oscar_colony,
        verbose=verbose,
        timestamp=True,
        # logger_name left as None -> configures the root logger
    )

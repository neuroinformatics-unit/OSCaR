from pathlib import Path

import fancylog

import oscar_colony


def init_logging(
    output_dir: str | Path = "tmp/logs/",
    filename: str | None = None,
    verbose: bool = True,
    log_to_console: bool = False,
) -> None:
    """Initialise logging for oscar_colony, via fancylog.

    Parameters
    ----------
    output_dir : str | Path
        Directory to write the log file to. Created if it doesn't exist.
    filename : str | None
        Name of the log file, without extension. If None, fancylog
        generates a name based on the current timestamp.
    verbose : bool
        Whether to log at DEBUG level (True) or INFO level (False).
    log_to_console : bool
        Whether to also print log output to the console.
    """
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    fancylog.start_logging(
        output_dir,
        oscar_colony,
        verbose=verbose,
        timestamp=True,
        filename=filename,
        log_to_console=log_to_console,
    )

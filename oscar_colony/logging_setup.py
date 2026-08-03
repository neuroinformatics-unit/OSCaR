from pathlib import Path

import fancylog

import oscar_colony


def init_logging(
    output_dir: str | Path = "tmp/logs/",
    filename: str | None = None,
    verbose: bool = True,
    file_log_level: str = "DEBUG",
    log_to_console: bool = False,
) -> None:
    """Initialise logging for oscar_colony, via fancylog.

    fancylog wraps the standard library `logging` module, additionally
    writing a header to the log file with the git commit, command-line
    arguments and Python version used for the run.

    Parameters
    ----------
    output_dir : str | Path
        Directory to write the log file to. Created if it doesn't exist.
    filename : str | None
        Name of the log file, without extension. If None, fancylog
        generates a name based on the current timestamp.
    verbose : bool
        Whether to log at DEBUG level (True) or INFO level (False).
    file_log_level : str
        Level of logging written to the log file, independent of the
        console log level set via `verbose`.
    log_to_console : bool
        Whether to also print log output to the console.
    """
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    fancylog.start_logging(
        output_dir,
        oscar_colony,
        verbose=verbose,
        file_log_level=file_log_level,
        timestamp=True,
        filename=filename,
        log_to_console=log_to_console,
    )

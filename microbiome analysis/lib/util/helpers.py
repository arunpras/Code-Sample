import os
import datetime
import click
import subprocess
import shutil
from enum import Enum
from contextlib import contextmanager


class MessageSeverity(Enum):
    MESSAGE = "MESSAGE"
    WARNING = "WARNING"
    HARD_ERROR = "HARD_ERROR"


def notify_app(msg, severity=MessageSeverity.MESSAGE, filepath=None):
    # TODO: consider using boto3 here instead
    if filepath:
        log_message(filepath, f"Sending a notification ({severity}): {msg}")
    if shutil.which("send_message"):
        subprocess.run(["send_message", f"--message \"{msg}\"", f"--severity {severity}"])
    else:
        click.echo("send_message command not found, skipping")


def log_message(filepath, msg):
    if not os.path.exists(filepath):
        with open(filepath, 'w') as f:
            fname = os.path.basename(filepath).replace(".log", ".py")
            f.write(f"{fname} log file\n\n")

    msg = f"{datetime.datetime.now().strftime('%Y-%m-%d:%H:%M:%S')} {msg}"
    click.echo(msg)
    with open(filepath, 'a') as f:
        f.write(f"{msg}\n")


@contextmanager
def log_start_end(filepath, name):
    log_message(filepath, f"{name} started.")
    yield
    log_message(filepath, f"{name} completed.")

#!/usr/bin/env python3
"""Execution-mode helpers for Ribo-seq Snakemake wrappers (native/conda/container)."""

__author__ = "Yangming Si"
__copyright__ = "Copyright 2026, Yangming Si"
__email__ = "siyangming1991@163.com"
__license__ = "MIT"

VALID_MODES = ("native", "conda", "container")


def build_container_command(image, bind, workdir, cmd):
    """Return argv for apptainer exec."""
    argv = [
        "apptainer",
        "exec",
        "--bind",
        bind,
        "--pwd",
        workdir,
        image,
    ]
    argv.extend(list(cmd))
    return argv


def container_run():
    return "apptainer exec --bind $(pwd):$(pwd) --pwd $(pwd) "


def _tool_cfg(config, tool_name):
    cfg = config.get(tool_name)
    if not isinstance(cfg, dict):
        raise ValueError(f"Missing tool config section: {tool_name}")
    return cfg


def _normalize_image(image):
    if not image:
        return ""
    if image.startswith(("docker://", "oras://", "library://", "http://", "https://", "shub://")):
        return image
    return f"docker://{image}"


def exec_wrapper_binary(config, tool_name, bin_key, default_bin):
    exec_mode = config.get("exec_mode", "conda")
    if exec_mode not in VALID_MODES:
        raise ValueError(
            f"Invalid exec_mode={exec_mode!r}; expected one of {VALID_MODES}"
        )
    tool = _tool_cfg(config, tool_name)

    if exec_mode == "container":
        image = _normalize_image(tool.get("container_image") or "")
        if not image:
            raise ValueError(
                f"Missing container_image under {tool_name} when exec_mode is container"
            )
        return f"{container_run()}{image} ", default_bin

    if exec_mode == "native":
        tool_bin = tool.get(bin_key) or ""
        if not tool_bin:
            raise ValueError(f"Missing {bin_key} in {tool_name} when exec_mode is native")
        return "", tool_bin

    return "", default_bin

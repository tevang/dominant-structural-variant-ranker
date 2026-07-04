#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import re
import subprocess

RUN_URL_RE = re.compile(
    r"^https://github\.com/(?P<owner>[^/]+)/(?P<repo>[^/]+)/actions/runs/(?P<run_id>\d+)(?:/.*)?$"
)


def _parse_run_target(value: str, repo_override: str | None) -> tuple[str, str]:
    match = RUN_URL_RE.match(value)
    if match:
        return f"{match.group('owner')}/{match.group('repo')}", match.group("run_id")

    if not value.isdigit():
        raise ValueError(
            "Run target must be a GitHub Actions run URL or a numeric run id."
        )
    if not repo_override:
        raise ValueError("Provide --repo when using a numeric run id.")
    return repo_override, value


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Inspect a failed GitHub Actions run and print failed-job logs.",
    )
    parser.add_argument(
        "run",
        help="GitHub Actions run URL or numeric run id.",
    )
    parser.add_argument(
        "--repo",
        help="Repository owner/name, required when passing a numeric run id.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)

    try:
        repo, run_id = _parse_run_target(args.run, args.repo)
    except ValueError as exc:
        parser.error(str(exc))

    gh = os.environ.get("GH", "gh")
    env = os.environ.copy()
    if "GH_TOKEN" not in env and "GITHUB_TOKEN" in env:
        env["GH_TOKEN"] = env["GITHUB_TOKEN"]

    command = [
        gh,
        "run",
        "view",
        run_id,
        "--repo",
        repo,
        "--log-failed",
    ]
    completed = subprocess.run(command, env=env)
    return completed.returncode


if __name__ == "__main__":
    raise SystemExit(main())

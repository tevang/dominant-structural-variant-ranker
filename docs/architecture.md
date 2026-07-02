# Architecture

DSVR is organized as a thin orchestration layer. Internal modules own typed
configuration, data models, file IO, subprocess execution, parsing, ranking, and
reporting. External chemistry engines remain outside the repository.


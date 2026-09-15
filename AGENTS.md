# AGENTS.md

- Follow `doc/StyleGuide.md` for all code you write or modify.
- Prioritize reusing existing capabilities:
    1. Search the `include/Physica` first. Reuse what you find there if it fits the requirement.
    2. Only write a manual implementation when nothing in the library meets the need. In that case, add a comment at the implementation site: `// TODO: <what Physica lacks that forced this manual implementation>`

## Constraints

- `find /` and any other full-disk or large-scale scanning commands must not be used:
    1. The location should preferably be obtained from build/link commands, build system caches, package managers, environment variables, etc.
    2. Only when all of the above have failed and the target directory is already known may a scope-limited search be performed within a specific directory.

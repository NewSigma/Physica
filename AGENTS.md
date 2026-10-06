# AGENTS.md

1. Follow `doc/StyleGuide.md` for all code you write or modify.
2. Prioritize reusing existing capabilities:
    - Search the `include/Physica` first. Reuse what you find there if it fits the requirement.
    - Only write a manual implementation when nothing in the library meets the need. In that case, add a comment at the implementation site: `// TODO: <what Physica lacks that forced this manual implementation>`

## Constraints

1. `find /` and any other full-disk or large-scale scanning commands must not be used:
    - The location should preferably be obtained from build system caches, package managers, IDE configurations, environment variables, etc.
    - Only when all of the above have failed and the target directory is already known may a scope-limited search be performed within a specific directory.

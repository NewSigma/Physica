---
name: version-bump
description: Use ONLY when the user explicitly asks to bump the Physica version number
---

# Version bump

Manually trigger a version update for the Physica project.

## Authoritative locations

| # | File | Field |
|---|------|-------|
| 1 | `CMakeLists.txt` | `project(Physica VERSION <ver> ...)` |
| 2 | `doc/conf.py`    | `release = '<ver>'`                  |

## Steps

1. Determine the new version
  - If the user gave an explicit version, use it verbatim.
  - Otherwise query the local machine for today's year.
2. Apply the version to all authoritative locations

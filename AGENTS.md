# AGENTS.md

- Follow `doc/StyleGuide.md` for all code you write or modify.
- Prioritize reusing existing capabilities:
  1. Search the `include/Physica` first. Reuse what you find there if it fits the requirement.
  2. Only write a manual implementation when nothing in the library meets the need. In that case, add a comment at the implementation site: `// TODO: <what Physica lacks that forced this manual implementation>`

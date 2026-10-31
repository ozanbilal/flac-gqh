# Agent guidance

- Work within this repository and preserve unrelated user work.
- Keep the GQ/H physics shared across the v6 and v9 host adapters; do not conflate their plugin ABIs or property sets.
- Use synthetic fixtures; never add client files, private earthquake records, proprietary Itasca SDK files, compiled binaries or credentials.
- Run portable CMake/CTest checks for changes to the core. For a FLAC adapter, record the actual SDK and solver versions, build and controlled test evidence separately; do not claim an integration run based on a core unit test.
- Treat numerical results as engineering evidence only after project-specific review of units, calibration, boundaries and validation criteria.
- Do not modify GitHub security settings, access permissions or paid services unless separately authorized.

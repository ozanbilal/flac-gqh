# Contributing

This project is licensed under the [MIT License](LICENSE).

For a proposed change, describe the problem, affected FLAC3D version and numerical conventions. Include small, synthetic tests and an independently derived numerical reference when applicable. Run the portable CMake/CTest checks and report results. Changes to a version-specific adapter also need a separate build and controlled test with the matching licensed SDK and solver.

Do not add proprietary Itasca SDK files or binaries, client data, non-redistributable motions, credentials, or private project logs. Do not upload existing project data or solver output without confirmed permission. Keep v6 and v9 ABI differences explicit and keep the constitutive physics shared unless a documented difference requires otherwise.

Report suspected security issues privately as described in [SECURITY.md](SECURITY.md); do not include secrets in an issue or pull request.

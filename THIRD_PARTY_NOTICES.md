# Third-party code

The AO metric diagnostics use the unmodified, header-only Eigen 3.4.0 library.
Its source is fetched from the versioned upstream archive and checked against
SHA-256 `8586084f71f9bde545ee7fa6d00288b264a2b7ac3607b974e54d13e7162c1c72`.

- Source: https://gitlab.com/libeigen/eigen/-/tree/3.4.0
- Archive: https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
- License: Mozilla Public License 2.0; relevant upstream notices are retained
  in `licenses/eigen-3.4.0/`. The `EIGEN_MPL2_ONLY` build definition prevents
  accidental inclusion of Eigen's optional LGPL components.
- COV does not modify Eigen. Eigen source and file notices remain available
  at the pinned source location and in CMake's fetched dependency directory.

GLFW and Dear ImGui retain their upstream source/license files in the respective
CMake dependency directories. Windows distribution packages must retain those
notices alongside this file and the COV license. NVIDIA driver/toolkit components
are separate installed dependencies; Gaussian is not distributed with COV.

#!/usr/bin/env bash
# Source this file before running the R experiments on a PACE (Georgia Tech) login shell:
#
#   source "scripts/local_scripts/setup_r_environment_pace.sh"
#   cd experiments_log
#   Rscript laplace.R 0.5 --cores=1 --trials=1 --burn=5
#
# Optional overrides before sourcing:
#   export R_MODULE_CANDIDATES="r/4.3.3-tidy r/4.3.3-bio r"
#   export SEQUENTIAL_FDP_R_LIBS_USER="/path/to/project-local/r_libs"

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
  echo "This script is meant to be sourced so it can update your current shell:"
  echo "  source \"$0\""
  exit 2
fi

_sequential_fdp_setup_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export SEQUENTIAL_FDP_AUDITING_ROOT="$(cd "${_sequential_fdp_setup_dir}/../.." && pwd)"

_sequential_fdp_has_command() {
  command -v "$1" >/dev/null 2>&1
}

_sequential_fdp_init_module_command() {
  if _sequential_fdp_has_command module; then
    return 0
  fi

  # PACE-specific Lmod initialization path.
  # shellcheck disable=SC1091
  [[ -r /usr/local/pace-apps/lmod/lmod/init/bash ]] && source /usr/local/pace-apps/lmod/lmod/init/bash

  # Fallback to common Lmod / Environment Modules locations.
  # shellcheck disable=SC1091
  [[ -r /etc/profile.d/modules.sh ]] && source /etc/profile.d/modules.sh
  # shellcheck disable=SC1091
  [[ -r /usr/share/Modules/init/bash ]] && source /usr/share/Modules/init/bash
  # shellcheck disable=SC1091
  [[ -r /usr/share/lmod/lmod/init/bash ]] && source /usr/share/lmod/lmod/init/bash

  _sequential_fdp_has_command module
}

_sequential_fdp_load_r_module() {
  if _sequential_fdp_has_command Rscript; then
    return 0
  fi

  if ! _sequential_fdp_init_module_command; then
    return 1
  fi

  # PACE R modules use lowercase "r/" with flavor suffixes (-tidy, -bio, -cuda, -geo).
  # Default preference order: tidy > bio > cuda > bare fallback.
  local candidate
  for candidate in ${R_MODULE_CANDIDATES:-r/4.3.3-tidy r/4.3.3-bio r/4.3.2-cuda r/4.2.1-tidy r/4.2.1-bio r}; do
    if module load "${candidate}" >/dev/null 2>&1; then
      echo "Loaded R module: ${candidate}"
      return 0
    fi
  done

  return 1
}

if ! _sequential_fdp_load_r_module; then
  echo "Rscript is not available yet."
  echo "If your cluster uses modules, try one of these first, then source this script again:"
  echo "  module spider r"
  echo "  export R_MODULE_CANDIDATES=\"r/<version>-<flavor> r\""
fi

# Wrap R and Rscript conditionally if they execute via apptainer.
# The host PACE modules aggressively configure GCC_ROOT/CC/CXX paths that leak
# into the apptainer container and override GCC's internal header searches,
# breaking package installation inside the container. Wait... wait,
# if we unset those variables, the container's GCC runs flawlessly.
if type Rscript 2>/dev/null | grep -q "apptainer"; then
  # Only wrap once if we source multiple times
  if ! type _orig_Rscript >/dev/null 2>&1; then
    eval "$(declare -f Rscript | sed '1s/^Rscript/_orig_Rscript/')"
    Rscript() {
      (
        unset CC CXX GCCROOT GCC_ROOT CPATH C_INCLUDE_PATH CPLUS_INCLUDE_PATH LIBRARY_PATH LD_LIBRARY_PATH
        _orig_Rscript "$@"
      )
    }
  fi
fi

if type R 2>/dev/null | grep -q "apptainer"; then
  if ! type _orig_R >/dev/null 2>&1; then
    eval "$(declare -f R | sed '1s/^R/_orig_R/')"
    R() {
      (
        unset CC CXX GCCROOT GCC_ROOT CPATH C_INCLUDE_PATH CPLUS_INCLUDE_PATH LIBRARY_PATH LD_LIBRARY_PATH
        _orig_R "$@"
      )
    }
  fi
fi

export SEQUENTIAL_FDP_R_LIBS_USER="${SEQUENTIAL_FDP_R_LIBS_USER:-${SEQUENTIAL_FDP_AUDITING_ROOT}/.r_libs}"
mkdir -p "${SEQUENTIAL_FDP_R_LIBS_USER}"
case ":${R_LIBS_USER:-}:" in
  *":${SEQUENTIAL_FDP_R_LIBS_USER}:"*) ;;
  *) export R_LIBS_USER="${SEQUENTIAL_FDP_R_LIBS_USER}${R_LIBS_USER:+:${R_LIBS_USER}}" ;;
esac

check_sequential_fdp_r_env() {
  if ! _sequential_fdp_has_command Rscript; then
    echo "Rscript was not found. Load an R module, then source this script again."
    return 1
  fi

  Rscript -e '
    pkgs <- c("KernSmooth", "rmutil", "pbmcapply")
    cat("R:", R.version.string, "\n")
    cat("Rscript:", Sys.which("Rscript"), "\n")
    cat("R_LIBS_USER:", Sys.getenv("R_LIBS_USER"), "\n")
    missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
    if (length(missing)) {
      cat("Missing packages:", paste(missing, collapse = ", "), "\n")
      quit(status = 1)
    }
    cat("Required packages are available:", paste(pkgs, collapse = ", "), "\n")
  '
}

install_sequential_fdp_r_packages() {
  if ! _sequential_fdp_has_command Rscript; then
    echo "Rscript was not found. Load an R module, then source this script again."
    return 1
  fi

  Rscript -e '
    dir.create(Sys.getenv("R_LIBS_USER"), recursive = TRUE, showWarnings = FALSE)
    pkgs <- c("KernSmooth", "rmutil", "pbmcapply")
    missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
    if (length(missing)) {
      install.packages(missing, lib = Sys.getenv("R_LIBS_USER"),
                       repos = "https://cloud.r-project.org")
    } else {
      cat("All required packages are already installed.\n")
    }
  '
}

run_laplace_smoke() {
  if ! check_sequential_fdp_r_env; then
    return 1
  fi

  (
    cd "${SEQUENTIAL_FDP_AUDITING_ROOT}/experiments_log" &&
      Rscript laplace.R "${1:-0.5}" --cores=1 --trials=1 --burn=5
  )
}

echo "Sequential FDP R environment configured (PACE)."
echo "Project root: ${SEQUENTIAL_FDP_AUDITING_ROOT}"
echo "R package library: ${SEQUENTIAL_FDP_R_LIBS_USER}"
echo "Check with: check_sequential_fdp_r_env"
echo "Optional package install: install_sequential_fdp_r_packages"
echo "Optional tiny laplace smoke test: run_laplace_smoke"

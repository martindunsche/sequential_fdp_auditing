#!/usr/bin/env bash
# Source this file before running the R experiments on a cluster login shell:
#
#   source "scripts/local scripts/setup_r_environment.sh"
#   cd experiments_log
#   Rscript laplace.R 0.5 --cores=1 --trials=1 --burn=5
#
# Optional overrides before sourcing:
#   export R_MODULE_CANDIDATES="R/4.3.2 R/4.2.3 R"
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

  # Common Environment Modules / Lmod initialization locations.
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

  local candidate
  for candidate in ${R_MODULE_CANDIDATES:-R r}; do
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
  echo "  module avail R"
  echo "  export R_MODULE_CANDIDATES=\"R/<version> R\""
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
    pkgs <- c("KernSmooth", "rmutil")
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

echo "Sequential FDP R environment configured."
echo "Project root: ${SEQUENTIAL_FDP_AUDITING_ROOT}"
echo "R package library: ${SEQUENTIAL_FDP_R_LIBS_USER}"
echo "Check with: check_sequential_fdp_r_env"
echo "Optional package install: install_sequential_fdp_r_packages"
echo "Optional tiny laplace smoke test: run_laplace_smoke"

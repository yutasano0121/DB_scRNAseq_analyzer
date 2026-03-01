#!/usr/bin/env bash
# install_launcher.sh — One-time setup for the scRNA-seq desktop shortcut.
#
# Run this script ONCE from a terminal:
#
#   bash /path/to/10Xapp/install_launcher.sh
#
# What it does:
#   1. Installs required system libraries (Linux: apt/dnf; macOS: Homebrew)
#   2. Creates a desktop shortcut (Linux only)
#
# After running this, double-click "scRNA-seq Analysis" on your Desktop
# to launch the app — no terminal or RStudio needed.
#
# Windows users: no system libraries are needed. Open RStudio and run
#   shiny::runApp("app.R") — R packages install automatically on first launch.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo ""
echo "=== scRNA-seq Analysis App — Setup ==="
echo ""

# ---------------------------------------------------------------------------
# 1. System library installation
# ---------------------------------------------------------------------------
# Several R packages (BPCells, sf/units, Cairo) need C/C++ libraries that
# must be installed via the OS package manager before R can compile them.
# On Windows these arrive pre-compiled in binary packages — nothing needed.

# Ask the user for confirmation before installing system packages.
# Usage: confirm_sys_install "apt" "pkg1 pkg2 ..."
# Returns 0 if the user agrees, 1 if they decline.
confirm_sys_install() {
  local manager="$1"
  local pkg_list="$2"

  echo ""
  echo "  The following system libraries are required to compile R packages"
  echo "  (BPCells, Cairo, sf/units, etc.) and will be installed via $manager:"
  echo ""
  for pkg in $pkg_list; do
    echo "    • $pkg"
  done
  echo ""
  printf "  Install these system libraries now? [y/N] "
  read -r answer
  case "$(echo "$answer" | tr '[:upper:]' '[:lower:]')" in
    y|yes) return 0 ;;
    *)
      echo ""
      echo "[sys] Skipped system library installation."
      echo "      Note: some R packages may fail to compile without these libraries."
      echo "      Re-run install_launcher.sh later to install them."
      echo ""
      return 1
      ;;
  esac
}

install_system_deps_debian() {
  local pkgs="libhdf5-dev libudunits2-dev libgdal-dev libgeos-dev libproj-dev \
libcairo2-dev libssl-dev libcurl4-openssl-dev libxml2-dev libfontconfig1-dev \
libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev"

  echo "[sys] Detected Debian/Ubuntu — system libraries needed via apt."
  if confirm_sys_install "apt" "$pkgs"; then
    sudo apt-get update -qq
    # shellcheck disable=SC2086
    sudo apt-get install -y $pkgs
    echo "[sys] System libraries installed."
  fi
}

install_system_deps_fedora() {
  local pkgs="hdf5-devel udunits2-devel gdal-devel geos-devel proj-devel \
cairo-devel openssl-devel libcurl-devel libxml2-devel fontconfig-devel \
harfbuzz-devel fribidi-devel freetype-devel libpng-devel libtiff-devel libjpeg-turbo-devel"

  echo "[sys] Detected Fedora/RHEL/CentOS — system libraries needed via dnf."
  if confirm_sys_install "dnf" "$pkgs"; then
    # shellcheck disable=SC2086
    sudo dnf install -y $pkgs
    echo "[sys] System libraries installed."
  fi
}

install_system_deps_macos() {
  local pkgs="hdf5 udunits gdal geos proj cairo openssl libxml2"

  echo "[sys] Detected macOS — system libraries needed via Homebrew."
  if ! command -v brew &>/dev/null; then
    echo ""
    echo "  Homebrew is not installed. Please install it first:"
    echo "    https://brew.sh"
    echo ""
    echo "  Then re-run this script."
    exit 1
  fi
  if confirm_sys_install "brew" "$pkgs"; then
    # shellcheck disable=SC2086
    brew install $pkgs
    echo "[sys] System libraries installed."
  fi
}

OS="$(uname -s)"

if [[ "$OS" == "Linux" ]]; then
  if command -v apt-get &>/dev/null; then
    install_system_deps_debian
  elif command -v dnf &>/dev/null; then
    install_system_deps_fedora
  elif command -v yum &>/dev/null; then
    # Older RHEL/CentOS — yum is a dnf alias on modern systems
    YUM_PKGS="hdf5-devel udunits2-devel gdal-devel geos-devel \
proj-devel cairo-devel openssl-devel libcurl-devel libxml2-devel"
    echo "[sys] Detected older RHEL/CentOS — system libraries needed via yum."
    if confirm_sys_install "yum" "$YUM_PKGS"; then
      # shellcheck disable=SC2086
      sudo yum install -y $YUM_PKGS
      echo "[sys] System libraries installed."
    fi
  else
    echo "[sys] WARNING: Could not detect package manager."
    echo "      Please install these libraries manually before continuing:"
    echo "        HDF5, udunits2, GDAL, GEOS, PROJ, Cairo, OpenSSL, libcurl, libxml2"
    echo ""
  fi
elif [[ "$OS" == "Darwin" ]]; then
  install_system_deps_macos
else
  echo "[sys] Windows detected or unknown OS — no system libraries needed."
fi

# ---------------------------------------------------------------------------
# 2. Validate launcher files
# ---------------------------------------------------------------------------

if [[ ! -f "$SCRIPT_DIR/launch.sh" ]]; then
  echo "ERROR: launch.sh not found in $SCRIPT_DIR"
  exit 1
fi

# ---------------------------------------------------------------------------
# 3. Desktop shortcut (Linux only)
# ---------------------------------------------------------------------------

if [[ "$OS" != "Linux" ]]; then
  echo ""
  if [[ "$OS" == "Darwin" ]]; then
    echo "macOS: no .desktop file needed."
    echo "To launch the app, run:  bash $SCRIPT_DIR/launch.sh"
    echo "Or open RStudio → Run App."
  fi
  echo ""
  echo "=== Setup complete! ==="
  echo "System libraries are installed. R packages will install on first launch."
  echo ""
  exit 0
fi

DESKTOP_FILE="$SCRIPT_DIR/10Xapp.desktop"
DESKTOP_DIR="$HOME/Desktop"
SHORTCUT="$DESKTOP_DIR/10Xapp.desktop"

if [[ ! -f "$DESKTOP_FILE" ]]; then
  echo "ERROR: 10Xapp.desktop not found in $SCRIPT_DIR"
  exit 1
fi

chmod +x "$SCRIPT_DIR/launch.sh"
echo "[1/3] Made launch.sh executable."

mkdir -p "$DESKTOP_DIR"

# Rewrite Exec= and Path= to the actual install location
sed -e "s|^Exec=.*|Exec=$SCRIPT_DIR/launch.sh|" \
    -e "s|^Path=.*|Path=$SCRIPT_DIR|" \
    "$DESKTOP_FILE" > "$SHORTCUT"
chmod +x "$SHORTCUT"
echo "[2/3] Desktop shortcut created at $SHORTCUT"

# Mark as trusted (required on Ubuntu/GNOME to enable double-click execution)
if command -v gio &>/dev/null; then
  gio set "$SHORTCUT" metadata::trusted true
  echo "[3/3] Marked shortcut as trusted."
elif command -v dbus-launch &>/dev/null; then
  dbus-launch gio set "$SHORTCUT" metadata::trusted true 2>/dev/null \
    && echo "[3/3] Marked shortcut as trusted." \
    || echo "[3/3] Could not mark as trusted — right-click the icon and choose 'Allow Launching'."
else
  echo "[3/3] Could not mark as trusted automatically."
  echo "      Right-click the icon on your Desktop and choose 'Allow Launching'."
fi

echo ""
echo "=== Setup complete! ==="
echo ""
echo "Double-click 'scRNA-seq Analysis' on your Desktop to start the app."
echo "The first launch will install required R packages (10–30 min, one time only)."
echo ""
echo "Tip: You can also search for 'scRNA-seq' in your application launcher."
echo ""

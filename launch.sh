#!/usr/bin/env bash
# launch.sh — Desktop launcher for the scRNA-seq Analysis app.
#
# Intended to be invoked by double-clicking the 10Xapp.desktop icon.
# No terminal window is shown during normal (packages-installed) use.
# On first launch, a terminal window opens for package installation,
# then closes automatically when done.
#
# One-time setup: run install_launcher.sh to create the Desktop shortcut.

set -euo pipefail

# ---------------------------------------------------------------------------
# 1. Resolve project root (the directory containing this script)
# ---------------------------------------------------------------------------

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOG_FILE="$SCRIPT_DIR/launch.log"

# ---------------------------------------------------------------------------
# 2. Locate Rscript
# ---------------------------------------------------------------------------

find_rscript() {
  # First: check if Rscript is already on PATH
  if command -v Rscript &>/dev/null; then
    command -v Rscript
    return 0
  fi
  # Second: check common install locations
  local candidates=(
    /usr/bin/Rscript
    /usr/local/bin/Rscript
    /usr/local/lib/R/bin/Rscript
    "$HOME/anaconda3/bin/Rscript"
    "$HOME/miniconda3/bin/Rscript"
    "$HOME/miniforge3/bin/Rscript"
    "$HOME/.local/bin/Rscript"
  )
  # Also glob POSIT/RStudio versioned installs: /opt/R/<version>/bin/Rscript
  for pattern in /opt/R/*/bin/Rscript; do
    candidates+=("$pattern")
  done

  for candidate in "${candidates[@]}"; do
    if [[ -x "$candidate" ]]; then
      echo "$candidate"
      return 0
    fi
  done
  return 1
}

RSCRIPT="$(find_rscript || true)"

if [[ -z "$RSCRIPT" ]]; then
  zenity --error \
    --title="scRNA-seq App — R Not Found" \
    --text="R (Rscript) could not be found on this system.\n\nPlease install R from:\nhttps://cran.r-project.org\n\nThen double-click the icon again." \
    --width=420 \
    2>/dev/null \
  || notify-send \
      --icon=dialog-error \
      "scRNA-seq App" \
      "R is not installed. Please install R from cran.r-project.org."
  exit 1
fi

cd "$SCRIPT_DIR"

# ---------------------------------------------------------------------------
# 3. Probe for missing packages (~1 second when all installed)
# ---------------------------------------------------------------------------
# Use if/else rather than capturing $? — commands inside an `if` condition
# are exempt from set -e, so a non-zero exit (packages missing) is safe here.

PACKAGES_MISSING=false
if ! "$RSCRIPT" --vanilla -e "
  source('R/utils.R')
  missing <- check_and_install_dependencies(auto_install = FALSE)
  quit(status = if (length(missing) == 0) 0 else 1)
" 2>/dev/null; then
  PACKAGES_MISSING=true
fi

# ---------------------------------------------------------------------------
# 4. First-run: show install UI if packages are missing
# ---------------------------------------------------------------------------

if [[ "$PACKAGES_MISSING" == "true" ]]; then
  # Inform user before opening install terminal
  zenity --info \
    --title="scRNA-seq App — First-Time Setup" \
    --text="Some required R packages are not yet installed.\n\nA terminal window will open showing the installation progress.\n\n\342\200\242 This only happens once.\n\342\200\242 It may take 10\342\200\22330 minutes.\n\342\200\242 Do not close the terminal window until it finishes." \
    --width=460 \
    2>/dev/null \
  || true  # Proceed even if zenity is unavailable

  # Build the install command to run inside the terminal
  INSTALL_CMD="
    echo '================================================'
    echo '  scRNA-seq App: Installing required packages'
    echo '================================================'
    echo ''
    cd '$SCRIPT_DIR'
    '$RSCRIPT' -e \"source('R/utils.R'); check_and_install_dependencies(auto_install = TRUE, prompt = FALSE)\"
    echo ''
    echo 'Setup complete. This window will close in 5 seconds...'
    sleep 5
  "

  # Open a terminal for the install; wait for it to finish before continuing.
  # Try gnome-terminal (with --wait), fall back to xterm (blocks by default).
  if command -v gnome-terminal &>/dev/null \
      && gnome-terminal --help 2>&1 | grep -q -- "--wait"; then
    gnome-terminal \
      --title="scRNA-seq App — Setup" \
      --wait \
      -- bash -c "$INSTALL_CMD"
  elif command -v xterm &>/dev/null; then
    xterm \
      -title "scRNA-seq App — Setup" \
      -fa "Monospace" -fs 11 \
      -e bash -c "$INSTALL_CMD"
  else
    # Last resort: run silently and log; user sees nothing but app opens when done
    bash -c "$INSTALL_CMD" >> "$LOG_FILE" 2>&1
  fi
fi

# ---------------------------------------------------------------------------
# 5. Launch the Shiny app (browser opens automatically)
# ---------------------------------------------------------------------------

# Rotate the log so it doesn't grow unboundedly across launches
if [[ -f "$LOG_FILE" ]]; then
  mv "$LOG_FILE" "${LOG_FILE}.prev" 2>/dev/null || true
fi

nohup "$RSCRIPT" -e "
  options(shiny.launch.browser = TRUE, shiny.port = 0)
  setwd('$SCRIPT_DIR')
  source('app.R')
" >> "$LOG_FILE" 2>&1 &
APP_PID=$!

# ---------------------------------------------------------------------------
# 6. Brief health check — catch immediate startup failures
# ---------------------------------------------------------------------------

sleep 4
if ! kill -0 "$APP_PID" 2>/dev/null; then
  zenity --error \
    --title="scRNA-seq App — Launch Failed" \
    --text="The app failed to start.\n\nPlease check the log file for details:\n$LOG_FILE" \
    --width=440 \
    2>/dev/null \
  || notify-send \
      --icon=dialog-error \
      "scRNA-seq App" \
      "Launch failed. See $LOG_FILE for details."
  exit 1
fi

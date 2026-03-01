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
# Helper: probe for missing packages
# ---------------------------------------------------------------------------
# Prints missing package names to stdout (one per line), exits 1 if any found.
# Log messages from R go to stderr and are suppressed by callers as needed.
probe_missing_packages() {
  "$RSCRIPT" --vanilla -e "
    source('R/utils.R')
    missing <- check_and_install_dependencies(auto_install = FALSE)
    if (length(missing) > 0) cat(paste(missing, collapse = '\n'))
    quit(status = if (length(missing) == 0) 0 else 1)
  " 2>/dev/null
}

# ---------------------------------------------------------------------------
# Helper: print fix hints for failed packages (terminal)
# ---------------------------------------------------------------------------
print_fix_hints() {
  echo "  Recommended fixes:"
  echo ""
  echo "  1. Missing system libraries (needed for BPCells, Cairo, sf/units):"
  echo "       bash \"$SCRIPT_DIR/install_launcher.sh\""
  echo ""
  echo "  2. Network error — check your connection, then re-run:"
  echo "       bash \"$SCRIPT_DIR/launch.sh\""
  echo ""
  echo "  3. Retry from inside the app:"
  echo "       Settings tab → \"Install Missing Packages\""
  echo ""
  echo "  4. See full error details:"
  echo "       $LOG_FILE"
}

# ---------------------------------------------------------------------------
# 3. Probe for missing packages (~1 second when all installed)
# ---------------------------------------------------------------------------

PACKAGES_MISSING=false
MISSING_LIST=""
if MISSING_LIST="$(probe_missing_packages)"; then
  : # all packages present — nothing to do
else
  PACKAGES_MISSING=true
fi

# ---------------------------------------------------------------------------
# 4. First-run: install missing packages if needed
# ---------------------------------------------------------------------------

if [[ "$PACKAGES_MISSING" == "true" ]]; then

  if [[ -t 0 ]]; then
    # -----------------------------------------------------------------------
    # Running from a terminal — show inline prompt (same experience as runApp)
    # -----------------------------------------------------------------------
    echo ""
    echo "------------------------------------------------------------"
    echo "The following R packages are not yet installed:"
    echo ""
    while IFS= read -r pkg; do
      [[ -n "$pkg" ]] && echo "  • $pkg"
    done <<< "$MISSING_LIST"
    echo ""
    echo "------------------------------------------------------------"
    printf "Install these packages now? [y/N] "
    read -r answer
    echo ""

    if [[ "$(echo "$answer" | tr '[:upper:]' '[:lower:]')" =~ ^(y|yes)$ ]]; then
      echo "Installing packages... (this may take 10–30 minutes on first run)"
      echo ""
      "$RSCRIPT" -e "
        source('R/utils.R')
        check_and_install_dependencies(auto_install = TRUE, prompt = FALSE)
      "
      echo ""

      # Re-probe: check whether anything failed to install
      STILL_MISSING=""
      if ! STILL_MISSING="$(probe_missing_packages)"; then
        echo "------------------------------------------------------------"
        echo "WARNING: The following packages failed to install:"
        echo ""
        while IFS= read -r pkg; do
          [[ -n "$pkg" ]] && echo "  • $pkg"
        done <<< "$STILL_MISSING"
        echo ""
        print_fix_hints
        echo "------------------------------------------------------------"
        echo ""
        printf "Launch the app anyway? [y/N] "
        read -r launch_anyway
        echo ""
        if [[ ! "$(echo "$launch_anyway" | tr '[:upper:]' '[:lower:]')" =~ ^(y|yes)$ ]]; then
          echo "Exiting. Fix the issues above, then re-run: bash launch.sh"
          exit 1
        fi
        echo "Launching with missing packages. Some analysis steps may not work."
      else
        echo "All packages installed successfully. Launching app..."
      fi

    else
      echo "Skipped. The app will launch, but modules requiring these packages"
      echo "may fail at runtime:"
      echo ""
      while IFS= read -r pkg; do
        [[ -n "$pkg" ]] && echo "  • $pkg"
      done <<< "$MISSING_LIST"
      echo ""
      print_fix_hints
    fi

  else
    # -----------------------------------------------------------------------
    # Running from desktop icon — use GUI dialogs (no terminal available)
    # -----------------------------------------------------------------------
    zenity --info \
      --title="scRNA-seq App — First-Time Setup" \
      --text="Some required R packages are not yet installed.\n\nA terminal window will open showing the installation progress.\n\n\342\200\242 This only happens once.\n\342\200\242 It may take 10\342\200\22330 minutes.\n\342\200\242 Do not close the terminal window until it finishes." \
      --width=460 \
      2>/dev/null \
    || true  # Proceed even if zenity is unavailable

    INSTALL_CMD="
      echo '================================================'
      echo '  scRNA-seq App: Installing required packages'
      echo '================================================'
      echo ''
      cd '$SCRIPT_DIR'
      '$RSCRIPT' -e \"source('R/utils.R'); check_and_install_dependencies(auto_install = TRUE, prompt = FALSE)\"
      echo ''
      echo 'Done. This window will close in 5 seconds...'
      sleep 5
    "

    # Open a terminal for the install; wait for it to finish before continuing.
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
      # Last resort: run silently and log
      bash -c "$INSTALL_CMD" >> "$LOG_FILE" 2>&1
    fi

    # Re-probe after desktop install: show a warning dialog if anything failed
    STILL_MISSING=""
    if ! STILL_MISSING="$(probe_missing_packages)"; then
      FAILED_NAMES="$(echo "$STILL_MISSING" | tr '\n' ' ')"
      zenity --warning \
        --title="scRNA-seq App — Partial Install" \
        --text="Some packages failed to install:\n\n${FAILED_NAMES}\n\nCommon fixes:\n\342\200\242 Run install_launcher.sh to add missing system libraries\n\342\200\242 Check your internet connection\n\342\200\242 In the app: Settings \342\206\222 Install Missing Packages\n\342\200\242 See details: $LOG_FILE\n\nThe app will still launch; affected steps may not work." \
        --width=500 \
        2>/dev/null \
      || true
    fi

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
  shiny::runApp('$SCRIPT_DIR')
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
